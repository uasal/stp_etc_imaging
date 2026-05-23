"""Target-list exposure-time handoff utilities.

This module is intentionally independent from budgie. It consumes a YAML target
list schema and computes ETC integration times per target/classification.

Targets YAML schema:

```yaml
observatory: UM            # "UM" or "STP"; per-target override allowed
snr_k: 3                   # default detection threshold
n_bands: 1                 # default number of independent bands
frame_exp_time_s: 60       # default per-frame exposure used by calc_int_time
classifications:
  exoplanet_host_stars:
    - name: HD 12345
      V_mag: 5.2
      desired_contrast: 1.0e-7
      sep_arcsec: 0.5
      spectral_type: G2V
      notes: ""
      snr_k: 5             # optional per-target override
      n_bands: 2           # optional per-target override
  debris_disks: [...]
  prime_targets: [...]
```

Semantics:
- `SNR_required = snr_k / desired_contrast` (post-processed contrast floor convention).
- `n_bands` assumes each band must independently reach the contrast floor, so
  `t_exp_total = n_bands * t_exp_per_band`.

Spectrum file resolution:
- Use explicit `pickles_support_dir` when provided.
- Else use `$UASAL_ARCHIVE/astr_obj_models/stars/pickles_models/dat_uvk`.
- Else fall back to current working directory behavior.
"""

from __future__ import annotations

import csv
import importlib.metadata
import os
from pathlib import Path
from typing import Any

import astropy.units as u
import yaml

from .ExposureTimeSNRCalculator import Observatory

_REQUIRED_TARGET_FIELDS = ("name", "V_mag", "desired_contrast")
_PICKLES_SUBDIR = "astr_obj_models/stars/pickles_models/dat_uvk"

_SPECTRAL_TYPE_TO_PICKLES = {
    "O5V": "pickles_uk_1.fits",
    "O9V": "pickles_uk_2.fits",
    "B0V": "pickles_uk_3.fits",
    "B1V": "pickles_uk_4.fits",
    "B3V": "pickles_uk_5.fits",
    "B8V": "pickles_uk_7.fits",
    "A0V": "pickles_uk_9.fits",
    "A2V": "pickles_uk_10.fits",
    "A5V": "pickles_uk_12.fits",
    "F0V": "pickles_uk_14.fits",
    "F5V": "pickles_uk_16.fits",
    "F8V": "pickles_uk_20.fits",
    "G0V": "pickles_uk_23.fits",
    "G2V": "pickles_uk_26.fits",
    "G5V": "pickles_uk_27.fits",
    "K0V": "pickles_uk_31.fits",
    "K2V": "pickles_uk_33.fits",
    "K5V": "pickles_uk_36.fits",
    "M0V": "pickles_uk_38.fits",
    "M2V": "pickles_uk_40.fits",
    "M4V": "pickles_uk_43.fits",
    "M5V": "pickles_uk_44.fits",
}


def contrast_to_snr(contrast: float, k: float = 3.0) -> float:
    """Return required SNR from desired post-processed contrast floor."""
    contrast = float(contrast)
    k = float(k)
    if contrast <= 0:
        raise ValueError("desired_contrast must be > 0")
    if k <= 0:
        raise ValueError("snr_k must be > 0")
    return k / contrast


def _normalize_target(
    target: dict[str, Any],
    *,
    classification: str,
    defaults: dict[str, Any],
) -> dict[str, Any]:
    merged = dict(defaults)
    merged.update(target)
    merged["classification"] = classification

    missing = [field for field in _REQUIRED_TARGET_FIELDS if field not in merged]
    if missing:
        raise ValueError(
            f"Target in classification '{classification}' is missing required fields: {', '.join(missing)}"
        )

    try:
        merged["V_mag"] = float(merged["V_mag"])
        merged["desired_contrast"] = float(merged["desired_contrast"])
        merged["snr_k"] = float(merged.get("snr_k", defaults["snr_k"]))
        merged["n_bands"] = int(merged.get("n_bands", defaults["n_bands"]))
        merged["frame_exp_time_s"] = float(
            merged.get("frame_exp_time_s", defaults["frame_exp_time_s"])
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Target '{merged.get('name', '<unknown>')}' has invalid numeric values"
        ) from exc

    if merged["desired_contrast"] <= 0:
        raise ValueError(f"Target '{merged['name']}' desired_contrast must be > 0")
    if merged["snr_k"] <= 0:
        raise ValueError(f"Target '{merged['name']}' snr_k must be > 0")
    if merged["n_bands"] < 1:
        raise ValueError(f"Target '{merged['name']}' n_bands must be >= 1")
    if merged["frame_exp_time_s"] <= 0:
        raise ValueError(f"Target '{merged['name']}' frame_exp_time_s must be > 0")

    merged["spectral_type"] = str(merged.get("spectral_type", "G2V")).strip().upper()
    merged["notes"] = str(merged.get("notes", ""))
    merged["observatory"] = str(merged.get("observatory", defaults["observatory"]))
    return merged


def load_targets(yaml_path) -> dict:
    """Load and validate target-list YAML schema and return normalized data."""
    yaml_path = Path(yaml_path)
    with yaml_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}

    if not isinstance(raw, dict):
        raise ValueError("Target YAML must be a mapping")

    if "classifications" not in raw:
        raise ValueError("Target YAML must define 'classifications'")
    if not isinstance(raw["classifications"], dict):
        raise ValueError("'classifications' must be a mapping of classification -> list[targets]")

    defaults = {
        "observatory": str(raw.get("observatory", "UM")),
        "snr_k": float(raw.get("snr_k", 3.0)),
        "n_bands": int(raw.get("n_bands", 1)),
        "frame_exp_time_s": float(raw.get("frame_exp_time_s", 60.0)),
    }

    if defaults["n_bands"] < 1:
        raise ValueError("Top-level n_bands must be >= 1")
    if defaults["snr_k"] <= 0:
        raise ValueError("Top-level snr_k must be > 0")
    if defaults["frame_exp_time_s"] <= 0:
        raise ValueError("Top-level frame_exp_time_s must be > 0")

    normalized_classifications: dict[str, list[dict[str, Any]]] = {}
    for classification, targets in raw["classifications"].items():
        if not isinstance(targets, list):
            raise ValueError(f"Classification '{classification}' must contain a list of targets")
        normalized_targets = [
            _normalize_target(target, classification=classification, defaults=defaults)
            for target in targets
        ]
        normalized_classifications[classification] = normalized_targets

    return {
        **defaults,
        "classifications": normalized_classifications,
    }


def _resolve_pickles_file(spectral_type: str) -> str:
    key = spectral_type.strip().upper()
    return _SPECTRAL_TYPE_TO_PICKLES.get(key, "pickles_uk_26.fits")


def _resolve_pickles_support_dir(explicit: str | None = None) -> str | None:
    """Return the directory containing pickles_uk_*.fits."""
    if explicit:
        return str(Path(explicit))
    archive = os.environ.get("UASAL_ARCHIVE")
    if archive:
        return str(Path(archive) / _PICKLES_SUBDIR)
    return None


def _background_for_target(observatory: Observatory, target: dict[str, Any], support_data_dir: str | None):
    if "background_file" in target:
        return target["background_file"], target.get("background_support_data_path", support_data_dir)

    data_observatory = getattr(observatory, "data_observatory", {}) or {}
    background_file = (
        data_observatory.get("astrophysics", {})
        .get("zodi", {})
        .get("profile")
    )
    if background_file is None:
        raise ValueError(
            f"Target '{target['name']}' is missing background_file and observatory zodi profile is unavailable"
        )
    support_path = getattr(observatory, "support_data_observatory", support_data_dir)
    return background_file, support_path


def _quantity_to_seconds(value: Any) -> float:
    if hasattr(value, "to"):
        try:
            return float(value.to(u.s).value)
        except (AttributeError, TypeError, ValueError, u.UnitConversionError):
            pass
    if hasattr(value, "value"):
        return float(value.value)
    return float(value)


def exposure_time_for_target(
    target: dict,
    observatory_name: str = "UM",
    n_bands: int = 1,
    snr_k: float = 3.0,
    frame_exp_time_s: float = 60.0,
    custom_toml_dir: str | None = None,
    support_data_dir: str | None = None,
    pickles_support_dir: str | None = None,
) -> dict:
    """Compute ETC integration time for one target and return result summary."""
    merged_target = {
        "observatory": observatory_name,
        "n_bands": n_bands,
        "snr_k": snr_k,
        "frame_exp_time_s": frame_exp_time_s,
        **target,
    }
    merged_target = _normalize_target(
        merged_target,
        classification=str(merged_target.get("classification", "unclassified")),
        defaults={
            "observatory": observatory_name,
            "n_bands": n_bands,
            "snr_k": snr_k,
            "frame_exp_time_s": frame_exp_time_s,
        },
    )

    observatory = Observatory(merged_target["observatory"])
    observatory.make_STP(custom_toml_dir=custom_toml_dir, support_data_dir=support_data_dir)

    spectral_type = merged_target["spectral_type"]
    source_pickles_file = _resolve_pickles_file(spectral_type)
    resolved_pickles_support_dir = _resolve_pickles_support_dir(pickles_support_dir)
    if resolved_pickles_support_dir is not None:
        expected_pickles_path = Path(resolved_pickles_support_dir) / source_pickles_file
        if not expected_pickles_path.exists():
            raise FileNotFoundError(
                f"Missing Pickles spectrum file at '{expected_pickles_path}'. "
                "Set UASAL_ARCHIVE to a uasal_archive clone or pass pickles_support_dir explicitly."
            )

    observatory.set_source(
        source_pickles_file=source_pickles_file,
        source_z=float(merged_target.get("source_z", 0.0)),
        support_data_path=resolved_pickles_support_dir,
    )

    background_file, background_support_path = _background_for_target(
        observatory, merged_target, support_data_dir
    )
    observatory.set_background(
        background_file=background_file,
        support_data_path=background_support_path,
    )

    observatory.make_observation(flux=merged_target["V_mag"], flux_units=u.ABmag)

    snr_required = contrast_to_snr(merged_target["desired_contrast"], merged_target["snr_k"])
    t_exp_per_band = observatory.calc_int_time(
        snr=snr_required,
        exp_time=merged_target["frame_exp_time_s"] * u.s,
    )
    t_exp_per_band_s = _quantity_to_seconds(t_exp_per_band)
    t_exp_total_s = merged_target["n_bands"] * t_exp_per_band_s

    notes = merged_target["notes"]
    if spectral_type not in _SPECTRAL_TYPE_TO_PICKLES:
        suffix = "Spectral type not in lookup; defaulted to G2V template."
        notes = f"{notes} {suffix}".strip()

    return {
        "name": merged_target["name"],
        "classification": merged_target["classification"],
        "V_mag": merged_target["V_mag"],
        "desired_contrast": merged_target["desired_contrast"],
        "snr_required": snr_required,
        "snr_k": merged_target["snr_k"],
        "n_bands": merged_target["n_bands"],
        "frame_exp_time_s": merged_target["frame_exp_time_s"],
        "t_exp_per_band_s": t_exp_per_band_s,
        "t_exp_total_s": t_exp_total_s,
        "observatory": merged_target["observatory"],
        "notes": notes.strip(),
    }


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "classification",
        "name",
        "observatory",
        "V_mag",
        "desired_contrast",
        "snr_k",
        "snr_required",
        "n_bands",
        "frame_exp_time_s",
        "t_exp_per_band_s",
        "t_exp_total_s",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _format_markdown_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "_No targets_\n"

    lines = [
        "| Target | V_mag | Desired contrast | SNR required | n_bands | t_exp_per_band [s] | t_exp_total [s] |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {name} | {V_mag:.3f} | {desired_contrast:.3e} | {snr_required:.3e} | {n_bands} | {t_exp_per_band_s:.3f} | {t_exp_total_s:.3f} |".format(
                **row
            )
        )
    return "\n".join(lines) + "\n"


def _write_markdown_report(summary: dict[str, Any], path: Path) -> None:
    lines = [
        "# Exposure Time Budget Report",
        "",
        "## Assumptions",
        "- Contrast-to-SNR convention: `SNR_required = snr_k / desired_contrast` (default `snr_k=3`).",
        "- `n_bands` semantics: each band must independently reach the desired contrast, so `t_exp_total = n_bands * t_exp_per_band`.",
        "",
        "## Run Metadata",
        f"- Observatory default: `{summary['defaults']['observatory']}`",
        f"- snr_k default: `{summary['defaults']['snr_k']}`",
        f"- n_bands default: `{summary['defaults']['n_bands']}`",
        f"- frame_exp_time_s default: `{summary['defaults']['frame_exp_time_s']}`",
        f"- stp_etc_imaging version: `{summary['software_versions']['stp_etc_imaging']}`",
        f"- astropy version: `{summary['software_versions']['astropy']}`",
        "",
        "## Totals",
        f"- Total targets: **{summary['totals']['n_targets']}**",
        f"- Total exposure time [s]: **{summary['totals']['t_exp_total_s']:.3f}**",
        "",
    ]

    for classification, section in summary["classifications"].items():
        lines.extend(
            [
                f"## Classification: `{classification}`",
                f"- Targets: **{section['totals']['n_targets']}**",
                f"- Total exposure time [s]: **{section['totals']['t_exp_total_s']:.3f}**",
                "",
                _format_markdown_table(section["targets"]),
            ]
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def run_target_list(
    yaml_path,
    output_dir,
    observatory_name: str = "UM",
    snr_k_default: float = 3.0,
    n_bands_default: int = 1,
    pickles_support_dir: str | None = None,
) -> dict:
    """Run ETC calculations for a target-list YAML and emit report artifacts."""
    loaded = load_targets(yaml_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    defaults = {
        "observatory": loaded.get("observatory", observatory_name),
        "snr_k": loaded.get("snr_k", snr_k_default),
        "n_bands": loaded.get("n_bands", n_bands_default),
        "frame_exp_time_s": loaded.get("frame_exp_time_s", 60.0),
    }

    class_summaries: dict[str, Any] = {}
    all_rows: list[dict[str, Any]] = []
    total_exp = 0.0

    for classification, targets in loaded["classifications"].items():
        results = []
        class_total = 0.0
        for target in targets:
            target_input = dict(target)
            target_input["classification"] = classification
            target_pickles_support_dir = target_input.get("pickles_support_dir", pickles_support_dir)
            result = exposure_time_for_target(
                target=target_input,
                observatory_name=target_input.get("observatory", defaults["observatory"]),
                n_bands=int(target_input.get("n_bands", defaults["n_bands"])),
                snr_k=float(target_input.get("snr_k", defaults["snr_k"])),
                frame_exp_time_s=float(target_input.get("frame_exp_time_s", defaults["frame_exp_time_s"])),
                pickles_support_dir=target_pickles_support_dir,
            )
            results.append(result)
            all_rows.append(result)
            class_total += result["t_exp_total_s"]
        class_summaries[classification] = {
            "targets": results,
            "totals": {
                "n_targets": len(results),
                "t_exp_total_s": class_total,
            },
        }
        total_exp += class_total

    summary = {
        "input_yaml": str(Path(yaml_path)),
        "defaults": defaults,
        "totals": {
            "n_targets": len(all_rows),
            "t_exp_total_s": total_exp,
        },
        "classifications": class_summaries,
        "software_versions": {
            "stp_etc_imaging": importlib.metadata.version("stp_etc_imaging"),
            "astropy": importlib.metadata.version("astropy"),
        },
    }

    _write_markdown_report(summary, output_dir / "exposure_time-report.md")
    _write_csv(all_rows, output_dir / "exposure_time-summary.csv")
    (output_dir / "exposure_time-summary.yaml").write_text(
        yaml.safe_dump(summary, sort_keys=False),
        encoding="utf-8",
    )

    return summary
