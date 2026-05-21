from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import astropy.units as u
import pytest
import yaml

from stp_etc_imaging.target_list import (
    contrast_to_snr,
    exposure_time_for_target,
    load_targets,
    run_target_list,
)


def _write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def _mock_observatory(mock_obs):
    mock_obs.make_STP.return_value = None
    mock_obs.set_source.return_value = None
    mock_obs.set_background.return_value = None
    mock_obs.make_observation.return_value = (1 * u.electron / u.s, 1 * u.electron / u.s)
    mock_obs.calc_int_time.return_value = 120 * u.s
    mock_obs.data_observatory = {"astrophysics": {"zodi": {"profile": "zodi_profile.fits"}}}
    mock_obs.support_data_observatory = "/support/"


def test_contrast_to_snr():
    assert contrast_to_snr(1e-7) == pytest.approx(3e7)
    assert contrast_to_snr(1e-7, k=5) == pytest.approx(5e7)


def test_load_targets_validation_and_overrides(tmp_path):
    good = {
        "observatory": "UM",
        "snr_k": 3,
        "n_bands": 1,
        "frame_exp_time_s": 60,
        "classifications": {
            "exoplanet_host_stars": [
                {
                    "name": "HD 12345",
                    "V_mag": 5.2,
                    "desired_contrast": 1e-7,
                    "snr_k": 5,
                }
            ]
        },
    }
    good_path = _write_yaml(tmp_path / "good.yaml", good)
    loaded = load_targets(good_path)
    target = loaded["classifications"]["exoplanet_host_stars"][0]
    assert target["snr_k"] == 5
    assert target["n_bands"] == 1
    assert target["observatory"] == "UM"

    bad = {
        "classifications": {
            "exoplanet_host_stars": [
                {
                    "name": "HD missing contrast",
                    "V_mag": 5.2,
                }
            ]
        }
    }
    bad_path = _write_yaml(tmp_path / "bad.yaml", bad)
    with pytest.raises(ValueError, match="missing required fields"):
        load_targets(bad_path)


def test_n_bands_scaling(tmp_path):
    target = {
        "name": "HD 12345",
        "classification": "exoplanet_host_stars",
        "V_mag": 5.2,
        "desired_contrast": 1e-7,
        "n_bands": 3,
        "snr_k": 3,
    }

    with patch("stp_etc_imaging.target_list.Observatory") as obs_cls:
        mock_obs = obs_cls.return_value
        _mock_observatory(mock_obs)

        result = exposure_time_for_target(target, observatory_name="UM")

    called_kwargs = mock_obs.calc_int_time.call_args.kwargs
    assert called_kwargs["snr"] == pytest.approx(3e7)
    assert result["t_exp_per_band_s"] == pytest.approx(120.0)
    assert result["t_exp_total_s"] == pytest.approx(360.0)


def test_run_target_list_outputs_and_totals(tmp_path):
    input_yaml = {
        "observatory": "UM",
        "snr_k": 3,
        "n_bands": 1,
        "frame_exp_time_s": 60,
        "classifications": {
            "exoplanet_host_stars": [
                {"name": "A1", "V_mag": 5.0, "desired_contrast": 1e-7},
                {"name": "A2", "V_mag": 6.0, "desired_contrast": 2e-7, "n_bands": 2},
            ],
            "debris_disks": [
                {"name": "B1", "V_mag": 7.0, "desired_contrast": 1e-6},
                {"name": "B2", "V_mag": 8.0, "desired_contrast": 2e-6},
            ],
            "prime_targets": [
                {"name": "C1", "V_mag": 9.0, "desired_contrast": 3e-7},
                {"name": "C2", "V_mag": 10.0, "desired_contrast": 4e-7},
            ],
        },
    }
    yaml_path = _write_yaml(tmp_path / "targets.yaml", input_yaml)

    with patch("stp_etc_imaging.target_list.Observatory") as obs_cls:
        mock_obs = obs_cls.return_value
        _mock_observatory(mock_obs)

        out_dir = tmp_path / "out"
        summary = run_target_list(yaml_path, out_dir)

    report_md = out_dir / "exposure_time-report.md"
    summary_csv = out_dir / "exposure_time-summary.csv"
    summary_yaml = out_dir / "exposure_time-summary.yaml"

    assert report_md.exists()
    assert summary_csv.exists()
    assert summary_yaml.exists()

    rows = summary_csv.read_text(encoding="utf-8").strip().splitlines()
    assert len(rows) == 7  # header + 6 targets

    # 5 single-band targets + 1 two-band target, each 120s per-band
    assert summary["totals"]["n_targets"] == 6
    assert summary["totals"]["t_exp_total_s"] == pytest.approx(120.0 * 7)


def test_budget_adapter_run_report_delegates(tmp_path):
    budgie = pytest.importorskip("budgie")

    from stp_etc_imaging.budget_adapter import ExposureTimeBudget

    assert issubclass(ExposureTimeBudget, budgie.Budget)

    budget = ExposureTimeBudget.__new__(ExposureTimeBudget)
    budget.budget_dir = str(tmp_path)
    budget.name = "targets.yaml"

    with patch("stp_etc_imaging.target_list.run_target_list", return_value={"ok": True}) as runner:
        result = budget.run_report(tmp_path / "reports")

    assert result == {"ok": True}
    called_yaml_path, called_output_dir = runner.call_args.args
    assert called_yaml_path == tmp_path / "targets.yaml"
    assert called_output_dir == tmp_path / "reports"
