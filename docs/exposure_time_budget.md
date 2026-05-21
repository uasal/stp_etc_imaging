# Exposure-time budget handoff

`stp_etc_imaging` provides a target-list runner and an optional budgie adapter so downstream projects can keep their own `targets.yaml` while exposure-time logic stays in this package.

## Targets YAML schema

```yaml
observatory: UM
snr_k: 3
n_bands: 1
frame_exp_time_s: 60
classifications:
  exoplanet_host_stars:
    - name: HD 12345
      V_mag: 5.2
      desired_contrast: 1.0e-7
      sep_arcsec: 0.5
      spectral_type: G2V
      notes: ""
      snr_k: 5
      n_bands: 2
  debris_disks: []
  prime_targets: []
```

`n_bands` semantics: each band must independently reach contrast, so `t_exp_total = n_bands * t_exp_per_band`.

Contrast convention: `SNR_required = snr_k / desired_contrast` with default `snr_k = 3`.

## Direct invocation

```python
from stp_etc_imaging.target_list import run_target_list

run_target_list("targets.yaml", "output/", observatory_name="UM")
```

Artifacts written in `output/`:
- `exposure_time-report.md`
- `exposure_time-summary.csv`
- `exposure_time-summary.yaml`

## Via budgie runner

After the budgie companion PR lands (`douglase/copilot/refine-forest-rendering`):

```bash
python scripts/run_report.py stp_etc_imaging.budget_adapter:ExposureTimeBudget targets.yaml
```
