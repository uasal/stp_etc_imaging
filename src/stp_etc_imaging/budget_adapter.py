"""Thin budgie adapter for exposure-time target lists."""

from __future__ import annotations

try:
    from budgie import Budget

    _HAS_BUDGIE = True
except ImportError:  # pragma: no cover - exercised in environments without budgie
    Budget = object
    _HAS_BUDGIE = False


class ExposureTimeBudget(Budget):
    """budgie.Budget adapter that delegates to stp_etc_imaging.target_list."""

    def __init__(self, name: str = "exposure_time.yaml"):
        if not _HAS_BUDGIE:
            raise ImportError(
                "budgie is required: pip install 'budgie @ git+https://github.com/uasal/budgie@develop'"
            )
        super().__init__(name)

    def run_report(self, output_dir):
        from pathlib import Path

        from .target_list import run_target_list

        yaml_path = Path(self.budget_dir) / self.name
        return run_target_list(yaml_path, output_dir)
