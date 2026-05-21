import importlib.metadata
from pathlib import Path
__version__ = importlib.metadata.version(__package__ or "config_stp")

from .target_list import (
    load_targets,
    contrast_to_snr,
    run_target_list,
    exposure_time_for_target,
)

__all__ = [
    "__version__",
    "load_targets",
    "contrast_to_snr",
    "run_target_list",
    "exposure_time_for_target",
]
