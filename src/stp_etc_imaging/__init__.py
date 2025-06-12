import importlib.metadata
from pathlib import Path
__version__ = importlib.metadata.version(__package__ or "config_stp")

__all__ = [ "__version__"]
