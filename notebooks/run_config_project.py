# set environment
from pathlib import Path
import sys

# Find project root based on this file's location
_THIS_FILE = Path(__file__).resolve()
FD_NBK     = _THIS_FILE.parent  # notebook folder
FD_PRJ     = FD_NBK.parent      # project root dir
FD_EXE     = FD_PRJ / "scripts" # script folder

# Put scripts on sys.path if not already there
spath = str(FD_EXE)
if spath not in sys.path:
    sys.path.insert(0, spath)

# Import the config
from config_project import *