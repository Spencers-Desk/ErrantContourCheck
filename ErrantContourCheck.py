# This script can be imported into RayStation's "Script Management" to act as an entry point for the main script.

import os, sys

# Ensure this script's directory and its 'src' subfolder are on sys.path
_BASE_DIR = os.path.dirname(os.path.abspath(__file__))
_SRC_DIR = os.path.join(_BASE_DIR, "src")
for _p in (_BASE_DIR, _SRC_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Robust import of main module
try:
    from archive import main as ecc_main  # package style (if src is a package)
except ImportError:
    try:
        import importlib
        ecc_main = importlib.import_module("main")  # direct module in added src path
    except Exception as e:
        raise ImportError("Could not import 'main' from src. Check that src/main.py exists.") from e