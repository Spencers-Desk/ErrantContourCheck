"""
ErrantCheck - Multi-ROI Contour Analysis Tool for RayStation

A comprehensive quality assurance tool for detecting and visualizing errant 
contours in RayStation treatment planning structures.

Modules:
    geometry_tools: Geometric calculations and contour analysis functions
    raystation_interface: RayStation data access and processing
    roi_selector_gui: Interactive ROI selection dialog
    main_analysis_gui: Comprehensive analysis visualization interface
    main: Main entry point and workflow orchestration

Usage:
    from errantCheck import main
    main.main()

Or run directly:
    python -m errantCheck.main
"""

__version__ = "1.0.0"
__author__ = "RayStation Tools Development"
__description__ = "Multi-ROI Contour Analysis Tool for RayStation"

# Import main components for easy access
from . import geometry_tools
from . import raystation_interface
from . import roi_selector_gui
from . import main_analysis_gui
from . import main

# Convenience imports
from .main import main as run_analysis

__all__ = [
    'geometry_tools',
    'raystation_interface', 
    'roi_selector_gui',
    'main_analysis_gui',
    'main',
    'run_analysis'
]
