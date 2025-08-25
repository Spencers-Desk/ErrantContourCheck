# ErrantContourCheck - Usage Guide

## Overview
ErrantContourCheck is a modular tool for analyzing ROI contours in RayStation to identify potential errors such as:
- Missing contours on slices (gaps)
- Multiple contours on the same slice
- Isolated contours 
- Very small contours (< 2mm²)
- Contours extending outside the External boundary

## Quick Start

### Method 1: Complete Workflow (Recommended)
```python
# In RayStation scripting environment:
exec(open(r'path\to\ErrantContourCheck\src\main.py').read())
```

### Method 2: Step-by-Step Usage
```python
# Import modules
from init import plan_data
from roi_selection_gui import get_roi_selection
from analysis_gui import create_analysis_gui

# Initialize and build plan data
plan_data.build()

# Select ROIs
selected_rois, external_roi = get_roi_selection(plan_data)

# Process data
plan_data.process_roi_data()

# Show analysis GUI
create_analysis_gui(plan_data)
```

## Module Structure

### utils.py
- **Calculation Functions**: All mathematical and statistical calculations
- **Geometric Utilities**: Point-in-polygon, contour segmentation, area calculation
- **Analysis Functions**: Error detection, statistics computation, parameter calculation
- **Color Logic**: Visualization color determination based on contour properties

### init.py
- **PlanData class**: Handles RayStation connection, ROI discovery, and data processing
- **plan_data instance**: Pre-configured global instance ready to use
- Key methods:
  - `build()`: Connect to RayStation and discover ROIs
  - `set_rois(selected_roi_names, external_roi_name)`: Set which ROIs to analyze
  - `process_roi_data()`: Process contour data and calculate statistics

### roi_selection_gui.py
- **get_roi_selection()**: Interactive GUI for selecting ROIs to analyze
- Validates ROI selections and prevents conflicts
- Returns selected ROIs and optional external ROI

### analysis_gui.py  
- **AnalysisGUI class**: Comprehensive visualization interface
- **create_analysis_gui()**: Convenience function to create and show GUI
- Features:
  - Contours per slice plots with error highlighting
  - Area distribution histograms
  - Interactive slice-by-slice viewers
  - Mouse and keyboard navigation
  - Error detection banners

### main.py
- Complete workflow orchestration
- Error handling and user feedback
- Single entry point for the entire analysis

## GUI Navigation

### Analysis Interface
- **Red bars/contours**: Issues detected (multiple contours, small areas, outside external)
- **Orange bars**: Isolated contours
- **Gray bars**: Missing contours (gaps)
- **Blue contours**: Normal contours

### Slice Navigation
- **Left Click**: Next slice (in clicked ROI viewer)
- **Right Click**: Previous slice (in clicked ROI viewer)  
- **Arrow Keys**: Navigate all ROI viewers simultaneously

## Requirements
- RayStation scripting environment
- Python 3.6+ (typical for RayStation)
- Libraries: matplotlib, tkinter (usually available in RayStation)

## Output Statistics
The tool provides:
- Number of contours per ROI
- Gap analysis (missing slices)
- Area statistics (total, average, min, max)
- Visual highlighting of potential issues

## Error Detection
- **Multiple contours per slice**: More than one contour on same z-level
- **Isolated contours**: Contours with no nearby neighbors (>2x slice thickness)
- **Small contours**: Area < 2.0 mm²
- **External violations**: Contours extending outside External ROI boundary
- **Missing contours**: Expected slices with no contours (gaps)

## File Organization
```
ErrantContourCheck/
├── src/
│   ├── init.py              # Data initialization and processing
│   ├── roi_selection_gui.py # ROI selection interface
│   ├── analysis_gui.py      # Main analysis visualization
│   ├── main.py              # Complete workflow entry point
├── tests/
│   ├── test_workflow.py     # Workflow demonstration test
│   └── test_utils.py        # Unit tests for utility functions
└── USAGE.md                 # This file
```

## Development Notes
- Built for Python 3.6 compatibility (RayStation standard)
- Modular design for maintainability and testing
- Uses matplotlib for visualization and tkinter for GUI
- Handles various RayStation geometry data formats
- Designed for public release with clean interfaces
