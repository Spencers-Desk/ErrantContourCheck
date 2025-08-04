# ErrantCheck - Multi-ROI Contour Analysis Tool

A comprehensive quality assurance tool for detecting and visualizing errant contours in RayStation treatment planning structures.

## Overview

ErrantCheck analyzes multiple Regions of Interest (ROIs) simultaneously to identify potential contouring errors including:

- **Gap Detection**: Missing contours between slices
- **Multiple Contour Detection**: Multiple contours on the same slice
- **Small Contour Detection**: Contours with area < 2 mm²
- **Isolation Detection**: Contours separated from neighbors by >2× slice thickness
- **Boundary Violations**: Contour portions extending outside External ROI

## File Structure

```
errantCheck/
├── __init__.py                 # Package initialization
├── main.py                     # Main entry point and workflow orchestration
├── run_analysis.py             # Simple launcher script
├── geometry_tools.py           # Geometric calculations and contour analysis
├── raystation_interface.py     # RayStation data access and processing
├── roi_selector_gui.py         # Interactive ROI selection dialog
├── main_analysis_gui.py        # Comprehensive analysis visualization
├── errant_contours.py          # Original monolithic script (backup)
└── errant_patients.txt         # Patient list file
```

## Module Descriptions

### geometry_tools.py
Contains core geometric calculation functions:
- `calculate_contour_area()`: Shoelace formula for contour area calculation
- `point_in_polygon()`: Ray casting algorithm for point-in-polygon testing
- `contour_outside_external()`: Check if contour extends outside external boundary
- `split_contour_by_external()`: Split contours into inside/outside segments

### raystation_interface.py
Handles all RayStation data interaction:
- `get_raystation_objects()`: Get patient, case, and exam objects
- `get_available_rois()`: Filter and list available target/external ROIs
- `process_roi_contours()`: Extract and process contour data
- `process_external_contours()`: Process external ROI for boundary checking
- `calculate_slice_analysis()`: Calculate global slice thickness and ranges
- `analyze_roi_statistics()`: Compute statistical analysis for each ROI

### roi_selector_gui.py
Provides the initial ROI selection interface:
- Fast, optimized ROI discovery and filtering
- Support for up to 4 target ROIs selection
- External ROI selection for boundary checking
- Input validation and error handling
- Performance timing and progress indicators

### main_analysis_gui.py
Creates the comprehensive analysis visualization:
- Multi-panel layout with contour plots, area distributions, and slice viewers
- Interactive slice navigation with mouse and keyboard controls
- Color-coded error highlighting (red for errors, orange for warnings)
- Pre-rendered slice data for instantaneous navigation
- Scrollable interface with error banners and statistics
- Real-time error detection and reporting

### main.py
Orchestrates the complete analysis workflow:
- Step-by-step process coordination
- Error handling and user feedback
- Integration of all modules
- Comprehensive logging and status reporting

## Usage

### Method 1: Direct Execution
```python
# In RayStation script console
exec(open(r'\\nasr200\rayStationScripts\dev\SRL\errantCheck\run_analysis.py').read())
```

### Method 2: Module Import
```python
# Add to Python path if needed
import sys
sys.path.append(r'\\nasr200\rayStationScripts\dev\SRL')

# Import and run
from errantCheck.main import main
main()
```

### Method 3: Package Import
```python
import errantCheck
errantCheck.run_analysis()
```

## Workflow

1. **RayStation Connection**: Automatically connects to current patient/case/exam
2. **ROI Selection**: User selects 1-4 target ROIs and an external ROI via GUI
3. **Data Processing**: Extracts and processes contour data for all selected ROIs
4. **Analysis**: Calculates statistics, detects errors, and prepares visualizations
5. **Visualization**: Launches comprehensive GUI with:
   - Contours per slice plots (color-coded for errors)
   - Area distribution histograms
   - Interactive slice viewers with navigation
   - Error detection banners
   - Detailed statistics panel

## Color Coding

- **Blue**: Normal contours
- **Red**: Problematic contours (small area <2mm² or outside external boundary)
- **Orange**: Isolated contours (separated by >2× slice thickness)
- **Gray**: Missing slices or external contour references
- **Light Gray**: Missing contour slices

## Requirements

- RayStation with Python scripting enabled
- Patient case loaded with contoured structures
- Target ROIs (PTV, CTV, GTV, Target types)
- External ROI (recommended for boundary checking)

## Performance Notes

- ROI discovery optimized for large case databases
- Slice data pre-rendered for instantaneous navigation
- Memory-efficient artist hide/show approach for slice viewing
- Progress indicators for long operations

## Error Detection Features

- **Gap Analysis**: Automatically detects missing slices based on uniform thickness
- **Multiplicity Check**: Flags slices with >1 contour
- **Size Validation**: Identifies unusually small contours
- **Isolation Detection**: Finds contours separated from main volume
- **Boundary Validation**: Checks for contours extending outside patient boundary
- **Statistical Analysis**: Comprehensive area and distribution analysis

## Output

The tool provides both visual and quantitative analysis:
- Interactive plots for visual inspection
- Detailed statistics per ROI
- Error summaries and counts
- Slice-by-slice navigation for detailed review
- Export-ready analysis data

## Troubleshooting

Common issues and solutions:
1. **Import Errors**: Ensure errantCheck directory is in Python path
2. **No ROIs Found**: Check that case has contoured Target-type structures
3. **Performance Issues**: Use External ROI selection to limit processing scope
4. **Display Issues**: Ensure matplotlib and tkinter are properly configured

## Development Notes

The modular structure facilitates:
- Easy maintenance and updates
- Individual component testing
- Feature additions and modifications
- Code reuse across other tools
- Clear separation of concerns

Each module has comprehensive docstrings and error handling for robust operation in the clinical environment.
