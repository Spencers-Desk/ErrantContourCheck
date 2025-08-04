#!/usr/bin/env python
"""
Multi-ROI Contour Analysis Tool for RayStation

This is the main entry point for the errant contour detection and analysis tool.
It orchestrates the workflow from ROI selection through data processing to 
visualization and analysis.

Usage:
    Run this script in RayStation with a patient case loaded.
    
Workflow:
    1. Get current RayStation objects (patient, case, exam)
    2. Present ROI selection GUI
    3. Process selected ROI contour data
    4. Process external ROI for boundary checking
    5. Analyze contour statistics and detect issues
    6. Launch comprehensive analysis GUI

Features:
    - Multi-ROI analysis (up to 4 ROIs simultaneously)
    - Gap detection between slices
    - Multiple contour detection per slice
    - Small contour detection (<2mm²)
    - Isolated contour detection
    - Boundary violation detection (outside External ROI)
    - Interactive slice-by-slice visualization
    - Comprehensive statistics and error reporting
"""

# Module imports
try:
    from .raystation_interface import (
        get_raystation_objects, 
        process_roi_contours, 
        process_external_contours,
        calculate_slice_analysis,
        analyze_roi_statistics
    )
    from .roi_selector_gui import get_roi_selection
    from .main_analysis_gui import create_analysis_gui
except ImportError:
    # Fallback for direct execution
    from raystation_interface import (
        get_raystation_objects, 
        process_roi_contours, 
        process_external_contours,
        calculate_slice_analysis,
        analyze_roi_statistics
    )
    from roi_selector_gui import get_roi_selection
    from main_analysis_gui import create_analysis_gui


def main():
    """
    Main function that orchestrates the complete analysis workflow.
    """
    print("="*60)
    print("Multi-ROI Contour Analysis Tool")
    print("="*60)
    
    # Step 1: Get RayStation objects
    print("\\n1. Connecting to RayStation...")
    try:
        patient, case, exam = get_raystation_objects()
        print(f"   Patient: {patient.Name}")
        print(f"   Case: {case.CaseName}")
        print(f"   Exam: {exam.Name}")
    except Exception as e:
        print(f"   Error: Could not connect to RayStation - {e}")
        return
    
    # Step 2: ROI Selection
    print("\\n2. Opening ROI selection dialog...")
    roi_names, external_roi_name = get_roi_selection(case, exam)
    
    if roi_names is None:
        print("   No ROIs selected. Exiting.")
        return
    
    print(f"   Selected ROIs: {roi_names}")
    print(f"   External ROI: {external_roi_name}")
    
    # Step 3: Process ROI contour data
    print("\\n3. Processing ROI contour data...")
    roi_data = process_roi_contours(case, exam, roi_names)
    
    if not roi_data:
        print("   No valid ROI data found. Exiting.")
        return
    
    # Step 4: Process External ROI
    print("\\n4. Processing External ROI for boundary checking...")
    # Get all z-values from target ROIs
    target_z_values = set()
    for roi_name, data in roi_data.items():
        target_z_values.update(data['z_values'])
    
    external_contours = process_external_contours(case, exam, external_roi_name, target_z_values)
    
    # Step 5: Calculate global analysis
    print("\\n5. Performing global analysis...")
    try:
        slice_thickness, min_z, max_z, all_z_values = calculate_slice_analysis(roi_data)
    except ValueError as e:
        print(f"   Error: {e}")
        return
    
    # Step 6: Analyze individual ROI statistics
    print("\\n6. Analyzing ROI statistics...")
    analyze_roi_statistics(roi_data, slice_thickness)
    
    # Step 7: Launch analysis GUI
    print("\\n7. Launching analysis GUI...")
    try:
        create_analysis_gui(roi_names, roi_data, external_contours, slice_thickness, all_z_values)
        print("   Analysis complete.")
    except Exception as e:
        print(f"   Error creating GUI: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
