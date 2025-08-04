"""
RayStation interface module for contour data extraction and processing.

This module handles:
- Getting current RayStation objects (patient, case, exam)
- Processing ROI contour data
- Extracting external contour boundaries
- Calculating slice thickness and z-ranges
"""

try:
    from connect import get_current
except ImportError:
    # Fallback for testing outside RayStation
    def get_current(obj_type):
        raise RuntimeError(f"Cannot get {obj_type} - not running in RayStation")

from collections import Counter

try:
    from .geometry_tools import calculate_contour_area
except ImportError:
    from geometry_tools import calculate_contour_area


def get_raystation_objects():
    """Get the core RayStation objects needed for analysis."""
    patient = get_current("Patient")
    case = get_current("Case")
    exam = get_current("Examination")
    return patient, case, exam


def get_available_rois(case, exam):
    """
    Get lists of available target and external ROIs with contours.
    
    Returns:
        tuple: (target_rois, external_rois) - lists of ROI names
    """
    # Use the RegionsOfInterest approach for better performance
    all_rois = case.PatientModel.RegionsOfInterest
    print(f"Found {len(all_rois)} total ROIs")
    
    # Filter for Target type ROIs and those with contours
    target_rois = []
    external_rois = []
    
    print("Filtering ROIs and checking for contours...")
    
    for i, roi in enumerate(all_rois):
        if i % 10 == 0:  # Progress indicator
            print(f"  Processing ROI {i+1}/{len(all_rois)}: {roi.Name}")
        
        roi_name = roi.Name
        roi_type = roi.Type
        
        # Check if ROI has contours more efficiently
        try:
            roi_geometry = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_name]
            # Just check if PrimaryShape exists, don't count contours yet
            has_contours = (hasattr(roi_geometry, 'PrimaryShape') and 
                           roi_geometry.PrimaryShape is not None)
        except:
            has_contours = False
        
        if has_contours:
            if roi_type == "External":
                external_rois.append(roi_name)
            elif roi_type in ["Ptv", "Ctv", "Gtv", "Target"]:  # Common target types
                target_rois.append(roi_name)
    
    print(f"Found {len(target_rois)} target ROIs, {len(external_rois)} external ROIs")
    
    # Sort the lists
    target_rois.sort()
    external_rois.sort()
    
    return target_rois, external_rois


def validate_roi_contours(case, exam, roi_names):
    """
    Validate that selected ROIs actually have contours.
    
    Returns:
        list: Valid ROI names that have contours
    """
    valid_selections = []
    for roi_name in roi_names:
        try:
            roi_geometry = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_name]
            if (hasattr(roi_geometry, 'PrimaryShape') and 
                roi_geometry.PrimaryShape is not None and
                hasattr(roi_geometry.PrimaryShape, 'Contours') and
                len(roi_geometry.PrimaryShape.Contours) > 0):
                valid_selections.append(roi_name)
                print(f"  {roi_name}: {len(roi_geometry.PrimaryShape.Contours)} contours")
            else:
                print(f"  {roi_name}: No contours found!")
        except Exception as e:
            print(f"  {roi_name}: Error checking contours - {e}")
    
    return valid_selections


def process_roi_contours(case, exam, roi_names):
    """
    Process contour data for selected ROIs.
    
    Returns:
        dict: ROI data dictionary with contour information
    """
    roi_data = {}
    
    for roi_name in roi_names:
        print(f"\nProcessing ROI: {roi_name}")
        
        try:
            geom = case.PatientModel.StructureSets[exam.Name].RoiGeometries[roi_name].PrimaryShape
            contours = geom.Contours
            
            number_of_contours = len(contours)
            print(f"  {number_of_contours} contours found")
            
            z_values = [None] * number_of_contours
            contour_areas = [None] * number_of_contours
            
            for index, contour in enumerate(contours):
                z_values[index] = round(contours[index][0]['z'], 2)
                contour_areas[index] = calculate_contour_area(contour)
            
            # Store data for this ROI
            roi_data[roi_name] = {
                'contours': contours,
                'number_of_contours': number_of_contours,
                'z_values': z_values,
                'contour_areas': contour_areas,
                'contours_per_slice': Counter(z_values)
            }
            
        except Exception as e:
            print(f"  Error processing {roi_name}: {e}")
            continue
    
    return roi_data


def process_external_contours(case, exam, external_roi_name, target_z_values):
    """
    Process external ROI contours for boundary checking.
    
    Args:
        case: RayStation case object
        exam: RayStation exam object
        external_roi_name: Name of external ROI
        target_z_values: Set of z-values from target ROIs
    
    Returns:
        dict: External contours by z-slice
    """
    external_contours = {}
    
    try:
        print(f"\nProcessing External ROI: {external_roi_name}")
        external_geom = case.PatientModel.StructureSets[exam.Name].RoiGeometries[external_roi_name].PrimaryShape
        
        external_contour_list = external_geom.Contours
        print(f"  {len(external_contour_list)} external contours found")
        
        # Process only external contours at z-levels we care about
        for contour in external_contour_list:
            z_val = round(contour[0]['z'], 2)
            if z_val in target_z_values:
                # Extract x, y coordinates for this external contour
                x_coords = []
                y_coords = []
                for point in contour:
                    if hasattr(point, 'x'):
                        x_coords.append(point.x)
                        y_coords.append(-point.y)  # Flip y-coordinate to match RayStation orientation
                    else:
                        x_coords.append(point['x'])
                        y_coords.append(-point['y'])  # Flip y-coordinate to match RayStation orientation
                
                if x_coords and y_coords:
                    # Store external contour coordinates for this z-slice
                    if z_val not in external_contours:
                        external_contours[z_val] = []
                    external_contours[z_val].append({
                        'x_coords': x_coords,
                        'y_coords': y_coords
                    })
        
        print(f"  External contours processed for {len(external_contours)} relevant z-slices")
        
    except Exception as e:
        print(f"  Warning: Could not process External ROI '{external_roi_name}': {e}")
        print("  Boundary checking will be disabled.")
    
    return external_contours


def calculate_slice_analysis(roi_data):
    """
    Calculate global slice thickness and analysis for all ROIs.
    
    Returns:
        tuple: (slice_thickness, min_z, max_z, all_z_values)
    """
    # Calculate global z-range and slice thickness from all ROIs
    all_z_values = []
    for roi_name, data in roi_data.items():
        all_z_values.extend(data['z_values'])
    
    if not all_z_values:
        raise ValueError("No contours found in any ROI!")
    
    # Sort z-values and calculate slice thickness
    ordered_z_values = sorted(all_z_values)
    z_differences = [round(a - b, 2) for a, b in zip(ordered_z_values[1:], ordered_z_values[:-1])]
    slice_thickness = min([diff for diff in z_differences if diff != 0])
    
    min_z = min(all_z_values)
    max_z = max(all_z_values)
    
    print(f"\nGlobal Analysis:")
    print(f"Slice thickness: {slice_thickness} mm")
    print(f"Z-range: {min_z:.1f} to {max_z:.1f} mm")
    
    return slice_thickness, min_z, max_z, all_z_values


def analyze_roi_statistics(roi_data, slice_thickness):
    """
    Analyze statistics for each ROI and add to roi_data.
    """
    for roi_name, data in roi_data.items():
        print(f"\n{roi_name} Analysis:")
        
        ordered_z = sorted(data['z_values'])
        unique_z = sorted(set(data['z_values']))
        z_diffs = [round(a - b, 2) for a, b in zip(ordered_z[1:], ordered_z[:-1])]
        gaps = sum(diff > slice_thickness for diff in z_diffs)
        
        total_area = sum(data['contour_areas'])
        avg_area = total_area / len(data['contour_areas']) if data['contour_areas'] else 0
        min_area = min(data['contour_areas']) if data['contour_areas'] else 0
        max_area = max(data['contour_areas']) if data['contour_areas'] else 0
        
        print(f"  Total contours: {data['number_of_contours']}")
        print(f"  Unique slices: {len(unique_z)}")
        print(f"  Gaps found: {gaps}")
        print(f"  Total area: {total_area:.2f} mm²")
        print(f"  Average area: {avg_area:.2f} mm²")
        print(f"  Min area: {min_area:.2f} mm²")
        print(f"  Max area: {max_area:.2f} mm²")
        
        # Store additional analysis
        data.update({
            'gaps': gaps,
            'total_area': total_area,
            'avg_area': avg_area,
            'min_area': min_area,
            'max_area': max_area
        })
