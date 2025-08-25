"""
Utility functions for ErrantContourCheck

This module contains all calculation and analysis utility functions including:
- Contour area calculations
- Point-in-polygon testing
- Contour segmentation
- Statistical analysis functions
- Geometric utilities
"""

from collections import Counter


def calculate_contour_area(contour):
    """
    Calculate the area of a contour using the Shoelace formula.
    
    Args:
        contour: A list of points with .x, .y, .z attributes (or dictionary access)
    
    Returns:
        float: Area in square units (typically mm²)
    """
    if len(contour) < 3:
        return 0.0
    
    # Extract x and y coordinates
    x_coords = []
    y_coords = []
    
    for point in contour:
        # Handle both dictionary-style and attribute-style access
        if hasattr(point, 'x'):
            x_coords.append(point.x)
            y_coords.append(point.y)
        else:
            x_coords.append(point['x'])
            y_coords.append(point['y'])
    
    # Shoelace formula
    n = len(x_coords)
    area = 0.0
    
    for i in range(n):
        j = (i + 1) % n  # Next vertex (wraps around to 0 for last vertex)
        area += x_coords[i] * y_coords[j]
        area -= x_coords[j] * y_coords[i]
    
    return abs(area) / 2.0


def point_in_polygon(x, y, polygon_x, polygon_y):
    """
    Check if a point is inside a polygon using ray casting algorithm.
    Returns True if point is inside, False otherwise.
    
    Args:
        x, y: Coordinates of the point to test
        polygon_x, polygon_y: Lists of x and y coordinates of polygon vertices
    
    Returns:
        bool: True if point is inside polygon, False otherwise
    """
    n = len(polygon_x)
    inside = False
    
    p1x, p1y = polygon_x[0], polygon_y[0]
    for i in range(1, n + 1):
        p2x, p2y = polygon_x[i % n], polygon_y[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    
    return inside


def split_contour_by_external(contour_x, contour_y, external_contours_at_z):
    """
    Split a contour into segments that are inside vs outside the external boundary.
    Returns a list of segments with 'inside' or 'outside' classification.
    
    Args:
        contour_x, contour_y: Lists of x and y coordinates of the contour
        external_contours_at_z: List of external contours at this z-level
    
    Returns:
        list: List of segments with 'x_coords', 'y_coords', and 'type' keys
    """
    if not external_contours_at_z:
        # No external contour, treat entire contour as inside
        return [{'x_coords': contour_x, 'y_coords': contour_y, 'type': 'inside'}]
    
    segments = []
    current_segment_x = []
    current_segment_y = []
    current_type = None
    
    for i in range(len(contour_x)):
        x, y = contour_x[i], contour_y[i]
        
        # Check if this point is inside any external contour
        inside_any_external = False
        for ext_contour in external_contours_at_z:
            if point_in_polygon(x, y, ext_contour['x_coords'], ext_contour['y_coords']):
                inside_any_external = True
                break
        
        point_type = 'inside' if inside_any_external else 'outside'
        
        # If this is the first point or type hasn't changed, add to current segment
        if current_type is None or current_type == point_type:
            current_segment_x.append(x)
            current_segment_y.append(y)
            current_type = point_type
        else:
            # Type changed, save current segment and start new one
            if len(current_segment_x) > 1:  # Only save segments with multiple points
                segments.append({
                    'x_coords': current_segment_x[:],
                    'y_coords': current_segment_y[:],
                    'type': current_type
                })
            
            # Start new segment with this point
            current_segment_x = [x]
            current_segment_y = [y]
            current_type = point_type
    
    # Save final segment
    if len(current_segment_x) > 1:
        segments.append({
            'x_coords': current_segment_x,
            'y_coords': current_segment_y,
            'type': current_type
        })
    
    return segments if segments else [{'x_coords': contour_x, 'y_coords': contour_y, 'type': 'inside'}]


def calculate_global_parameters(roi_data):
    """
    Calculate global slice thickness and z-range parameters.
    
    Args:
        roi_data: Dictionary of ROI data with z_values for each ROI
    
    Returns:
        dict: Dictionary with 'all_z_values', 'slice_thickness', 'min_z', 'max_z'
    """
    all_z_values = []
    for roi_name, data in roi_data.items():
        all_z_values.extend(data['z_values'])
    
    if not all_z_values:
        return {
            'all_z_values': [],
            'slice_thickness': 0,
            'min_z': 0,
            'max_z': 0
        }
    
    # Sort z-values and calculate slice thickness
    ordered_z_values = sorted(all_z_values)
    z_differences = [round(a - b, 2) for a, b in zip(ordered_z_values[1:], ordered_z_values[:-1])]
    slice_thickness = min([diff for diff in z_differences if diff != 0]) if z_differences else 0
    
    return {
        'all_z_values': all_z_values,
        'slice_thickness': slice_thickness,
        'min_z': min(all_z_values),
        'max_z': max(all_z_values)
    }


def calculate_roi_statistics(roi_data, slice_thickness):
    """
    Calculate comprehensive statistics for each ROI.
    
    Args:
        roi_data: Dictionary of ROI data to update with statistics
        slice_thickness: Slice thickness in mm for gap analysis
    
    Returns:
        dict: Updated roi_data with statistics added
    """
    updated_data = {}
    
    for roi_name, data in roi_data.items():
        # Copy existing data
        updated_roi_data = data.copy()
        
        # Calculate statistics
        ordered_z = sorted(data['z_values'])
        unique_z = sorted(set(data['z_values']))
        z_diffs = [round(a - b, 2) for a, b in zip(ordered_z[1:], ordered_z[:-1])]
        gaps = sum(diff > slice_thickness for diff in z_diffs) if slice_thickness > 0 else 0
        
        total_area = sum(data['contour_areas'])
        avg_area = total_area / len(data['contour_areas']) if data['contour_areas'] else 0
        min_area = min(data['contour_areas']) if data['contour_areas'] else 0
        max_area = max(data['contour_areas']) if data['contour_areas'] else 0
        
        # Add statistics to data
        updated_roi_data.update({
            'gaps': gaps,
            'total_area': total_area,
            'avg_area': avg_area,
            'min_area': min_area,
            'max_area': max_area
        })
        
        updated_data[roi_name] = updated_roi_data
    
    return updated_data


def calculate_max_frequency(roi_data):
    """
    Calculate maximum frequency for area distribution plots.
    
    Args:
        roi_data: Dictionary of ROI data containing contour_areas
    
    Returns:
        int: Maximum frequency across all ROIs
    """
    max_frequency = 0
    for roi_name, data in roi_data.items():
        if 'contour_areas' in data and data['contour_areas']:
            areas = data['contour_areas']
            for area in areas:
                bin_count = sum(1 for a in areas if int(a) == int(area))
                max_frequency = max(max_frequency, bin_count)
    return max_frequency


def detect_global_issues(roi_data, expected_z_values, slice_thickness, external_contours=None):
    """
    Detect global issues across all ROIs for error reporting.
    
    Args:
        roi_data: Dictionary of ROI data
        expected_z_values: List of expected z-values for gap detection
        slice_thickness: Slice thickness for isolation analysis
        external_contours: Optional external contour data for boundary checking
    
    Returns:
        dict: Dictionary with boolean flags for different issue types
    """
    global_has_multiple_contours = False
    global_has_isolated_contours = False
    global_has_missing_contours = False
    global_has_external_violations = False
    
    isolation_threshold = 2 * slice_thickness if slice_thickness > 0 else 5.0
    
    for roi_name, data in roi_data.items():
        # Check for multiple contours and missing slices
        counts = []
        for expected_z in expected_z_values:
            counts.append(data['contours_per_slice'].get(expected_z, 0))
        
        if any(count > 1 for count in counts):
            global_has_multiple_contours = True
        if any(count == 0 for count in counts):
            global_has_missing_contours = True
        
        # Check for isolated contours
        for i, z in enumerate(expected_z_values):
            if counts[i] > 0:
                min_distance_to_neighbor = float('inf')
                for j, other_z in enumerate(expected_z_values):
                    if i != j and counts[j] > 0:
                        distance = abs(z - other_z)
                        min_distance_to_neighbor = min(min_distance_to_neighbor, distance)
                
                if min_distance_to_neighbor > isolation_threshold:
                    global_has_isolated_contours = True
                    break
        
        # Check for external violations (simplified check)
        if external_contours:
            global_has_external_violations = True  # Simplified - would need detailed analysis
    
    return {
        'multiple_contours': global_has_multiple_contours,
        'isolated_contours': global_has_isolated_contours,
        'missing_contours': global_has_missing_contours,
        'external_violations': global_has_external_violations
    }


def find_isolated_slices(expected_z_values, counts, slice_thickness):
    """
    Find slices with isolated contours based on distance to neighbors.
    
    Args:
        expected_z_values: List of z-values
        counts: List of contour counts for each z-value
        slice_thickness: Slice thickness for isolation threshold
    
    Returns:
        set: Set of z-values that are isolated
    """
    isolation_threshold = 2 * slice_thickness if slice_thickness > 0 else 5.0
    isolated_slices = set()
    
    for i, z in enumerate(expected_z_values):
        if counts[i] > 0:
            min_distance_to_neighbor = float('inf')
            for j, other_z in enumerate(expected_z_values):
                if i != j and counts[j] > 0:
                    distance = abs(z - other_z)
                    min_distance_to_neighbor = min(min_distance_to_neighbor, distance)
            
            if min_distance_to_neighbor > isolation_threshold:
                isolated_slices.add(z)
    
    return isolated_slices


def process_contour_coordinates(contour):
    """
    Process contour coordinates from RayStation format to lists.
    
    Args:
        contour: RayStation contour object with points
    
    Returns:
        tuple: (x_coords, y_coords) lists
    """
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
        x_coords.append(x_coords[0])  # Close contour
        y_coords.append(y_coords[0])
    
    return x_coords, y_coords


def determine_contour_color(area, segment_type, is_external=False):
    """
    Determine the appropriate color for a contour segment.
    
    Args:
        area: Contour area in mm²
        segment_type: 'inside', 'outside', or 'external'
        is_external: Whether this is an external contour
    
    Returns:
        tuple: (color, alpha) for visualization
    """
    if is_external:
        return 'gray', 0.5
    elif area < 2.0:
        return 'red', 0.8  # Small area always red
    elif segment_type == 'outside':
        return 'red', 0.8  # Outside external always red
    else:
        return 'blue', 0.6  # Normal contour
