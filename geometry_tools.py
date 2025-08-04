"""
Geometry calculation and contour analysis tools for errant contour detection.

This module contains functions for:
- Calculating contour areas using the Shoelace formula
- Point-in-polygon testing using ray casting
- Contour boundary checking against external ROIs
- Splitting contours into inside/outside segments
"""

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


def contour_outside_external(contour_x, contour_y, external_contours_at_z):
    """
    Check if any part of a contour is outside all external contours at this z-level.
    Returns True if contour has points outside external boundary.
    """
    if not external_contours_at_z:
        return False  # No external contour to check against
    
    # Check several points along the contour (sample every few points for efficiency)
    sample_indices = range(0, len(contour_x), max(1, len(contour_x) // 10))
    
    for i in sample_indices:
        x, y = contour_x[i], contour_y[i]
        inside_any_external = False
        
        # Check if this point is inside any external contour
        for ext_contour in external_contours_at_z:
            if point_in_polygon(x, y, ext_contour['x_coords'], ext_contour['y_coords']):
                inside_any_external = True
                break
        
        # If this point is outside all external contours, the contour is problematic
        if not inside_any_external:
            return True
    
    return False


def split_contour_by_external(contour_x, contour_y, external_contours_at_z):
    """
    Split a contour into segments that are inside vs outside the external boundary.
    Returns a list of segments with 'inside' or 'outside' classification.
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
