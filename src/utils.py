from config import debug  # Use relative import for config

def ensure_roi_has_contours(plan_data, roi_name):
    if debug: print(f"ensure_roi_has_contours called for {roi_name}")
    """
    Ensure the ROI has a .Contours attribute. If not, create a copy and convert to contour representation.
    Returns the name of the ROI with contours.
    """
    try:
        roi_geom = plan_data.case.PatientModel.StructureSets[plan_data.exam.Name].RoiGeometries[roi_name]
        primary_shape = getattr(roi_geom, 'PrimaryShape', None)
        if primary_shape is not None and hasattr(primary_shape, 'Contours'):
            if debug: print(f"ROI {roi_name} already has contours")
            return roi_name  # Already has contours
    except Exception as e:
        if debug: print(f"Exception in ensure_roi_has_contours geometry check: {e}")

    # If we reach here, need to create a copy and convert to contour
    import random
    import string
    new_roi_name = f"{roi_name}_contour_{''.join(random.choices(string.ascii_uppercase + string.digits, k=4))}"
    if debug: print(f"Creating copy ROI for {roi_name}")
    try:
        plan_data.case.PatientModel.CreateRoi(Name=new_roi_name, Color="Yellow", Type="Organ", TissueName=None, RbeCellTypeName=None, RoiMaterial=None)
        plan_data.case.PatientModel.RegionsOfInterest[new_roi_name].CreateAlgebraGeometry(
            Examination=plan_data.exam,
            Algorithm="Contours",
            ExpressionA={ 'Operation': "Union", 'SourceRoiNames': [roi_name], 'MarginSettings': { 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 } },
            ExpressionB={ 'Operation': "Union", 'SourceRoiNames': [], 'MarginSettings': { 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 } },
            ResultOperation="None",
            ResultMarginSettings={ 'Type': "Expand", 'Superior': 0, 'Inferior': 0, 'Anterior': 0, 'Posterior': 0, 'Right': 0, 'Left': 0 })
        plan_data.case.PatientModel.StructureSets[plan_data.exam.Name].RoiGeometries[new_roi_name].SetRepresentation(Representation="Contours")
        if debug: print(f"Successfully created contour ROI: {new_roi_name}")
        # Ensure temp_rois exists
        if not hasattr(plan_data, 'temp_rois'):
            plan_data.temp_rois = []
        plan_data.temp_rois.append(new_roi_name)  # Track for later deletion
        if debug: print(f"Added {new_roi_name} to temp_rois: {plan_data.temp_rois}")
        return new_roi_name
    except Exception as e:
        if debug: print(f"Failed to create contour ROI for '{roi_name}': {e}")
        return None

def point_in_polygon(x, y, polygon_x, polygon_y):
    """
    Check if a point is inside a polygon using ray casting algorithm.
    Returns True if point is inside, False otherwise.
    """
    n = len(polygon_x)
    inside = False
    if n < 3:
        return False
    p1x, p1y = polygon_x[0], polygon_y[0]
    for i in range(1, n + 1):
        p2x, p2y = polygon_x[i % n], polygon_y[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    else:
                        xinters = p1x
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
    sample_step = max(1, len(contour_x) // 10)
    sample_indices = range(0, len(contour_x), sample_step)
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

def calculate_contour_area(contour):
    """
    Calculate the area of a contour using the Shoelace formula.
    Args:
        contour: A list of points with .x, .y, .z attributes (or dictionary access)
    Returns:
        float: Area in square units (typically mm²)
    """
    # Impossible to have area > 0 with less than 3 points
    if len(contour) < 3:
        return 0.0
    # Extract x and y coordinates
    x_coords = []
    y_coords = []
    for point in contour:
        # Handle both dictionary-style and attribute-style access
        if hasattr(point, 'x') and hasattr(point, 'y'):
            x_coords.append(point.x)
            y_coords.append(point.y)
        elif isinstance(point, dict) and 'x' in point and 'y' in point:
            x_coords.append(point['x'])
            y_coords.append(point['y'])
        else:
            raise ValueError("Contour point must have 'x' and 'y' attributes or keys")
    # Shoelace formula
    n = len(x_coords)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n  # Next vertex (wraps around to 0 for last vertex)
        area += x_coords[i] * y_coords[j]
        area -= x_coords[j] * y_coords[i]
    return abs(area) / 2.0

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
    if not segments:
        return [{'x_coords': contour_x, 'y_coords': contour_y, 'type': 'inside'}]
    return segments