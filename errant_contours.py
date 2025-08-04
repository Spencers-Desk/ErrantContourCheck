from connect import get_current
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
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

# Outline
## Check for gaps
## Check for >1 contours on a single slice

# Assumptions
## Uniform slice thickness
## true tumor shape is contiguous in all slice - cannot account for bilateral nodal chains for instance

# 1. Grab the core RayStation objects
patient = get_current("Patient")
case    = get_current("Case")
exam    = get_current("Examination")

# 2. Get available ROIs and create selection GUI
def get_roi_selection():
    """Create a GUI for selecting ROIs to analyze."""
    import tkinter as tk
    from tkinter import ttk, messagebox
    import time
    
    print("Starting ROI selection GUI...")
    start_time = time.time()
    
    # Get all ROIs from the case more efficiently (similar to api.py approach)
    try:
        print("Getting ROI list...")
        roi_start = time.time()
        
        # Use the RegionsOfInterest approach like in api.py for better performance
        all_rois = case.PatientModel.RegionsOfInterest
        print(f"Found {len(all_rois)} total ROIs in {time.time() - roi_start:.2f}s")
        
        # Filter for Target type ROIs and those with contours
        target_rois = []
        external_rois = []
        
        print("Filtering ROIs and checking for contours...")
        filter_start = time.time()
        
        for i, roi in enumerate(all_rois):
            if i % 10 == 0:  # Progress indicator
                print(f"  Processing ROI {i+1}/{len(all_rois)}: {roi.Name}")
            
            roi_name = roi.Name
            roi_type = roi.Type
            
            # Check if ROI has contours more efficiently
            # Skip the slow contour checking for now - just check if geometry exists
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
        
        print(f"ROI filtering completed in {time.time() - filter_start:.2f}s")
        print(f"Found {len(target_rois)} target ROIs, {len(external_rois)} external ROIs")
        
        # Sort the lists
        target_rois.sort()
        external_rois.sort()
        
    except Exception as e:
        print(f"Error during ROI processing: {e}")
        messagebox.showerror("Error", f"Could not retrieve ROI list: {e}")
        return None, None
    
    if not target_rois:
        messagebox.showwarning("No Target ROIs", "No Target-type ROIs with geometry found!")
        return None, None
    
    print(f"Starting GUI creation after {time.time() - start_time:.2f}s total...")
    gui_start = time.time()
    
    # Create selection window
    root = tk.Tk()
    root.title("Select ROIs for Contour Analysis")
    root.geometry("500x400")
    root.resizable(False, False)
    
    # Variables to store selections
    selected_rois = []
    selected_external = tk.StringVar()
    
    # Main frame
    main_frame = tk.Frame(root, padx=20, pady=20)
    main_frame.pack(fill=tk.BOTH, expand=True)
    
    # Title
    title_label = tk.Label(main_frame, text="Multi-ROI Contour Analysis", 
                          font=('Arial', 16, 'bold'))
    title_label.pack(pady=(0, 20))
    
    # Instructions
    instr_label = tk.Label(main_frame, 
                          text="Select 1-4 Target ROIs to analyze for errant contours:",
                          font=('Arial', 10))
    instr_label.pack(pady=(0, 15))
    
    # ROI selection frame
    roi_frame = tk.LabelFrame(main_frame, text="Target ROIs", font=('Arial', 10, 'bold'))
    roi_frame.pack(fill=tk.X, pady=(0, 15))
    
    # Create 4 Combobox widgets for ROI selection
    roi_vars = []
    roi_combos = []
    
    for i in range(4):
        row_frame = tk.Frame(roi_frame)
        row_frame.pack(fill=tk.X, padx=10, pady=5)
        
        label = tk.Label(row_frame, text=f"ROI {i+1}:", width=8, anchor='w')
        label.pack(side=tk.LEFT)
        
        var = tk.StringVar()
        roi_vars.append(var)
        
        combo = ttk.Combobox(row_frame, textvariable=var, values=[""] + target_rois,
                            state="readonly", width=40)
        combo.pack(side=tk.LEFT, padx=(10, 0))
        roi_combos.append(combo)
    
    # External ROI selection frame
    ext_frame = tk.LabelFrame(main_frame, text="External ROI (for boundary checking)", 
                             font=('Arial', 10, 'bold'))
    ext_frame.pack(fill=tk.X, pady=(0, 20))
    
    ext_row_frame = tk.Frame(ext_frame)
    ext_row_frame.pack(fill=tk.X, padx=10, pady=5)
    
    ext_label = tk.Label(ext_row_frame, text="External:", width=8, anchor='w')
    ext_label.pack(side=tk.LEFT)
    
    # Pre-select "External" if available
    default_external = "External" if "External" in external_rois else ""
    if default_external:
        selected_external.set(default_external)
    
    ext_combo = ttk.Combobox(ext_row_frame, textvariable=selected_external, 
                            values=[""] + external_rois,
                            state="readonly", width=40)
    ext_combo.pack(side=tk.LEFT, padx=(10, 0))
    
    # Button frame
    button_frame = tk.Frame(main_frame)
    button_frame.pack(fill=tk.X, pady=(10, 0))
    
    def validate_and_proceed():
        print("Validating selections...")
        # Get selected ROIs (remove empty selections)
        selected = [var.get() for var in roi_vars if var.get().strip()]
        
        if not selected:
            messagebox.showerror("Error", "Please select at least one ROI!")
            return
        
        # Check for duplicates
        if len(selected) != len(set(selected)):
            messagebox.showerror("Error", "Please select different ROIs (no duplicates)!")
            return
        
        print(f"Selected ROIs: {selected}")
        
        # Now do the more expensive contour count check only for selected ROIs
        print("Verifying selected ROIs have contours...")
        valid_selections = []
        for roi_name in selected:
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
        
        if not valid_selections:
            messagebox.showerror("Error", "None of the selected ROIs have contours!")
            return
        
        if len(valid_selections) != len(selected):
            missing = set(selected) - set(valid_selections)
            messagebox.showwarning("Warning", 
                f"These ROIs have no contours and will be skipped: {', '.join(missing)}")
        
        # Store selections
        selected_rois.extend(valid_selections)
        
        # Close window
        root.quit()
        root.destroy()
    
    def cancel_selection():
        root.quit()
        root.destroy()
    
    # Buttons
    tk.Button(button_frame, text="Cancel", command=cancel_selection,
              width=15, bg='lightcoral').pack(side=tk.RIGHT, padx=(10, 0))
    
    tk.Button(button_frame, text="Analyze ROIs", command=validate_and_proceed,
              width=15, bg='lightgreen').pack(side=tk.RIGHT)
    
    # Summary info
    info_text = f"Available Target ROIs: {len(target_rois)}\n"
    info_text += f"Available External ROIs: {len(external_rois)}"
    
    info_label = tk.Label(main_frame, text=info_text, font=('Arial', 9), 
                         fg='gray', justify=tk.LEFT)
    info_label.pack(side=tk.BOTTOM, anchor='w', pady=(10, 0))
    
    print(f"GUI creation completed in {time.time() - gui_start:.2f}s")
    print(f"Total time to show GUI: {time.time() - start_time:.2f}s")
    
    # Start the GUI
    root.mainloop()
    
    # Return selections
    if selected_rois:
        external_roi = selected_external.get() if selected_external.get() else "External"
        return selected_rois, external_roi
    else:
        return None, None

# Get user selections
print("Opening ROI selection dialog...")
roi_names, external_roi_name = get_roi_selection()

if roi_names is None:
    print("No ROIs selected. Exiting.")
    exit()

print(f"Selected ROIs: {roi_names}")
print(f"External ROI: {external_roi_name}")

# 3. Process each ROI and collect data
roi_data = {}  # Dictionary to store data for each ROI

for roi_name in roi_names:
    print(f"\nProcessing ROI: {roi_name}")
    
    # Get the geometric‐contour container for this ROI
    try:
        geom = case.PatientModel.\
                 StructureSets[exam.Name].\
                 RoiGeometries[roi_name]. \
                 PrimaryShape
        
        contours = geom.Contours
        # contours is a list of loops; each loop is itself a list of points with .X, .Y, .Z
        
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

# Process External ROI for boundary checking
external_contours = {}  # Dictionary to store external contours by z-slice
try:
    print(f"\nProcessing External ROI: {external_roi_name}")
    external_geom = case.PatientModel.\
             StructureSets[exam.Name].\
             RoiGeometries[external_roi_name]. \
             PrimaryShape
    
    external_contour_list = external_geom.Contours
    print(f"  {len(external_contour_list)} external contours found")
    
    # Get all z-values from our target ROIs to know which external slices we need
    target_z_values = set()
    for roi_name, data in roi_data.items():
        target_z_values.update(data['z_values'])
    
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

# Calculate global z-range and slice thickness from all ROIs
all_z_values = []
for roi_name, data in roi_data.items():
    all_z_values.extend(data['z_values'])

if not all_z_values:
    print("No contours found in any ROI!")
    exit()

# Sort z-values and calculate slice thickness
ordered_z_values = sorted(all_z_values)
z_differences = [round(a - b, 2) for a, b in zip(ordered_z_values[1:], ordered_z_values[:-1])]
slice_thickness = min([diff for diff in z_differences if diff != 0])

print(f"\nGlobal Analysis:")
print(f"Slice thickness: {slice_thickness} mm")
print(f"Z-range: {min(all_z_values):.1f} to {max(all_z_values):.1f} mm")

# Print analysis for each ROI
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


# Create GUI with visualizations for multiple ROIs
def create_gui():
    root = tk.Tk()
    roi_names_str = ", ".join(roi_names)
    root.title(f"Multi-ROI Contour Analysis - {roi_names_str}")
    root.geometry("1600x1000")
    
    # Filter roi_names to only include ROIs that have data
    valid_roi_names = [roi for roi in roi_names if roi in roi_data]
    num_rois = len(valid_roi_names)
    
    if num_rois == 0:
        print("No valid ROI data found!")
        return
    
    # Create figure with 3 rows of subplots (contours per slice, area distribution, slice viewer)
    # Adjust subplot heights: row 1 normal, row 2 shorter (2/3 height), row 3 normal
    fig = plt.figure(figsize=(16, 10))
    
    # Define custom subplot grid with height ratios [1, 0.67, 1]
    gs = fig.add_gridspec(3, num_rois, height_ratios=[1, 0.67, 1])
    
    # Calculate global z-range for consistent y-axis
    min_z = min(all_z_values)
    max_z = max(all_z_values)
    
    # Create complete range of z-values at slice thickness intervals
    expected_z_values = []
    current_z = min_z
    while current_z <= max_z + slice_thickness/2:
        expected_z_values.append(round(current_z, 2))
        current_z += slice_thickness
    
    y_positions = list(range(len(expected_z_values)))
    y_labels = [f"{z:.1f}" for z in expected_z_values]
    
    # Track global issues
    global_has_multiple_contours = False
    global_has_isolated_contours = False
    global_has_missing_contours = False
    global_has_external_violations = False
    
    # Row 1: Contours per slice plots
    contour_axes = []
    for roi_idx, roi_name in enumerate(valid_roi_names):
        # Create subplot using gridspec (row 0, column roi_idx)
        ax = fig.add_subplot(gs[0, roi_idx])
        contour_axes.append(ax)
        
        data = roi_data[roi_name]
        
        # Create counts array with 0 for missing slices
        counts = []
        for expected_z in expected_z_values:
            counts.append(data['contours_per_slice'].get(expected_z, 0))
        
        # Create bars for this ROI (no gaps between bars)
        bars = ax.barh(y_positions, counts, height=1.0, 
                      color='skyblue', alpha=0.7, 
                      edgecolor='navy', linewidth=0.5)
        
        # Detect isolated contours for this ROI
        isolation_threshold = 2 * slice_thickness
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
        
        # Color bars based on issues
        has_multiple_contours = False
        has_isolated_contours = len(isolated_slices) > 0
        has_missing_contours = any(count == 0 for count in counts)
        
        for i, count in enumerate(counts):
            z_value = expected_z_values[i]
            if count == 0:  # Missing slice
                bars[i].set_color('lightgray')
                bars[i].set_alpha(0.3)
            elif count > 1:  # Multiple contours on same slice
                bars[i].set_color('red')
                bars[i].set_alpha(0.8)
                has_multiple_contours = True
            elif z_value in isolated_slices:  # Isolated contour
                bars[i].set_color('orange')
                bars[i].set_alpha(0.8)
        
        # Update global flags
        global_has_multiple_contours |= has_multiple_contours
        global_has_isolated_contours |= has_isolated_contours
        global_has_missing_contours |= has_missing_contours
        
        # Remove the value labels on bars to reduce clutter
        
        # Set up axes for this subplot
        ax.set_title(f'{roi_name}\nContours per Slice', fontsize=10)
        ax.grid(axis='x', alpha=0.3)
        
        # Only show y-axis labels on the leftmost plot (every 3rd slice)
        if roi_idx == 0:
            # Show every 3rd slice to reduce crowding
            sparse_positions = y_positions[::3]  # Every 3rd position
            sparse_labels = [y_labels[i] for i in range(0, len(y_labels), 3)]
            ax.set_yticks(sparse_positions)
            ax.set_yticklabels(sparse_labels)
            ax.set_ylabel('Slice Z-position (mm)')
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])
        
        # Set x-axis
        max_count = max(counts) if counts else 1
        ax.set_xticks(range(max_count + 1))
        ax.set_xlim(0, max_count + 0.5)
        ax.set_xlabel('# Contours')
    
    # Row 2: Area distribution plots (shorter height)
    area_axes = []
    for roi_idx, roi_name in enumerate(valid_roi_names):
        # Create subplot using gridspec (row 1, column roi_idx) - this will be 2/3 height
        ax = fig.add_subplot(gs[1, roi_idx])
        area_axes.append(ax)
        
        data = roi_data[roi_name]
        contour_areas = data['contour_areas']
        
        if contour_areas:
            # Create histogram
            max_area = max(contour_areas)
            bins = range(0, int(max_area) + 2, 1)
            
            n, bins_used, patches = ax.hist(contour_areas, bins=bins, 
                     alpha=0.7, color='lightgreen',
                     edgecolor='darkgreen', linewidth=0.5)
            
            # Color bars red for areas < 2mm²
            for i, patch in enumerate(patches):
                if bins_used[i] < 2:
                    patch.set_color('red')
                    patch.set_alpha(0.8)
        
        # Set up axes for this subplot
        ax.set_title(f'{roi_name}\nArea Distribution', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xlabel('Contour Area (mm²)')
        
        # Only show y-axis labels on the leftmost plot
        if roi_idx == 0:
            ax.set_ylabel('Frequency')
        else:
            ax.set_yticklabels([])
    
    # Align y-axes for contour plots (top row)
    for ax in contour_axes:
        ax.set_ylim(-0.5, len(expected_z_values) - 0.5)
    
    # Align y-axes for area plots (bottom row) - find max frequency
    max_frequency = 0
    for roi_name in valid_roi_names:
        if roi_name in roi_data and roi_data[roi_name]['contour_areas']:
            areas = roi_data[roi_name]['contour_areas']
            max_area = max(areas) if areas else 0
            # Simple frequency count for each bin
            for area in areas:
                bin_count = sum(1 for a in areas if int(a) == int(area))
                max_frequency = max(max_frequency, bin_count)
    
    for ax in area_axes:
        ax.set_ylim(0, max_frequency + 1)
    
    # Row 3: Interactive slice viewer plots
    slice_axes = []
    current_slice_indices = {}  # Track current slice for each ROI
    prerendered_slices = {}  # Store pre-rendered slice data for fast switching
    slice_artists = {}  # Store matplotlib artists for blitting
    slice_backgrounds = {}  # Store axis backgrounds for blitting
    
    for roi_idx, roi_name in enumerate(valid_roi_names):
        # Create subplot using gridspec (row 2, column roi_idx)
        ax = fig.add_subplot(gs[2, roi_idx])
        slice_axes.append(ax)
        
        # Initialize to show the first slice (slice 1)
        data = roi_data[roi_name]
        available_slices = sorted(data['contours_per_slice'].keys())
        if available_slices:
            current_slice_indices[roi_name] = 0  # Start at first slice (slice 1)
        else:
            current_slice_indices[roi_name] = 0
        
        # Pre-process all slices for this ROI to speed up navigation
        prerendered_slices[roi_name] = {}
        
        # Calculate axis limits once for consistency and speed
        all_x_coords = []
        all_y_coords = []
        max_contours_per_slice = 0  # Track max contours to pre-allocate artists
        
        for slice_idx, z_slice in enumerate(available_slices):
            # Find contours at this z-level
            contours_at_z = []
            contour_areas_at_z = []
            
            for i, z_val in enumerate(data['z_values']):
                if abs(z_val - z_slice) < 0.01:  # Account for floating point precision
                    contours_at_z.append(data['contours'][i])
                    contour_areas_at_z.append(data['contour_areas'][i])
            
            # Pre-process contour coordinates and collect bounds
            processed_contours = []
            for contour, area in zip(contours_at_z, contour_areas_at_z):
                # Extract x, y coordinates
                x_coords = []
                y_coords = []
                
                for point in contour:
                    if hasattr(point, 'x'):
                        x_coords.append(point.x)
                        y_coords.append(-point.y)  # Flip y-coordinate to match RayStation orientation
                    else:
                        x_coords.append(point['x'])
                        y_coords.append(-point['y'])  # Flip y-coordinate to match RayStation orientation
                
                # Close the contour
                if x_coords and y_coords:
                    x_coords.append(x_coords[0])
                    y_coords.append(y_coords[0])
                    
                    # Collect coordinates for bounds calculation
                    all_x_coords.extend(x_coords)
                    all_y_coords.extend(y_coords)
                    
                    # Split contour into inside/outside segments for partial highlighting
                    segments = split_contour_by_external(
                        x_coords, y_coords, external_contours.get(z_slice, []))
                    
                    # Check if any part is outside external boundary
                    outside_external = any(seg['type'] == 'outside' for seg in segments)
                    
                    # Add each segment as a separate contour element
                    for segment in segments:
                        # Determine color based on area and segment type
                        if area < 2.0:
                            color = 'red'  # Small area always red
                            alpha = 0.8
                        elif segment['type'] == 'outside':
                            color = 'red'  # Outside external always red
                            alpha = 0.8
                        else:
                            color = 'blue'  # Normal contour
                            alpha = 0.6
                        
                        processed_contours.append({
                            'x_coords': segment['x_coords'],
                            'y_coords': segment['y_coords'],
                            'area': area,
                            'segment_type': segment['type'],
                            'outside_external': outside_external,
                            'color': color,
                            'alpha': alpha
                        })
            
            # Add external contour outlines for visual reference (don't scale to them)
            external_contours_at_z = external_contours.get(z_slice, [])
            for ext_contour in external_contours_at_z:
                # Add external contour as a reference (gray line, no fill)
                processed_contours.append({
                    'x_coords': ext_contour['x_coords'] + [ext_contour['x_coords'][0]],  # Close contour
                    'y_coords': ext_contour['y_coords'] + [ext_contour['y_coords'][0]],
                    'area': 0,  # External contour doesn't count as ROI area
                    'segment_type': 'external',
                    'outside_external': False,
                    'color': 'gray',
                    'alpha': 0.5,
                    'is_external': True  # Flag to identify external contours
                })
            
            # Track maximum contours per slice for artist pre-allocation
            max_contours_per_slice = max(max_contours_per_slice, len(processed_contours))
            
            # Store slice info with external boundary violations
            num_contours = len(contours_at_z)
            small_contours = sum(1 for area in contour_areas_at_z if area < 2.0)
            # Count unique contours that have external violations (not individual segments)
            contours_with_violations = set()
            for contour in processed_contours:
                if contour.get('outside_external', False) and not contour.get('is_external', False):
                    # Use area as a simple way to group segments from the same original contour
                    contours_with_violations.add(contour['area'])
            external_violations = len(contours_with_violations)
            
            # Track global external violations
            if external_violations > 0:
                global_has_external_violations = True
            
            slice_info = f'Slice {slice_idx + 1}/{len(available_slices)}'
            if small_contours > 0 or external_violations > 0:
                details = []
                if small_contours > 0:
                    details.append(f'{small_contours} small')
                if external_violations > 0:
                    details.append(f'{external_violations} outside')
                slice_info += f'\n{num_contours} contours ({", ".join(details)})'
            else:
                slice_info += f'\n{num_contours} contours'
            
            prerendered_slices[roi_name][slice_idx] = {
                'contours': processed_contours,
                'title': f'{roi_name}\n{slice_info}',
                'z_value': z_slice
            }
        
        # Set up axis limits and properties once (for consistency and speed)
        if all_x_coords and all_y_coords:
            margin = 5  # mm margin
            ax.set_xlim(min(all_x_coords) - margin, max(all_x_coords) + margin)
            ax.set_ylim(min(all_y_coords) - margin, max(all_y_coords) + margin)
        
        # Set up axis properties
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Remove x and y axes for cleaner display
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        
        # Set title above the plot (not inside)
        ax.set_title('', fontsize=9)  # Will be updated by update_slice_display
        
        # Pre-allocate matplotlib artists for ALL slices (hide/show approach)
        # This is more efficient than updating artist data - we just toggle visibility
        slice_artists[roi_name] = {
            'slices': {},  # Store artists for each slice
            'ax': ax
        }
        
        # Create artists for each slice
        for slice_idx in range(len(available_slices)):
            if slice_idx in prerendered_slices[roi_name]:
                slice_data = prerendered_slices[roi_name][slice_idx]
                contours = slice_data['contours']
                
                # Create artists for this specific slice
                slice_line_artists = []
                slice_fill_artists = []
                
                for contour_data in contours:
                    # Create line artist
                    line, = ax.plot(contour_data['x_coords'], contour_data['y_coords'], 
                                  color=contour_data['color'], 
                                  linewidth=3 if contour_data.get('is_external', False) else 2,
                                  alpha=contour_data['alpha'], visible=False)
                    slice_line_artists.append(line)
                    
                    # Create fill artist (only for non-external contours)
                    if not contour_data.get('is_external', False):
                        fill = ax.fill(contour_data['x_coords'], contour_data['y_coords'], 
                                     color=contour_data['color'], alpha=0.2, visible=False)[0]
                        slice_fill_artists.append(fill)
                    else:
                        # For external contours, add a dummy invisible fill to maintain indexing
                        fill = ax.fill([], [], alpha=0, visible=False)[0]
                        slice_fill_artists.append(fill)
                
                # Store artists for this slice
                slice_artists[roi_name]['slices'][slice_idx] = {
                    'line_artists': slice_line_artists,
                    'fill_artists': slice_fill_artists,
                    'title': slice_data['title']
                }
        
        # Set proper axis limits to make contours fill the plot area
        if all_x_coords and all_y_coords:
            margin_x = (max(all_x_coords) - min(all_x_coords)) * 0.05  # 5% margin
            margin_y = (max(all_y_coords) - min(all_y_coords)) * 0.05
            ax.set_xlim(min(all_x_coords) - margin_x, max(all_x_coords) + margin_x)
            ax.set_ylim(min(all_y_coords) - margin_y, max(all_y_coords) + margin_y)
        else:
            # Fallback for empty ROIs
            ax.set_xlim(-50, 50)
            ax.set_ylim(-50, 50)
    
    # Function to update slice display (optimized hide/show approach)
    def update_slice_display(roi_name, slice_index):
        artists = slice_artists[roi_name]
        ax = artists['ax']
        
        # Hide all artists for this ROI first
        for slice_idx, slice_artists_data in artists['slices'].items():
            for line_artist in slice_artists_data['line_artists']:
                line_artist.set_visible(False)
            for fill_artist in slice_artists_data['fill_artists']:
                fill_artist.set_visible(False)
        
        # Show artists for the current slice and update title
        if slice_index in artists['slices']:
            slice_artists_data = artists['slices'][slice_index]
            for line_artist in slice_artists_data['line_artists']:
                line_artist.set_visible(True)
            for fill_artist in slice_artists_data['fill_artists']:
                fill_artist.set_visible(True)
            
            # Update title above the plot
            ax.set_title(slice_artists_data['title'], fontsize=9)
        else:
            # No data for this slice
            ax.set_title(f'{roi_name}\nNo Data', fontsize=9)
    
    # Initialize all slice displays
    for roi_name in valid_roi_names:
        update_slice_display(roi_name, current_slice_indices[roi_name])
    
    # Initial canvas refresh to ensure proper rendering
    fig.canvas.draw_idle()
    
    # Mouse click handler for slice navigation with boundary checking
    def on_click(event):
        if event.inaxes in slice_axes:
            roi_idx = slice_axes.index(event.inaxes)
            roi_name = valid_roi_names[roi_idx]
            
            data = roi_data[roi_name]
            available_slices = sorted(data['contours_per_slice'].keys())
            if not available_slices:
                return
            
            current_idx = current_slice_indices[roi_name]
            max_idx = len(available_slices) - 1
            
            # Left click: next slice (with boundary check)
            if event.button == 1:  # Left click
                if current_idx < max_idx:  # Only advance if not at last slice
                    current_slice_indices[roi_name] = current_idx + 1
                    update_slice_display(roi_name, current_slice_indices[roi_name])
                    # Only redraw the specific subplot region for speed
                    roi_idx = valid_roi_names.index(roi_name)
                    slice_axes[roi_idx].figure.canvas.draw_idle()
            
            # Right click: previous slice (with boundary check)
            elif event.button == 3:  # Right click
                if current_idx > 0:  # Only go back if not at first slice
                    current_slice_indices[roi_name] = current_idx - 1
                    update_slice_display(roi_name, current_slice_indices[roi_name])
                    # Only redraw the specific subplot region for speed
                    roi_idx = valid_roi_names.index(roi_name)
                    slice_axes[roi_idx].figure.canvas.draw_idle()
    
    # Keyboard handler for slice navigation with boundary checking
    def on_key(event):
        if event.key in ['up', 'down', 'left', 'right']:
            # Apply to all ROIs simultaneously with boundary checking
            for roi_name in valid_roi_names:
                available_slices = sorted(roi_data[roi_name]['contours_per_slice'].keys())
                if not available_slices:
                    continue
                
                current_idx = current_slice_indices[roi_name]
                max_idx = len(available_slices) - 1
                
                if event.key in ['up', 'right']:
                    if current_idx < max_idx:  # Only advance if not at last slice
                        current_slice_indices[roi_name] = current_idx + 1
                elif event.key in ['down', 'left']:
                    if current_idx > 0:  # Only go back if not at first slice
                        current_slice_indices[roi_name] = current_idx - 1
                
                update_slice_display(roi_name, current_slice_indices[roi_name])
            
            # Use draw_idle() for better performance
            fig.canvas.draw_idle()
    
    # Connect event handlers (click and keyboard only)
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    plt.tight_layout()
    
    # Create a scrollable frame for the entire content
    # Main canvas for scrolling
    main_canvas = tk.Canvas(root)
    scrollbar = tk.Scrollbar(root, orient="vertical", command=main_canvas.yview)
    scrollable_frame = tk.Frame(main_canvas)
    
    scrollable_frame.bind(
        "<Configure>",
        lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
    )
    
    main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    main_canvas.configure(yscrollcommand=scrollbar.set)
    
    # Add mouse wheel scrolling support
    def on_mousewheel(event):
        main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
    
    # Bind mouse wheel to canvas (Windows)
    main_canvas.bind("<MouseWheel>", on_mousewheel)
    # For Linux
    main_canvas.bind("<Button-4>", lambda e: main_canvas.yview_scroll(-1, "units"))
    main_canvas.bind("<Button-5>", lambda e: main_canvas.yview_scroll(1, "units"))
    
    # Make sure the canvas can receive focus for scrolling
    main_canvas.focus_set()
    
    # Pack the canvas and scrollbar
    main_canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")
    
    # Embed plot in the scrollable frame
    canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
    canvas.draw()
    
    # Add error detection banner at the top if issues found
    if global_has_multiple_contours or global_has_isolated_contours or global_has_missing_contours or global_has_external_violations:
        error_text = "⚠️ Possible Errant Contours Detected"
        details = []
        if global_has_missing_contours:
            details.append("Gaps found")
        if global_has_isolated_contours:
            details.append("Isolated contours (orange)")
        if global_has_multiple_contours:
            details.append("Multiple contours per slice (red)")
        if global_has_external_violations:
            details.append("Contours outside External (red)")
        
        if details:
            error_text += f" ({', '.join(details)})"
        
        error_frame = tk.Frame(scrollable_frame, bg='red', relief=tk.RAISED, bd=2)
        error_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(5, 5))
        
        error_label = tk.Label(error_frame, text=error_text, bg='red', fg='white', 
                              font=('Arial', 12, 'bold'))
        error_label.pack(pady=5)
    
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(5, 5))
    
    # Re-connect event handlers after canvas is embedded
    canvas.mpl_connect('button_press_event', on_click)
    canvas.mpl_connect('key_press_event', on_key)
    
    # Make sure the matplotlib canvas can receive focus for keyboard events
    canvas.get_tk_widget().focus_set()
    
    # Create statistics text frame at bottom with separate rows for each ROI
    stats_frame = tk.Frame(scrollable_frame, bg='wheat', relief=tk.RAISED, bd=2)
    stats_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
    
    # Build statistics text with separate rows for each ROI
    for roi_idx, roi_name in enumerate(roi_names):
        if roi_name in roi_data:
            data = roi_data[roi_name]
            stats_text = f"{roi_name}: {data['number_of_contours']} contours | "
            stats_text += f"{data['gaps']} gaps | "
            stats_text += f"Total Area: {data['total_area']:.0f}mm² | "
            stats_text += f"Avg Area: {data['avg_area']:.0f}mm² | "
            stats_text += f"Min: {data['min_area']:.0f}mm² | "
            stats_text += f"Max: {data['max_area']:.0f}mm²"
            
            stats_label = tk.Label(stats_frame, text=stats_text, bg='wheat', font=('Arial', 9))
            stats_label.pack(pady=2)
    
    # Add slice thickness info and navigation instructions
    slice_info = tk.Label(stats_frame, text=f"Slice Thickness: {slice_thickness}mm", 
                         bg='wheat', font=('Arial', 9, 'bold'))
    slice_info.pack(pady=2)
    
    # Add navigation instructions
    nav_info = tk.Label(stats_frame, 
                       text="Slice Navigation: Left Click = Next | Right Click = Previous | Arrow Keys = Navigate All ROIs", 
                       bg='wheat', font=('Arial', 8, 'italic'))
    nav_info.pack(pady=1)
    
    root.mainloop()

# Launch GUI
create_gui()