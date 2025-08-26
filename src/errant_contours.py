from connect import get_current
from config import debug
import config
import utils
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from collections import Counter

# Import split_contour_by_external from utils
from utils import split_contour_by_external

if debug: print("Script started")  # Top-level entry

class PlanData:
    def __init__(self):
        self.patient = None
        self.case = None
        self.exam = None

        self.roi_list = []
        self.ptv_list = []
        self.external = None

        self.slice_thickness = None

        self.temp_rois = []

        self.build()

    def build(self):
        """Build the plan data structure."""
        self.patient = get_current("Patient")
        self.case = get_current("Case")
        self.exam = get_current("Examination")

        self.roi_list = [roi for roi in self.case.PatientModel.RegionsOfInterest]
        # Only include ROIs with Type == "Ptv"
        self.ptv_list = [roi for roi in self.roi_list if getattr(roi, "Type", None) == "Ptv"]
        # External ROI name (if available)
        self.external = getattr(config, "ExternalName", None)

class RoiData:
    def __init__(self, roi, plan_data):
        self.Name = roi.Name
        self.contours = None
        self.number_of_contours = 0
        self.z_values = []
        self.contour_areas = []
        self.contours_per_slice = Counter()
        self.gaps = 0
        self.total_area = 0.0
        self.avg_area = 0.0
        self.min_area = 0.0
        self.max_area = 0.0
        self.build(roi, plan_data)

    def build(self, roi, plan_data):
        self.contours = getattr(plan_data.case.PatientModel.StructureSets[plan_data.exam.Name].RoiGeometries[roi.Name].PrimaryShape, "Contours", None)

def get_roi_selection(plan_data):
    if debug: print("Entering get_roi_selection")
    import tkinter as tk
    from tkinter import ttk, messagebox
    import time

    if debug: print("Starting ROI selection GUI...")
    start_time = time.time()

    # Get viable PTVs for user to select
    viable_rois = []
    try:
        if debug: print("Getting ROI list...")
        roi_start = time.time()
        if debug: print(f"Found {len(plan_data.roi_list)} total ROIs in {time.time() - roi_start:.2f}s")
        if debug: print("Filtering ROIs and checking for contours...")
        filter_start = time.time()
        for i, roi in enumerate(plan_data.ptv_list):
            if i % 10 == 0:
                if debug: print(f"  Processing ROI {i+1}/{len(plan_data.ptv_list)}: {roi.Name}")
            try:
                roi_geometry = plan_data.case.PatientModel.StructureSets[plan_data.exam.Name].RoiGeometries[roi.Name]
                has_contours = (hasattr(roi_geometry, 'PrimaryShape') and 
                               roi_geometry.PrimaryShape is not None)
            except:
                has_contours = False
            if has_contours:
                viable_rois.append(roi.Name)
        if debug: print(f"ROI filtering completed in {time.time() - filter_start:.2f}s")
        if debug: print(f"Found {len(viable_rois)} target ROIs")
        viable_rois.sort()
    except Exception as e:
        if debug: print(f"Error during ROI processing: {e}")
        messagebox.showerror("Error", f"Could not retrieve ROI list: {e}")
        return None
    if not viable_rois:
        messagebox.showwarning("No Target ROIs", "No Target-type ROIs with geometry found!")
        return None
    if debug: print(f"Starting GUI creation after {time.time() - start_time:.2f}s total...")
    gui_start = time.time()

    # Create selection window
    root = tk.Tk()
    root.title("Select ROIs for Contour Analysis")
    root.geometry("500x400")
    root.resizable(False, False)

    # Variables to store selections
    selected_rois = []

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
        combo = ttk.Combobox(row_frame, textvariable=var, values=[""] + viable_rois,
                            state="readonly", width=40)
        combo.pack(side=tk.LEFT, padx=(10, 0))
        roi_combos.append(combo)

    # Button frame
    button_frame = tk.Frame(main_frame)
    button_frame.pack(fill=tk.X, pady=(10, 0))

    def validate_and_proceed():
        if debug: print("validate_and_proceed called")
        if debug: print("Validating selections...")
        selected = [var.get() for var in roi_vars if var.get().strip()]
        if not selected:
            messagebox.showerror("Error", "Please select at least one ROI!")
            return
        if len(selected) != len(set(selected)):
            messagebox.showerror("Error", "Please select different ROIs (no duplicates)!")
            return
        if debug: print(f"Selected ROIs: {selected}")
        selected_rois.extend(selected)
        if debug: print(f"Selected ROIs after validation: {selected_rois}")
        if debug: print("Closing ROI selection window")
        root.quit()
        root.destroy()

    def cancel_selection():
        if debug: print("ROI selection cancelled")
        root.quit()
        root.destroy()

    # Buttons
    tk.Button(button_frame, text="Cancel", command=cancel_selection,
              width=15, bg='lightcoral').pack(side=tk.RIGHT, padx=(10, 0))
    tk.Button(button_frame, text="Analyze ROIs", command=validate_and_proceed,
              width=15, bg='lightgreen').pack(side=tk.RIGHT)

    # Summary info
    info_text = f"Available Target ROIs: {len(viable_rois)}\n"
    info_label = tk.Label(main_frame, text=info_text, font=('Arial', 9), 
                         fg='gray', justify=tk.LEFT)
    info_label.pack(side=tk.BOTTOM, anchor='w', pady=(10, 0))

    if debug: print(f"GUI creation completed in {time.time() - gui_start:.2f}s")
    if debug: print(f"Total time to show GUI: {time.time() - start_time:.2f}s")

    # Start the GUI
    if debug: print("ROI selection GUI ready, starting mainloop")
    root.mainloop()
    if debug: print("ROI selection GUI closed")

    final_rois = []

    # Return selections
    if selected_rois:
        if debug: print(f"Ensuring contour attribute exists for selected ROIs")
        for roi in selected_rois:
            new_roi_name = utils.ensure_roi_has_contours(plan_data, roi)
            if new_roi_name is None:
                final_rois.append(roi)
            elif new_roi_name != roi:
                final_rois.append(new_roi_name)
            else:
                final_rois.append(roi)
        if debug: print(f"Returning selected_rois: {final_rois}")
        return final_rois
    else:
        return None

def analysis(plan_data, roi_data):
    if debug: print("Starting ROI analysis loop")
    for roi_name, roi_obj in roi_data.items():
        if debug: print(f"Processing ROI: {roi_name}")
        try:
            contours = roi_obj.contours
            number_of_contours = len(contours) if contours else 0
            if debug: print(f"  {number_of_contours} contours found")

            z_values = []
            contour_areas = []

            for contour in contours or []:
                z_values.append(round(contour[0]['z'], 2))
                contour_areas.append(utils.calculate_contour_area(contour))

            roi_obj.number_of_contours = number_of_contours
            roi_obj.z_values = z_values
            roi_obj.contour_areas = contour_areas
            roi_obj.contours_per_slice = Counter(z_values)

            if debug: print(f"Finished processing ROI: {roi_name}")
        except Exception as e:
            if debug: print(f"  Error processing {roi_name}: {e}")
            continue

    if debug: print("Finished ROI analysis loop")

    if debug: print("Processing External ROI for boundary checking")
    external_contours = {}
    try:
        external_roi_name = plan_data.external
        if not external_roi_name or external_roi_name not in [roi.Name for roi in plan_data.roi_list]:
            print(f"WARNING: External ROI '{external_roi_name}' not found in case. External contours will not be shown.")
        else:
            if debug: print(f"\nProcessing External ROI: {external_roi_name}")
            external_geom = plan_data.case.PatientModel.\
                    StructureSets[plan_data.exam.Name].\
                    RoiGeometries[external_roi_name]. \
                    PrimaryShape

            external_contour_list = getattr(external_geom, "Contours", [])
            if debug: print(f"  {len(external_contour_list)} external contours found")

            for contour in external_contour_list:
                z_val = round(contour[0]['z'], 2)
                x_coords = []
                y_coords = []
                for point in contour:
                    if hasattr(point, 'x'):
                        x_coords.append(point.x)
                        y_coords.append(-point.y)
                    else:
                        x_coords.append(point['x'])
                        y_coords.append(-point['y'])
                if x_coords and y_coords:
                    if z_val not in external_contours:
                        external_contours[z_val] = []
                    external_contours[z_val].append({
                        'x_coords': x_coords,
                        'y_coords': y_coords
                    })
                    if debug: print(f"Added external contour for z={z_val}: {len(x_coords)} points")
                else:
                    if debug: print(f"WARNING: External contour for z={z_val} has empty coordinates!")
            if debug: print(f"  External contours processed for {len(external_contours)} z-slices")
            if debug:
                for z, contours in external_contours.items():
                    print(f"External z={z}: {len(contours)} contour(s)")
    except Exception as e:
        print(f"  Warning: Could not process External ROI '{plan_data.external}': {e}")
        print("  Boundary checking will be disabled.")

    all_z_values = []
    for roi_obj in roi_data.values():
        all_z_values.extend(roi_obj.z_values)

    if not all_z_values:
        print("No contours found in any ROI!")
        exit()

    ordered_z_values = sorted(all_z_values)
    z_differences = [round(a - b, 2) for a, b in zip(ordered_z_values[1:], ordered_z_values[:-1])]
    slice_thickness = min([diff for diff in z_differences if diff != 0]) if z_differences else 0

    print(f"\nGlobal Analysis:")
    print(f"Slice thickness: {slice_thickness} mm")
    print(f"Z-range: {min(all_z_values):.1f} to {max(all_z_values):.1f} mm")

    for roi_name, roi_obj in roi_data.items():
        if debug: print(f"\n{roi_name} Analysis:")
        
        ordered_z = sorted(roi_obj.z_values)
        unique_z = sorted(set(roi_obj.z_values))
        z_diffs = [round(a - b, 2) for a, b in zip(ordered_z[1:], ordered_z[:-1])]
        gaps = sum(diff > slice_thickness for diff in z_diffs)
        
        total_area = sum(roi_obj.contour_areas)
        avg_area = total_area / len(roi_obj.contour_areas) if roi_obj.contour_areas else 0
        min_area = min(roi_obj.contour_areas) if roi_obj.contour_areas else 0
        max_area = max(roi_obj.contour_areas) if roi_obj.contour_areas else 0
        
        if debug: print(f"  Total contours: {roi_obj.number_of_contours}")
        if debug: print(f"  Unique slices: {len(unique_z)}")
        if debug: print(f"  Gaps found: {gaps}")
        if debug: print(f"  Total area: {total_area:.2f} mm²")
        if debug: print(f"  Average area: {avg_area:.2f} mm²")
        if debug: print(f"  Min area: {min_area:.2f} mm²")
        if debug: print(f"  Max area: {max_area:.2f} mm²")
        
        roi_obj.gaps = gaps
        roi_obj.total_area = total_area
        roi_obj.avg_area = avg_area
        roi_obj.min_area = min_area
        roi_obj.max_area = max_area

    # Return external_contours and slice_thickness for GUI
    return external_contours, slice_thickness, all_z_values

# Create GUI with visualizations for multiple ROIs
def create_gui(roi_data, external_contours, slice_thickness, all_z_values):
    if debug: print("Entered create_gui")
    if debug: print("Creating Tk root window")
    root = tk.Tk()
    roi_names_str = ", ".join(roi_data.keys())
    if debug: print(f"Setting window title: Multi-ROI Contour Analysis - {roi_names_str}")
    root.title(f"Multi-ROI Contour Analysis - {roi_names_str}")
    if debug: print("Setting window geometry")
    root.geometry("1600x1000")
    
    # ADD THIS HANDLER TO ENSURE MAINLOOP EXITS ON WINDOW CLOSE
    def on_window_close():
        if debug: print("Window close event triggered")
        root.quit()
        root.destroy()
    root.protocol("WM_DELETE_WINDOW", on_window_close)

    if debug: print("Filtering valid_roi_names")
    valid_roi_names = [roi for roi in roi_data.keys()]
    num_rois = len(valid_roi_names)
    if debug: print(f"Number of valid ROIs: {num_rois}")
    
    if num_rois == 0:
        if debug: print("No valid ROI data found! Returning from create_gui")
        return
    
    if debug: print("Creating matplotlib figure")
    fig = plt.figure(figsize=(16, 10))
    if debug: print("Creating gridspec for subplots")
    gs = fig.add_gridspec(3, num_rois, height_ratios=[1, 0.67, 1])
    
    if debug: print("Calculating global z-range")
    min_z = min(all_z_values)
    max_z = max(all_z_values)
    if debug: print(f"min_z={min_z}, max_z={max_z}")
    
    if debug: print("Building expected_z_values")
    expected_z_values = []
    current_z = min_z
    while current_z <= max_z + slice_thickness/2:
        expected_z_values.append(round(current_z, 2))
        current_z += slice_thickness
    if debug: print(f"expected_z_values: {expected_z_values}")
    
    y_positions = list(range(len(expected_z_values)))
    y_labels = [f"{z:.1f}" for z in expected_z_values]
    
    if debug: print("Initializing global issue flags")
    global_has_multiple_contours = False
    global_has_isolated_contours = False
    global_has_missing_contours = False
    global_has_external_violations = False
    
    if debug: print("Creating contour axes")
    contour_axes = []
    for roi_idx, roi_name in enumerate(valid_roi_names):
        if debug: print(f"Creating subplot for contour axes: {roi_name}")
        ax = fig.add_subplot(gs[0, roi_idx])
        contour_axes.append(ax)
        
        data = roi_data[roi_name]
        
        # Create counts array with 0 for missing slices
        counts = []
        for expected_z in expected_z_values:
            counts.append(data.contours_per_slice.get(expected_z, 0))
        
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
    if debug: print("Creating area axes")
    area_axes = []
    for roi_idx, roi_name in enumerate(valid_roi_names):
        if debug: print(f"Creating subplot for area axes: {roi_name}")
        ax = fig.add_subplot(gs[1, roi_idx])
        area_axes.append(ax)
        
        data = roi_data[roi_name]
        contour_areas = data.contour_areas
        
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
        if roi_name in roi_data and roi_data[roi_name].contour_areas:
            areas = roi_data[roi_name].contour_areas
            max_area = max(areas) if areas else 0
            # Simple frequency count for each bin
            for area in areas:
                bin_count = sum(1 for a in areas if int(a) == int(area))
                max_frequency = max(max_frequency, bin_count)
    
    for ax in area_axes:
        ax.set_ylim(0, max_frequency + 1)
    
    if debug: print("Creating slice axes")
    slice_axes = []
    current_slice_indices = {}
    prerendered_slices = {}
    slice_artists = {}
    slice_backgrounds = {}
    
    for roi_idx, roi_name in enumerate(valid_roi_names):
        if debug: print(f"Creating subplot for slice axes: {roi_name}")
        ax = fig.add_subplot(gs[2, roi_idx])
        slice_axes.append(ax)
        
        # Initialize to show the first slice (slice 1)
        data = roi_data[roi_name]
        available_slices = sorted(data.contours_per_slice.keys())
        if debug: print(f"Available slices for {roi_name}: {available_slices}")
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
            
            for i, z_val in enumerate(data.z_values):
                if abs(z_val - z_slice) < 0.01:  # Account for floating point precision
                    contours_at_z.append(data.contours[i])
                    contour_areas_at_z.append(data.contour_areas[i])
            
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
            if debug: print(f"Slice {slice_idx} z={z_slice}: {len(external_contours_at_z)} external contour(s)")
            for ext_contour in external_contours_at_z:
                if not ext_contour['x_coords'] or not ext_contour['y_coords']:
                    if debug: print(f"WARNING: External contour for z={z_slice} is empty!")
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
        if debug: print(f"update_slice_display called for {roi_name}, slice_index={slice_index}")
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
    
    if debug: print("Initializing all slice displays")
    for roi_name in valid_roi_names:
        if debug: print(f"Initializing slice display for {roi_name}")
        update_slice_display(roi_name, current_slice_indices[roi_name])
    
    if debug: print("Drawing initial canvas")
    fig.canvas.draw_idle()
    
    if debug: print("Setting up event handlers")
    def on_click(event):
        if debug: print(f"on_click event: {event}")
        if event.inaxes in slice_axes:
            roi_idx = slice_axes.index(event.inaxes)
            roi_name = valid_roi_names[roi_idx]
            
            data = roi_data[roi_name]
            available_slices = sorted(data.contours_per_slice.keys())
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
    
    def on_key(event):
        if debug: print(f"on_key event: {event}")
        if event.key in ['up', 'down', 'left', 'right']:
            # Apply to all ROIs simultaneously with boundary checking
            for roi_name in valid_roi_names:
                available_slices = sorted(roi_data[roi_name].contours_per_slice.keys())
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
    
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    if debug: print("Calling plt.tight_layout()")
    plt.tight_layout()
    
    if debug: print("Creating scrollable frame")
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
    
    if debug: print("Setting up mouse wheel scrolling")
    # Add mouse wheel scrolling support
    def on_mousewheel(event):
        if debug: print(f"on_mousewheel event: {event}")
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
    
    if debug: print("Embedding matplotlib figure in scrollable frame")
    canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
    canvas.draw()
    
    if debug: print("Checking for error detection banner")
    if global_has_multiple_contours or global_has_isolated_contours or global_has_missing_contours or global_has_external_violations:
        if debug: print("Displaying error detection banner")
        # Add error detection banner at the top if issues found
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
    
    if debug: print("Packing canvas widget")
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(5, 5))
    
    if debug: print("Reconnecting event handlers after embedding canvas")
    canvas.mpl_connect('button_press_event', on_click)
    canvas.mpl_connect('key_press_event', on_key)
    canvas.get_tk_widget().focus_set()
    
    if debug: print("Creating statistics frame")
    stats_frame = tk.Frame(scrollable_frame, bg='wheat', relief=tk.RAISED, bd=2)
    stats_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
    for roi_idx, roi_name in enumerate(valid_roi_names):
        if roi_name in roi_data:
            if debug: print(f"Adding stats for {roi_name}")
            data = roi_data[roi_name]
            stats_text = f"{roi_name}: {data.number_of_contours} contours | "
            stats_text += f"{data.gaps} gaps | "
            stats_text += f"Total Area: {data.total_area:.0f}mm² | "
            stats_text += f"Avg Area: {data.avg_area:.0f}mm² | "
            stats_text += f"Min: {data.min_area:.0f}mm² | "
            stats_text += f"Max: {data.max_area:.0f}mm²"
            
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
    
    if debug: print("About to start GUI mainloop")
    root.mainloop()
    if debug: print("Exited GUI mainloop")
    if debug: print("Exiting create_gui")

def main():
    plan_data = PlanData()
    roi_data = {}
    if debug: print("Opening ROI selection dialog...")
    selected_rois = get_roi_selection(plan_data)
    if debug: print(f"ROI selection returned: selected_rois={selected_rois}")

    if not selected_rois:
        if debug: print("No ROIs selected. Exiting.")
        exit()

    if debug: print(f"Selected ROIs: {selected_rois}")

    for roi_name in selected_rois:
        # Try to find ROI object in ptv_list, then roi_list, else create dummy
        roi_obj = next((roi for roi in plan_data.ptv_list if roi.Name == roi_name), None)
        if not roi_obj:
            roi_obj = next((roi for roi in plan_data.roi_list if roi.Name == roi_name), None)
        if not roi_obj:
            # Create a dummy object with Name attribute for RoiData to fetch geometry
            class DummyRoi:
                pass
            roi_obj = DummyRoi()
            roi_obj.Name = roi_name
        roi_data[roi_name] = RoiData(roi_obj, plan_data)

    external_contours, slice_thickness, all_z_values = analysis(plan_data, roi_data)

    create_gui(roi_data, external_contours, slice_thickness, all_z_values)

main()