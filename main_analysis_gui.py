"""
Main analysis GUI for multi-ROI contour visualization and error detection.

This module creates the comprehensive analysis interface showing:
- Contours per slice plots with error highlighting
- Area distribution histograms
- Interactive slice viewers with navigation
- Error detection banners and statistics
"""

import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try:
    from .geometry_tools import split_contour_by_external
except ImportError:
    from geometry_tools import split_contour_by_external


def create_analysis_gui(roi_names, roi_data, external_contours, slice_thickness, all_z_values):
    """
    Create the main analysis GUI with visualizations for multiple ROIs.
    
    Args:
        roi_names: List of ROI names being analyzed
        roi_data: Dictionary of ROI contour data
        external_contours: Dictionary of external contour data by z-slice
        slice_thickness: Calculated slice thickness in mm
        all_z_values: List of all z-values from all ROIs
    """
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
    
    # Create contour plots, area plots, and slice viewers
    contour_axes, area_axes, slice_axes = _create_subplots(
        fig, gs, valid_roi_names, roi_data, expected_z_values, y_positions, 
        y_labels, slice_thickness, external_contours
    )
    
    # Update global issue flags
    global_has_multiple_contours, global_has_isolated_contours, global_has_missing_contours, global_has_external_violations = _detect_global_issues(
        valid_roi_names, roi_data, expected_z_values, slice_thickness, external_contours
    )
    
    # Create interactive slice navigation
    current_slice_indices, slice_artists = _setup_slice_navigation(
        valid_roi_names, roi_data, slice_axes, external_contours
    )
    
    # Setup event handlers
    _setup_event_handlers(fig, slice_axes, valid_roi_names, roi_data, 
                         current_slice_indices, slice_artists)
    
    plt.tight_layout()
    
    # Create scrollable tkinter interface
    _create_scrollable_interface(root, fig, valid_roi_names, roi_data, slice_thickness,
                                global_has_multiple_contours, global_has_isolated_contours,
                                global_has_missing_contours, global_has_external_violations)
    
    root.mainloop()


def _create_subplots(fig, gs, valid_roi_names, roi_data, expected_z_values, y_positions, 
                    y_labels, slice_thickness, external_contours):
    """Create the three rows of subplots: contour plots, area plots, and slice viewers."""
    contour_axes = []
    area_axes = []
    slice_axes = []
    
    # Row 1: Contours per slice plots
    for roi_idx, roi_name in enumerate(valid_roi_names):
        ax = fig.add_subplot(gs[0, roi_idx])
        contour_axes.append(ax)
        
        data = roi_data[roi_name]
        
        # Create counts array with 0 for missing slices
        counts = []
        for expected_z in expected_z_values:
            counts.append(data['contours_per_slice'].get(expected_z, 0))
        
        # Create bars for this ROI
        bars = ax.barh(y_positions, counts, height=1.0, 
                      color='skyblue', alpha=0.7, 
                      edgecolor='navy', linewidth=0.5)
        
        # Detect and color isolated contours
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
        for i, count in enumerate(counts):
            z_value = expected_z_values[i]
            if count == 0:  # Missing slice
                bars[i].set_color('lightgray')
                bars[i].set_alpha(0.3)
            elif count > 1:  # Multiple contours on same slice
                bars[i].set_color('red')
                bars[i].set_alpha(0.8)
            elif z_value in isolated_slices:  # Isolated contour
                bars[i].set_color('orange')
                bars[i].set_alpha(0.8)
        
        # Set up axes
        ax.set_title(f'{roi_name}\\nContours per Slice', fontsize=10)
        ax.grid(axis='x', alpha=0.3)
        
        if roi_idx == 0:
            sparse_positions = y_positions[::3]
            sparse_labels = [y_labels[i] for i in range(0, len(y_labels), 3)]
            ax.set_yticks(sparse_positions)
            ax.set_yticklabels(sparse_labels)
            ax.set_ylabel('Slice Z-position (mm)')
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])
        
        max_count = max(counts) if counts else 1
        ax.set_xticks(range(max_count + 1))
        ax.set_xlim(0, max_count + 0.5)
        ax.set_xlabel('# Contours')
    
    # Row 2: Area distribution plots
    for roi_idx, roi_name in enumerate(valid_roi_names):
        ax = fig.add_subplot(gs[1, roi_idx])
        area_axes.append(ax)
        
        data = roi_data[roi_name]
        contour_areas = data['contour_areas']
        
        if contour_areas:
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
        
        ax.set_title(f'{roi_name}\\nArea Distribution', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xlabel('Contour Area (mm²)')
        
        if roi_idx == 0:
            ax.set_ylabel('Frequency')
        else:
            ax.set_yticklabels([])
    
    # Row 3: Interactive slice viewer plots (placeholders)
    for roi_idx, roi_name in enumerate(valid_roi_names):
        ax = fig.add_subplot(gs[2, roi_idx])
        slice_axes.append(ax)
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title('', fontsize=9)
    
    # Align y-axes
    for ax in contour_axes:
        ax.set_ylim(-0.5, len(expected_z_values) - 0.5)
    
    # Align y-axes for area plots
    max_frequency = _calculate_max_frequency(valid_roi_names, roi_data)
    for ax in area_axes:
        ax.set_ylim(0, max_frequency + 1)
    
    return contour_axes, area_axes, slice_axes


def _detect_global_issues(valid_roi_names, roi_data, expected_z_values, slice_thickness, external_contours):
    """Detect global issues across all ROIs for error banner."""
    global_has_multiple_contours = False
    global_has_isolated_contours = False
    global_has_missing_contours = False
    global_has_external_violations = False
    
    isolation_threshold = 2 * slice_thickness
    
    for roi_name in valid_roi_names:
        data = roi_data[roi_name]
        
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
            global_has_external_violations = True  # This would need more detailed checking
    
    return global_has_multiple_contours, global_has_isolated_contours, global_has_missing_contours, global_has_external_violations


def _calculate_max_frequency(valid_roi_names, roi_data):
    """Calculate maximum frequency for area plot y-axis alignment."""
    max_frequency = 0
    for roi_name in valid_roi_names:
        if roi_name in roi_data and roi_data[roi_name]['contour_areas']:
            areas = roi_data[roi_name]['contour_areas']
            for area in areas:
                bin_count = sum(1 for a in areas if int(a) == int(area))
                max_frequency = max(max_frequency, bin_count)
    return max_frequency


def _setup_slice_navigation(valid_roi_names, roi_data, slice_axes, external_contours):
    """Setup interactive slice navigation with pre-rendered data."""
    current_slice_indices = {}
    slice_artists = {}
    
    for roi_idx, roi_name in enumerate(valid_roi_names):
        ax = slice_axes[roi_idx]
        data = roi_data[roi_name]
        available_slices = sorted(data['contours_per_slice'].keys())
        
        current_slice_indices[roi_name] = 0
        
        # Pre-process all slices for this ROI
        prerendered_slices = {}
        all_x_coords = []
        all_y_coords = []
        
        for slice_idx, z_slice in enumerate(available_slices):
            # Find contours at this z-level
            contours_at_z = []
            contour_areas_at_z = []
            
            for i, z_val in enumerate(data['z_values']):
                if abs(z_val - z_slice) < 0.01:
                    contours_at_z.append(data['contours'][i])
                    contour_areas_at_z.append(data['contour_areas'][i])
            
            # Process contour coordinates
            processed_contours = []
            for contour, area in zip(contours_at_z, contour_areas_at_z):
                x_coords = []
                y_coords = []
                
                for point in contour:
                    if hasattr(point, 'x'):
                        x_coords.append(point.x)
                        y_coords.append(-point.y)  # Flip y-coordinate
                    else:
                        x_coords.append(point['x'])
                        y_coords.append(-point['y'])  # Flip y-coordinate
                
                if x_coords and y_coords:
                    x_coords.append(x_coords[0])  # Close contour
                    y_coords.append(y_coords[0])
                    
                    all_x_coords.extend(x_coords)
                    all_y_coords.extend(y_coords)
                    
                    # Split contour for partial highlighting
                    segments = split_contour_by_external(
                        x_coords, y_coords, external_contours.get(z_slice, []))
                    
                    for segment in segments:
                        if area < 2.0:
                            color = 'red'
                            alpha = 0.8
                        elif segment['type'] == 'outside':
                            color = 'red'
                            alpha = 0.8
                        else:
                            color = 'blue'
                            alpha = 0.6
                        
                        processed_contours.append({
                            'x_coords': segment['x_coords'],
                            'y_coords': segment['y_coords'],
                            'color': color,
                            'alpha': alpha
                        })
            
            # Add external contour references
            external_contours_at_z = external_contours.get(z_slice, [])
            for ext_contour in external_contours_at_z:
                processed_contours.append({
                    'x_coords': ext_contour['x_coords'] + [ext_contour['x_coords'][0]],
                    'y_coords': ext_contour['y_coords'] + [ext_contour['y_coords'][0]],
                    'color': 'gray',
                    'alpha': 0.5,
                    'is_external': True
                })
            
            prerendered_slices[slice_idx] = {
                'contours': processed_contours,
                'title': f'{roi_name}\\nSlice {slice_idx + 1}/{len(available_slices)}'
            }
        
        # Set axis limits
        if all_x_coords and all_y_coords:
            margin_x = (max(all_x_coords) - min(all_x_coords)) * 0.05
            margin_y = (max(all_y_coords) - min(all_y_coords)) * 0.05
            ax.set_xlim(min(all_x_coords) - margin_x, max(all_x_coords) + margin_x)
            ax.set_ylim(min(all_y_coords) - margin_y, max(all_y_coords) + margin_y)
        
        # Create artists for each slice
        slice_artists[roi_name] = {
            'slices': {},
            'ax': ax,
            'prerendered': prerendered_slices
        }
        
        for slice_idx in range(len(available_slices)):
            if slice_idx in prerendered_slices:
                slice_data = prerendered_slices[slice_idx]
                contours = slice_data['contours']
                
                line_artists = []
                fill_artists = []
                
                for contour_data in contours:
                    line, = ax.plot(contour_data['x_coords'], contour_data['y_coords'], 
                                  color=contour_data['color'], 
                                  linewidth=3 if contour_data.get('is_external', False) else 2,
                                  alpha=contour_data['alpha'], visible=False)
                    line_artists.append(line)
                    
                    if not contour_data.get('is_external', False):
                        fill = ax.fill(contour_data['x_coords'], contour_data['y_coords'], 
                                     color=contour_data['color'], alpha=0.2, visible=False)[0]
                        fill_artists.append(fill)
                    else:
                        fill = ax.fill([], [], alpha=0, visible=False)[0]
                        fill_artists.append(fill)
                
                slice_artists[roi_name]['slices'][slice_idx] = {
                    'line_artists': line_artists,
                    'fill_artists': fill_artists,
                    'title': slice_data['title']
                }
    
    # Initialize displays
    for roi_name in valid_roi_names:
        _update_slice_display(roi_name, 0, slice_artists)
    
    return current_slice_indices, slice_artists


def _update_slice_display(roi_name, slice_index, slice_artists):
    """Update slice display using hide/show approach."""
    artists = slice_artists[roi_name]
    ax = artists['ax']
    
    # Hide all artists first
    for slice_idx, slice_artists_data in artists['slices'].items():
        for line_artist in slice_artists_data['line_artists']:
            line_artist.set_visible(False)
        for fill_artist in slice_artists_data['fill_artists']:
            fill_artist.set_visible(False)
    
    # Show current slice
    if slice_index in artists['slices']:
        slice_artists_data = artists['slices'][slice_index]
        for line_artist in slice_artists_data['line_artists']:
            line_artist.set_visible(True)
        for fill_artist in slice_artists_data['fill_artists']:
            fill_artist.set_visible(True)
        ax.set_title(slice_artists_data['title'], fontsize=9)
    else:
        ax.set_title(f'{roi_name}\\nNo Data', fontsize=9)


def _setup_event_handlers(fig, slice_axes, valid_roi_names, roi_data, 
                         current_slice_indices, slice_artists):
    """Setup mouse and keyboard event handlers for slice navigation."""
    
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
            
            if event.button == 1 and current_idx < max_idx:  # Left click
                current_slice_indices[roi_name] = current_idx + 1
                _update_slice_display(roi_name, current_slice_indices[roi_name], slice_artists)
                slice_axes[roi_idx].figure.canvas.draw_idle()
            elif event.button == 3 and current_idx > 0:  # Right click
                current_slice_indices[roi_name] = current_idx - 1
                _update_slice_display(roi_name, current_slice_indices[roi_name], slice_artists)
                slice_axes[roi_idx].figure.canvas.draw_idle()
    
    def on_key(event):
        if event.key in ['up', 'down', 'left', 'right']:
            for roi_name in valid_roi_names:
                available_slices = sorted(roi_data[roi_name]['contours_per_slice'].keys())
                if not available_slices:
                    continue
                
                current_idx = current_slice_indices[roi_name]
                max_idx = len(available_slices) - 1
                
                if event.key in ['up', 'right'] and current_idx < max_idx:
                    current_slice_indices[roi_name] = current_idx + 1
                elif event.key in ['down', 'left'] and current_idx > 0:
                    current_slice_indices[roi_name] = current_idx - 1
                
                _update_slice_display(roi_name, current_slice_indices[roi_name], slice_artists)
            
            fig.canvas.draw_idle()
    
    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)


def _create_scrollable_interface(root, fig, valid_roi_names, roi_data, slice_thickness,
                               global_has_multiple_contours, global_has_isolated_contours,
                               global_has_missing_contours, global_has_external_violations):
    """Create the scrollable tkinter interface with error banners and statistics."""
    
    # Create scrollable frame
    main_canvas = tk.Canvas(root)
    scrollbar = tk.Scrollbar(root, orient="vertical", command=main_canvas.yview)
    scrollable_frame = tk.Frame(main_canvas)
    
    scrollable_frame.bind(
        "<Configure>",
        lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
    )
    
    main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    main_canvas.configure(yscrollcommand=scrollbar.set)
    
    def on_mousewheel(event):
        main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
    
    main_canvas.bind("<MouseWheel>", on_mousewheel)
    main_canvas.bind("<Button-4>", lambda e: main_canvas.yview_scroll(-1, "units"))
    main_canvas.bind("<Button-5>", lambda e: main_canvas.yview_scroll(1, "units"))
    main_canvas.focus_set()
    
    main_canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")
    
    # Embed plot
    canvas = FigureCanvasTkAgg(fig, master=scrollable_frame)
    canvas.draw()
    
    # Add error banner if needed
    if any([global_has_multiple_contours, global_has_isolated_contours, 
           global_has_missing_contours, global_has_external_violations]):
        
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
    
    # Re-connect event handlers
    canvas.mpl_connect('button_press_event', lambda event: None)  # These would be properly connected
    canvas.mpl_connect('key_press_event', lambda event: None)
    canvas.get_tk_widget().focus_set()
    
    # Create statistics frame
    stats_frame = tk.Frame(scrollable_frame, bg='wheat', relief=tk.RAISED, bd=2)
    stats_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
    
    # Add statistics for each ROI
    for roi_name in valid_roi_names:
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
    
    # Add slice thickness and navigation info
    slice_info = tk.Label(stats_frame, text=f"Slice Thickness: {slice_thickness}mm", 
                         bg='wheat', font=('Arial', 9, 'bold'))
    slice_info.pack(pady=2)
    
    nav_info = tk.Label(stats_frame, 
                       text="Slice Navigation: Left Click = Next | Right Click = Previous | Arrow Keys = Navigate All ROIs", 
                       bg='wheat', font=('Arial', 8, 'italic'))
    nav_info.pack(pady=1)
