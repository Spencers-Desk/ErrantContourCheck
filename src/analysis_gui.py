"""
Main Analysis GUI for ErrantContourCheck

This module provides the comprehensive analysis visualization interface that displays:
- Contours per slice plots with error highlighting
- Area distribution histograms  
- Interactive slice viewers with navigation
- Error detection banners and statistics
- Scrollable interface with mouse/keyboard controls
"""

import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from archive.utils import (
    split_contour_by_external, 
    calculate_max_frequency,
    detect_global_issues,
    find_isolated_slices,
    process_contour_coordinates,
    determine_contour_color
)

class AnalysisGUI:
    def __init__(self, plan_data):
        """
        Initialize the Analysis GUI.
        
        Args:
            plan_data: PlanData object with processed ROI data
        """
        self.plan_data = plan_data
        self.root = None
        self.fig = None
        self.canvas = None
        
        # Check that plan_data is properly processed
        if not self.plan_data.data_processed:
            raise ValueError("PlanData must have processed ROI data before creating Analysis GUI")
        
        # GUI state variables
        self.current_slice_indices = {}
        self.slice_artists = {}
        self.prerendered_slices = {}
        
        # Get data references for convenience
        self.roi_names = self.plan_data.selected_roi_names
        self.roi_data = self.plan_data.roi_data
        self.external_contours = self.plan_data.external_contours
        self.slice_thickness = self.plan_data.slice_thickness
        self.all_z_values = self.plan_data.all_z_values
    
        self.slice_axes = []          # NEW: store slice axes for later event hookup
        self.valid_roi_names_for_events = []
        self._closed = False
    
    def _on_close(self):
        """Handle window close: set stop flag and destroy root."""
        try:
            self.plan_data.stop_requested = True
            print("Stop requested from Analysis GUI")  # Add logging for troubleshooting
        except Exception as e:
            print(f"Error setting stop flag: {e}")
        self._closed = True
        if self.root:
            try:
                self.root.destroy()
            except Exception as e:
                print(f"Error destroying root: {e}")
    
    def create_gui(self):
        """Create and display the main analysis GUI."""
        print("Creating Analysis GUI...")
        
        # Filter roi_names to only include ROIs that have data
        valid_roi_names = [roi for roi in self.roi_names if roi in self.roi_data]
        num_rois = len(valid_roi_names)
        
        if num_rois == 0:
            print("No valid ROI data found!")
            return False  # Return False instead of None
        
        # Create main window
        self.root = tk.Tk()
        roi_names_str = ", ".join(self.roi_names)
        self.root.title(f"Multi-ROI Contour Analysis - {roi_names_str}")
        self.root.geometry("1600x1000")
        
        # Create matplotlib figure
        self._create_matplotlib_figure(valid_roi_names, num_rois)
        
        # Create scrollable tkinter interface
        self._create_scrollable_interface(valid_roi_names)
        
        # Ensure stop flag set when user closes window
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        
        # Show the GUI
        self.root.mainloop()
        
        # After mainloop returns (window closed) - more robust handling
        print(f"GUI mainloop ended, closed={self._closed}, stop_requested={self.plan_data.stop_requested}")
        self.plan_data.stop_requested = True  # Always set stop_requested after GUI closes
        return self._closed  # Return closed state

    def _create_matplotlib_figure(self, valid_roi_names, num_rois):
        """Create the matplotlib figure with all subplots."""
        # Create figure with 3 rows of subplots (contours per slice, area distribution, slice viewer)
        self.fig = plt.figure(figsize=(16, 10))
        
        # Define custom subplot grid with height ratios [1, 0.67, 1]
        gs = self.fig.add_gridspec(3, num_rois, height_ratios=[1, 0.67, 1])
        
        # Calculate global z-range for consistent y-axis
        min_z = min(self.all_z_values)
        max_z = max(self.all_z_values)
        
        # Create complete range of z-values at slice thickness intervals
        expected_z_values = []
        current_z = min_z
        while current_z <= max_z + self.slice_thickness/2:
            expected_z_values.append(round(current_z, 2))
            current_z += self.slice_thickness
        
        y_positions = list(range(len(expected_z_values)))
        y_labels = [f"{z:.1f}" for z in expected_z_values]
        
        # Track global issues for error banner
        global_issues = detect_global_issues(
            {roi: self.roi_data[roi] for roi in valid_roi_names}, 
            expected_z_values, 
            self.slice_thickness, 
            self.external_contours
        )
        
        # Create all three rows of subplots
        contour_axes, area_axes, slice_axes = self._create_all_subplots(
            gs, valid_roi_names, expected_z_values, y_positions, y_labels
        )
        
        # Setup interactive slice navigation
        self._setup_slice_navigation(valid_roi_names, slice_axes)
        
        # DEFERRED: event handlers must be connected only AFTER FigureCanvasTkAgg is created
        # self._setup_event_handlers(slice_axes, valid_roi_names)
        
        # Store for later connection once TkAgg canvas exists
        self.slice_axes = slice_axes
        self.valid_roi_names_for_events = valid_roi_names
    
        plt.tight_layout()
        
        # Store global issues for error banner
        self.global_issues = global_issues
    

    
    def _create_all_subplots(self, gs, valid_roi_names, expected_z_values, y_positions, y_labels):
        """Create all three rows of subplots."""
        contour_axes = []
        area_axes = []
        slice_axes = []
        
        # Row 1: Contours per slice plots
        for roi_idx, roi_name in enumerate(valid_roi_names):
            ax = self.fig.add_subplot(gs[0, roi_idx])
            contour_axes.append(ax)
            
            self._create_contour_plot(ax, roi_name, expected_z_values, y_positions, y_labels, roi_idx)
        
        # Row 2: Area distribution plots
        for roi_idx, roi_name in enumerate(valid_roi_names):
            ax = self.fig.add_subplot(gs[1, roi_idx])
            area_axes.append(ax)
            
            self._create_area_plot(ax, roi_name, roi_idx)
        
        # Row 3: Interactive slice viewer plots
        for roi_idx, roi_name in enumerate(valid_roi_names):
            ax = self.fig.add_subplot(gs[2, roi_idx])
            slice_axes.append(ax)
            
            self._setup_slice_plot(ax, roi_name)
        
        # Align y-axes
        for ax in contour_axes:
            ax.set_ylim(-0.5, len(expected_z_values) - 0.5)
        
        # Align y-axes for area plots
        max_frequency = calculate_max_frequency({roi: self.roi_data[roi] for roi in valid_roi_names})
        for ax in area_axes:
            ax.set_ylim(0, max_frequency + 1)
        
        return contour_axes, area_axes, slice_axes
    
    def _create_contour_plot(self, ax, roi_name, expected_z_values, y_positions, y_labels, roi_idx):
        """Create contours per slice plot for one ROI."""
        data = self.roi_data[roi_name]
        
        # Create counts array with 0 for missing slices
        counts = []
        for expected_z in expected_z_values:
            counts.append(data['contours_per_slice'].get(expected_z, 0))
        
        # Create bars for this ROI
        bars = ax.barh(y_positions, counts, height=1.0, 
                      color='skyblue', alpha=0.7, 
                      edgecolor='navy', linewidth=0.5)
        
        # Detect isolated contours
        isolated_slices = find_isolated_slices(expected_z_values, counts, self.slice_thickness)
        
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
        ax.set_title(f'{roi_name}\nContours per Slice', fontsize=10)
        ax.grid(axis='x', alpha=0.3)
        
        if roi_idx == 0:
            # Show every 3rd slice to reduce crowding
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
    
    def _create_area_plot(self, ax, roi_name, roi_idx):
        """Create area distribution plot for one ROI."""
        data = self.roi_data[roi_name]
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
        
        ax.set_title(f'{roi_name}\nArea Distribution', fontsize=10)
        ax.grid(axis='y', alpha=0.3)
        ax.set_xlabel('Contour Area (mm²)')
        
        if roi_idx == 0:
            ax.set_ylabel('Frequency')
        else:
            ax.set_yticklabels([])
    
    def _setup_slice_plot(self, ax, roi_name):
        """Setup slice viewer plot for one ROI."""
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Remove x and y axes for cleaner display
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title('', fontsize=9)  # Will be updated by navigation
    

    
    def _setup_slice_navigation(self, valid_roi_names, slice_axes):
        """Setup interactive slice navigation with pre-rendered data."""
        for roi_idx, roi_name in enumerate(valid_roi_names):
            ax = slice_axes[roi_idx]
            data = self.roi_data[roi_name]
            available_slices = sorted(data['contours_per_slice'].keys())
            
            self.current_slice_indices[roi_name] = 0
            
            # Pre-process all slices for this ROI
            self._prerender_slices(roi_name, available_slices, ax)
            
            # Initialize to show first slice
            self._update_slice_display(roi_name, 0)
    
    def _prerender_slices(self, roi_name, available_slices, ax):
        """Pre-render all slice data for fast navigation."""
        data = self.roi_data[roi_name]
        self.prerendered_slices[roi_name] = {}
        all_x_coords = []
        all_y_coords = []
        
        for slice_idx, z_slice in enumerate(available_slices):
            # Find contours at this z-level
            contours_at_z = []
            contour_areas_at_z = []
            
            for i, z_val in enumerate(data['z_values']):
                if abs(z_val - z_slice) < 0.01:  # Account for floating point precision
                    contours_at_z.append(data['contours'][i])
                    contour_areas_at_z.append(data['contour_areas'][i])
            
            # Process contour coordinates
            processed_contours = []
            for contour, area in zip(contours_at_z, contour_areas_at_z):
                x_coords, y_coords = process_contour_coordinates(contour)
                
                if x_coords and y_coords:
                    all_x_coords.extend(x_coords)
                    all_y_coords.extend(y_coords)
                    
                    # Split contour for partial highlighting
                    segments = split_contour_by_external(
                        x_coords, y_coords, self.external_contours.get(z_slice, []))
                    
                    # Check if any part is outside external boundary
                    outside_external = any(seg['type'] == 'outside' for seg in segments)
                    
                    for segment in segments:
                        # Determine color based on area and segment type
                        color, alpha = determine_contour_color(area, segment['type'])
                        
                        processed_contours.append({
                            'x_coords': segment['x_coords'],
                            'y_coords': segment['y_coords'],
                            'area': area,
                            'segment_type': segment['type'],
                            'outside_external': outside_external,
                            'color': color,
                            'alpha': alpha
                        })
            
            # Add external contour outlines for visual reference
            external_contours_at_z = self.external_contours.get(z_slice, [])
            for ext_contour in external_contours_at_z:
                color, alpha = determine_contour_color(0, 'external', is_external=True)
                processed_contours.append({
                    'x_coords': ext_contour['x_coords'] + [ext_contour['x_coords'][0]],
                    'y_coords': ext_contour['y_coords'] + [ext_contour['y_coords'][0]],
                    'area': 0,
                    'segment_type': 'external',
                    'outside_external': False,
                    'color': color,
                    'alpha': alpha,
                    'is_external': True
                })
            
            # Store slice info
            num_contours = len(contours_at_z)
            small_contours = sum(1 for area in contour_areas_at_z if area < 2.0)
            
            # Count external violations
            contours_with_violations = set()
            for contour in processed_contours:
                if contour.get('outside_external', False) and not contour.get('is_external', False):
                    contours_with_violations.add(contour['area'])
            external_violations = len(contours_with_violations)
            
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
            self.prerendered_slices[roi_name][slice_idx] = {
                'contours': processed_contours,
                'title': f'{roi_name}\n{slice_info}',
                'z_value': z_slice
            }
        
        # Set axis limits
        if all_x_coords and all_y_coords:
            margin_x = (max(all_x_coords) - min(all_x_coords)) * 0.05
            margin_y = (max(all_y_coords) - min(all_y_coords)) * 0.05
            ax.set_xlim(min(all_x_coords) - margin_x, max(all_x_coords) + margin_x)
            ax.set_ylim(min(all_y_coords) - margin_y, max(all_y_coords) + margin_y)
        
        # Create artists for each slice
        self._create_slice_artists(roi_name, available_slices, ax)
    
    def _create_slice_artists(self, roi_name, available_slices, ax):
        """Create matplotlib artists for all slices."""
        self.slice_artists[roi_name] = {
            'slices': {},
            'ax': ax
        }
        
        for slice_idx in range(len(available_slices)):
            if slice_idx in self.prerendered_slices[roi_name]:
                slice_data = self.prerendered_slices[roi_name][slice_idx]
                contours = slice_data['contours']
                
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
                        # For external contours, add a dummy invisible fill
                        fill = ax.fill([], [], alpha=0, visible=False)[0]
                        slice_fill_artists.append(fill)
                
                self.slice_artists[roi_name]['slices'][slice_idx] = {
                    'line_artists': slice_line_artists,
                    'fill_artists': slice_fill_artists,
                    'title': slice_data['title']
                }
    
    def _update_slice_display(self, roi_name, slice_index):
        """Update slice display using hide/show approach."""
        artists = self.slice_artists[roi_name]
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
            ax.set_title(f'{roi_name}\nNo Data', fontsize=9)
    
    def _setup_event_handlers(self, slice_axes, valid_roi_names):
        """Setup mouse and keyboard event handlers."""
        def on_click(event):
            # NEW: ensure canvas regains focus so arrow keys work after any click
            if hasattr(self, 'canvas') and self.canvas:
                try:
                    self.canvas.get_tk_widget().focus_set()
                except Exception:
                    pass
            if event.inaxes in slice_axes:
                roi_idx = slice_axes.index(event.inaxes)
                roi_name = valid_roi_names[roi_idx]
                
                data = self.roi_data[roi_name]
                available_slices = sorted(data['contours_per_slice'].keys())
                if not available_slices:
                    return
                
                current_idx = self.current_slice_indices[roi_name]
                max_idx = len(available_slices) - 1
                
                # Left click: next slice
                if event.button == 1 and current_idx < max_idx:
                    self.current_slice_indices[roi_name] = current_idx + 1
                    self._update_slice_display(roi_name, self.current_slice_indices[roi_name])
                    slice_axes[roi_idx].figure.canvas.draw_idle()
                
                # Right click: previous slice
                elif event.button == 3 and current_idx > 0:
                    self.current_slice_indices[roi_name] = current_idx - 1
                    self._update_slice_display(roi_name, self.current_slice_indices[roi_name])
                    slice_axes[roi_idx].figure.canvas.draw_idle()
        
        def on_key(event):
            if event.key in ['up', 'down', 'left', 'right']:
                # Apply to all ROIs simultaneously
                for roi_name in valid_roi_names:
                    available_slices = sorted(self.roi_data[roi_name]['contours_per_slice'].keys())
                    if not available_slices:
                        continue
                    
                    current_idx = self.current_slice_indices[roi_name]
                    max_idx = len(available_slices) - 1
                    
                    if event.key in ['up', 'right'] and current_idx < max_idx:
                        self.current_slice_indices[roi_name] = current_idx + 1
                    elif event.key in ['down', 'left'] and current_idx > 0:
                        self.current_slice_indices[roi_name] = current_idx - 1
                    
                    self._update_slice_display(roi_name, self.current_slice_indices[roi_name])
                
                self.fig.canvas.draw_idle()
        
        self.fig.canvas.mpl_connect('button_press_event', on_click)
        self.fig.canvas.mpl_connect('key_press_event', on_key)
    
    def _create_scrollable_interface(self, valid_roi_names):
        """Create the scrollable tkinter interface."""
        # Create scrollable frame
        main_canvas = tk.Canvas(self.root)
        scrollbar = tk.Scrollbar(self.root, orient="vertical", command=main_canvas.yview)
        scrollable_frame = tk.Frame(main_canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        )
        
        main_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        main_canvas.configure(yscrollcommand=scrollbar.set)
        
        # Add mouse wheel scrolling
        def on_mousewheel(event):
            main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        main_canvas.bind("<MouseWheel>", on_mousewheel)
        main_canvas.bind("<Button-4>", lambda e: main_canvas.yview_scroll(-1, "units"))
        main_canvas.bind("<Button-5>", lambda e: main_canvas.yview_scroll(1, "units"))
        main_canvas.focus_set()
        
        main_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Embed plot
        self.canvas = FigureCanvasTkAgg(self.fig, master=scrollable_frame)
        self.canvas.draw()
        # NEW: bind focus on any click inside canvas widget (ensures key events)
        self.canvas.get_tk_widget().bind("<Button-1>", lambda e: e.widget.focus_set())
        
        # Add error banner if needed
        self._add_error_banner(scrollable_frame)
        
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True, pady=(5, 5))
        
        # REMOVED two no-op event connections that prevented reliable key/click handling
        # (previously added lambdas swallowing events)
        # self.canvas.mpl_connect('button_press_event', lambda event: None)
        # self.canvas.mpl_connect('key_press_event', lambda event: None)
        
        # Ensure initial focus so arrow keys work immediately
        try:
            self.canvas.get_tk_widget().focus_set()
        except Exception:
            pass
        # NEW: now that TkAgg canvas is active, (re)attach matplotlib event handlers
        if self.slice_axes and self.valid_roi_names_for_events:
            self._setup_event_handlers(self.slice_axes, self.valid_roi_names_for_events)

        # Add statistics frame
        self._add_statistics_frame(scrollable_frame, valid_roi_names)
    
    def _add_error_banner(self, parent_frame):
        """Add error detection banner if issues found."""
        issues = self.global_issues
        if any(issues.values()):
            error_text = "⚠️ Possible Errant Contours Detected"
            details = []
            if issues['missing_contours']:
                details.append("Gaps found")
            if issues['isolated_contours']:
                details.append("Isolated contours (orange)")
            if issues['multiple_contours']:
                details.append("Multiple contours per slice (red)")
            if issues['external_violations']:
                details.append("Contours outside External (red)")
            
            if details:
                error_text += f" ({', '.join(details)})"
            
            error_frame = tk.Frame(parent_frame, bg='red', relief=tk.RAISED, bd=2)
            error_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=(5, 5))
            
            error_label = tk.Label(error_frame, text=error_text, bg='red', fg='white', 
                                  font=('Arial', 12, 'bold'))
            error_label.pack(pady=5)
    
    def _add_statistics_frame(self, parent_frame, valid_roi_names):
        """Add statistics frame at bottom."""
        stats_frame = tk.Frame(parent_frame, bg='wheat', relief=tk.RAISED, bd=2)
        stats_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        
        # Add close button at the top of stats frame
        button_frame = tk.Frame(stats_frame, bg='wheat')
        button_frame.pack(pady=5)
        
        close_button = tk.Button(button_frame, text="Close", 
                                command=self._on_close,
                                bg='red', fg='white', 
                                font=('Arial', 12, 'bold'),
                                padx=20, pady=5)
        close_button.pack()
        
        # Add statistics for each ROI
        for roi_name in valid_roi_names:
            if roi_name in self.roi_data:
                data = self.roi_data[roi_name]
                stats_text = f"{roi_name}: {data['number_of_contours']} contours | "
                stats_text += f"{data['gaps']} gaps | "
                stats_text += f"Total Area: {data['total_area']:.0f}mm² | "
                stats_text += f"Avg Area: {data['avg_area']:.0f}mm² | "
                stats_text += f"Min: {data['min_area']:.0f}mm² | "
                stats_text += f"Max: {data['max_area']:.0f}mm²"
                
                stats_label = tk.Label(stats_frame, text=stats_text, bg='wheat', font=('Arial', 9))
                stats_label.pack(pady=2)
        
        # Add slice thickness and navigation info
        slice_info = tk.Label(stats_frame, text=f"Slice Thickness: {self.slice_thickness}mm", 
                             bg='wheat', font=('Arial', 9, 'bold'))
        slice_info.pack(pady=2)
        
        nav_info = tk.Label(stats_frame, 
                           text="Slice Navigation: Left Click = Next | Right Click = Previous | Arrow Keys = Navigate All ROIs", 
                           bg='wheat', font=('Arial', 8, 'italic'))
        nav_info.pack(pady=1)


def create_analysis_gui(plan_data):
    """
    Convenience function to create and show the analysis GUI.
    
    Args:
        plan_data: PlanData object with processed ROI data
    
    Returns:
        bool: True if GUI was closed normally, False if error occurred
    """
    try:
        gui = AnalysisGUI(plan_data)
        result = gui.create_gui()
        # Make sure stop is requested regardless of how GUI was closed
        plan_data.stop_requested = True
        print(f"Analysis GUI closed with result: {result}")
        return True
    except Exception as e:
        print(f"Error creating Analysis GUI: {e}")
        import traceback
        traceback.print_exc()
        plan_data.stop_requested = True  # Set stop flag even on error
        return False