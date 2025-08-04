"""
ROI Selection GUI for multi-ROI contour analysis.

This module provides the initial dialog for selecting target ROIs and external ROI
for analysis. It includes validation and performance optimization.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import time


def get_roi_selection(case, exam):
    """
    Create a GUI for selecting ROIs to analyze.
    
    Args:
        case: RayStation case object
        exam: RayStation exam object
    
    Returns:
        tuple: (selected_rois, external_roi_name) or (None, None) if cancelled
    """
    try:
        from .raystation_interface import get_available_rois, validate_roi_contours
    except ImportError:
        from raystation_interface import get_available_rois, validate_roi_contours
    
    print("Starting ROI selection GUI...")
    start_time = time.time()
    
    # Get all available ROIs
    try:
        target_rois, external_rois = get_available_rois(case, exam)
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
        valid_selections = validate_roi_contours(case, exam, selected)
        
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
