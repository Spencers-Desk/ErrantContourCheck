"""
ROI Selection GUI for ErrantContourCheck

This module provides an interactive GUI for selecting target ROIs and external ROI
for errant contour analysis. It integrates with the PlanData object to get available
ROIs and validate selections.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import time


class ROISelectionGUI:
    def __init__(self, plan_data):
        """
        Initialize the ROI Selection GUI.
        
        Args:
            plan_data: PlanData object that has been initialized with build()
        """
        self.plan_data = plan_data
        self.selected_rois = []
        self.selected_external = None
        self.root = None
        
        # Check that plan_data is properly initialized
        if not self.plan_data.initialized or not self.plan_data.rois_discovered:
            raise ValueError("PlanData must be initialized and have ROIs discovered before creating GUI")
    
    def show_selection_dialog(self):
        """
        Show the ROI selection dialog and return the user's selections.
        
        Returns:
            tuple: (selected_roi_names, external_roi_name) or (None, None) if cancelled
        """
        print("Starting ROI selection GUI...")
        start_time = time.time()
        
        # Get available ROIs from plan_data
        target_rois = self.plan_data.available_target_rois
        external_rois = self.plan_data.available_external_rois
        
        if not target_rois:
            messagebox.showwarning("No Target ROIs", "No Target-type ROIs with geometry found!")
            return None, None
        
        print(f"Starting GUI creation after {time.time() - start_time:.2f}s total...")
        gui_start = time.time()
        
        # Create selection window
        self.root = tk.Tk()
        self.root.title("Select ROIs for Contour Analysis")
        self.root.geometry("500x400")
        self.root.resizable(False, False)
        
        # Variables to store selections
        self.selected_rois = []
        selected_external_var = tk.StringVar()
        
        # Main frame
        main_frame = tk.Frame(self.root, padx=20, pady=20)
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
        selected_external_var.set(default_external)
        
        ext_combo = ttk.Combobox(ext_row_frame, textvariable=selected_external_var, 
                                values=[""] + external_rois,
                                state="readonly", width=40)
        ext_combo.pack(side=tk.LEFT, padx=(10, 0))
        
        # Button frame
        button_frame = tk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))
        
        def validate_and_proceed():
            """Validate selections and close dialog."""
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
            
            # Get external ROI selection (may be empty string)
            external_selection = selected_external_var.get().strip()
            print(f"Selected External ROI: '{external_selection}'")
            
            # Use plan_data to validate ROI selections
            if self.plan_data.set_selected_rois(selected, external_selection if external_selection else None):
                # Store the validated selections
                self.selected_rois = self.plan_data.selected_roi_names
                self.selected_external = self.plan_data.external_roi_name
                
                print(f"Final validated external ROI: {self.selected_external}")
                
                # Close window
                self.root.quit()
                self.root.destroy()
            else:
                messagebox.showerror("Error", "None of the selected ROIs have valid contours!")
        
        def cancel_selection():
            """Cancel selection and close dialog."""
            self.selected_rois = []
            self.selected_external = None
            self.root.quit()
            self.root.destroy()
        
        # Buttons
        tk.Button(button_frame, text="Cancel", command=cancel_selection,
                  width=15, bg='red').pack(side=tk.RIGHT, padx=(10, 0))
        
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
        self.root.mainloop()
        
        # Return selections
        if self.selected_rois:
            return self.selected_rois, self.selected_external
        else:
            return None, None


def get_roi_selection(plan_data):
    """
    Convenience function to show ROI selection dialog.
    
    Args:
        plan_data: Initialized PlanData object
        
    Returns:
        tuple: (selected_roi_names, external_roi_name) or (None, None) if cancelled
    """
    try:
        gui = ROISelectionGUI(plan_data)
        return gui.show_selection_dialog()
    except Exception as e:
        print(f"Error creating ROI selection GUI: {e}")
        messagebox.showerror("Error", f"Could not create ROI selection dialog: {e}")
        return None, None


# Example usage:
if __name__ == "__main__":
    # This would be used for testing the GUI independently
    print("ROI Selection GUI module")
    print("This module should be imported and used with a PlanData object")