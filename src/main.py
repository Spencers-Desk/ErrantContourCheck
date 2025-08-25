"""
Main entry point for ErrantContourCheck analysis tool.

This script provides the complete workflow:
1. Initialize RayStation connection and discover ROIs
2. Present ROI selection dialog
3. Process selected ROI data
4. Display comprehensive analysis GUI with visualization
"""

import sys

def main():
    """
    Main entry point for the ErrantContourCheck tool.
    """
    try:
        print("="*70)
        print("               ErrantContourCheck Analysis Tool")
        print("="*70)
        
        # Import modules
        print("\nLoading modules...")
        from init import plan_data  # (on_gui_closed already registered inside init)
        from roi_selection_GUI import get_roi_selection
        from analysis_gui import create_analysis_gui
        
        # Step 1: Initialize connection and build plan data
        print("\n1. Connecting to RayStation and discovering ROIs...")
        if not plan_data.build():
            print("[ERROR] Failed to initialize plan data!")
            print("   Make sure RayStation is running and a plan is loaded.")
            return False
        
        print(f"[OK] Found {len(plan_data.available_target_rois)} ROIs in current plan")
        
        # Step 2: ROI selection
        print("\n2. Opening ROI selection dialog...")
        selected_rois, external_roi = get_roi_selection(plan_data)
        
        if selected_rois is None:
            print("[ERROR] No ROIs selected. Analysis cancelled.")
            return False
        
        print(f"[OK] Selected {len(selected_rois)} ROI(s): {', '.join(selected_rois)}")
        if external_roi:
            print(f"   External ROI: {external_roi}")
        
        # Step 3: Process ROI data
        print("\n3. Processing ROI contour data...")
        if not plan_data.process_roi_data():
            # Abort early if GUI stop was requested mid-processing
            if plan_data.stop_requested:
                print("[INFO] Processing aborted due to GUI close request.")
                return False
            print("[ERROR] Failed to process ROI data!")
            return False
        
        print("[OK] ROI data processed successfully")
        
        create_analysis_gui(plan_data)  # This should block until GUI is closed
        
        print("\n[OK] Analysis complete.")
        return True
    
    except ImportError as e:
        print(f"\n[ERROR] Import error: {e}")
        print("   Make sure all required modules are available")
        print("   This tool requires RayStation scripting environment")
        return False
    

if __name__ == "__main__":

    main()
