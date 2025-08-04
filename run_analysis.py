#!/usr/bin/env python
"""
ErrantCheck Launcher

Simple launcher script for the Multi-ROI Contour Analysis Tool.
This script can be run directly in RayStation to start the analysis.

Usage in RayStation:
    1. Load a patient case with contoured structures
    2. Run this script
    3. Follow the GUI prompts to select ROIs
    4. Review the analysis results
"""

# Import and run the main analysis
try:
    # Try package import first
    try:
        from errantCheck.main import main
    except ImportError:
        # Fallback to direct import
        from main import main
    
    if __name__ == "__main__":
        main()
        
except ImportError as e:
    print("Error importing errantCheck module:")
    print(f"  {e}")
    print("")
    print("Make sure the errantCheck directory is in your Python path.")
    print("You may need to add it using:")
    print("  import sys")
    print("  sys.path.append(r'\\\\nasr200\\rayStationScripts\\dev\\SRL')")
    print("")
    print("Or run directly from the errantCheck directory:")
    print("  exec(open('run_analysis.py').read())")
    
except Exception as e:
    print("Error running errant contour analysis:")
    print(f"  {e}")
    import traceback
    traceback.print_exc()
