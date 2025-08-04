#!/usr/bin/env python
"""
Test script to verify module imports work correctly.
Run this to test the errantCheck package structure.
"""

def test_imports():
    """Test that all modules can be imported successfully."""
    print("Testing errantCheck module imports...")
    
    try:
        print("  Testing geometry_tools...")
        import geometry_tools
        print("    ✓ geometry_tools imported successfully")
        
        print("  Testing raystation_interface...")
        import raystation_interface
        print("    ✓ raystation_interface imported successfully")
        
        print("  Testing roi_selector_gui...")
        import roi_selector_gui
        print("    ✓ roi_selector_gui imported successfully")
        
        print("  Testing main_analysis_gui...")
        import main_analysis_gui
        print("    ✓ main_analysis_gui imported successfully")
        
        print("  Testing main...")
        import main
        print("    ✓ main imported successfully")
        
        print("\n✓ All modules imported successfully!")
        print("The errantCheck package is properly structured.")
        
        return True
        
    except ImportError as e:
        print(f"\n✗ Import error: {e}")
        return False
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        return False

def test_functions():
    """Test that key functions are accessible."""
    print("\nTesting key functions...")
    
    try:
        from geometry_tools import calculate_contour_area, point_in_polygon
        print("  ✓ Geometry functions accessible")
        
        # Note: RayStation functions will fail outside RayStation, but import should work
        from raystation_interface import get_available_rois
        print("  ✓ RayStation interface functions accessible")
        
        from roi_selector_gui import get_roi_selection
        print("  ✓ ROI selector GUI function accessible")
        
        from main_analysis_gui import create_analysis_gui
        print("  ✓ Main analysis GUI function accessible")
        
        from main import main
        print("  ✓ Main entry point accessible")
        
        print("\n✓ All key functions are accessible!")
        return True
        
    except Exception as e:
        print(f"\n✗ Function access error: {e}")
        return False

if __name__ == "__main__":
    print("="*50)
    print("ErrantCheck Package Structure Test")
    print("="*50)
    
    success1 = test_imports()
    success2 = test_functions()
    
    print("\n" + "="*50)
    if success1 and success2:
        print("✓ ALL TESTS PASSED")
        print("The errantCheck package is ready for use!")
    else:
        print("✗ SOME TESTS FAILED")
        print("Please check the error messages above.")
    print("="*50)
