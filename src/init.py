from connect import get_current
from collections import Counter
import time
import os  # added
import threading, sys  # added
from archive.utils import calculate_contour_area, calculate_global_parameters, calculate_roi_statistics

# Global stop event (GUI close sets this)
STOP_EVENT = threading.Event()

class PlanData:
    def __init__(self):
        # RayStation objects
        self.patient = None
        self.case = None
        self.exam = None
        
        # ROI lists
        self.available_target_rois = []
        self.available_external_rois = []
        
        # Selected ROIs for analysis
        self.selected_roi_names = []
        self.external_roi_name = "External"
        
        # Processed ROI data
        self.roi_data = {}
        self.external_contours = {}
        
        # Analysis parameters
        self.slice_thickness = None
        self.all_z_values = []
        self.min_z = None
        self.max_z = None
        
        # Status flags
        self.initialized = False
        self.rois_discovered = False
        self.data_processed = False
        self.stop_requested = False  # <- added earlier
        self.has_any_gaps = False          # <- new aggregate flags
        self.has_any_outside = False       # <- new aggregate flags
        # Tolerance / criteria for external boundary evaluation
        self.edge_tolerance = 0.25  # mm: distance to polygon edge considered "on" (inside)
        self.inside_ratio_threshold = 0.95  # proportion of contour points that must be inside
        # Internal snapshot of issue fields before stats
        self._pre_stats_issue_fields = {}
        
    def reset(self):
        """Reset discovered and processed data while keeping current RS object handles."""
        self.available_target_rois = []
        self.available_external_rois = []
        self.selected_roi_names = []
        self.external_roi_name = "External"
        self.roi_data = {}
        self.external_contours = {}
        self.slice_thickness = None
        self.all_z_values = []
        self.min_z = None
        self.max_z = None
        self.rois_discovered = False
        self.data_processed = False
        self.stop_requested = False  # reset stop flag
        self.has_any_gaps = False
        self.has_any_outside = False
        self.edge_tolerance = 0.25
        self.inside_ratio_threshold = 0.95
        print("PlanData state has been reset.")
        
    def build(self):
        """Initialize the RayStation connection and discover available ROIs."""
        print("Initializing PlanData...")
        
        try:
            # Get RayStation objects
            self.patient = get_current("Patient")
            self.case = get_current("Case")
            self.exam = get_current("Examination")
            
            print(f"Connected to Patient: {self.patient.Name}")
            print(f"Case: {self.case.CaseName}")
            print(f"Exam: {self.exam.Name}")
            
            self.initialized = True
            
            # Discover available ROIs
            self._discover_rois()
            
            return True
            
        except Exception as e:
            print(f"Error initializing PlanData: {type(e).__name__}: {e}")
            return False
    
    def _discover_rois(self):
        """Discover and categorize available ROIs."""
        print("Discovering available ROIs...")
        start_time = time.time()
        
        try:
            # Get all ROIs from the case
            all_rois = self.case.PatientModel.RegionsOfInterest
            print(f"Found {len(all_rois)} total ROIs")
            
            # Filter for Target type ROIs and those with contours
            target_rois = []
            external_rois = []
            
            for i, roi in enumerate(all_rois):
                self.abort_if_requested()
                if i % 10 == 0:  # Progress indicator
                    print(f"  Processing ROI {i+1}/{len(all_rois)}: {roi.Name}")
                
                roi_name = roi.Name
                roi_type = roi.Type
                
                # Check if ROI has contours
                try:
                    roi_geometry = self.case.PatientModel.StructureSets[self.exam.Name].RoiGeometries[roi_name]
                    has_contours = (hasattr(roi_geometry, 'PrimaryShape') and 
                                   roi_geometry.PrimaryShape is not None)
                except Exception:
                    has_contours = False
                
                if has_contours:
                    if roi_type == "External":
                        external_rois.append(roi_name)
                    elif roi_type in ["Ptv", "Ctv", "Gtv", "Target"]:  # Common target types
                        target_rois.append(roi_name)
            
            # Sort the lists
            self.available_target_rois = sorted(target_rois)
            self.available_external_rois = sorted(external_rois)
            
            print(f"ROI discovery completed in {time.time() - start_time:.2f}s")
            print(f"Found {len(self.available_target_rois)} target ROIs, {len(self.available_external_rois)} external ROIs")
            
            self.rois_discovered = True
            return True
            
        except Exception as e:
            print(f"Error discovering ROIs: {type(e).__name__}: {e}")
            return False
    
    def set_selected_rois(self, roi_names, external_roi_name="External"):
        """
        Set the ROIs to be analyzed.
        
        Args:
            roi_names: List of target ROI names to analyze
            external_roi_name: Name of external ROI for boundary checking
        """
        # Validate that ROIs exist and have contours
        valid_rois = []
        
        for roi_name in roi_names:
            if roi_name in self.available_target_rois:
                try:
                    contours, used_temp = self._get_roi_contours(roi_name)
                    if contours:
                        print(f"  {roi_name}: {len(contours)} contours{' (via temp copy)' if used_temp else ''}")
                        valid_rois.append(roi_name)
                    else:
                        print(f"  {roi_name}: No contours found (even after conversion).")
                except Exception as e:
                    print(f"  {roi_name}: Error checking contours - {e}")
            else:
                print(f"  {roi_name}: Not in available target ROIs")
        
        if not valid_rois:
            print("No valid ROIs selected!")
            return False
        
        self.selected_roi_names = valid_rois
        self.external_roi_name = external_roi_name
        
        print(f"Selected ROIs: {self.selected_roi_names}")
        print(f"External ROI: {self.external_roi_name}")
        
        return True
    
    def process_roi_data(self):
        """Process contour data for selected ROIs."""
        if not self.initialized:
            print("Cannot process: PlanData not initialized (call build()).")
            return False
        if not self.rois_discovered:
            print("Cannot process: ROIs not yet discovered.")
            return False
        if not self.selected_roi_names:
            print("No ROIs selected for processing!")
            return False
        
        print("Processing ROI contour data...")
        self.roi_data = {}
        
        # Process each selected ROI
        for roi_name in self.selected_roi_names:
            self.abort_if_requested()
            print(f"\nProcessing ROI: {roi_name}")
            try:
                contours, used_temp = self._get_roi_contours(roi_name)
                if not contours:
                    print(f"  Error: No contours available for {roi_name}")
                    continue
                if used_temp:
                    print("  (Contours obtained via temporary conversion)")
                
                number_of_contours = len(contours)
                print(f"  {number_of_contours} contours found")
                
                z_values = [None] * number_of_contours
                contour_areas = [None] * number_of_contours
                
                for index, contour in enumerate(contours):
                    self.abort_if_requested()
                    # Robust z-value extraction (handles dict or object points)
                    try:
                        first_point = contour[0]
                        if isinstance(first_point, dict):
                            z_raw = first_point.get('z')
                        else:
                            z_raw = getattr(first_point, 'z', None)
                        z_values[index] = round(float(z_raw), 2) if z_raw is not None else None
                    except Exception as z_err:
                        print(f"    Warning: Failed to read z for contour {index}: {type(z_err).__name__}: {z_err}")
                        z_values[index] = None
                    contour_areas[index] = calculate_contour_area(contour)
                
                # Store data for this ROI
                filtered_z = [z for z in z_values if z is not None]
                self.roi_data[roi_name] = {
                    'contours': contours,
                    'number_of_contours': number_of_contours,
                    'z_values_raw': z_values,              # keep raw (may include None)
                    'z_values': filtered_z,                # cleaned list used for analysis
                    'contour_areas': contour_areas,
                    'contours_per_slice': Counter(filtered_z),
                    'outside_external': 0,                 # default until evaluated
                    'external_check_performed': False,
                    'gaps': []  # placeholder; will be replaced by statistics function
                }
            except Exception as e:
                print(f"  Error processing {roi_name}: {type(e).__name__}: {e}")
                continue
        
        # Process external ROI for boundary checking
        self.abort_if_requested()
        self._process_external_roi()
        self.abort_if_requested()
        self._evaluate_external_bounds()
        self.abort_if_requested()
        
        # Calculate global analysis parameters
        self._calculate_global_parameters()
        self.abort_if_requested()
        # Calculate statistics for each ROI
        self._calculate_roi_statistics()
        
        self.data_processed = True
        print("ROI data processing complete!")
        return True
    
    def _process_external_roi(self):
        """Process external ROI contours for boundary checking."""
        self.external_contours = {}
        
        try:
            print(f"\nProcessing External ROI: {self.external_roi_name}")
            external_geom = self.case.PatientModel.StructureSets[self.exam.Name].RoiGeometries[self.external_roi_name].PrimaryShape
            
            external_contour_list = external_geom.Contours
            print(f"  {len(external_contour_list)} external contours found")
            
            # Get all z-values from target ROIs
            target_z_values = set()
            for roi_name, data in self.roi_data.items():
                target_z_values.update(data['z_values'])
            
            # Process only external contours at z-levels we care about
            for contour in external_contour_list:
                self.abort_if_requested()
                try:
                    p0 = contour[0]
                    z_val = round((p0['z'] if isinstance(p0, dict) else p0.z), 2)
                except Exception:
                    continue
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
                        if z_val not in self.external_contours:
                            self.external_contours[z_val] = []
                        self.external_contours[z_val].append({
                            'x_coords': x_coords,
                            'y_coords': y_coords
                        })
            
            print(f"  External contours processed for {len(self.external_contours)} relevant z-slices")
            
        except Exception as e:
            print(f"  Warning: Could not process External ROI '{self.external_roi_name}': {type(e).__name__}: {e}")
            print("  Boundary checking will be disabled.")
    
    def _evaluate_external_bounds(self):
        """
        Tolerance & orientation-robust test:
        - Accept contour as inside if at least inside_ratio_threshold of points lie inside (or on edge)
          of ANY external polygon.
        - Tries both y-flipped and non-flipped orientations to avoid sign convention mismatch.
        """
        if not self.external_contours:
            for roi_name, data in self.roi_data.items():
                data['outside_external'] = 0
                data['external_check_performed'] = False
            return

        tol = self.edge_tolerance
        ratio_req = self.inside_ratio_threshold

        def point_on_segment(px, py, qx, qy, x, y):
            # Distance from point to segment (px,py)-(qx,qy) within tolerance?
            vx, vy = qx - px, qy - py
            wx, wy = x - px, y - py
            seg_len2 = vx*vx + vy*vy
            if seg_len2 == 0:
                # Degenerate segment
                dx, dy = x - px, y - py
                return (dx*dx + dy*dy) ** 0.5 <= tol
            t = max(0.0, min(1.0, (wx*vx + wy*vy) / seg_len2))
            cx, cy = px + t*vx, py + t*vy
            dx, dy = x - cx, y - cy
            return (dx*dx + dy*dy) ** 0.5 <= tol

        def point_in_poly_or_edge(x, y, px, py):
            # Ray casting + edge proximity
            inside = False
            n = len(px)
            for i in range(n):
                j = (i - 1) % n
                xi, yi = px[i], py[i]
                xj, yj = px[j], py[j]
                # Edge proximity
                if point_on_segment(xj, yj, xi, yi, x, y):
                    return True
                intersect = ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / ( (yj - yi) if (yj - yi) != 0 else 1e-12 ) + xi)
                if intersect:
                    inside = not inside
            return inside

        def inside_ratio(points, polygons):
            # Return best inside ratio across polygons
            best = 0.0
            for poly in polygons:
                px = poly['x_coords']
                py = poly['y_coords']
                # Quick bbox reject (expanded by tolerance)
                minx, maxx = min(px)-tol, max(px)+tol
                miny, maxy = min(py)-tol, max(py)+tol
                pts_in_bbox = [p for p in points if (minx <= p[0] <= maxx and miny <= p[1] <= maxy)]
                if not pts_in_bbox:
                    continue
                inside_count = 0
                for x, y in pts_in_bbox:
                    if point_in_poly_or_edge(x, y, px, py):
                        inside_count += 1
                ratio = inside_count / len(points)
                if ratio > best:
                    best = ratio
                if best >= ratio_req:
                    break
            return best

        for roi_name, data in self.roi_data.items():
            self.abort_if_requested()
            outside_count = 0
            contours = data['contours']
            z_vals = data['z_values_raw']
            for idx, contour in enumerate(contours):
                self.abort_if_requested()
                z = z_vals[idx]
                if z is None or z not in self.external_contours:
                    continue
                # Gather both orientation sets of points
                pts_flip = []
                pts_raw = []
                try:
                    for pt in contour:
                        if hasattr(pt, 'x'):
                            pts_raw.append((pt.x, pt.y))
                            pts_flip.append((pt.x, -pt.y))
                        else:
                            pts_raw.append((pt['x'], pt['y']))
                            pts_flip.append((pt['x'], -pt['y']))
                except Exception:
                    continue
                if not pts_raw:
                    continue
                ext_polys = self.external_contours[z]

                # Test with stored external polygons (already y-flipped)
                ratio_flip = inside_ratio(pts_flip, ext_polys)

                # Build unflipped versions of external for orientation fallback
                ext_polys_unflipped = []
                for poly in ext_polys:
                    ext_polys_unflipped.append({
                        'x_coords': poly['x_coords'],
                        'y_coords': [-yy for yy in poly['y_coords']]
                    })
                ratio_raw = inside_ratio(pts_raw, ext_polys_unflipped)

                best_ratio = max(ratio_flip, ratio_raw)
                if best_ratio < ratio_req:
                    outside_count += 1
            data['outside_external'] = outside_count
            data['external_check_performed'] = True

    def _calculate_global_parameters(self):
        """Calculate global slice thickness and z-range parameters."""
        self.abort_if_requested()
        global_params = calculate_global_parameters(self.roi_data)
        self.abort_if_requested()
        
        self.all_z_values = global_params['all_z_values']
        self.slice_thickness = global_params['slice_thickness']
        self.min_z = global_params['min_z']
        self.max_z = global_params['max_z']
        
        if not self.all_z_values:
            print("No contours found in any ROI!")
            return
        if self.slice_thickness is None or self.slice_thickness <= 0:
            print("Warning: Unreliable slice thickness computed.")
        
        print(f"\nGlobal Analysis:")
        print(f"Slice thickness: {self.slice_thickness} mm")
        print(f"Z-range: {self.min_z:.1f} to {self.max_z:.1f} mm")
    
    def _calculate_roi_statistics(self):
        """Calculate statistics for each ROI and add to roi_data."""
        self.abort_if_requested()
        
        # Preserve issue-related fields that calculate_roi_statistics might drop
        self._pre_stats_issue_fields = {
            roi: {
                'outside_external': data.get('outside_external', 0),
                'external_check_performed': data.get('external_check_performed', False)
            } for roi, data in self.roi_data.items()
        }
        self.roi_data = calculate_roi_statistics(self.roi_data, self.slice_thickness)
        # Re-attach preserved fields & normalize
        self._preserve_issue_fields()
        # Post-process to normalize gap/outside flags
        self.has_any_gaps = False
        self.has_any_outside = False
        for roi_name, data in self.roi_data.items():
            gaps_val = data.get('gaps', [])
            if isinstance(gaps_val, (list, tuple, set)):
                gap_count = len(gaps_val)
            else:
                try:
                    gap_count = int(gaps_val)
                    if gap_count == 0:
                        gaps_val = []
                except Exception:
                    gap_count = 0
            data['gap_count'] = gap_count
            data['has_gaps'] = gap_count > 0
            if data['has_gaps']:
                self.has_any_gaps = True
            if data.get('outside_external', 0) > 0:
                self.has_any_outside = True
        # Print statistics for user feedback
        for roi_name, data in self.roi_data.items():
            print(f"\n{roi_name} Analysis:")
            print(f"  Total contours: {data['number_of_contours']}")
            print(f"  Unique slices: {len(set(data['z_values']))}")
            print(f"  Gaps found: {data['gap_count']}")
            print(f"  Total area: {data['total_area']:.2f} mm²")
            print(f"  Average area: {data['avg_area']:.2f} mm²")
            print(f"  Min area: {data['min_area']:.2f} mm²")
            print(f"  Max area: {data['max_area']:.2f} mm²")
        print("\nIssue Summary:")
        print(f"  Any gaps: {self.has_any_gaps}")
        print(f"  Any contours outside external: {self.has_any_outside}")
        if self.has_any_outside:
            offenders = [r for r,d in self.roi_data.items() if d.get('outside_external',0) > 0]
            print(f"  ROIs with outside contours: {offenders}")

    def _preserve_issue_fields(self):
        """Re-attach outside/external flags after stats and ensure defaults."""
        for roi_name, data in self.roi_data.items():
            preserved = self._pre_stats_issue_fields.get(roi_name, {})
            data.setdefault('outside_external', preserved.get('outside_external', 0))
            data.setdefault('external_check_performed', preserved.get('external_check_performed', False))
            # Defensive normalization
            if not isinstance(data['outside_external'], int):
                try:
                    data['outside_external'] = int(data['outside_external'])
                except Exception:
                    data['outside_external'] = 0
            if not isinstance(data['external_check_performed'], bool):
                data['external_check_performed'] = bool(data['external_check_performed'])

    def finalize_issue_flags(self):
        """Optional: call after processing if GUI wants to be extra sure fields exist."""
        self._preserve_issue_fields()

    def get_banner_flags(self):
        """
        Explicit banner booleans for GUI:
          show_gap_banner: True if any gaps
          show_outside_banner: True if any outside contours
        """
        return {
            'processed': self.data_processed,
            'show_gap_banner': self.has_any_gaps,
            'show_outside_banner': self.has_any_outside
        }

    def get_summary(self):
        """Get a summary of the current state."""
        summary = {
            'initialized': self.initialized,
            'rois_discovered': self.rois_discovered,
            'data_processed': self.data_processed,
            'available_target_rois': len(self.available_target_rois),
            'available_external_rois': len(self.available_external_rois),
            'selected_rois': len(self.selected_roi_names),
            'processed_rois': len(self.roi_data)
        }
        
        if self.data_processed:
            summary.update({
                'slice_thickness': self.slice_thickness,
                'z_range': f"{self.min_z:.1f} to {self.max_z:.1f} mm",
                'total_slices': len(self.all_z_values),
                'has_any_gaps': self.has_any_gaps,
                'has_any_outside': self.has_any_outside
            })
        
        return summary

    def get_issue_status(self):
        """
        Returns a concise dict for GUI logic to decide banner state.
        Red banner should only show if has_any_gaps or has_any_outside is True.
        """
        return {
            'processed': self.data_processed,
            'has_any_gaps': self.has_any_gaps,
            'has_any_outside': self.has_any_outside,
            'selected_rois': list(self.selected_roi_names)
        }

    def _get_roi_contours(self, roi_name):
        """
        Attempt to retrieve contours for an ROI. If the ROI's primary shape
        doesn't have contours available, create a temporary copy, convert
        to 'Contours' representation, extract contours, then delete the temp ROI.
        
        Returns: (contours_list or None, used_temp_copy: bool)
        """
        try:
            geom = self.case.PatientModel.StructureSets[self.exam.Name].RoiGeometries[roi_name]
            
            # Check if primary shape exists and has contours
            if (hasattr(geom, 'PrimaryShape') and geom.PrimaryShape and
                hasattr(geom.PrimaryShape, 'Contours') and geom.PrimaryShape.Contours and
                len(geom.PrimaryShape.Contours) > 0):
                return geom.PrimaryShape.Contours, False
        except Exception:
            pass  # Fall through to temp conversion
        
        # Primary shape doesn't have contours - need temporary conversion
        temp_name = f"__TMP__{roi_name}_{int(time.time()*1000)}"
        model = self.case.PatientModel
        struct_set = model.StructureSets[self.exam.Name]
        contours = None
        
        try:
            # Create temporary ROI
            model.CreateRoi(Name=temp_name, Color="255,0,0", Type="Undefined")
            self.abort_if_requested()
            
            # Copy geometry from source ROI
            model.CopyRoiGeometries(
                SourceRoiNames=[roi_name],
                TargetRoiNames=[temp_name],
                ExaminationNames=[self.exam.Name]
            )
            self.abort_if_requested()
            
            # Get temporary ROI geometry and convert to contours representation
            temp_geom = struct_set.RoiGeometries[temp_name]
            try:
                temp_geom.SetRepresentation(Representation='Contours')
                self.abort_if_requested()
            except Exception as rep_err:
                raise RuntimeError(f"SetRepresentation to Contours failed: {rep_err}")
            
            # Extract contours from the converted temporary ROI
            if (hasattr(temp_geom, 'PrimaryShape') and temp_geom.PrimaryShape and
                hasattr(temp_geom.PrimaryShape, 'Contours') and temp_geom.PrimaryShape.Contours):
                contours = temp_geom.PrimaryShape.Contours
            else:
                raise RuntimeError("No contours available even after conversion to Contours representation")
                
        except Exception as e:
            print(f"  Warning: Temporary contour extraction failed for {roi_name}: {type(e).__name__}: {e}")
        finally:
            # Clean up temporary ROI
            try:
                model.RegionsOfInterest[temp_name].DeleteRoi()
            except Exception:
                pass
        
        return (contours if contours else None), True

# Create global instance
plan_data = PlanData()