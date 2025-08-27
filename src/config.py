# General debug flag
debug = False

# Name of the External ROI (used for boundary checking)
ExternalName = "External"

# Area threshold for flagging small contours (mm²)
contour_area_threshold = 1  # Used for coloring small contours red

# Gap threshold (number of slice thicknesses to consider a gap)
gap_threshold = 1

# Isolation threshold (multiple of slice thickness to consider a contour isolated)
isolation_threshold_factor = 1  # Used for orange bars

# Maximum number of ROIs user can select in GUI
max_roi_selection = 4

# GUI window sizes
roi_selection_gui_size = "500x400"
analysis_gui_size = "1600x1000"

# Color settings for plots
color_normal = "blue"
color_small = "red"
color_isolated = "orange"
color_external = "gray"
color_missing = "lightgray"
color_multiple = "red"

# Minimum area for a contour to be considered valid (mm²)
min_valid_contour_area = 2.0

# Margin (mm) for plot axis limits
plot_margin_mm = 5

# Histogram bin size for area plots
area_histogram_bin_size = 1

# Font sizes for GUI labels/titles
font_title = ('Arial', 16, 'bold')
font_label = ('Arial', 10)
font_stats = ('Arial', 9)
font_error_banner = ('Arial', 12, 'bold')
font_nav_info = ('Arial', 8, 'italic')

# Colors for GUI backgrounds and error banners
gui_bg_color = 'wheat'
error_banner_color = 'red'
error_banner_fg = 'white'

# Whether to flip y-coordinates in plots
flip_y_coords = True