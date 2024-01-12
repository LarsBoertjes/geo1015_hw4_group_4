from cropping_code import get_points
from vegetation_filter import filter_vegetation
from create_canopy_height import main as create_raster_canopy

input_las_file = "../data/69BZ2_13.LAZ"
input_dtm_file = "dtm.asc"
output_las_cropped_file = "cropped_file.las"
output_csv_cropped_file = "cropped_file.csv"
output_vegetation_las = "classified_vegetation.las"
dtm_canopy_height = "canopy_height.tiff"

# Our boundaries
min_x = 187465.5
max_x = min_x + 500
min_y = 315228.5
max_y = min_y + 500

boundary = [min_x, max_x, min_y, max_y]
pd_pcl = get_points(boundary, input_las_file, output_las_cropped_file, output_csv_cropped_file)

filter_vegetation(output_las_cropped_file, output_csv_cropped_file, output_vegetation_las)
create_raster_canopy(output_vegetation_las, input_dtm_file, dtm_canopy_height)