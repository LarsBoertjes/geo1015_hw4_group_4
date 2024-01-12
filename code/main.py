# Make sure this file is the AHN4 tile
input_las_file = "69BZ2_13.LAZ"
# Make sure this file is the dtm created by the C++ code
input_dtm_file = "dtm.asc"

# Our boundaries to crop the AHN4 tile
min_x = 187465.5
min_y = 315228.5
max_x = min_x + 500
max_y = min_y + 500

# Imported functions from python files
from cropping_code import get_points
from vegetation_filter import filter_vegetation
from create_canopy_height import main as create_canopy_height_function
from create_chm import subtract_rasters as create_chm_function

# Output files, possible to change the name of the output files if you want to
output_las_cropped_file = "cropped_file.las"
output_csv_cropped_file = "cropped_file.csv"
output_vegetation_las = "classified_vegetation.las"
canopy_height_tiff = "canopy_height.tiff"
output_chm_file = "chm.tiff"

boundary = [min_x, max_x, min_y, max_y]
print('Cropping the LAZ file to a smaller LAS and csv file...')
pd_pcl = get_points(boundary, input_las_file, output_las_cropped_file, output_csv_cropped_file)
print(f'Cropped data exported in LAS file ({output_las_cropped_file}) and CSV file ({output_csv_cropped_file})')
print('Extracting the vegetation points from AHN4...')
filter_vegetation(output_las_cropped_file, output_csv_cropped_file, output_vegetation_las)
print(f'Vegetation points saved in a LAS file ({output_vegetation_las})')
print('Makes a grid of the vegetation LAS file...')
create_canopy_height_function(output_vegetation_las, input_dtm_file, canopy_height_tiff)
print(f'GRID file is created (canopy_height_tiff)')
print('Canopy height tiff is substracted form the DTM...')
create_chm_function(canopy_height_tiff, input_dtm_file, output_chm_file)
print(f'CHM is created ({output_chm_file})')