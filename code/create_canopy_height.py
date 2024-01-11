import sys
import rasterio
import numpy as np
import laspy
from rasterio.transform import from_origin

def read_geospatial_file(file_path):
    file_extension = file_path.split('.')[-1].lower()

    if file_extension == 'asc':
        with rasterio.open(file_path) as src:
            data = src.read(1)
    elif file_extension == 'las':
        las_file = laspy.file.File(file_path, mode='r')
        data = las_file.points
    else:
        raise ValueError(f"Unsupported file type: {file_extension}. Supported types are ASC and LAS.")

    return data

def process_vegetation_points(las, dtm_data, min_bounding_box_x, min_bounding_box_y):
    for point in las.points:
        x = point.array['X'] / 100
        y = point.array['Y'] / 100
        z = point.array['Z'] / 100
        diff_x = (x - min_bounding_box_x) * 0.9999999
        diff_y = (y - min_bounding_box_y) * 0.9999999
        x_index = int(diff_x / 0.5)
        y_index = 999 - int(diff_y / 0.5)
        if dtm_data[y_index, x_index] < z:
            dtm_data[y_index, x_index] = z

def write_output_tiff(output_file, dtm_data, min_bounding_box_x, max_bounding_box_y):
    top_left_x = min_bounding_box_x
    top_left_y = max_bounding_box_y
    pixel_size_x = 0.5
    pixel_size_y = 0.5
    transform = from_origin(top_left_x, top_left_y, pixel_size_x, pixel_size_y)

    with rasterio.open(output_file, 'w',
                       driver='GTiff',
                       height=dtm_data.shape[0],
                       width=dtm_data.shape[1],
                       count=1,
                       dtype=dtm_data.dtype,
                       crs=rasterio.crs.CRS.from_epsg(7415),
                       transform=transform) as dst:
        dst.write(dtm_data, 1)

def main(pointcloud_file, dtm_file, output_file):
    print(f"Reading DTM file: {dtm_file}")
    dtm_data = read_geospatial_file(dtm_file)

    with laspy.open(pointcloud_file, mode='r') as fh:
        las = fh.read()
        print(f'Points from data: {len(las.points)}')

        min_bounding_box_x = min(las['x'])
        min_bounding_box_y = min(las['y'])

    print("rasterizing vegetation pointcloud...")
    process_vegetation_points(las, dtm_data, min_bounding_box_x, min_bounding_box_y)

    print("Writing output TIFF file...")
    write_output_tiff(output_file, dtm_data, min_bounding_box_x, min_bounding_box_y)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Arguments to be passed: pointcloud_file.las dtm_file.asc output_file.tiff")
    else:
        pointcloud_file = sys.argv[1]
        dtm_file = sys.argv[2]
        output_file = sys.argv[3]
        main(pointcloud_file, dtm_file, output_file)
        print(f'Canopy height raster written to: {output_file}')
