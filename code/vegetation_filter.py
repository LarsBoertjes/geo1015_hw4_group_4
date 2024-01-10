import sys
import laspy
import numpy as np

def filter_vegetation(input_file, output_file):
    las = laspy.read(input_file)

    unclassified_points = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
    unclassified_points.points = las.points[las.classification == 1]

    points_with_one_return = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
    points_with_one_return.points = las.points[las.number_of_returns == 1]

    # use set operations on coordinates
    common_points_indices = np.setdiff1d(np.where(las.classification == 1)[0], np.where(las.number_of_returns == 1)[0])

    # Create a new LAS file with the common points
    result = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
    result.points = las.points[common_points_indices]

    result.write(output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vegetation_filter.py input.las output.las")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_vegetation(input_file, output_file)

    print(f'filtered pointcloud written to {output_file}')
