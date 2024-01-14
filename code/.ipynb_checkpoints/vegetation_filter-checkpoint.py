import sys
import laspy
import numpy as np
import pandas as pd

def filter_vegetation(input_file_las, input_file_csv, output_file):
    las = laspy.read(input_file_las)
    pd_pcl = pd.read_csv(input_file_csv)
    
    unclassified_points = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
    unclassified_points.points = las.points[las.classification == 1]

    pd_classification_1 = pd_pcl[pd_pcl['Classification']==1]
    pd_classification_1.loc[:, 'NDVI'] = pd_classification_1.apply(lambda row: (row['Infrared'] - row['Red']) / (row['Infrared'] + row['Red']), axis=1)
    pd_classified_NDVI = pd_classification_1[(pd_classification_1['NDVI'] > 0.2) | (pd_classification_1['NDVI'] > 1)]
    index_NDVI = list(pd_classified_NDVI.index)
    
    # use set operations on coordinates
    common_points_indices = np.setdiff1d(np.where(las.classification == 1)[0], np.where(las.number_of_returns == 1)[0])

    final_index = set(list(list(common_points_indices) + index_NDVI))
    # Create a new LAS file with the common points
    result = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
    result.points = las.points[list(final_index)]

    result.write(output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: !python vegetation_filter.py input.las input.csv output.las")
        sys.exit(1)

    input_file_las = sys.argv[1]
    input_file_csv = sys.argv[2]
    output_file = sys.argv[3]

    filter_vegetation(input_file_las, input_file_csv, output_file)

    print(f'filtered pointcloud written to {output_file}')
