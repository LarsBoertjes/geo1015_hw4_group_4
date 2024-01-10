import laspy
import numpy as np

las = laspy.read('data\cropped_file.las')

unclassified_points = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
unclassified_points.points = las.points[las.classification == 1]

points_with_one_return = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
points_with_one_return.points = las.points[las.number_of_returns == 1]

# Extract point positions as a NumPy array
# points1 = np.vstack((unclassified_points.x, unclassified_points.y, unclassified_points.z)).T
# points2 = np.vstack((points_with_one_return.x, points_with_one_return.y, points_with_one_return.z)).T

#use set operations on coordinates
common_points_indices = np.setdiff1d(np.where(las.classification == 1)[0], np.where(las.number_of_returns == 1)[0])

# Create a new LAS file with the common points
result = laspy.create(point_format=las.header.point_format, file_version=las.header.version)
result.points = las.points[common_points_indices]

result.write('unclassified_points_with_one_return.las')




# classification = point.array['classification']
# in array the following can be extracted (according to point format 0 in documentation:
# intensity, return_number, number_of_returns, scan_direction_flag, edge_of_flight_line, classification,
# synthetic, key_point, withheld, scan_angle_rank, user_data, point_source_id

