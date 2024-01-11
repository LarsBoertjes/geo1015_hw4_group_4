## This files returns a:
# pandas dataframe in a csv file
# LAS file where every point contains all attributes except the attributes infrared, deviation, amplitude, reflectance

import pandas as pd
from pyntcloud import PyntCloud
import pdal
import json

input_file = '../data/69BZ2_13.LAZ'
# Our boundaries
min_x = 187465.5
max_x = min_x + 500
min_y = 315228.5
max_y = min_y + 500

def get_points(bounds, file):
    # make pipelines to input 
    _pipeline = {
        "pipeline": [
            file,  # input las file
            {
                "type": "filters.crop",  # The crop filter removes points that fall outside or inside a cropping bounding box 
                "bounds": f"([{bounds[0]},{bounds[1]}],[{bounds[2]},{bounds[3]}])"  # input parameters for the filter, in the case the boundaries for our window
            },
            {
                "type": "writers.las", 
                "filename": "cropped_file.las" # Output file name without the dimensions infrared, deviation, amplitude, reflectance
            }
        ]
    }
    pipeline = pdal.Pipeline(json.dumps(_pipeline))  # make pipeline to arrays to stage
    pipeline.execute()  # execute the pipeline
    assert pipeline.arrays  # validate if the pipeline works, it's expected that there will be some arrays in the pipeline result.
    point_cloud = pipeline.arrays[0]  # extract the point cloud data from the first column 
    pd_pcl = pd.DataFrame(point_cloud)  # convert the point_cloud arrays to Panda DataFrame
    attribute_names = list(point_cloud.dtype.names)
    pd_pcl.columns = attribute_names  # the column names for the Panda DataFrame   
    return pd_pcl  # create a PyntCloud object from the Pandas DataFrame and return it.

def main():
    boundary = [min_x, max_x, min_y, max_y]
    pd_pcl = get_points(boundary, input_file)
    pd_pcl.to_csv("cropped_file.csv")

if __name__ == '__main__':
    main()