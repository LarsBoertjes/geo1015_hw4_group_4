## This files returns a:
# pandas dataframe in a csv file
# LAS file where every point contains all attributes except the attributes infrared, deviation, amplitude, reflectance

import pandas as pd
import pdal
import json
import sys

def get_points(bounds, input_las_file, output_las_file, output_csv_file):
    # make pipelines to input 
    _pipeline = {
        "pipeline": [
            input_las_file,  # input las file
            {
                "type": "filters.crop",  # The crop filter removes points that fall outside or inside a cropping bounding box 
                "bounds": f"([{bounds[0]},{bounds[1]}],[{bounds[2]},{bounds[3]}])"  # input parameters for the filter, in the case the boundaries for our window
            },
            {
                "type": "writers.las", 
                "filename": output_las_file # Output file name
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
    pd_pcl.to_csv(output_csv_file)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: !python cropping_code.py input.las min_x max_x min_y max_y output.las output.csv")
        print("Extra information: [min_x, max_x, min_y, max_y] are the boundaries of your cropped file.")
        sys.exit(1)

    input_file_las = sys.argv[1]
    bounds = [sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]]
    output_las_file = sys.argv[6]
    output_csv_file = sys.argv[7]

    get_points(bounds, input_file_las, output_las_file, output_csv_file)

    print(f'cropped LAS file written to {output_las_file} {output_csv_file}')