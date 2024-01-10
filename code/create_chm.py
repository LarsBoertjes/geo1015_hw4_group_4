import argparse
import rasterio
from rasterio.transform import from_origin
import numpy as np

def subtract_rasters(canopy_height, dtm, output_path):
    # Read the first raster
    with rasterio.open(canopy_height) as src1:
        data1 = src1.read(1)
        transform1 = src1.transform

    # Read the second raster
    with rasterio.open(dtm) as src2:
        data2 = src2.read(1)
        transform2 = src2.transform

    # Check if the rasters have the same shape and transformation
    if data1.shape != data2.shape or transform1 != transform2:
        raise ValueError("Input rasters must have the same shape and transformation.")

    # Perform the subtraction
    result_data = data1 - data2

    # Write the result to a new raster file
    with rasterio.open(
        output_path,
        'w',
        driver='GTiff',
        height=result_data.shape[0],
        width=result_data.shape[1],
        count=1,
        dtype=result_data.dtype,
        crs=src1.crs,
        transform=transform1,
    ) as dst:
        dst.write(result_data, 1)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Create CHM from canopy height and DTM.")
    parser.add_argument("dtm", help="Path to the DTM raster file (.asc)")
    parser.add_argument("canopy_height", help="Path to the canopy_height raster file (.tiff)")
    parser.add_argument("output", help="Path to the output raster file.")
    args = parser.parse_args()

    # Call the subtract_rasters function with the provided arguments
    subtract_rasters(args.canopy_height, args.dtm, args.output)
    print(f'created CHM raster as: {args.output}')

if __name__ == "__main__":
    main()

