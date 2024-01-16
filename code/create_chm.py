import argparse
import rasterio
import sys
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

    tolerance = 1e-8
    
    # Check if the rasters have the same shape and transformation
    if data1.shape != data2.shape or not np.allclose(transform1, transform2, atol=tolerance):
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

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: !python create_chm.py canopy_height.tiff dtm_file.asc output_file.tiff")
    else:
        canopy_height = sys.argv[1]
        dtm_file = sys.argv[2]
        output_file = sys.argv[3]
        subtract_rasters(canopy_height, dtm_file, output_file)
