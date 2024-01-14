*** INSTALLING PACKAGES ***
In able to run the C++ code, make sure the following libraries are installed in your C++ environment:
- @Lars

In able to run the Python code, make sure the following libraries are installed in your Python environment:
- pandas
- pdal
- json
- sys
- laspy
- numpy
- argparse
- rasterio

--------------------------------------------------------------------------
*** CREATING CHM.TIFF FILE ***
STEP 1.
C++ code explaining - @Lars, die moet namelijk als eerst

STEP 2.
Go to main.py and check the input files and boundaries:
- input_las_file (line 2): make sure you have the required AHN4 tile in the "code" folder.
- input_dtm_file (line 4): make sure you have completed step 1 above. This will in fact provide the required dtm (dtm.asc) of your AHN4 tile in the "code" folder.
- min_x, min_x, max_x, max_y (lines 7-10): Change the x and y values to get the desired boundaries for your cropped file. Change min_x and min_y to the minimum x and minimum y values for your cropped file. Maximum x and maximum y are automatically generated since it is a 500x500m grid.

Now run the main.py file, this will create five files in the "code" folder. If you have not changed the names of the output files in lines 19-23, the following files are created:
"cropped_file.las": the cropped area as LAS file. 
"cropped_file.csv": the cropped area exported as CSV.
"classified_vegetation.las": all points that are classified (by our algorithm) as vegetation.
"canopy_height.tiff": the DSM of vegetation.
"chm.tiff": the canopy heigt model (CHM) of your area.

--------------------------------------------------------------------------
*** RUNNING INDIVIDUAL STEPS OF THE ASSIGNMENT *** 
If you want to check one of the functions (the assignments steps seperately), please read the following. Below it is possible to see which file belongs to which step of the assignment. You can run the seperate python files in the console line.

- Step 2: cropping_code.py crops the AHN4 tile into a 500mX500M area. Type in the console line:
!python cropping_code.py input.las min_x max_x min_y max_y output.las output.csv
Where min_x and max_x are the minimum and maximum x values of the cropped area, respectively. And min_y and max_y are the minimum and maximum y values of the cropped area, respectively.
Input is the AHN4 tile and the boundaries. The output is an LAS file and CSV file that contain all points in that area defined by the boundaries. The LAS file does not contain the attributes infrared, deviation, amplitude, reflectance. The CSV does.

- Step 4: vegetation_filter.py returns only the points that are classified by this algorithm as vegetation. Type in the console line:
!python vegetation_filter.py input.las input.csv output.las
The function needs as input as LAS file and CSV file. The output is a LAS file with all the points that are classified as vegetation.

- Step 4: create_canopy_height.py creates a grid that contains the vegetation on a DTM. Type in the console line:
!python pointcloud_file.las dtm_file.asc output_file.tiff
Input is the LAS file with all the points that are classified as vegetation and the DTM. Output is a canopy height tiff file.

- Step 5: create_chm.py creates the CHM. Type in the console line:
!python create_chm.py canopy_height.tiff dtm_file.asc output_file.tiff
Input is the DSM of vegetation tiff and DTM file. Output is a CHM tiff file.