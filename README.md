# geo1015_hw4_group_4
Lars Boertjes - 4704541,
Noah Alting - 4828968,
Marieke van Arnhem - 4918738

*** DOCUMENTATION ***

report editor: https://www.overleaf.com/3335222756vbxkvhxqjdzb#79d9ec

todo's: https://docs.google.com/document/d/1QMNfYO7wIZpxRiS44uYSKETZf_nUEOBe1tBgOgcwHUM/edit

--------------------------------------------------------------------------
*** INSTALLING PACKAGES ***

Being able to run the C++ code, make sure the following libraries are installed in your C++ environment:
- CGAL
- lasreader
- laswriter
- DatasetASC writer/reader

Being able to run the Python code, make sure the following libraries are installed in your Python environment:
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
Go to main.cpp and check the input files and boundaries:
- read_lasfile(path, thin) (line 48) : make sure you have the required AHN4 tile on your system and specify the path as the first parameter for the read function.
- specify your min_x, min_y, max_x, max_y in lines 51-54.
- If you want to change the GFTIN parameters you can do so in lines 77 & 78.
- Specify the name and output path of your ASCII output file in line 126.
- Make an run main.cpp file
- Also make sure to use the DataSetASC.cpp file we uploaded in this repository.  

STEP 2.
Go to main.py and check the input files and boundaries:
- input_las_file (line 2): make sure you have the required AHN4 tile in the "code" folder and changed the name in line 2 to your AHN4 tile name.
- input_dtm_file (line 4): make sure you have completed step 1 above. This will in fact provide the required dtm (dtm.asc) of your AHN4 tile in the "code" folder.
- min_x, min_y, max_x, max_y (lines 7-10): Change the x and y values to get the desired boundaries for your cropped file. Change min_x and min_y to the minimum x and minimum y values for your cropped file. Maximum x and maximum y are automatically generated since it is a 500x500m grid.

Now run the main.py file, this will create five files in the "code" folder. If you have not changed the names of the output files in lines 19-23, the following files are created:
- "cropped_file.las": the cropped area as LAS file. 
- "cropped_file.csv": the cropped area exported as CSV.
- "classified_vegetation.las": all points that are classified (by our algorithm) as vegetation.
- "canopy_height.tiff": the DSM of vegetation.
- "chm.tiff": the canopy heigt model (CHM) of your area.

--------------------------------------------------------------------------
*** RUNNING INDIVIDUAL STEPS OF THE ASSIGNMENT *** 

If you want to check one of the functions (the assignments steps seperately), please read the following. Below it is possible to see which file belongs to which step of the assignment. You can run the seperate python files in the console line.

- Step 2: cropping_code.py crops the AHN4 tile into a 500mX500M area. Input is the AHN4 tile and the boundaries. The output is an LAS file and CSV file that contain all points in that area defined by the boundaries. The LAS file does not contain the attributes infrared, deviation, amplitude, reflectance. The CSV does. Type in the console line:

!python cropping_code.py input.las min_x max_x min_y max_y output.las output.csv

Where min_x and max_x are the minimum and maximum x values of the cropped area, respectively. And min_y and max_y are the minimum and maximum y values of the cropped area, respectively.

- Step 4: vegetation_filter.py returns only the points that are classified by this algorithm as vegetation. The function needs as input as LAS file and CSV file. The output is a LAS file with all the points that are classified as vegetation. Type in the console line:

!python vegetation_filter.py input.las input.csv output.las

- Step 4: create_canopy_height.py creates a grid that contains the vegetation on a DTM. Input is the LAS file with all the points that are classified as vegetation and the DTM. Output is a canopy height tiff file. Type in the console line:

!python pointcloud_file.las dtm_file.asc output_file.tiff

- Step 5: create_chm.py creates the CHM. Input is the DSM of vegetation tiff and DTM file. Output is a CHM tiff file. Type in the console line:

!python create_chm.py canopy_height.tiff dtm_file.asc output_file.tiff
