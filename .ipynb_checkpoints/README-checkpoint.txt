STEP 1.
C++ code explaining - Lars, die moet namelijk als eerst

STEP 2.
Go to main.py and check the input files and boundaries:
- input_las_file: make sure to drag the needed AHN4 tile to the folder "code".
- input_dtm_file: make sure you completed the steps above to create a dtm.tif file in the folder "code".
- change the boundaries to the desired x- and y-value. Change min_x and min_y to the minimum x and minimum y values for your cropped file. Maximum x and maximum y are automatically generated since it is a 500x500m grid.
Now run the main.py file, this will create five files in the "code" folder: output_las_cropped_file, output_csv_cropped_file, output_vegetation_las, canopy_height_tiff, output_chm_file.
The output_chm_file is the CHM file.

---------------------
If you want to check one of the functions, please read the following. Below it is possible to see which file belongs to which step of the assignment.
In total four functions are imported from four python files. Description per python file:
- Step 2: cropping_code.py crops the AHN4 tile into a 500mX500M area. Input is the AHN4 tile and the boundaries. The output is an LAS file and CSV file that contain all points in that area defined by the boundaries. The LAS file does not contain the attributes infrared, deviation, amplitude, reflectance. The CSV does.
- Step 4: vegetation_filter.py returns only the points that are classified by this algorithm as vegetation. The function needs as input as LAS file and CSV file. The output is a LAS file with all the points that are classified as vegetation.
- Step 4: create_canopy_height.py creates a grid that contains the vegetation on a DTM. Input is the LAS file with all the points that are classified as vegetation and the DTM. Output is a canopy height tiff file.
- Step 5: create_chm.py creates the CHM. Input is the canopy height tiff and DTM tiff file. Output is a CHM tiff file.