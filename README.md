# geo1015_hw4_group_4
Processing point clouds

report editor: https://www.overleaf.com/3335222756vbxkvhxqjdzb#79d9ec

todo's: https://docs.google.com/document/d/1QMNfYO7wIZpxRiS44uYSKETZf_nUEOBe1tBgOgcwHUM/edit


order of code execution:

input                                   code                            output
raw_data.las                    -> create_dtm.cpp               -> dtm.tiff
raw_data.las                    -> cropping_code.py             -> attributes.csv + pointcloud.las
attribute.csv + pointcloud.las  -> vegetation_filter.py         -> vegetationpoints.las
vegetationpoints.las            -> create_canopy_height.py      -> canopy.tiff
canopy.tiff + dtm.tiff          -> create_chm.py                -> chm.tiff


