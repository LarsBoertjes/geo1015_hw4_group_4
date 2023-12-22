
## Help with this code

### Replacement

  - **startinpy**: the [DT of CGAL](https://doc.cgal.org/latest/Triangulation_2/index.html#Chapter_2D_Triangulations) is used, and the concepts are the same almost as in startinpy (infinite, adjacency, etc.). You need to undertand the data structure used also: https://doc.cgal.org/latest/TDS_2/index.html#Chapter_2D_Triangulation_Data_Structure
  - **laspy**: LAStools is used, and it is a small drama to make it work with LAZ so I copied the necessary code from the original code (https://github.com/LAStools/LAStools/tree/master/LASlib) and it compiles fine (some warning are thrown but you can ignore those)
  - **numpy**: you can just use a 2D list: `std::vector<std::vector<double>>` and I added a class `DatasetASC` that mimics the behaviour of rasterio (kinda)
  - **rasterio**: using GDAL with C++ is doable but installation under Windows is a it rocky I think, so I think it's better if you use the [Esri ASCII raster format (ASC)](https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm). Yes it's from Esri but it's good! I give an example in the `./data/` folder, and the `DatasetASC` has a write() function that works.

### How to install CGAL 

If macOS: `brew install cgal` + `brew install eigen`

If Windows, then I guess using WSL is best here, see instructions there: https://tudelft3d.github.io/geogeeks/cpp/wslclion/

Although I heard 2 days ago that [Anaconda has a CGAL installation](https://anaconda.org/conda-forge/cgal) and one colleague said it works fine (which is incredible but I haven't had time to test it yet).

### How to run the code

```
    mkdir build
    cd build
    cmake ..
    make -j8
    ./hw04
```

If you use CLion, you can just open the folder and the `CMakeLists.txt` will be detected and all should work.

### Useful stuff

  - https://doc.cgal.org/latest/Manual/tuto_gis.html
  - the best docs for C++ functions and data structures is that one IMO: https://cplusplus.com/reference/vector/vector/
  - OBJ specs: https://en.wikipedia.org/wiki/Object_file