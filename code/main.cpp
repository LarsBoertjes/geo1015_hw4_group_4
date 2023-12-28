//-- demo code to help you kickstart hw04
//-- Hugo Ledoux <h.ledoux@tudelft.nl>
//-- 2023-12-20

#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
#include <string>


//-- we're using CGAL DT instead of startinpy
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

//-- definition for using CGAL, no need to grasp these
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<double, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;
typedef K::Point_2  Point2;
typedef K::Point_3  Point3;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Face_handle Face_handle;

// -- LAS reading and writing
#include <lasreader.hpp>
#include <laswriter.hpp>


#include "DatasetASC.h"

//-- forward declarations of all functions
std::vector<Point3> read_lasfile(std::string filename, int thin = 1);
void write_obj(Delaunay &dt);
void write_lasfile(const char filename[256], std::vector<Point3> points);

int main(int argc, char** argv)
{
  // read pointcloud from input file
  std::vector<Point3> lsPts = read_lasfile("../data/69BZ2_13.LAZ", 1000);
  double shift_x = 186980.00;
  double shift_y = 314980.00;

  // specify the bounding box of 500x500 meters
  double min_x = 187465.5;
  double max_x = min_x + 500;
  double min_y = 315228.5;
  double max_y = min_y + 500;

  std::cout << "Bounding box EPSG:28992 " << std::endl << "--------------" << std::endl << "x-range: " << min_x << " to " << max_x << std::endl << "y-range: " << min_y << " to " << max_y << "\n\n";
  std::cout << "Bounding box (local CloudCompare) " << std::endl << "--------------" << std::endl << "x-range: " << min_x - shift_x << " to " << max_x - shift_x << std::endl << "y-range: " << min_y - shift_y << " to " << max_y - shift_y << "\n\n";

  // extract points within bounding box
  std::vector<Point3> filteredPts;
  for (auto pt : lsPts) {
    if (pt.x() >= min_x && pt.x() <= max_x && pt.y() >= min_y && pt.y() <= max_y) {
      filteredPts.push_back(pt);
    }
  }

  // write filteredPts to LAS file, so I can check it using cloudcompare
  std::string ofilename = "cropped.las";
  //write_lasfile(ofilename.c_str(), filteredPts);

  // lines to check if it did actually crop the points
  std::cout << "Number of points " << lsPts.size() << std::endl;
  std::cout << "Number of points after cropping " << filteredPts.size() << std::endl;


}

void write_obj(Delaunay& dt) {
  std::map<Vertex_handle,int> vids;
  std::ofstream sout("mydt.obj");
  sout << std::setprecision(5) << std::fixed;
  int count = 1;
  for (Vertex_handle vh : dt.finite_vertex_handles()) {
    vids[vh] = count;
    count++;
    sout << "v " << vh->point().x() << " " << vh->point().y() << " " << vh->info() << std::endl;
  }
  for (Face_handle fh : dt.finite_face_handles()) {
    sout << "f " << vids.find(fh->vertex(0))->second << " " << vids.find(fh->vertex(1))->second << " " << vids.find(fh->vertex(2))->second << std::endl;
  }
  sout.close();
  std::cout << "OBJ file written to: 'mydt.obj'" << std::endl;
}

std::vector<Point3> read_lasfile(std::string filename, int thin) {
  /*
  Function to read points from a LAS/LAZ file

  Inputs:
    filename:   the filename to read the LAS file to

  Returns:
    a std::vector<Point3> with the points from the LAS file
  */
  LASreadOpener lasreadopener;
  lasreadopener.set_file_name(filename.c_str());
  LASreader* lasreader = lasreadopener.open();
  if (!lasreader){
    std::cerr << "cannot read las file: " << filename << "\n";
    exit(1);
  }
  //-- store each point3 in a CGAL Point_3 object
  //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
  std::vector<Point3> points;
  int count = 0;
  while (lasreader->read_point()) {
    // auto c = int(lasreader->point.get_classification());
    // std::cout << c << std::endl;
    // auto i = int(lasreader->point.get_intensity());
    if (count % thin == 0)  {
      points.push_back(
        Point3(
          lasreader->point.get_x(),
          lasreader->point.get_y(),
          lasreader->point.get_z()
        )
      );
    }
    count++;
  }
  lasreader->close();
  delete lasreader;
  return points;
}


