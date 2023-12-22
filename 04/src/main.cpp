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

int main(int argc, char** argv)
{
  //-- read pointcloud from input file
  std::vector<Point3> lsPts = read_lasfile("../data/crop.laz", 1);
  // std::vector<Point3> lsPts = read_lasfile("../data/crop.laz", 100);
  std::cout << "Number of points " << lsPts.size() << std::endl;
  std::cout << lsPts[2] << std::endl;

  Delaunay dt;
  Delaunay::Vertex_handle vh;
  for (auto pt : lsPts) {
    vh  = dt.insert(Point2(pt.x(), pt.y()));
    vh->info() = pt.z();
  }

  std::cout << "DT #vertices: " << dt.number_of_vertices() << std::endl;
  std::cout << "DT #triangles: " << dt.number_of_faces() << std::endl;

  //-- save the DT to an OBJ file called "mydt.obj"
  write_obj(dt);

  //-- find triangle containing a point q
  Point2 q(84656.0, 447059.1);
  Face_handle fh = dt.locate(q); //-- you get a Face: https://doc.cgal.org/latest/TDS_2/classTriangulationDataStructure__2_1_1Face.html
  std::cout << fh->vertex(0)->point() << std::endl;
  std::cout << fh->vertex(1)->point() << std::endl;
  std::cout << fh->vertex(2)->point() << std::endl;
  std::cout << "area triangle: " << CGAL::area(fh->vertex(0)->point(), 
                                               fh->vertex(1)->point(), 
                                               fh->vertex(2)->point()) << std::endl;


  //--fetch its neighbour (there are 3)
  Face_handle n0 = fh->neighbor(0);
  std::cout << dt.is_infinite(n0) << std::endl;
  std::cout << n0->vertex(0)->point() << std::endl;
  std::cout << n0->vertex(0)->info() << std::endl;
  //-- equals a rasterio dataset
  DatasetASC d = DatasetASC(3, 4, 100.0, 100.0, 10.0, -9999.0);
  
  //-- you could modify it
  d.data[1][2] = 666.666;
  //-- row-col => x-y
  double x, y;
  d.rc2xy(2, 0, x, y);
  std::cout << "(" << x << ", " << y << ")" << std::endl;
  //-- x-y => row-col
  int r, c;
  d.xy2rc(138, 122, r, c);
  std::cout << "(" << r << ", " << c << ")" << std::endl;

  //-- write the array to a ASC file 
  d.write("out.asc"); 

  //-- we're done, return 0 to say all went fine
  return 0;
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
