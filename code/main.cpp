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
std::vector<Point3> extractLowestPtsInGrid (const std::vector<Point3>& points, int gridRows, int gridCols, double min_x, double max_x, double min_y, double max_y);

int main(int argc, char** argv)
{
  // read pointcloud from input file
  std::vector<Point3> lsPts = read_lasfile("../data/69BZ2_13.LAZ", 1000);

  // specify the bounding box of 500x500 meters
  double min_x = 187465.5;
  double max_x = min_x + 500;
  double min_y = 315228.5;
  double max_y = min_y + 500;

  std::cout << "Bounding box EPSG:28992 " << std::endl << "--------------" << std::endl << "x-range: " << min_x << " to " << max_x << std::endl << "y-range: " << min_y << " to " << max_y << "\n\n";
  

  // extract points within bounding box
  std::vector<Point3> croppedPts;
  for (auto pt : lsPts) {
    if (pt.x() >= min_x && pt.x() <= max_x && pt.y() >= min_y && pt.y() <= max_y) {
      croppedPts.push_back(pt);
    }
  }

  // lines to check if it did actually crop the points
  std::cout << "Number of points " << lsPts.size() << std::endl;
  std::cout << "Number of points after cropping " << croppedPts.size() << std::endl;

  // Step 0 of Ground filtering with TIN refinement: create array with locally lowest points in 2D grid
  // -- The size of the largest building is around 23x20 meters so we will use a grid resolution of almost 30x30m (500x500 / 17 rows, columns = 29.4 x 29.4)
  // -- Divide the 500x500m into a grid with 17 columns and rows and per cell push the cell with the lowest pt.z() value to groundPts
  std::vector<Point3> groundPts = extractLowestPtsInGrid(croppedPts, 17, 17, min_x, max_x, min_y, max_y);

  std::cout << "Number of points in groundPts " << groundPts.size() << std::endl;
  // this should be 17 * 17 = 289

  // Step 1 of Ground filtering with TIN refinement: Construction of a rudimentary initial TIN
  Delaunay dt;
  Delaunay::Vertex_handle vh;
  for (auto pt : groundPts) {
    vh = dt.insert(Point2(pt.x(), pt.y()));
    vh->info() = pt.z();
  }

  std::cout << "DT #vertices: " << dt.number_of_vertices() << std::endl;
  std::cout << "DT #triangles: " << dt.number_of_faces() << std::endl;

  // Step 2: Computation of two geometric properties for each point that is not already labelled as ground
  // push all new ground points to groundPts vector
  double d_max = 5;
  double alpha_max = 30;

  for (auto pt : croppedPts) {
    // check if pt is not already in groundPts
    if (std::find(groundPts.begin(), groundPts.end(), pt) == groundPts.end()) {
      // find the triangle containing the point
      Point2 pt_xy(pt.x(), pt.y());
      Face_handle fh = dt.locate(pt_xy);

      // to do: check if distance between pt.z() and plane fh is not more than d_max
      //double d = dt.
      
      for (int i = 0; i < 3; i++) {
        double max_alpha = 0.0;

        if (dt.is_infinite(fh->vertex(i))) {
          // do not do computations using infinite vertices
          std::cout << fh->vertex(i)->point() << " is infinite vertex" << std::endl;
        } else {
          std::cout << fh->vertex(i)->point() << std::endl;
          std::cout << "z van het punt: " << fh->vertex(i)->point() << std::endl;
          // to do: compute angle between 2 points by comparing the horizontal difference and vertical difference
          // to do: store angle in max_alpha if angle > max_alpha

        }
      }

      std::cout << "--------------------" << std::endl;

    }
  }
  

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

std::vector<Point3> extractLowestPtsInGrid (const std::vector<Point3>& points, int gridRows, int gridCols, double min_x, double max_x, double min_y, double max_y) {
  double cellWidth = (max_x - min_x) / gridCols;
  double cellHeight = (max_y - min_y) / gridRows;

  // create a 2D vector to store the lowest point in each cell
  std::vector<std::vector<Point3>> lowestPoints(gridRows, std::vector<Point3>(gridCols, Point3()));

  // initialize lowest z values to positve infinity
  std::vector<std::vector<double>> lowestZValues(gridRows, std::vector<double>(gridCols, std::numeric_limits<double>::infinity()));

  // iterate through the points and update the lowest point in each cell
    for (const auto& pt : points) {
        int col = static_cast<int>((pt.x() - min_x) / cellWidth);
        int row = static_cast<int>((pt.y() - min_y) / cellHeight);

        if (col >= 0 && col < gridCols && row >= 0 && row < gridRows) {
            if (pt.z() < lowestZValues[row][col]) {
                lowestZValues[row][col] = pt.z();
                lowestPoints[row][col] = pt;
            }
        }
    }

    // flatten the 2D vector into a 1D vector of initial groundPoints
    std::vector<Point3> result;
    for (int row = 0; row < gridRows; ++row) {
      for (int col = 0; col < gridCols; ++col) {
        result.push_back(lowestPoints[row][col]);
      }
    }

    return result;
}
