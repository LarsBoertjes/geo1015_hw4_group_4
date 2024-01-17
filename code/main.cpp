//-- demo code to help you kickstart hw04
//-- Hugo Ledoux <h.ledoux@tudelft.nl>
//-- 2023-12-20

#include <iostream>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;
#include <string>
#include <cmath>


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
std::vector<Point3> extractLowestPtsInGrid (const std::vector<Point3>& points, int gridRows, int gridCols, double min_x, double max_x, double min_y, double max_y);
double barycentricInterpolation(double x, double y, std::vector<double>& x_coord, std::vector<double>& y_coord, std::vector<double>& z_vals);
double maxAngleDegrees(double x, double y, double z, std::vector<double>& x_coord, std::vector<double>& y_coord, std::vector<double>& z_vals);
void laplaceInterpolation(DatasetASC &dtm, const std::vector<Point3> &groundPts, Delaunay dt);
double computeDistance(float Ax, float Ay, float Bx, float By);
bool isValidCell(int row, int col, int numRows, int numCols);

int main(int argc, char** argv)
{
  // read pointcloud from input file
  std::vector<Point3> lsPts = read_lasfile("../data/69BZ2_13.LAZ", 1);

  // specify the bounding box of 500x500 meters
  double min_x = 187465.5;
  double max_x = min_x + 500;
  double min_y = 315228.5;
  double max_y = min_y + 500;

  // extract points within bounding box + margin
  std::vector<Point3> croppedPts;
  for (auto pt : lsPts) {
    if (pt.x() >= min_x - 25 && pt.x() <= max_x + 25 && pt.y() >= min_y - 25 && pt.y() <= max_y + 25) {
      croppedPts.push_back(pt);
    }
  }

  // Step 0 of Ground filtering with TIN refinement: create array with locally lowest points in 2D grid
  std::vector<Point3> groundPts = extractLowestPtsInGrid(croppedPts, 20, 20, min_x - 25, max_x + 25, min_y - 25, max_y + 25);

  // Step 1 of Ground filtering with TIN refinement: Construction of a rudimentary initial TIN
  Delaunay dt;
  Delaunay::Vertex_handle vh;
  for (auto pt : groundPts) {
    vh = dt.insert(Point2(pt.x(), pt.y()));
    vh->info() = pt.z();
  }

  // Step 2: Computation of two geometric properties for each point that is not already labelled as ground
  // push all new ground points to groundPts vector
  double d_max = 1.5;
  double alpha_max = 7;

  for (int i = 0; i < croppedPts.size(); ++i) {
    // find the triangle containing the point
    Point2 pt_xy(croppedPts[i].x(), croppedPts[i].y());
    Face_handle fh = dt.locate(pt_xy);

    std::vector<double> x_coord;
    std::vector<double> y_coord;
    std::vector<double> z_vals;

    for (int i = 0; i < 3; i++) {
      // if one of the vertices is of the big triangle don't add the vertex to list
      if (dt.is_infinite(fh->vertex(i))) {
        continue;
      } else {
        x_coord.push_back(fh->vertex(i)->point().x());
        y_coord.push_back(fh->vertex(i)->point().y());
        z_vals.push_back(fh->vertex(i)->info());
      }
    }

    double heightPlaneAtP;

    if (x_coord.size() == 3) {
      heightPlaneAtP = barycentricInterpolation(croppedPts[i].x(), croppedPts[i].y(), x_coord, y_coord, z_vals);
    } else if (x_coord.size() == 2) {
      // else do a weighted average of the 2 finite points
    } heightPlaneAtP = (z_vals[0] + z_vals[1]) / 2;

    double heightDifference = abs(heightPlaneAtP - croppedPts[i].z());
    double maxAngle = maxAngleDegrees(croppedPts[i].x(), croppedPts[i].y(), croppedPts[i].z(), x_coord, y_coord, z_vals);

    // check if point is within groundPts criteria
    if (maxAngle < alpha_max && heightDifference < d_max) {
      groundPts.push_back(croppedPts[i]);
      vh = dt.insert(Point2(croppedPts[i].x(), croppedPts[i].y()));
      vh->info() = croppedPts[i].z();
    }
    
  }

  // creating empty grid and filling it using laplacaInterpolation
  DatasetASC d = DatasetASC(1000, 1000, min_x, min_y, 0.5, -9999.0);

  laplaceInterpolation(d, groundPts, dt);

  // Write the gridded DTM to an ASC file
  d.write("../data/dtm.asc");
    
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

std::vector<Point3> extractLowestPtsInGrid (const std::vector<Point3>& points, int gridRows, int gridCols, double min_x, double max_x, double min_y, double max_y) {
  double cellWidth = (max_x - min_x) / gridCols;
  double cellHeight = (max_y - min_y) / gridRows;

  // create a 2D vector to store the lowest point in each cell
  std::vector<std::vector<Point3>> lowestPoints(gridRows, std::vector<Point3>(gridCols, Point3()));

  // initialize lowest z values to positive infinity
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

double barycentricInterpolation(double x, double y, std::vector<double>& x_coord, 
                                                    std::vector<double>& y_coord, 
                                                    std::vector<double>& z_vals) {
    
    double x1 = x_coord[0];
    double y1 = y_coord[0];
    double z1 = z_vals[0];

    double x2 = x_coord[1];
    double y2 = y_coord[1];
    double z2 = z_vals[1];

    double x3 = x_coord[2];
    double y3 = y_coord[2];
    double z3 = z_vals[2];

    double detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
    double alpha = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
    double beta = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
    double gamma = 1.0 - alpha - beta;

    double interpolatedValue = alpha * z1 + beta * z2 + gamma * z3;

    return interpolatedValue;
}

double maxAngleDegrees(double x, double y, double z, std::vector<double>& x_coord,
                                                 std::vector<double>& y_coord,
                                                 std::vector<double>& z_vals) {

    double maxAngle = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < x_coord.size(); ++i) {
      double xy_diff = std::sqrt(std::pow(x - x_coord[i], 2) + std::pow(y - y_coord[i], 2));
      double z_diff = z - z_vals[i];

      // calculate angle between two points using the xy_diff and z_diff
      double angle = atan(z_diff / xy_diff);

      if (angle > maxAngle) {
        maxAngle = angle;
      }
    }

    return maxAngle * (180.0 / M_PI); 
}

void laplaceInterpolation(DatasetASC &dtm, const std::vector<Point3> &groundPts, Delaunay dt) {
  for (int row = 0; row < dtm._nrows; ++row) {
    for (int col = 0; col < dtm._ncols; ++col) {

      // get coordinates of the center of the cells
      double cell_center_x = dtm._xllcorner + (dtm._cellsize * col) + (dtm._cellsize / 2.0);
      double cell_center_y = dtm._yllcorner + (dtm._cellsize * row) + (dtm._cellsize / 2.0);
      Point2 cell_center(cell_center_x, cell_center_y);  
      
      // insert vertex_cell_center into DT and iterate over incident vertices
      Delaunay::Vertex_handle vertex_cell_center = dt.insert(cell_center);
      Delaunay::Vertex_circulator vc = dt.incident_vertices(vertex_cell_center);
      Delaunay::Vertex_circulator done = vc;

      // circulate over incident faces of vertex_cell_center
      Delaunay::Face_circulator fc = dt.incident_faces(vertex_cell_center);
      Delaunay::Face_circulator done_face = fc;

      int incidentVertexNum = 0;
      int incidentFaceNum = 0;          
          
      double sum_weighted_height = 0;
      double total_weight = 0;

      double edgeVx;
      double xpi;

      if (vc != 0 && fc != 0) {
        do {
          // iterate over the faces from face_circulator fc 
          Delaunay::Vertex_handle vh1 = fc->vertex(0);
          Delaunay::Vertex_handle vh2 = fc->vertex(1);
          Delaunay::Vertex_handle vh3 = fc->vertex(2);

          std::vector<Delaunay::Vertex_handle> secondTri;

          for (int i = 0; i < 3; ++i) {
            std::vector<Point2> faceVertices = {};
            Delaunay::Vertex_handle uh1 = fc->neighbor(i)->vertex(0);
            Delaunay::Vertex_handle uh2 = fc->neighbor(i)->vertex(1);
            Delaunay::Vertex_handle uh3 = fc->neighbor(i)->vertex(2);
            faceVertices.push_back(uh1->point());
            faceVertices.push_back(uh2->point());
            faceVertices.push_back(uh3->point());

            // check if cell center & vc->point() are in faceVertices
            if (std::find(faceVertices.begin(), faceVertices.end(), cell_center) != faceVertices.end() &&
                std::find(faceVertices.begin(), faceVertices.end(), vc->point()) != faceVertices.end()) {
                    secondTri.push_back(uh1);
                    secondTri.push_back(uh2);
                    secondTri.push_back(uh3);
                  } 
                } 

            Point2 u1 = secondTri[0]->point();
            Point2 u2 = secondTri[1]->point();
            Point2 u3 = secondTri[2]->point();

            Point2 circumcenter_a = circumcenter(vh1->point(), vh2->point(), vh3->point());
            Point2 circumcenter_b = circumcenter(u1, u2, u3);

            edgeVx = computeDistance(circumcenter_a.x(), circumcenter_a.y(), circumcenter_b.x(), circumcenter_b.y());
            xpi = computeDistance(vertex_cell_center->point().x(), vertex_cell_center->point().y(), vc->point().x(), vc->point().y());

            double weight_i = edgeVx / xpi;
            double weighted_height = weight_i * vc->info();

            // handling of nodata points
            if ((!isnan(edgeVx) && edgeVx > 0) && (!isnan(xpi) && xpi > 0)) {
              sum_weighted_height += weighted_height;
              total_weight += weight_i;
            }

            ++vc; 
            ++fc;
            } while (vc != done);
          }

          double interpolated_z = sum_weighted_height / total_weight;

          // remove vertex_cell_center from DT and proceed
          dt.remove(vertex_cell_center);

          // add interpolated z value to the empty grid
          dtm.data[dtm._nrows - row - 1][col] = interpolated_z;
        }
      }
    }


double computeDistance(float Ax, float Ay, float Bx, float By) {
  double x_diff = abs(Ax - Bx);
  double y_diff = abs(Ay - By);

  double distance = sqrt((x_diff * x_diff) + (y_diff * y_diff));

  return distance;
}

bool isValidCell(int row, int col, int numRows, int numCols) {
  return (row >= 0 && row < numRows && col >= 0 && col < numCols);
}