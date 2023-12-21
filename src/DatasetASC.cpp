
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "DatasetASC.h"


DatasetASC::DatasetASC(int nrows, int ncols, double xllcorner, double yllcorner, double cellsize, double nodata) {
  _ncols = ncols;
  _nrows = nrows;
  _xllcorner = xllcorner;
  _yllcorner = yllcorner;
  _cellsize = cellsize;
  _nodata = nodata;
  
  //-- read the data in a 2D array (here a vector or vector)
  //-- read each row/line
  for (int i = 0; i < _nrows; i++) {
    std::vector<double> onerow(_ncols, 0.0);
    data.push_back(onerow);
  }
}



void DatasetASC::write(std::string ofile) 
{
  std::ofstream sout(ofile);
  //-- write to the file
  sout << "ncols" << " " << _ncols << std::endl;
  sout << "nrows" << " " << _nrows << std::endl;
  sout << std::setprecision(3) << std::fixed << "xllcorner" << " " << _xllcorner << std::endl;
  sout << std::setprecision(3) << std::fixed << "yllcorner" << " " << _yllcorner << std::endl;
  sout << "cellsize" << " " << _cellsize << std::endl;
  // nodata = 99999;
  sout << "NODATA_value" << " " << _nodata << std::endl;
  for (int i = 0; i < _nrows; i++) {
    for (int j = 0; j < _ncols; j++) {
      // if (data[i][j] > (_nodata - 1))
      //   sout << _nodata << " ";  
      // else
        sout << data[i][j] << " ";  
    }
    sout << std::endl;
    
  }
  // Close the file
  sout.close();
  std::cout << "ASC file written to: " << ofile << std::endl;

}



//-- returns centre of the pixel from (row, col)
//-- origin is top-left, 0-based index
//-- false if outside
bool DatasetASC::rc2xy(int row, int col, double& x, double& y)
{
  if ( (row < 0) || (row > (_nrows - 1)) ) {
    return false;
  }
  if ( (col < 0) || (col > (_ncols - 1)) ) {
    return false;
  }
  x = _xllcorner + (col * _cellsize) + (_cellsize / 2);
  y = _yllcorner + ((_nrows - 1 - row) * _cellsize) + (_cellsize / 2);
  return true;
}


//-- returns the row/col in the variables passed by reference
//-- if xy inside the ds, then true; otherwise false
bool DatasetASC::xy2rc(double x, double y, int& row, int& col){
  if ( (x < _xllcorner) || (x > (_xllcorner + (_ncols * _cellsize))) ) {
    return false;
  }
  if ( (y < _yllcorner) || (y > (_yllcorner + (_nrows * _cellsize))) ) {
    return false;
  }
  col = (x - _xllcorner) / _cellsize;
  row = _nrows - ((y - _yllcorner) / _cellsize);
  return true;
}


