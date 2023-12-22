
#include <vector>

class DatasetASC
{
  public:
    int _nrows;
    int _ncols;
    double _xllcorner;
    double _yllcorner;
    double _cellsize;
    double _nodata;
    std::vector<std::vector<double>> data;

    DatasetASC(int ncols, int nrows, double xllcorner, double yllcorner, double cellsize, double nodata);

    void write(std::string ofile);
    bool rc2xy(int row, int col, double& x, double& y);
    bool xy2rc(double x, double y, int& row, int& col);

}; 