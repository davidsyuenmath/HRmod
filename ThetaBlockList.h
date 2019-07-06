/***************************************
* This class encodes Borcherds Product data
***************************************/
// include guard
#ifndef __ThetaBlockList_H_INCLUDED__
#define __ThetaBlockList_H_INCLUDED__

//=================================
// forward declared dependencies, if any
//class Foo;

//=================================
// included dependencies
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>
//#include <gmp.h>
//#include <gmpxx.h>
//#include <mpfr.h>
//#include <mpc.h>
#include <flint.h>
#include <flintxx.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_mat.h>
#include <fmpz_vec.h>
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "QZSeriesWH.h"

using namespace std;
//=================================

class ThetaBlockList {
  private:

  public:
    int num, arraysize;
    int weight, index;
    vector<int>* d;

    ThetaBlockList(int weight, int index, int ell);
    ThetaBlockList(string directory, int weight, int index, int ell);
    ThetaBlockList(string directory, int weight, int index, int ell, int withDenom);
    ~ThetaBlockList(); //destructor
    void destroy();
    void init(string directory, int weight, int index, int ell);
    void initWD(string directory, int weight, int index, int ell, int withDenom);
    void insertThetaBlock(vector<int> c);
    void increaseSize();

};

#endif // __ThetaBlockList_H_INCLUDED__
