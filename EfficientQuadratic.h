/***************************************
* An efficient way to evaluate a quadratic in MySeriesTruncMod
***************************************/
// include guard
#ifndef __EfficientQuadratic_H_INCLUDED__
#define __EfficientQuadratic_H_INCLUDED__

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
#include <flint.h>
#include <flintxx.h>
#include <fmpz.h>
#include <fmpz_mod_poly.h>
#include "MySeriesTruncMod.h"


using namespace std;
//=================================

class EfficientQuadratic {
  private:
    int numTB;
    int ** qc;  ///square matrix storing the coefficients of the quadratic combination
    fmpz_t tmp;

  public:
    EfficientQuadratic(vector<int>coeff, vector<int>TBi, vector<int>TBi2);
    ~EfficientQuadratic();
    MySeriesTruncMod* evaluate(MySeriesTruncMod**qs);
};

#endif // __QZSeriesWH_H_INCLUDED__
