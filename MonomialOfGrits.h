/***************************************
* An efficient way to evaluate a Monomial in MySeriesTruncMod
***************************************/
// include guard
#ifndef __MonomialOfGrits_H_INCLUDED__
#define __MonomialOfGrits_H_INCLUDED__

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

class MonomialOfGrits {
  private:
    int maxTBnum;
    vector<int> ind; 
    vector<int> exponent;

  public:
    MonomialOfGrits( vector<int>  ind0, vector<int> exponent0);
    ~MonomialOfGrits();
    MySeriesTruncMod* evaluate(MySeriesTruncMod**qs);
    MySeriesTruncMod* evaluate(MySeriesTruncMod**qs, int upto);
    void printRepresentation();
};

#endif // __QZSeriesWH_H_INCLUDED__
