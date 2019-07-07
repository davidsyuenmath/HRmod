/***************************************
* An efficient way to evaluate a multivariable polynomial in MySeriesTruncMod
***************************************/
// include guard
#ifndef __EfficientMultiPoly_H_INCLUDED__
#define __EfficientMultiPoly_H_INCLUDED__

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

class EfficientMultiPoly {
  private:
    int maxTBnum;
    int oneExisting; //Any one TB index in use.
    fmpz_t* coeffs; vector<int>* ind; int indSize;
    int doRecurse;
    EfficientMultiPoly *polyA, *polyB; int splittingIndex;
    int contains(vector<int>v, int a);
    vector<int> copyWithDeleteOne(vector<int>v, int a);

  public:
    EfficientMultiPoly(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0
    );
    EfficientMultiPoly(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0, int oneExisting0
    );
    void init(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0, int oneExisting0
    );
    ~EfficientMultiPoly();
    MySeriesTruncMod* evaluate(MySeriesTruncMod**qs);
    void makeEfficient();
    MySeriesTruncMod* evaluateNoFrills(MySeriesTruncMod**qs);
    void printRepresentation();
};

#endif // __QZSeriesWH_H_INCLUDED__
