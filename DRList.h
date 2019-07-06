/***************************************
* This class encodes a Laurent polynomial as
* z^(minexp) times a FLINT polynomial p
***************************************/
// include guard
#ifndef __DRList_H_INCLUDED__
#define __DRList_H_INCLUDED__

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

#include "Diagnostics.h"

using namespace std;
//=================================

class DRList {
  private:
    int minDcalculated;

  public:
    vector<int> DList, RList;
    int minD;

    DRList();
    ~DRList(); //destructor
    void clear();
    int alreadyOnList(int d, int r);
    void insert(int d, int r);
    void insert(vector<int> d, vector<int> r);
    vector<int> getDList();
    vector<int> getRList();
    int getLength();
    int getMinD();
    void saveToFile(ofstream&f, string varname);
    void print();

};

#endif // __DRList_H_INCLUDED__
