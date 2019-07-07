/***************************************
* Simple struct containing specifications of a JacobiForm including tweaks in future
***************************************/
// include guard
#ifndef __JacobiFormData_H_INCLUDED__
#define __JacobiFormData_H_INCLUDED__

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
#include "QZSeriesWH.h"


using namespace std;
//=================================

class JacobiFormData {
  private:
    int jfType; //0= Theta block with no character, possible with Hecke ops
    vector<int> d;
    int tee; //power of eta
    int weight;
    int index;

  public:
    JacobiFormData(
        int tee0, vector<int> d0 
    );
    JacobiFormData(
        int jfType0, int weight0, int index0, int tee0, vector<int> d0 
    );
    ~JacobiFormData();
    QZSeriesWH* makeJF(int uptoN);
};

#endif // __QZSeriesWH_H_INCLUDED__
