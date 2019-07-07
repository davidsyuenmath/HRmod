/***************************************
* This class
***************************************/
// include guard
#ifndef __RepnumGroup_H_INCLUDED__
#define __RepnumGroup_H_INCLUDED__

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
#include "MySeriesTruncMod.h"


using namespace std;
//=================================

class RepnumGroup {
  private:
    slong level;
    int findIndex(slong sa, slong sb, slong sc);
    int mindet;

  public:
    vector<slong> a, b, c;
    vector<MyRepnumList*> rpn;


    RepnumGroup(int level0, int mindet0);
    RepnumGroup(int level0);
    ~RepnumGroup();
    void clearGroup();
    MyRepnumList* getRepnumList(slong sa, slong sb, slong sc, int uptoDot);
    void makeRepnumList(slong sa, slong sb, slong sc, int uptoDot);
};

#endif // __QZSeriesWH_H_INCLUDED__
