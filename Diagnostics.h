/***************************************
* This class encodes a Laurent polynomial as
* z^(minexp) times a FLINT polynomial p
***************************************/
// include guard
#ifndef __Diagnostics_H_INCLUDED__
#define __Diagnostics_H_INCLUDED__

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


using namespace std;
//=================================

class Diagnostics {
  private:


  public:
    static ofstream f;
    static int checkSymmetry;
    static int saveDiagnosticFile;
    static int debugMode;
    static int skipMathematica;
    static int useNegativeSingularOnly;
    static fmpz_t HumVecEntry00;
    static void init();
    static void init(string fname);
    static void printVector(vector<int> c);
    static string vectorToString(vector<int> c);

};

#endif // __Diagnostics_H_INCLUDED__
