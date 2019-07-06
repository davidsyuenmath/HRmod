/***************************************
* This class encodes a Laurent polynomial as
* z^(minexp) times a FLINT polynomial p
***************************************/
// include guard
#ifndef __AXBsolver_H_INCLUDED__
#define __AXBsolver_H_INCLUDED__

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


using namespace std;
//=================================

class AXBsolver {
  private:
    fmpz_mat_t A;
    fmpz_t * b;

  public:
    int nrows, ncols;
    long long** sol;
    fmpz_mat_t solM;

    int numSolutions;
    AXBsolver(fmpz_mat_t A0, fmpz_t* b0, int nrows0, int ncols0);
    AXBsolver(fmpz_mat_t A0, fmpz_mat_t b0, int nrows0, int ncols0);
    AXBsolver(int** A0, int* b0, int nrows0, int ncols0);
    ~AXBsolver(); //destructor
    void destroy();
    void solve();
};

#endif // __AXBsolver_H_INCLUDED__
