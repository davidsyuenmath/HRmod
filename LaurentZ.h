/***************************************
* This class encodes a Laurent polynomial as
* z^(minexp) times a FLINT polynomial p
* Feb 4, 2018.  Now with read/save capabilities
***************************************/
// include guard
#ifndef __LAURENTZ_H_INCLUDED__
#define __LAURENTZ_H_INCLUDED__

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

class LaurentZ {
  private:
    fmpz_t tmp;

  public:
    fmpz_poly_t p;
    int minexp;
    LaurentZ(fmpz_t c);
    LaurentZ();
    LaurentZ(int c);
    LaurentZ(int c, int minexp0);
    LaurentZ(vector<int> coeffs);
    LaurentZ(vector<int> coeffs, int minexp0);
    ~LaurentZ(); //destructor
    void init(vector<int> coeffs, int minexp0);
    void normalize();
    void shiftExponents(int shift);
    bool isZero();
    string getString(string var);
    string getString();
    LaurentZ* copy();
    static LaurentZ* multiply(LaurentZ* a, LaurentZ *b);
    void multiplyWith(LaurentZ* b);
    static LaurentZ* add(LaurentZ *a, LaurentZ* b);
    void addWith(LaurentZ* b);
    void subtractWith(LaurentZ* b);
    static LaurentZ* addScalarMultiple(LaurentZ *a, fmpz_t m, LaurentZ* b);
    void addScalarMultipleWith(fmpz_t m, LaurentZ* b);
    int divideWith(LaurentZ* b); //returns 1 if success and no remainder
    void setRemainder(fmpz_poly_t r, LaurentZ* b); //b must be "doubly monic".
    static LaurentZ* scalarMultiply(LaurentZ* b, fmpz_t m);
    void scalarMultiply(fmpz_t m);
    static void multiply(LaurentZ* res, LaurentZ* a, LaurentZ *b);
    int isUnit();
    static LaurentZ* BinomialZnMinusOne(int d);
    static LaurentZ* BinomialZnPlusZMinusn(int d);
    static LaurentZ* BTB(vector<int> d);
    static LaurentZ* BabyThetaBlock(vector<int> d);
    static LaurentZ* Germ(int weight, vector<int> d);
    bool divideByIntegerWithCheck(fmpz_t c); //returns 1 if divide exactly
    void divideByIntegerWithoutCheck(fmpz_t c); //assumes division will be exact
    void getCoeff(fmpz_t coeff, int n);
    static LaurentZ* expand(LaurentZ* b, int L);
    void expand(int L);
    static LaurentZ* shrink(LaurentZ* b, int L);
    void shrink(int L);
    void addToTerm(fmpz_t c, int n);
    int getMaxExponent();
    int equals(LaurentZ* b);
    int checkSymmetry();
    void negate();
    void saveToFile(FILE *f);
    LaurentZ(FILE *f);

};

#endif // __LAURENTZ_H_INCLUDED__
