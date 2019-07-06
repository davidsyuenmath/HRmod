/***************************************
* This class encodes a truncated series
* whose coefficients are Laurent polynomial LaurentZ
* sum_{i=0}^trunc (LaurentZ coefficients) * q^i + O(q^(trunc+1))
* 9/8/2017 now with sigma, E4, E6.
***************************************/
// include guard
#ifndef __QZSERIESZ_H_INCLUDED__
#define __QZSERIESZ_H_INCLUDED__

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
#include <fmpz_poly.h>
#include "LaurentZ.h"

using namespace std;
//=================================

class QZSeriesZ {
  private:
    void resizeArray(int newsize);
    static QZSeriesZ* OnePlusSkip(int trunc0, int skip);

  public:
    int coeffArraySize; //trunc+1 at the beginning
    LaurentZ** coeff;
    int trunc;
    QZSeriesZ(); //Nothing initialized
    QZSeriesZ(int trunc0); //initialize to zero series up to q^trunc0
    QZSeriesZ(int trunc0, vector<int> coeffs);
    QZSeriesZ(int trunc0, int topOne, int topCoeff); //1 + topCoeff*q^topOne
    static QZSeriesZ* EtaFunctionWithoutFractionalExponent(int trunc0);
    static QZSeriesZ* InverseEtaFunctionWithoutFractionalExponent(int trunc0);
    //static QZSeriesZ* ThetaFunctionWithoutFractionalExponentWithoutBTB(int trunc0, int d);
    static QZSeriesZ* E4(int trunc0);
    static QZSeriesZ* E6(int trunc0);
    ~QZSeriesZ(); //destructor
    void destroy();
    void shiftExponents(int shift); //this could resize cpeffArraySize
    bool isZero();
    string getString(string qvar, string zvar);
    string getString();
    QZSeriesZ* copy();
    static QZSeriesZ* multiply(QZSeriesZ* a, QZSeriesZ *b);
    void multiplyWith(QZSeriesZ* b);
    static QZSeriesZ* add(QZSeriesZ *a, QZSeriesZ* b);
    void addWith(QZSeriesZ* b);
    void subtractWith(QZSeriesZ* b);
    static QZSeriesZ* addScalarMultiple(QZSeriesZ *a, fmpz_t m, QZSeriesZ* b);
    void addScalarMultipleWith(fmpz_t m, QZSeriesZ* b);
    static QZSeriesZ* divideByMonic(QZSeriesZ *a, QZSeriesZ* b); //b must be monic
    //void divideByMonic(QZSeriesZ* b); //returns 1 if success and no remainder
    int divideByLaurentZ(LaurentZ* b);
    static QZSeriesZ* multiplyByLaurentZ(QZSeriesZ *a, LaurentZ* b); //b must be monic
    void multiplyByLaurentZ(LaurentZ* b);
    bool divideByIntegerWithCheck(fmpz_t c); //returns 1 if divide exactly
    void divideByIntegerWithoutCheck(fmpz_t c); //assumes division will be exact
    void negate();
    void multiplyByScalar(fmpz_t m);
    static void sigma(fmpz_t result, int k, int n);
    void saveToFile(FILE *f);
    QZSeriesZ(FILE *f);

};

#endif // __QZSERIESZ_H_INCLUDED__
