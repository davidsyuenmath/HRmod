/***************************************
* This class encodes a truncated series
* whose coefficients are Laurent polynomial LaurentZ
* sum_{i=0}^trunc (LaurentZ coefficients) * q^i + O(q^(trunc+1))
***************************************/
// include guard
#ifndef __QZWSeriesQ_H_INCLUDED__
#define __QZWSeriesQ_H_INCLUDED__

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
#include "QZSeriesWH.h"


using namespace std;
//=================================

class QZWSeriesQ {
  private:
    void resizeArray(int newsize);

  public:
    int coeffArraySize; //trunc+1 at the beginning
    fmpz_t denom;
    QZSeriesWH** coeff;
    int trunc, defaultQTrunc;
    QZWSeriesQ(); //Nothing initialized
    QZWSeriesQ(int wTrunc, int qTrunc); //initialize to zero series up to w^wTrunc0 with q^qTrunc
    ~QZWSeriesQ(); //destructor
    void destroy();
    void shiftExponents(int shift); //this could resize cpeffArraySize
    string getString(string qvar, string zvar, string wvar);
    string getString();
    QZWSeriesQ* copy();
    static QZWSeriesQ* multiply(QZWSeriesQ* a, QZWSeriesQ *b);
    void multiplyWith(QZWSeriesQ* b);
    static QZWSeriesQ* add(QZWSeriesQ *a, QZWSeriesQ* b);
    void addWith(QZWSeriesQ* b);
    void subtractWith(QZWSeriesQ* b);
    static QZWSeriesQ* addScalarMultiple(QZWSeriesQ *a, fmpz_t m, QZWSeriesQ* b);
    void addScalarMultipleWith(fmpz_t m, QZWSeriesQ* b);
    void multiplyByQZSeriesWH(QZSeriesWH* b);
    bool divideByDenominator(); //returns 1 if divided exactly
    void divideByInteger(fmpz_t c);
    void negate();
    void multiplyByScalar(fmpz_t m);
    static QZWSeriesQ* makeBorchLift(QZSeriesWH* psi, QZSeriesWH* tb, int level, int leadingJC, int numJC);
    static QZWSeriesQ* makeGritLift(QZSeriesWH* phi, int level, int weight, int numJC);
    QZWSeriesQ* exp(int terms);
    int getFC(fmpz_t fc, int n, int r, int m, int level);
    void saveToFile(ofstream &f, string varname);



};

#endif // __QZWSeriesQ_H_INCLUDED__
