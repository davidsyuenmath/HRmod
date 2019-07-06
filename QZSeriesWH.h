/***************************************
* This class encodes a truncated weakly holomorphic Jacobi form series
* q^(frac24/24) * q^(qord) * QZSeriesZ
* 12/27/2017: added methods getA(), getB(), getC() associated to weight 0 forms
* 5/6/2018: Now with methods for getting Grit and GritTweak coefficients.
***************************************/
// include guard
#ifndef __QZSeriesWH_H_INCLUDED__
#define __QZSeriesWH_H_INCLUDED__

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
#include <fmpz_mat.h>
#include <fmpz_poly.h>
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "DRList.h"
#include "ThetaBlockList.h"

using namespace std;
//=================================

class QZSeriesWH {
  private:
    void normalizePartially();
    static int signum(int x);

  public:
    static fmpz_mat_t diagMat;
    static int printGG;
    static int etaPowerArraySize, thetaArraySize, E4PowerArraySize, E6PowerArraySize;
    static int *thetaPowerArraySize;
    static QZSeriesWH** etaPowers, **E4Powers, **E6Powers;
    static int* etaTruncations, *E4Truncations, *E6Truncations;
    static QZSeriesWH*** thetaPowers;
    static int ** thetaTruncations;
    static void initPowers();
    static QZSeriesWH* ThetaFunctionDPower(int trunc0, int d, int pow);
    static QZSeriesWH* ThetaBlockEfficient(int trunc0, int weight, vector<int> d);
    static QZSeriesWH* ThetaBlockEfficient(int trunc0, int weight, vector<int> d, int allowTweak);
    static QZSeriesWH* ThetaBlockInefficient(int trunc0, int weight, vector<int> d);

    int qExpAdjustOver24, zExpAdjustOver2, qMinExp;
    QZSeriesZ *qzs;
    QZSeriesWH(); //does not initialize qzs!
    QZSeriesWH(int trunc0); //initialize to zero series up to q^trunc0
    QZSeriesWH(int trunc0, int initialCoeff); //initialize to zero series up to q^trunc0
    QZSeriesWH(int trunc0, vector<int> coeffs);
    QZSeriesWH(QZSeriesWH* a);
    static QZSeriesWH* EtaFunction(int trunc0);
    static QZSeriesWH* InverseEtaFunction(int trunc0);
    static QZSeriesWH* EtaFunctionPower(int trunc0, int pow);
    static QZSeriesWH* ThetaFunctionD(int trunc0, int d);
    ~QZSeriesWH(); //destructor
    void destroy();
    bool isZero();
    void normalize();
    void shiftExponents(int shift);
    string getStringNoNormalize(string qvar, string zvar);
    string getStringNoNormalize();
    string getString(string qvar, string zvar);
    string getString();
    QZSeriesWH* copy();
    static QZSeriesWH* multiply(QZSeriesWH* a, QZSeriesWH *b);
    void multiplyWith(QZSeriesWH* b);
    static QZSeriesWH* add(QZSeriesWH *a, QZSeriesWH* b);
    void addWith(QZSeriesWH* b);
    void subtractWith(QZSeriesWH* b);
    static QZSeriesWH* addScalarMultiple(QZSeriesWH *a, fmpz_t m, QZSeriesWH* b);
    void addScalarMultipleWith(fmpz_t m, QZSeriesWH* b);
    static QZSeriesWH* divideByMonic(QZSeriesWH *a, QZSeriesWH* b); //b must be monic
    void divideByMonic( QZSeriesWH* b); //b must be monic
    int divideByLaurentZ(LaurentZ* b);
    static QZSeriesWH* multiplyByLaurentZ(QZSeriesWH *a, LaurentZ* b); //b must be monic
    void multiplyByLaurentZ(LaurentZ* b);
    void negate();


    int getMaxExponent();
    int getMinExponent();
    int divideBy(QZSeriesWH* b); //returns 1 if successful
    static QZSeriesWH* linearCombo(QZSeriesWH** qq, fmpz_t* c, int len);
    int divideByIntegerWithCheck(fmpz_t c); //returns 1 if divide exactly
    void divideByIntegerWithoutCheck(fmpz_t c); //assumes division will be exact
    static QZSeriesWH* rationalLinearCombo(QZSeriesWH** qq, fmpq_t* c, int len);
    static QZSeriesWH* ThetaBlockWithOps(int trunc0, int weight, vector<int> d, vector<int> heckeOps);
    static QZSeriesWH* ThetaBlockWithOpsjfFormat(int targetIndex, int trunc0, int weight, vector<int> d, vector<int> heckeOps);
    static QZSeriesWH* ThetaBlockProductjfFormat(int trunc0, vector<int> d);
    static QZSeriesWH* ThetaBlock(int trunc0, int weight, vector<int> d);
    static QZSeriesWH* ThetaBlockTweak(int trunc0, int weight, vector<int> d);
    static QZSeriesWH* E4(int trunc0);
    static QZSeriesWH* E4Power(int trunc0, int pow);
    static QZSeriesWH* E6(int trunc0);
    static QZSeriesWH* E6Power(int trunc0, int pow);
    static QZSeriesWH* ThetaBlockWithOps(int trunc0, int weight, vector<int> b);
    void applyUp(int index,int weight,int L);
    static int myHeckeDownMaxq(int targetm, int n,int h);
    static int myHeckeDownMaxqInverse(int targetm, int n,int h);
    void applyDown(int index, int weight, int L);
    void applyDoubleUp(int L);
    void getFC_noreduce(fmpz_t fc, slong n, slong r);
    void getFC(fmpz_t fc, slong n, slong r, slong index);
    void getFC_Dr(fmpz_t fc, slong D, slong r, slong index);
    static long long mysqrtfloor(long long x);
    LaurentZ* getCoeffCopy(int n);
    int getCoeffMaxr(int n);
    void saveToFile(ofstream &f, string varname);
    int checkSymmetry();
    void putIntoDRList(DRList* dr, int index, int uptoN);
    void getHumbertMultiplicity(fmpz_t h, int minD, int d, int r, int index);
    int hasNonzeroNewDR(DRList* dr, int index, int uptoN);
    void makeHumbertVector(fmpz_mat_t M, DRList* dr, int index, int colIndex);
    void makeHumbertVector(fmpz_mat_t M, DRList* dr, int index);
    void truncate(int newtrunc);
    void makeSingularVector(fmpz_mat_t M, DRList* dr, int index, int colIndex);
    void makeSingularVector(fmpz_mat_t M, DRList* dr, int index);
    double getOrd(int index);
    int getSymmetry(int index);
    void multiplyByScalar(fmpz_t m);

    void getABC(int&A, int&B, int&C);
    void saveToFile(FILE *f);
    QZSeriesWH(FILE *f);
    //static QZSeriesWH* makeTBEcombo(ThetaBlockList *tbList, fmpz_t* combo, fmpz_t denom, int maxdegree);;
    static QZSeriesWH* makeTBEcombo(ThetaBlockList *tbList, fmpz_t* combo, fmpz_t denom, int maxdegree);;

    static vector<slong> getDivisors(slong a);
    static slong GCDL(slong a, slong b);
    void getFCGrit(fmpz_t fc, slong n, slong r, slong m, slong weight, slong level);
    void getFCGritTweak(fmpz_t fc, slong nx2, slong rx2, slong mx2, slong weight, slong level);
    void getFCx2(fmpz_t fc, slong nx2, slong rx2, slong index);

};

#endif // __QZSeriesWH_H_INCLUDED__
