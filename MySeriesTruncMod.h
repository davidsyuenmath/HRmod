/***************************************
* This class encodes a truncated q series mod modn + O(q)^trunc
* That is, the top exponent is q^(trunc-1) = q^(upto)
* 4/27/2018 Now allows q^(integer/2) powers with int inSqrtQ
* Now has static int checkSameTrunc which when 0 just truncates both operands to smaller trunc
***************************************/
// include guard
#ifndef __MySeriesTruncMod_H_INCLUDED__
#define __MySeriesTruncMod_H_INCLUDED__

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


using namespace std;
//=================================

class MySeriesTruncMod {
  private:

  public:
    fmpz_mod_poly_t p;
    fmpz_t modn;
    slong trunc;
    static int ulongMaxSetUp;
    static fmpz_t ulongmax;
    static int checkSameTrunc;
    int inSqrtQ;

    MySeriesTruncMod(fmpz_t modn0, int upto);
    MySeriesTruncMod(int modn0, int upto);
    MySeriesTruncMod(fmpz_t modn0, int upto, int inSqrtQ);
    MySeriesTruncMod(int modn0, int upto, int intSqrtQ);
    MySeriesTruncMod(MySeriesTruncMod*f);
    static void checkOperands(MySeriesTruncMod *f, MySeriesTruncMod *g);
    static void checkSameinSqrtQ(MySeriesTruncMod *f, MySeriesTruncMod *g);
    static MySeriesTruncMod* monicBinomial(fmpz_t modn0, int upto, fmpz_t oneCoeff, int e);  //Initializes to 1 + c*x^e
    int setMonicBinomialBPSpecial(fmpz_t oneCoeff, int e);
    static MySeriesTruncMod* monomial(fmpz_t modn0, int upto, fmpz_t oneCoeff, int e);  //Initializes to  c*x^e

    ~MySeriesTruncMod();
    static MySeriesTruncMod* add(MySeriesTruncMod* f, MySeriesTruncMod* g);
    void addBy(MySeriesTruncMod* g);
    void addAndAdjustTrunc(MySeriesTruncMod* g);
    static MySeriesTruncMod* scalarMultiplity(MySeriesTruncMod* g, fmpz_t a);
    void scalarMultiplityBy(fmpz_t a);
    static MySeriesTruncMod* multiply(MySeriesTruncMod* f, MySeriesTruncMod* g);
    void multiplyBy(MySeriesTruncMod* g);
    static MySeriesTruncMod* pow(MySeriesTruncMod* f, slong e);
    static MySeriesTruncMod* pow(MySeriesTruncMod* f, fmpz_t e);
    void multiplyByPower(MySeriesTruncMod* f, slong e);
    void multiplyByPower(MySeriesTruncMod* f, fmpz_t e);
    static void setUp();

    void multiplyByPowerOfQandChangeTrunc(int e);
    MySeriesTruncMod*  expandExpBy(int factor, int newupto);
    MySeriesTruncMod*  contractExpBy(int factor);

    void getCoeff(fmpz_t c, int e);
    void setCoeff(fmpz_t c, int e);
    void setCoeff(ulong c, int e);
    void updateCoeffAdd(fmpz_t c, int e);
    void printstdout();
    void printstdoutWithTruncationOrder();
    void printstdout(const char *x);
    void printFile(FILE *outf);
    void printFile(FILE *outf, const char *x);

    void clearSqrts(); //convert inSqrtQ=1 to 0.
    static MySeriesTruncMod* divide(MySeriesTruncMod* f, MySeriesTruncMod* g);
    void divideBy(MySeriesTruncMod* g);
    MySeriesTruncMod* makeCopy();
    int getVanishingOrder();
    void shift(int e);
};

int MySeriesTruncMod::ulongMaxSetUp=0;
fmpz_t MySeriesTruncMod::ulongmax;

#endif // __QZSeriesWH_H_INCLUDED__
