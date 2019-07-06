/***************************************
* This class stores a list of (n,r,m) that dot with (a,b,c) to <= uptoDot
* and related info and does one restriction given "character" (BAa, BAb, BAc).
* 5/6/2018: Now does restriction of Grit and Tweaks
***************************************/
// include guard
#ifndef __MyRepnumList_H_INCLUDED__
#define __MyRepnumList_H_INCLUDED__

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
#include "SigObjectOldStyle.h"


using namespace std;
//=================================

class MyRepnumList {
  private:


  public:
    long long a, b, c, lev, mindet, delta, uptoDot;
    long long BPA,BPB,BPC;
    int maxq, maxqEff;
    fmpz_t singExp;
    vector<long long> n, r, m, reducedmn, reducedr, trList;
    vector<long long> nmodL, rmodL, mmodL, nEff, rEff, mEff, trEff, reducedmnEff, reducedrEff;
    long long L;
    int LrootPowersLen; fmpz_t* LrootPowers, *negPowers;
    fmpz_t *jfCoeffs; int jfCoeffsSize;

    MyRepnumList(int a0, int b0, int c0, int lev0, int uptoDot0, int mindet0);
    void initBasics(int a0, int b0, int c0, int lev0, int uptoDot0, int mindet0);
    ~MyRepnumList();
    static long long myfloor(double x);
    static long long myceiling(double x);
    static void myreduce(long long n, long long r, long long p, long long &newn, long long &newr);
    void prepEfficientSpecial(int L0, QZSeriesWH *jf, int Lroot, fmpz_t modn);
    void prepEfficientSpecialGeneralPowersLen(int L0, QZSeriesWH *jf, int Lroot, fmpz_t modnm, int powersLen);
    MySeriesTruncMod* makeRestriction();
    MySeriesTruncMod* oneBPCosetRestrictionPrepped(long long BAa, long long BAb, long long BAc, long long Ldenom, fmpz_t modn, int oneLroot, int uptoN);
    MySeriesTruncMod* oneBPCosetRestrictionPrepped(long long BAa, long long BAb, long long BAc, long long Ldenom, fmpz_t modn, fmpz_t oneLroot, int uptoN0);
    void fprintfDiagnostics(FILE* F, int uptoN);
    int hasDesiredDot(long long targetDot);

    MySeriesTruncMod* restrictGrit(QZSeriesWH* jf, long long BAa, long long BAb, long long BAc,
                                   long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoN, int weight);
    void printAbstractRestrictionToFile(FILE* OUTFILE, SigObjectOldStyle* sigObj,long long BAa, long long BAb, long long BAc,
                                   long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoN, int weight);

    MySeriesTruncMod* restrictGritNoCharacter(QZSeriesWH* jf,
                                    long long Ldenom, fmpz_t modn, int uptoN, int weight);
    MySeriesTruncMod* restrictGritTweak(QZSeriesWH* jf, long long BAa, long long BAb, long long BAc,
                                        long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoN, int weight);
    int getMinDotGrit(QZSeriesWH* jf, int weight);
};

#endif // __QZSeriesWH_H_INCLUDED__
