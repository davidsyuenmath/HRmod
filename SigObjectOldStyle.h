/***************************************
* SigObjectOldStyle encodes paramodular matrix indices in a standard order
***************************************/
// include guard
#ifndef __SigObjectOldStyle_H_INCLUDED__
#define __SigObjectOldStyle_H_INCLUDED__

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

using namespace std;
//=================================

class SigObjectOldStyle
{
  private:
    long long p, p2;
    int* Ninv, *Ngcd;
    int* u; int* w; int* z;
    vector<int*> vset;
    int* divisorToIndex;
    vector<int*> vsetv2;
    vector<int> vsetv1;
    int verbose;
    void printOUTarray(int* v, int n);
    static int MyGCD0(int a, int b);
    void MyExtendedGCD(int a, int b, int &c, int &d);
    void initRest();
    void initFirst();
    long long QFormEval(long long a, long long b, long long c,
     long long x, long long y);

    int myfloor(double x);
    long long myDet3(int a, int b, int c);

    int myceiling(double x);

    void mycopy(int a[][2], int b[][2]);
    void swap(int &x, int &y);
    long long mysqrtfloor(long long  x);
    int alreadyInList(int** x, int n, int *y, int len);

    int mycomp(int*x, int*y, int n);
    int eqvecQ(int*x, int*y, int n);
    void initSigObjectOldStyle(int p0, int d0);

  public:
    int uptoDet;
    int SIGLEN, MAXDET;
    int** sigs;  //starts at 1 !!
    int** startend;
    SigObjectOldStyle(int p, int d, string filename);
    SigObjectOldStyle(int p, int d0);
    SigObjectOldStyle(int p0, int d0, int verbose0);
    void destroy();
    void P1Reduce(int* v);
    void ReduceUPart(int a, int b, int c, int* v);
    void PX2Reduce(int a, int b, int c, int x, int y,  int* v5);
    void PX2ReduceLong(long long a, long long b, long long c, int x, int y,  int* v5);
    int getIndexFromSig(int *v);
    int getIndexFromSig(long long a, long long b, long long c);
    static int MyGCD(int a, int b);
    int quickDet3(int *a);
    int getIndexFromDet(int d);
    int getDetFromIndex(int in);
    int twinIndex(int in);
    void doStandardUNRnoverify(int A, int B, int C, int v1, int v2, long long* unr);
    void mult22L(long long a[][2], long long b[][2], long long c[][2]);
    void setMat22L(long long a[][2], int x, int y, int z, int w);
    void setTranspose22L(long long b[][2], long long a[][2]);

    void ReduceUPart(int a, int b, int c, int* v, long long &pm);
    void PX2ReduceLong(long long a, long long b, long long c, int x, int y,  int* v5, long long &pm);
    int getIndexFromSig(long long a, long long b, long long c, long long &pm);
    int twinIndex(int in, long long &pm);


};   //This semicolon is essential!


#endif // __SigObjectOldStyle_H_INCLUDED__
