/***************************************
* This class encodes a set of weakly holomorphic Jacobi forms
* with routines for finding nullspace, basis over integers, etc.
* 9/7/2017: add capability to read a JBasis-wt-index.jf file
* 1/9/2018: finally programmed general myLinearSolveOverIntegers
***************************************/
// include guard
#ifndef __JacobiFormBasis_H_INCLUDED__
#define __JacobiFormBasis_H_INCLUDED__

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
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "QZSeriesWH.h"
#include "AXBsolver.h"
#include "ThetaBlockList.h"
#include "Diagnostics.h"

using namespace std;
//=================================

class JacobiFormBasis {
  private:

  public:
    QZSeriesWH** jfArray;
    vector<int>* tbArray;
    int num, maxDegree, weight, index;
    vector<int> nList, maxrList, nsingList,rsingList, DList, rList;
    int maxrMade, DrMade;

    JacobiFormBasis(); //does nothing other than to make num=0
    JacobiFormBasis(int index0, int trunc); //weight 0 initialization
    JacobiFormBasis(int weight, int index, int trunc);
    JacobiFormBasis(int weight, int index, int trunc, fmpz_t* combo);
    JacobiFormBasis(string directory, int weight, int index, int trunc);
    JacobiFormBasis(string directory, int weight, int index, int trunc, fmpz_t* combo);
    ~JacobiFormBasis(); //destructor
    void destroy();
    static QZSeriesWH** readJBasisFile(int weight, int index, int trunc, int &num);
    static QZSeriesWH** readJBasisFile(string directory, int weight, int index, int trunc, int &num);
    static QZSeriesWH** readjfJBasisFile(string directory, int weight, int index, int trunc, int &num);
    static QZSeriesWH** readJBasisFile(int weight, int index, int trunc, int &num, fmpz_t* combo);
    static QZSeriesWH** readJBasisFile(string directory, int weight, int index, int trunc, int &num, fmpz_t* combo);
    static QZSeriesWH** readjfJBasisFile(string directory, int weight, int index, int trunc, int &num, fmpz_t* combo);
    void weight0AppendMorejfv2overjf(ThetaBlockList* tbList);
    void divideBy(QZSeriesWH* jf, int weight0, int index0);
    static void mySmithNormalForm(fmpz_mat_t U, fmpz_mat_t A, fmpz_mat_t V, fmpz_mat_t M); //returns M=UAV
    static void myInvert(fmpz_mat_t A);
    static void mySetIdentityMatrix(fmpz_mat_t A);
    static int myDiagonalForm(fmpz_mat_t A);
    void makeMaxrList();
    void makeDrList();
    QZSeriesWH* makeLinearCombo(fmpz_mat_t M, int row, fmpz_t denom);
    QZSeriesWH* makeLinearCombo(fmpz_mat_t M, int row);
    QZSeriesWH* makeLinearComboColumn(fmpz_mat_t M, int col);
    QZSeriesWH* makeLinearCombo(fmpz_t* combo, fmpz_t denom);
    QZSeriesWH* makeLinearCombo(fmpz_t* combo);
    void saveToFile(ofstream &f, string varname);
    static void saveMatToFile(ofstream &f, fmpz_mat_t m, string varname);
    static void printMat(fmpz_mat_t m, string varname);
    int getMinExponent();
    JacobiFormBasis* getWeightZeroIntegerBasis();
    int getMaxr(int n);
    JacobiFormBasis* getNullSpace(int targetqExponent);
    JacobiFormBasis* getNullSpaceUpto(int targetqExponent);
    JacobiFormBasis* getNullSpace(int targetqExponentFrom, int targetqExponentTo);
    JacobiFormBasis* getNullSpaceAll();

    static fmpz_t* myLinearSolve(fmpz_mat_t M, fmpz_t* b, fmpz_t denom);
    static fmpz_t* myLinearSolveOverIntegersOLD(fmpz_mat_t M, fmpz_t* b);
    static fmpz_t* myLinearSolveOverIntegers(fmpz_mat_t M, fmpz_t* b);
    static fmpz_t* myLinearSolve(fmpz_mat_t M, fmpz_mat_t b, fmpz_t denom);


    QZSeriesWH* getOneSolution(LaurentZ* targetLaurentPoly, int targetqExponent); //Returns NULL if no solution
    int checkSymmetry();
    void putIntoDRList(DRList* dr, int uptoN);
    void makeHumbertMatrix(fmpz_mat_t M, DRList* dr);
    JacobiFormBasis* getGoodJFs(fmpz_mat_t HumMatrix,fmpz_mat_t HumVec, QZSeriesWH* combinedSol);
    int allNonnegative(fmpz_mat_t M);

    JacobiFormBasis* subspaceMultipleOf(LaurentZ* pz);
    fmpz_t* matchSingularPart(QZSeriesWH* jfweak,  fmpz_t denom);
    void makeSingularMatrix(fmpz_mat_t M, DRList* dr);

    static int myAttemptDivide(fmpz_t x, fmpz_t a, fmpz_t b);  //Solves bx=a (x=a/b when b!=0)

};

#endif // __JacobiFormBasis_H_INCLUDED__
