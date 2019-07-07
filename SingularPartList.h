/***************************************
* Efficient making a Humbert matrix without saving JFs
***************************************/
// include guard
#ifndef __SingularPartList_H_INCLUDED__
#define __SingularPartList_H_INCLUDED__

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

#include <flint.h>
#include <flintxx.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_mat.h>
#include <fmpz_vec.h>
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "QZSeriesWH.h"
#include "ThetaBlockList.h"
#include "Diagnostics.h"
#include "DRList.h"



using namespace std;
//=================================

class SingularPartList {
  private:
    DRList * dr;

  public:
    vector<int> rowLen;
    fmpz_t** mat;
    int index;
    vector<int>* tbArray;
    int matSize, maxRowLen;
    fmpz_mat_t fmat; int fmatMade;

    SingularPartList(int index0);
    ~SingularPartList(); //destructor
    int getNumRows();
    void appendRow(QZSeriesWH* jf, int index, int uptoN);
    void appendRow(QZSeriesWH* jf, int index, int uptoN, vector<int> tb);
    void appendRows(ThetaBlockList* tbList, int index, int uptoN);
    void appendRows(ThetaBlockList* tbList, ThetaBlockList *tbList2, int index, int uptoN);
    void increaseMatSize();
    void makeMat(fmpz_mat_t &m, int beginRowNumber, int numRows);
    void checkfmat();
    void clearfmat();
    fmpz_t* getVec(int rowNumber);
    void saveToFile(ofstream&f, string varname);
    fmpz_t* matchSingularPart(QZSeriesWH *jf, fmpz_t denom);
};

#endif // __DRList_H_INCLUDED__
