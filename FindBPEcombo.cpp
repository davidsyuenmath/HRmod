#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <gmp.h>
#include <gmpxx.h>
#include <mpfr.h>
#include <flint.h>
#include <flintxx.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_mat.h>
#include <fmpz_vec.h>

#include "Diagnostics.h"
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "QZSeriesWH.h"
#include "AXBsolver.h"
#include "JacobiFormBasis.h"
#include "DRList.h"
//#include "BPdata.h"
//#include "FindBPunified.h"
#include "ThetaBlockList.h"
//#include "BPdataList.h"
#include "SigObjectOldStyle.h"
//#include "BPcheckCuspiness.h"
#include "SingularPartList.h"
//#include "WeightTwelveBases.h"
//#include "QZWSeriesQ.h"


//INCLUDING cpp files because I couldn't figure out how to compile them separately in IDE.
#include "Diagnostics.cpp"
#include "LaurentZ.cpp"
#include "QZSeriesZ.cpp"
#include "QZSeriesWH.cpp"
#include "AXBsolver.cpp"
#include "JacobiFormBasis.cpp"
#include "DRList.cpp"
//#include "BPdata.cpp"
//#include "FindBPunified.cpp"
#include "ThetaBlockList.cpp"
//#include "BPdataList.cpp"
#include "SigObjectOldStyle.cpp"
//#include "BPcheckCuspiness.cpp"
#include "SingularPartList.cpp"
//#include "WeightTwelveBases.cpp"
//#include "QZWSeriesQ.cpp"
//#include "MySeriesTruncMod.cpp"
//#include "MyRepnumList.cpp"


using namespace std;

  const int DIAGNOSE1=1; ///Whether to print DIAGNOSE1 messages PROGRESS MESSAGES.
  const int DIAGNOSE2=1; ///Whether to run test code.
  const int DIAGNOSE3=0; ///Whether to print DIAGNOSE3 messages.
  const int DIAGNOSE4=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE5=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE6=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE7=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE8=0; ///Whether to print DIAGNOSE4 messages.

  string infileName;
  int level, weight, uptoN;
  string saveFileName; stringstream saveFileNameStream; int saveFileQ;
  int maxqorder, maxd;
  string outfileRoot;
  int useFile;
  //string whjfFilename;
  char whjfStr[255];

  //CODE 1: sum/Delta^r
  //CODE 2: sum of BPE
  //CODE 3: (a*phi|V2 + c*thetaBuddy)/phi/b
  vector<int> TBphi, TBTheta; int tVal, ell;


int main(int argc, char* argv[])
{

    if(argc<2){
        cout << "Usage requires one or more arguments: parametersfilename,  [uptoBPEweight] [saveFileName].\n";
        cout<<"saveFileName output-parametersfilename-uptoBPEweight.ma is used by default.\n";
        cout<<"Info about this machine:\n";
        printf("sizeof(slong) = %d;\n", sizeof(slong));
        printf("sizeof(int) = %d;\n", sizeof(int));
        printf("sizeof(long long) = %d;\n", sizeof(long long));

    return 1;
    }
    int uptoBPEweight = 2;
    infileName = argv[1];
    if(argc>2){
     uptoBPEweight = atoi(argv[2]);
    }

    if(argc>3){
     saveFileNameStream<<argv[3];
     saveFileNameStream>>saveFileName;
     saveFileQ=1;
    }else{
     saveFileNameStream<<"output-"<<infileName<<"-"<<uptoBPEweight<<".ma";
     saveFileNameStream>>saveFileName;
     saveFileQ=1;
    }
    cout<<"Will save output to "<<saveFileName;

    ofstream outf;
    outf.open(saveFileName.c_str());
    if(outf.fail()){cout<<"File write open failure: "<<saveFileName<<"\n";exit(1);}



    FILE* PFILE, *WFILE;
    cout<<"Reading parameters from "<<infileName<<" ...\n";
    PFILE = fopen(infileName.c_str(),"r");
    if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
    int goalWeight, goalIndex, code;
    int x, i, j, len, ell;
    fscanf(PFILE, "%d", &goalIndex);
    uptoN=2+goalIndex/4;
    cout<<"goalIndex="<<goalIndex
        <<";uptoN="<<uptoN<<".\n";

    int numWJFs;
    fscanf(PFILE, "%d", &numWJFs);
    cout<<"Will read "<<numWJFs<<" WHJF files.\n";
    QZSeriesWH** whjfArray;
    whjfArray = new QZSeriesWH*[numWJFs];
    for(i=0;i<numWJFs;i++){
      fscanf(PFILE, "%s", whjfStr);
      WFILE = fopen(whjfStr,"r");
      if(!WFILE){cout<<" no such WHJF file.\n"; exit( 1); }
      cout<<"Reading Jacobi forms from "<<whjfStr<<".\n";
      whjfArray[i] = new QZSeriesWH(WFILE);
      fclose(WFILE);
    }

    JacobiFormBasis * wtzerojfs;
    wtzerojfs = new JacobiFormBasis(goalIndex, uptoN);
    stringstream sss; 
    ThetaBlockList *tbList; int k;
    for(k=2; k<=uptoBPEweight; k++){
      tbList=new ThetaBlockList(".", k, goalIndex, 12-k, 0);
      wtzerojfs->weight0AppendMorejfv2overjf(tbList);
      delete tbList;
    }
    fmpz_t denom;
    fmpz_init(denom);
    fmpz_t* combo;
    int comboLen = wtzerojfs->num;
    for(i=0;i<numWJFs;i++){
        cout<<"Working on JF # "<<i<<"\n";
        outf<<"entryNum="<<i<<";\n";

//NOTE JacobiFormBasis new function here. 20200601

        combo = wtzerojfs->matchSingularPartSaveToFile(whjfArray[i],denom, outf);
        if(combo!=NULL){
          cout<<"PROOF FOUND!\n";
          cout<<"proof = {";
          for(j=0;j<comboLen;j++){
              if(j>0){cout<<",";}
              cout<< fmpz_get_str(NULL,10,combo[j]);
          }
          cout<<"}/";
          cout<< fmpz_get_str(NULL,10,denom);
          cout<<";\n";
        }else{
          cout<<"no proof.\n";
        }
        outf<<"doSomething;\n\n";
    }

  return 0;
}


