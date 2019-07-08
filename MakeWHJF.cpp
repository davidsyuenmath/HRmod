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
//#include "SingularPartList.cpp"
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
        cout << "Usage requires one or more arguments: filename,  [uptoN], [saveFileName].\n";
        cout<<"uptoN by default is read from the file. Use 0 for default.\n";
        cout<<"saveFileName from file is used by default.\n";
        cout<<"Info about this machine:\n";
        printf("sizeof(slong) = %lu;\n", sizeof(slong));
        printf("sizeof(int) = %lu;\n", sizeof(int));
        printf("sizeof(long long) = %lu;\n", sizeof(long long));

    return 1;
    }
    infileName = argv[1];
    uptoN=0;
    if(argc>2){
     uptoN = atoi(argv[2]);
    }
    if(argc>3){
     saveFileNameStream<<argv[6];
     saveFileNameStream>>saveFileName;
     saveFileQ=1;
    }


    FILE* PFILE;
    cout<<"Reading parameters from "<<infileName<<" ...\n";
    PFILE = fopen(infileName.c_str(),"r");
    if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
    int goalWeight, goalIndex, code;
    int x, i, len, ell;
    fscanf(PFILE, "%s", whjfStr);
    if(!saveFileQ){
        stringstream ssff;
        ssff<<whjfStr;
        ssff>>outfileRoot;
    }
    fscanf(PFILE, "%d", &goalWeight);
    fscanf(PFILE, "%d", &goalIndex);
    fscanf(PFILE, "%d", &x); if(uptoN==0){uptoN=x;}
    fscanf(PFILE, "%d", &code);
    cout<<"goalWeight="<<goalWeight<<";goalIndex="<<goalIndex
        <<";uptoN="<<uptoN<<"; code="<<code<<".\n";
    TBphi.clear(); TBTheta.clear();
    QZSeriesWH* jf;
    int r, k, N, num;

    if((code==1)||(code==2)){
        fscanf(PFILE, "%d", &r);
        fscanf(PFILE, "%d", &k);
        fscanf(PFILE, "%d", &N);
        fscanf(PFILE, "%d", &num);
        fmpz_t* combo = new fmpz_t[num];
        fmpz_t denom;
        fmpz_init(denom);
        fmpz_fread(PFILE, denom);
        for(i=0;i<num;i++){
            fmpz_init(combo[i]);
            fmpz_fread(PFILE, combo[i]);
        }
        cout<<"Parameters read. Now making weight "<<k<<", index "<<N<<" Jacobi cusp form basis...\n";
        stringstream ss;
        ss<<"JBasis-wt"<<k;
        string dir;
        ss>>dir;
        JacobiFormBasis* jfb = new JacobiFormBasis(dir,k,N,uptoN+r,combo);
        if(jfb->num!=num){
            cout<<"JFBasis has a different num than num in parameters file.\n";
        }

        cout<<"Basis made. Now making rational linear combination.\n";
        jf = jfb->makeLinearCombo(combo, denom);
        if(r>0){
            cout<<"Combination made. Now dividing by Delta^"<<r<<".\n";
            QZSeriesWH* deltaPower = QZSeriesWH::EtaFunctionPower(uptoN+r, 24*r);
            cout<<"Delta power made, now dividing.\n";
            jf->divideBy(deltaPower);
            cout<<"Done dividing.\n";
            delete deltaPower;
        }
        delete jfb;
    }

    if(code==2){ //Assumes BPE has nu=1 or 2.
        int coeff;
        fscanf(PFILE, "%d", &k);
        fscanf(PFILE, "%d", &ell);
        TBTheta.clear();
        for(i=0;i<ell;i++){
            fscanf(PFILE, "%d", &x);
            TBTheta.push_back(x);
        }
        fscanf(PFILE, "%d", &coeff);
        QZSeriesWH *jfpsi, *jfphi;
        int jfTrunc=uptoN+2;
        if(DIAGNOSE1)cout<<"Making TBphi (with target maxqorder = "<<uptoN<<")...";
        jfphi = QZSeriesWH::ThetaBlock(jfTrunc*2, k, TBTheta);
        jfpsi = jfphi->copy();
        if(DIAGNOSE1)cout<<"Applying V2...";
        jfpsi->applyUp(goalIndex, k, 2);
        if(DIAGNOSE1)cout<<"Divide by TBphi...\n";
        jfpsi->divideBy(jfphi);
        delete jfphi;
        if(coeff==1){
            jf->addWith(jfpsi);
        }else if(coeff==-1){
            jf->subtractWith(jfpsi);
        }else{
            cout<<"CODE 2 with BPE coeff ("<<coeff<<") not 1 or -1 not supported. Abort.\n"; exit(1);
        }
        delete jfpsi;
    }


    if(code==3){
        fscanf(PFILE, "%d", &k);
        fscanf(PFILE, "%d", &N);
        fscanf(PFILE, "%d", &ell);
        fscanf(PFILE, "%d", &num);
        fmpz_t* combo = new fmpz_t[num];
        fmpz_t denom;
        fmpz_init(denom);
        fmpz_fread(PFILE, denom);
        for(i=0;i<num;i++){
            fmpz_init(combo[i]);
            fmpz_fread(PFILE, combo[i]);
        }
        cout<<"Parameters read. Now reading weight "<<k<<", index "<<N<<" theta blocks, making BPE combos...\n";
        stringstream ss7; string dir7;
        ss7<<"AllTB-wt"<<k;   //for now, let's just do TB without denom
        ss7>>dir7;
        ThetaBlockList* tbList=new ThetaBlockList(dir7, k, N, ell, 0);
        jf=QZSeriesWH::makeTBEcombo(tbList, combo, denom, uptoN);
    }


    if(code==4){

    }
    fclose(PFILE);



            //fprintf(OUTFILE, "L=%d;\n", L);

    FILE* OUTFILE;
    string f;
    stringstream ffs;
    if(saveFileQ){
        f=saveFileName;
    }else{
        if(1){
            ffs<<outfileRoot<<"-"<<uptoN<<".txt";
        }
        //ffs<<outfileRoot<<"-"<<modn<<"-"<<Lroot<<"-"<<beginHecke<<"-"<<endHecke<<"-"<<uptoN<<".ma";
        ffs>>f;
    }
    cout<<"Output file will be "<<f<<"\n";
    OUTFILE = fopen(f.c_str(), "w");
    if(OUTFILE==0){cout<<"File write open failure: "<<f<<"\n";exit(1);}

    jf->saveToFile(OUTFILE);
    fclose(OUTFILE);

    return 0;
}


