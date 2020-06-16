/*****
 * 20190120 improvement of 
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke/(sa*sc-sb*sb); //20190120
 *
 * ****/
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
#include <mpc.h>
#include <flint.h>
#include <flintxx.h>
#include <fmpz.h>
#include <fmpz_poly.h>
#include <fmpz_mat.h>
#include <fmpz_vec.h>
#include <fmpz_mod_poly.h>

#include "Diagnostics.h"
#include "LaurentZ.h"
#include "QZSeriesZ.h"
#include "QZSeriesWH.h"
#include "AXBsolver.h"
#include "JacobiFormBasis.h"
#include "DRList.h"
#include "SigObjectOldStyle.h"
//#include "SingularPartList.h"
#include "QZWSeriesQ.h"
#include "MySeriesTruncMod.h"
#include "MyRepnumList.h"
#include "RepnumGroup.h"
#include "EfficientQuadratic.h"
#include "EfficientMultiPoly.h"
#include "MonomialOfGrits.h"
#include "JacobiFormData.h"

//INCLUDING cpp files because I couldn't figure out how to compile them separately in IDE.
#include "Diagnostics.cpp"
#include "LaurentZ.cpp"
#include "QZSeriesZ.cpp"
#include "QZSeriesWH.cpp"
#include "AXBsolver.cpp"
#include "JacobiFormBasis.cpp"
#include "DRList.cpp"
#include "SigObjectOldStyle.cpp"
#include "QZWSeriesQ.cpp"
#include "MySeriesTruncMod.cpp"
#include "MyRepnumList.cpp"
#include "RepnumGroup.cpp"
#include "EfficientQuadratic.cpp"
#include "EfficientMultiPoly.cpp"
#include "MonomialOfGrits.cpp"
#include "JacobiFormData.cpp"


using namespace std;

  const int VERBOSITY=1; ///Whether to print DIAGNOSE4 messages.
  const int SMRTDIAGNOSE=0; ///Whether to print DIAGNOSE4 messages.
  const int PRINTTYPE=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE1=1; ///Whether to print DIAGNOSE1 messages PROGRESS MESSAGES.
  const int DIAGNOSE2=1; ///Whether to run test code.
  const int DIAGNOSE3=0; ///Whether to print DIAGNOSE3 messages.
  const int DIAGNOSE4=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE5=0; ///Whether to print DIAGNOSE4 messages.
  const int DIAGNOSE6=0; ///Whether to print DIAGNOSE4 messages.

  string GQfileName, infileName;
  int level, weight, uptoN, L, mindet;
  char HRfileCode[255];
  string saveFileName; stringstream saveFileNameStream; int saveFileQ;
  int maxqorder, maxd;
  vector<int> TBphi, TBTheta; int tVal, ell, ell2, ell3,tVal2, tVal3;
  vector<int>* spanningBlocks; int numBlocks;
  vector<int>* spanningBlocksWt4; int numBlocksWt4;
  vector<int>* spanningBlocksTw; int numBlocksTw;
    vector<int> numerCoeff, numerIndex1, numerIndex2,
            numerCoeffWt4, numerIndexWt4,
            numerCoeffT2, numerIndex1T2, numerIndex2T2,
            numerCoeffTw, numerIndex1Tw, numerIndex2Tw,
            denomCoeff, denomIndex;
    int doQuad, doQuadT2, doWt4, doTw;
    QZSeriesWH *jfwt4, *jftmp, *jfDenom;
    QZSeriesWH **jfGritsWt4, **jfTweaks;
    MySeriesTruncMod *jfwt4qs,*jfDenomqs, **jfGritsqs, ***jfGritsT2qs,**jfTweaksqs, *numerqs;
    EfficientQuadratic* quadGrits, *quadTweaks, *quadGritsT2;
    fmpz_t*detFactorT2;

  MonomialOfGrits* theGritMonNumer, *theGritMonDenom;
  int denomIsLinear;
  QZSeriesWH **jfGrits;
  QZSeriesWH **allJFs;
  MySeriesTruncMod **allQs;
  int numJFs;
  JacobiFormData** jfData;
  vector<int> usedJF, allWeights;

  int sa, sb, sc; //the restriction matrix
  vector<int> heckeArray, modnArray, LrootArray;
  fmpz_t* modnfmpzArray;
  int hecke, Lroot; int numHecke;
  fmpz_t modn, tmpfmpz;
  fmpz_t* LrootPowers, *LrootPowersT2;
  string outfileRoot;
  int useFile;
  //string whjfFilename;
  char whjfStr[255];
    FILE* OUTFILE;


void readParametersFilefmpzHR();
void readParametersFileGritMon();
void myfmpz_powm(fmpz_t ans, fmpz_t g, slong e, fmpz_t m);
vector<long long> RL4GH(long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, long long c, int k);
vector<long long> RL3GH(long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, int k);
vector<long long> RL1GH(long long L, long long sa, long long sb, long long sc, int k);
vector<long long> RL2GH(long long L, long long sa, long long sb, long long sc, long long a, int k);

vector<long long>* RLT2List(int ttype, long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, long long c);
MySeriesTruncMod* doOneCoset(vector<long long>* ghList, long long L, int targetDot, RepnumGroup* RNG);
MySeriesTruncMod* doOneCosetOLD(vector<long long>* ghList, long long Ldenom, int targetDot, RepnumGroup* RNG);
void printFileAbstract(FILE* OUTFILE, SigObjectOldStyle* sigObj, vector<long long>* ghList,
                       long long Ldenom, int targetDot, RepnumGroup* RNG);
void freeAllQs(MySeriesTruncMod** jfqs, vector<int> used);
void getAllQs(MySeriesTruncMod** jfqs, vector<int> used, QZSeriesWH **jfGrits,
MyRepnumList* rpL,
 long long BAa,long long BAb,long long BAc,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,int weight);
void getAllQsWithDesiredPrecision(MySeriesTruncMod** jfqs, vector<int> used,QZSeriesWH **jfGrits,
 MyRepnumList* &rpL,
 long long BAa,long long BAb,long long BAc,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,vector<int> allWeights, int desiredPrecision, RepnumGroup* RNG,
 long long sa, long long sb, long long sc);

void myfmpz_powm(fmpz_t ans, fmpz_t g, slong e, fmpz_t m){
    if(e>=0){
        fmpz_powm_ui(ans, g, (ulong)e, m);
    }else{
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_invmod(tmp, g, m);
        fmpz_powm_ui(ans, tmp, (ulong)(-e), m);
        fmpz_clear(tmp);
    }
}


int main(int argc, char* argv[])
{
    int doAll=1;
    int beginCtr=0, endCtr=0;
    int printEvery=1000;
    MySeriesTruncMod::checkSameTrunc=0;
    if(0){//testcode
        //FILE* f = fopen("testdata-3.txt", "r");
        //QZSeriesWH* qqq=new QZSeriesWH(f);
        //fclose(f);
        //cout<<qqq->getString()<<"\n";
        //qqq->multiplyWith(qqq);
        //FILE* g = fopen("testdata-4.txt", "w");
        //qqq->saveToFile(g);
        MySeriesTruncMod* s=new MySeriesTruncMod(17, 6, 0);
        MySeriesTruncMod* t=new MySeriesTruncMod(17, 7, 0);
        s->setCoeff(5,2);s->setCoeff(7,4);s->setCoeff(3,5);
        t->setCoeff(3,1);t->setCoeff(9,4);t->setCoeff(1,6);
        s->printstdout();cout<<"\n";
        //s->shift(3);
        //s->printstdout();cout<<"\n";
        t->printstdout();cout<<"\n";
        //t->shift(-1);
        //t->printstdout();cout<<"\n";
        MySeriesTruncMod* u=MySeriesTruncMod::divide(s,t);
        u->printstdout();cout<<"\n";
        cout<<"\n"<<MySeriesTruncMod::checkSameTrunc<<"\n";
        switch(doAll){
            case 0: cout<<"0\n"; break;
            case 1: cout<<"1\n"; break;
            default: cout<<"default\n Abort.\n";exit(1);
        }
        return(0);
    }

    if(argc<3){
        cout << "Usage requires two or more arguments: GritMon-file HR-file,  [uptoN], [beginCtr], [endCtr], [printEvery], [saveFileName].\n";
        cout<<"uptoN by default is read from the HR file. Use 0 for default.\n";
        cout<<"Runs cosets numbered beginCtr to endCtr if set. Set beginCtr=0 and endCtr=0 if you want to do all cosets.\n";
        cout<<"Output progress every printEvery. Use 0 for default.\n";
        cout<<"saveFileName from file is used by default.\n";
        cout<<"Info about this machine:\n";
        printf("sizeof(slong) = %d;\n", sizeof(slong));
        printf("sizeof(int) = %d;\n", sizeof(int));
        printf("sizeof(long long) = %d;\n", sizeof(long long));

    return 1;
    }
    GQfileName = argv[1];
    infileName = argv[2];
    uptoN=0;
    if(argc>3){
     uptoN = atoi(argv[3]);
    }
    if(argc>5){
     beginCtr = atoi(argv[4]);
     endCtr = atoi(argv[5]);
     if((beginCtr>0)&&(endCtr>0)){doAll=0;}
    }
    if(argc>6){
     printEvery = atoi(argv[6]);
     if(printEvery==0){printEvery=1000;}
    }
    if(argc>7){
     saveFileNameStream<<argv[7];
     saveFileNameStream>>saveFileName;
     saveFileQ=1;
    }

    fmpz_init(modn);    fmpz_init(tmpfmpz);
    readParametersFilefmpzHR();
    readParametersFileGritMon();


    string f;
    stringstream ffs;
    ffs<<"output-"<<GQfileName<<"-"<<infileName;
    ffs>>outfileRoot;
    ffs.clear();
    if(saveFileQ){
        f=saveFileName;
    }else{
        if(doAll){
            ffs<<outfileRoot<<"-"<<uptoN<<".ma";
        }else{
            ffs<<outfileRoot<<"-"<<uptoN<<"-"<<beginCtr<<"-"<<endCtr<<".ma";
        }
        //ffs<<outfileRoot<<"-"<<modn<<"-"<<Lroot<<"-"<<beginHecke<<"-"<<endHecke<<"-"<<uptoN<<".ma";
        ffs>>f;
    }
    cout<<"Output file will be "<<f<<"\n";
    OUTFILE = fopen(f.c_str(), "w");
    if(OUTFILE==0){cout<<"File write open failure: "<<f<<"\n";exit(1);}

    //OUTFILE.open(f.c_str());
    //if(OUTFILE.fail()){cout<<"File write open failure: "<<f<<"\n";exit(1);}


    int H, a, b, c, i, j, k;
    int maxHecke=heckeArray[0];
    for(H=1;H<numHecke;H++){
        maxHecke=max(maxHecke, heckeArray[H]);
    }

    int maxUptoDet;
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke; //See notebook for theory
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke/(sa*sc-sb*sb); //20190120
    SigObjectOldStyle* sigObj;
    if(SMRTDIAGNOSE){sigObj = new SigObjectOldStyle(level,maxUptoDet);}  //if diagnostics
    else{sigObj = NULL;} //new SigObjectOldStyle(level,maxUptoDet);

    maxqorder = (level*level+4*maxUptoDet)/(4*level) + uptoN; //20190104 try
    maxqorder = (level*level+4*maxUptoDet)/(4*level) + 200; //20190104 try
    maxqorder = (level*level+4*maxUptoDet)/(4*level);
    maxqorder = (level*level+4*maxUptoDet)/(4*level) + 10; //20200601 call this GQHR200
    cout<<"maxUptoDet = "<<maxUptoDet<<".\n";
    cout<<"maxqorder (for expanding TBs) = "<<maxqorder<<".\n";

    jfGrits=new QZSeriesWH*[numJFs];
    //jfGrits=new QZSeriesWH*[numBlocks];
    //jfGritsWt4=new QZSeriesWH*[numBlocksWt4];
    //jfTweaks=new QZSeriesWH*[numBlocksTw];

    int jfTrunc=maxqorder+2;    //Do theory later
    for(i=0;i<numJFs;i++){
            cout<<"Making JF #"<<i<<"\n";
            jfGrits[i]=jfData[i]->makeJF(jfTrunc);
    }

    //cout<<"HERE B1.\n";
    allQs=new MySeriesTruncMod*[numJFs];
    for(i=0;i<numJFs;i++)allQs[i]=NULL;
    //cout<<"HERE B1.\n";

    detFactorT2=new fmpz_t[16]; //For 1 to 15 reall
    for(i=0;i<16;i++){fmpz_init(detFactorT2[i]);}
    fmpz_set_si(tmpfmpz,2);
        for(i=1;i<16;i++){
            if(i==1){myfmpz_powm(detFactorT2[i],tmpfmpz,2*2*weight-3,modn);}
            else if((2<=i)&&(i<=7)){
                    myfmpz_powm(detFactorT2[i],tmpfmpz,2*weight-3,modn);
            }
            else if((8<=i)&&(i<=15)){
                    myfmpz_powm(detFactorT2[i],tmpfmpz,-3,modn);
            }
        }
    //cout<<"HERE B3.\n";
    RepnumGroup* RNG=new RepnumGroup(level,1);  //mindet is 1 because all Jacobi cusp forms.
    //cout<<"HERE B4.\n";
    MySeriesTruncMod* totqs, *oneqs, *partialsum;
    MySeriesTruncMod* phifqs=NULL, *origqs=NULL;
    fmpz_t detAfactor, fmpzL;
    fmpz_init(detAfactor);  fmpz_init(fmpzL);
    int ctr, totctr, globalctr;
    vector<long long>* ghL;
    //cout<<"HERE B5.\n";
    MyRepnumList* RL;
    int useShortcut, upperBound;



    for(H=0;H<numHecke;H++){
        hecke=heckeArray[H];
        fmpz_set(modn,modnfmpzArray[H]);
        Lroot=LrootArray[H];
        L=hecke; fmpz_set_si(fmpzL, L);
    //cout<<"HERE B6.\n";
        LrootPowers=new fmpz_t[L];
        for(i=0;i<L;i++){fmpz_init(LrootPowers[i]);}
        LrootPowersT2=new fmpz_t[2*L];
        for(i=0;i<2*L;i++){fmpz_init(LrootPowersT2[i]);}

        cout<<"Temporarily setting LrootPowers and LrootPowersT2 for L=1.\n";
        fmpz_set_si(LrootPowers[0],1);
        fmpz_set_si(LrootPowersT2[0],1);
        fmpz_set_si(LrootPowersT2[1],-1);

        cout<<"Getting ghL to make oriqqs... \n";
        ghL=RLT2List(4, 1, sa, sb, sc, 0, 0, 0);
        cout<<"Making origqs... \n";

        origqs = doOneCoset(ghL,1, uptoN, RNG);

        cout<<"origqs=";origqs->printstdout(); cout<<"\n";
        fprintf(OUTFILE,"origqs = (");origqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");
        fprintf(OUTFILE,"origqsVan=%d; origqsLen=%d; origqsTrunc=%d;\n", origqs->van, origqs->len, origqs->getTrunc());
        if(SMRTDIAGNOSE){fprintf(OUTFILE,"origqsAbstract=");printFileAbstract(OUTFILE,sigObj,ghL,1,uptoN,RNG);}
        delete[] ghL;

//cout<<"END TEST\n"; return 0;

        fmpz_set_si(tmpfmpz,Lroot);
        for(i=0;i<L;i++){
            myfmpz_powm(LrootPowers[i], tmpfmpz, 2*i, modn);  //Notice the 2*i so as to be consistent with LrootPowersT2
        }
    //cout<<"HERE B7.\n";
        fmpz_set_si(tmpfmpz,-Lroot);
        for(i=0;i<2*L;i++){
            myfmpz_powm(LrootPowersT2[i], tmpfmpz, i, modn);
        }
    //cout<<"HERE B8.\n";




        cout<<"**************** L="<<L<<" ***********************\n";
        totctr=1+1+L+L*L; ctr=0; globalctr=0;
        totqs=new MySeriesTruncMod(modn, uptoN*L);

        //TYPE 1
        cout<<"Doing type 1 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<".\n";
        globalctr++;
        if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
            ctr++;
            if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                cout<<"\n";
            }
            oneqs = origqs->expandExpBy(L*L, uptoN*L);
            myfmpz_powm(detAfactor, fmpzL, weight, modn);
            //cout<<"DIAGNOSTICS A: detAfactor="<<fmpz_get_si(detAfactor)<<".\n";
            oneqs->scalarMultiplityBy(detAfactor);
            totqs->addBy(oneqs);
            if(PRINTTYPE){fprintf(OUTFILE,"TYPE1 = (");oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}

            delete oneqs;
        }
        //TYPE 2
        cout<<"Doing type 2 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        if(sa%L==0){cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;}
        else{cout<<"Shortcut available.\n";useShortcut=1;upperBound=1;}
        for(a=0;a<upperBound;a++){
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    //cout<<"\n";
                    cout<<" [ type 2 a="<<a<<"]\n";
                }
                //oneqs=RL2->oneBPCosetRestrictionPrepped(a,0,0, L, modn, Lroot, uptoN*L);
                //detAfactor is 1 for Type 2 (multiply everything by L^(weight-3) later
                ghL=RLT2List(2, L, sa, sb, sc, a, 0, 0);
                oneqs = doOneCoset(ghL,L, uptoN, RNG);
                //cout<<"HERE Q1\n";
                if(useShortcut)oneqs->scalarMultiplityBy(fmpzL);
                totqs->addBy(oneqs);
                //cout<<"HERE Q2\n";
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE2X%d=(",a);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                if(SMRTDIAGNOSE){fprintf(OUTFILE,"TYPE2Abstract%d=",a);printFileAbstract(OUTFILE,sigObj,ghL,L,uptoN,RNG);}
                delete[] ghL;
                delete oneqs;
            }
       }

        //TYPE 3
        cout<<"Doing type 3 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=0;a<L;a++){
         RL=RNG->getRepnumList(L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa,L*uptoN);
         if(RL->hasDesiredDot(L*uptoN)){

//cout<<"TEST CODE ONLY\n";
//delete RL;
//cout<<"HERE F7\n";
//return 0;

         if(VERBOSITY>4)cout<<"RNG->getRepnumList("<<L*L*sa<<","<< L*(sb+a*sa)<<","<< sc+2*a*sb+a*a*sa<<","<<L*uptoN<<");\n";
          if((sc+2*a*sb+a*a*sa)%L==0){cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;}
          else{cout<<"Shortcut available.\n";useShortcut=1;upperBound=1;}
          for(b=0;b<upperBound;b++){
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 3 a="<<a<<",b="<<b<<"]\n";
                    cout.flush();
                }
                //oneqs=RL3[a]->oneBPCosetRestrictionPrepped(0,0,b, L, modn, Lroot, uptoN*L);
                //detAfactor is 1 for Type 3
                //cout<<"DIAGNOSTICS A: type 3: ";oneqs->printstdout();cout<<".\n";
                ghL=RLT2List(3, L, sa, sb, sc, a, b, 0);
                oneqs = doOneCoset(ghL,L, uptoN, RNG);
                //cout<<"HERE Q1\n";
                if(useShortcut)oneqs->scalarMultiplityBy(fmpzL);
                totqs->addBy(oneqs);
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE3X%dX%X = (",a,b);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                if(SMRTDIAGNOSE){fprintf(OUTFILE,"TYPE3Abstract%dX%d=",a,b);printFileAbstract(OUTFILE,sigObj,ghL,L,uptoN,RNG);}
                delete oneqs;
                delete[] ghL;
           }
          }
          }
                //cout<<"HERE F6\n";
          //delete RL;
          //if(useShortcut)delete RL;
"SOLUTION IS TO JUST ALLOW MEMORY LEAK?????\n";
"WE can track down the bug later.";
"ACTUALLY YOU MUST NOT DELETE MyRepnumList's that were gotten from RepnumGroup!!!";
                //cout<<"HERE F7\n";
       }

        //TYPE 4
        cout<<"Doing type 4 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        partialsum=new MySeriesTruncMod(modn, uptoN*L);
        if(sa%L==0){cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;}
        else{cout<<"Shortcut available.\n";useShortcut=1;upperBound=1;}
        for(a=0;a<upperBound;a++)
            for(b=-((int)(L/2));b<-((int)(L/2))+L;b++)
            for(c=-((int)(L/2));c<-((int)(L/2))+L;c++){
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<". partialsum=";partialsum->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<") ("
                        <<(100*(globalctr-beginCtr+1))/(endCtr-beginCtr+1)<<"%)";}
                    cout<<" [ type 4 a="<<a<<",b="<<b<<",c="<<c<<"]\n";
                }
                //oneqs=RL4->oneBPCosetRestrictionPrepped(a,b,c, L, modn, Lroot, uptoN*L);
                ghL=RLT2List(4, L, sa, sb, sc, a, b, c);
                oneqs = doOneCoset(ghL,L, uptoN, RNG);
                //cout<<"HERE Q1\n";
                partialsum->addBy(oneqs);
                if((beginCtr>0)&&(endCtr==beginCtr)){
                    fprintf(OUTFILE,"abc={%d,%d,%d};\n",a,b,c);
                }
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE4X%dX%dX%d = (",(a+L)%L,(b+L)%L,(c+L)%L);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                if(SMRTDIAGNOSE){fprintf(OUTFILE,"TYPE4Abstract%dX%dX%d=",(a+L)%L,(b+L)%L,(c+L)%L);printFileAbstract(OUTFILE,sigObj,ghL,L,uptoN,RNG);}
                delete[] ghL;
                delete oneqs;
            }
        }
        if(useShortcut)partialsum->scalarMultiplityBy(fmpzL);
        myfmpz_powm(detAfactor, fmpzL, -weight, modn);
        partialsum->scalarMultiplityBy(detAfactor);
        totqs->addBy(partialsum);
        if(PRINTTYPE){fprintf(OUTFILE,"TYPE4partialsum = (");partialsum->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
        delete partialsum;

        myfmpz_powm(detAfactor, fmpzL, weight-3, modn);
        totqs->scalarMultiplityBy(detAfactor);



        cout<<"Made totqs = ";totqs->printstdout();cout<<"\n";

        fprintf(OUTFILE, "L=%d;\n", L);
        fprintf(OUTFILE, "doAll=%d;\n", doAll);
        fprintf(OUTFILE, "beginCtr=%d;\n", beginCtr);
        fprintf(OUTFILE, "endCtr=%d;\n", endCtr);
        fprintf(OUTFILE, "modn=");
        fmpz_fprint(OUTFILE,modn);
        fprintf(OUTFILE, ";\n");
        fprintf(OUTFILE, "Lroot=%d;\n", Lroot);
        fprintf(OUTFILE, "sizeofslong = %d;\n", sizeof(slong));
        fprintf(OUTFILE, "sizeofint = %d;\n", sizeof(int));
        fprintf(OUTFILE, "sizeoflonglong = %d;\n", sizeof(long long));

        fprintf(OUTFILE,"totqs = (");totqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");
        fprintf(OUTFILE,"totqsVan=%d; totqsLen=%d; totqsTrunc=%d;\n", totqs->van, totqs->len, totqs->getTrunc());
        fprintf(OUTFILE, "domeSomething;\nseparator=nothing;\n\n");
        fflush(OUTFILE);



        delete totqs;
    }
return 0;

    MyRepnumList* RL4;
    MyRepnumList *RL2;
    MyRepnumList **RL3;
    int *useRL3;


    for(H=0;H<numHecke;H++){
        hecke=heckeArray[H];
        //fmpz_set_si(modn,modnArray[H]);

        fmpz_set(modn,modnfmpzArray[H]);

        Lroot=LrootArray[H];
        L=hecke;
        cout<<"**************** L="<<L<<" ***********************\n";


        //totctr=1+L+L*L+L*L*L; ctr=0; globalctr=0;
        totctr=1+1+L+L*L; ctr=0; globalctr=0;





        for(a=0;a<L;a++){
                //cout<<"HEREA.";
            delete RL3[a];
                //cout<<"HEREB.";
        }
        delete[] RL3; //cout<<"HEREC.";
        delete[] useRL3;
        delete RL4; //cout<<"HERED.";
        delete RL2; //cout<<"HEREE.";
        delete totqs; totqs=NULL;
        delete phifqs; phifqs=NULL;
    }






    fmpz_clear(modn); fmpz_clear(detAfactor); fmpz_clear(fmpzL); fmpz_clear(tmpfmpz);
    for(i=0;i<numHecke;i++){
            fmpz_clear(modnfmpzArray[i]);
    }
    delete[] modnfmpzArray;

    return 0;
}

void readParametersFileGritMon(){
  FILE* PFILE;
  int x;
  cout<<"Reading GritMon parameters from "<<GQfileName<<" ...\n";
  PFILE = fopen(GQfileName.c_str(),"r");
  if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
  fscanf(PFILE, "%d", &x);
  if((weight>0)&&(weight!=x)){
    cout<<"weight "<<weight<<" does not matching existing weight. Abort.\n";
    exit(1);
  }
  weight=x; //Check against existing weight??
  fscanf(PFILE, "%d", &x);
  if((level>0)&&(level!=x)){
    cout<<"level "<<level<<" does not matching existing level. Abort.\n";
    exit(1);
  }
  level=x;
  fscanf(PFILE, "%d", &x);
  numJFs=x;
  cout<<"Will read "<<numJFs<<" JFs.\n";
  jfData=new JacobiFormData*[numJFs];
  usedJF.clear(); allWeights.clear();
  int jfType, wt, lev, tee,len,i,j, numH;
  vector<int> ds;
  for(i=0;i<numJFs;i++){
    usedJF.push_back(0);
    fscanf(PFILE, "%d", &x);
    //cout<<x<<" ";
    if(x!=i){
      cout<<"Numbering of jf's not consecutive. Aborting to be safe.\n";
      exit(0);
    }
    fscanf(PFILE, "%d", &jfType); 
    fscanf(PFILE, "%d", &lev); 
    fscanf(PFILE, "%d", &wt); 
    fscanf(PFILE, "%d", &tee); 
    fscanf(PFILE, "%d", &len); 
    allWeights.push_back(wt);
    ds.clear();
    fscanf(PFILE, "%d", &numH); 
    ds.push_back(numH);
    for(j=0;j<len+2*numH;j++){
      fscanf(PFILE, "%d", &x);
      ds.push_back(x);
    }
    jfData[i]=new JacobiFormData(jfType,wt,lev,tee,ds);
  }
  int indSize, indx, expx;
  vector<int> ind, exponent;

  fscanf(PFILE, "%d", &indSize); 
  if(VERBOSITY>1){cout<<"indSize="<<indSize<<"\n";}
  ind.clear(); exponent.clear();
  for(i=0;i<indSize;i++){
    fscanf(PFILE, "%d", &indx); 
    fscanf(PFILE, "%d", &expx); 
    ind.push_back(indx);exponent.push_back(expx);
    usedJF[indx]=1;
    if(VERBOSITY>1){cout<<indx<<","<<expx<<"\n";}
  }
  theGritMonNumer = new MonomialOfGrits(ind, exponent);
  if(VERBOSITY>3){cout<<"[A] theGritMonNumer=";theGritMonNumer->printRepresentation();cout<<"\n";}

  fscanf(PFILE, "%d", &indSize);
  if(VERBOSITY>1){cout<<"indSize="<<indSize<<"\n";}
  ind.clear(); exponent.clear();
  for(i=0;i<indSize;i++){
    fscanf(PFILE, "%d", &indx);
    fscanf(PFILE, "%d", &expx);
    ind.push_back(indx);exponent.push_back(expx);
    usedJF[indx]=1;
    if(VERBOSITY>1){cout<<indx<<","<<expx<<"\n";}
  }
  theGritMonDenom = new MonomialOfGrits(ind, exponent);
  if(VERBOSITY>3){cout<<"[A] theGritMonDenom=";theGritMonDenom->printRepresentation();cout<<"\n";}

}

void readParametersFilefmpzHR(){
  FILE* PFILE;
  cout<<"Reading Hecke Restriction parameters from "<<infileName<<" ...\n";
  PFILE = fopen(infileName.c_str(),"r");
  if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
  int x, i, j, len;
  fscanf(PFILE, "%d", &weight);
  fscanf(PFILE, "%d", &level);



    fscanf(PFILE, "%d", &x);  sa=x;
    fscanf(PFILE, "%d", &x);  sb=x;
    fscanf(PFILE, "%d", &x);  sc=x;
    fscanf(PFILE, "%d", &x);  if(uptoN==0){uptoN=x;}

    //fscanf(PFILE, "%s", HRfileCode);

    fscanf(PFILE, "%d", &x);  numHecke=x;


    cout<<"weight = "<<weight<<".\n";
    cout<<"level = "<<level<<".\n";

  cout<<"restriction matrix = "<<sa<<","<<sb<<","<<sc<<".\n";

  cout<<"Desired up to dot product, aka uptoN = "<<uptoN<<".\n";

  heckeArray.clear(); modnArray.clear(); LrootArray.clear();
  cout<<"Will read "<<numHecke<<" sets of (hecke, modulus, Lroot):\n";
  modnfmpzArray=new fmpz_t[numHecke];
  for(i=0;i<numHecke;i++){
        fscanf(PFILE, "%d", &x);  heckeArray.push_back(x);
        cout<<"Hecke = "<<x<<", ";
        fmpz_fread(PFILE, modn);
        fmpz_init_set(modnfmpzArray[i],modn);
        cout<<"modulus = "<<fmpz_get_si(modn)<<", ";
        fscanf(PFILE, "%d", &x);   LrootArray.push_back(x);
        cout<<"Lroot = "<<x<<".\n";
  }
  fclose(PFILE);
}

vector<long long> RL1GH(long long L, long long sa, long long sb, long long sc, int k){
    vector<long long>ans;
    switch(k){
        case 0: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
        case 1: ans.push_back(4*L*L*sa);
            ans.push_back(4*L*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 2: ans.push_back(L*L*sa);
            ans.push_back(2*L*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 3: ans.push_back(L*L*sa);
            ans.push_back(2*L*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 4: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 5: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 6: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(L*sa + L*sb));
            ans.push_back(L*(L*sa + 2*L*sb + L*sc));
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 7: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(L*sa + L*sb));
            ans.push_back(L*(L*sa + 2*L*sb + L*sc));
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 8: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 9: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 10: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(L);
            ans.push_back(0);
            break;
        case 11: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(0);
            ans.push_back(L);
            ans.push_back(L);
            break;
        case 12: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 13: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 14: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(L);
            ans.push_back(L);
            ans.push_back(0);
            break;
        case 15: ans.push_back(L*L*sa);
            ans.push_back(L*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(L);
            ans.push_back(L);
            ans.push_back(L);
            break;
        default: cout<<"ERROR dfas3w8dsf. Abort.\n";exit(1);
    };
    return ans;
}
vector<long long> RL2GH(long long L, long long sa, long long sb, long long sc, long long a, int k){
    vector<long long>ans;
    switch(k){
        case 0: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a);
            ans.push_back(0);
            ans.push_back(0);
        case 1: ans.push_back(4*sa);
            ans.push_back(4*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(4*a);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 2: ans.push_back(sa);
            ans.push_back(2*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(a);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 3: ans.push_back(sa);
            ans.push_back(2*L*sb);
            ans.push_back(4*L*L*sc);
            ans.push_back(a + L);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 4: ans.push_back(4*sa);
            ans.push_back(2*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(4*a);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 5: ans.push_back(4*sa);
            ans.push_back(2*L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(4*a);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 6: ans.push_back(4*sa);
            ans.push_back(2*(sa + L*sb));
            ans.push_back(sa + 2*L*sb + L*L*sc);
            ans.push_back(4*a);
            ans.push_back(2*a);
            ans.push_back(a);
            break;
        case 7: ans.push_back(4*sa);
            ans.push_back(2*(sa + L*sb));
            ans.push_back(sa + 2*L*sb + L*L*sc);
            ans.push_back(4*a);
            ans.push_back(2*a);
            ans.push_back(a + L);
            break;
        case 8: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 9: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 10: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a);
            ans.push_back(L);
            ans.push_back(0);
            break;
        case 11: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a);
            ans.push_back(L);
            ans.push_back(L);
            break;
        case 12: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a + L);
            ans.push_back(0);
            ans.push_back(0);
            break;
        case 13: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a + L);
            ans.push_back(0);
            ans.push_back(L);
            break;
        case 14: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a + L);
            ans.push_back(L);
            ans.push_back(0);
            break;
        case 15: ans.push_back(sa);
            ans.push_back(L*sb);
            ans.push_back(L*L*sc);
            ans.push_back(a + L);
            ans.push_back(L);
            ans.push_back(L);
            break;
        default: cout<<"ERROR dfas3w8dsf. Abort.\n";exit(1);
    };
    return ans;
}
vector<long long> RL3GH(long long L, long long sa, long long sb, long long sc, long long a, long long b, int k){
    vector<long long>ans;
    switch(k){
        case 0: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b);
        case 1: ans.push_back(4*L*L*sa);
            ans.push_back(4*L*(a*sa + sb));
            ans.push_back(4*(a*a*sa + 2*a*sb + sc));
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(4*b);
            break;
        case 2: ans.push_back(L*L*sa);
            ans.push_back(2*L*(a*sa + sb));
            ans.push_back(4*(a*a*sa + 2*a*sb + sc));
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(4*b);
            break;
        case 3: ans.push_back(L*L*sa);
            ans.push_back(2*L*(a*sa + sb));
            ans.push_back(4*(a*a*sa + 2*a*sb + sc));
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(4*b);
            break;
        case 4: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b);
            break;
        case 5: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b + L);
            break;
        case 6: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(a*sa + L*sa + sb));
            ans.push_back(a*a*sa + 2*a*L*sa + L*L*sa + 2*a*sb + 2*L*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b);
            break;
        case 7: ans.push_back(4*L*L*sa);
            ans.push_back(2*L*(a*sa + L*sa + sb));
            ans.push_back(a*a*sa + 2*a*L*sa + L*L*sa + 2*a*sb + 2*L*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b + L);
            break;
        case 8: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b);
            break;
        case 9: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(0);
            ans.push_back(b + L);
            break;
        case 10: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(L);
            ans.push_back(b);
            break;
        case 11: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(0);
            ans.push_back(L);
            ans.push_back(b + L);
            break;
        case 12: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(b);
            break;
        case 13: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(L);
            ans.push_back(0);
            ans.push_back(b + L);
            break;
        case 14: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(L);
            ans.push_back(L);
            ans.push_back(b);
            break;
        case 15: ans.push_back(L*L*sa);
            ans.push_back(L*(a*sa + sb));
            ans.push_back(a*a*sa + 2*a*sb + sc);
            ans.push_back(L);
            ans.push_back(L);
            ans.push_back(b + L);
            break;
        default: cout<<"ERROR dfas3w8dsf. Abort.\n";exit(1);
    };
    return ans;
}
vector<long long> RL4GH(long long L, long long sa, long long sb, long long sc, long long a, long long b, long long c, int k){
    vector<long long>ans;
    switch(k){
        case 0: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a);
            ans.push_back(b);
            ans.push_back(c);
        case 1: ans.push_back(4*sa);
            ans.push_back(4*sb);
            ans.push_back(4*sc);
            ans.push_back(4*a);
            ans.push_back(4*b);
            ans.push_back(4*c);
            break;
        case 2: ans.push_back(sa);
            ans.push_back(2*sb);
            ans.push_back(4*sc);
            ans.push_back(a);
            ans.push_back(2*b);
            ans.push_back(4*c);
            break;
        case 3: ans.push_back(sa);
            ans.push_back(2*sb);
            ans.push_back(4*sc);
            ans.push_back(a + L);
            ans.push_back(2*b);
            ans.push_back(4*c);
            break;
        case 4: ans.push_back(4*sa);
            ans.push_back(2*sb);
            ans.push_back(sc);
            ans.push_back(4*a);
            ans.push_back(2*b);
            ans.push_back(c);
            break;
        case 5: ans.push_back(4*sa);
            ans.push_back(2*sb);
            ans.push_back(sc);
            ans.push_back(4*a);
            ans.push_back(2*b);
            ans.push_back(c + L);
            break;
        case 6: ans.push_back(4*sa);
            ans.push_back(2*(sa + sb));
            ans.push_back(sa + 2*sb + sc);
            ans.push_back(4*a);
            ans.push_back(2*(a + b));
            ans.push_back(a + 2*b + c);
            break;
        case 7: ans.push_back(4*sa);
            ans.push_back(2*(sa + sb));
            ans.push_back(sa + 2*sb + sc);
            ans.push_back(4*a);
            ans.push_back(2*(a + b));
            ans.push_back(a + 2*b + c + L);
            break;
        case 8: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a);
            ans.push_back(b);
            ans.push_back(c);
            break;
        case 9: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a);
            ans.push_back(b);
            ans.push_back(c + L);
            break;
        case 10: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a);
            ans.push_back(b + L);
            ans.push_back(c);
            break;
        case 11: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a);
            ans.push_back(b + L);
            ans.push_back(c + L);
            break;
        case 12: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a + L);
            ans.push_back(b);
            ans.push_back(c);
            break;
        case 13: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a + L);
            ans.push_back(b);
            ans.push_back(c + L);
            break;
        case 14: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a + L);
            ans.push_back(b + L);
            ans.push_back(c);
            break;
        case 15: ans.push_back(sa);
            ans.push_back(sb);
            ans.push_back(sc);
            ans.push_back(a + L);
            ans.push_back(b + L);
            ans.push_back(c + L);
            break;
        default: cout<<"ERROR dfas3w8dsf. Abort.\n";exit(1);
    };
    return ans;
}



vector<long long>* RLT2List(int ttype, long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, long long c){
    vector<long long>*ans = new vector<long long>[16];
    int k;
    //cout<<"HERE A.\n";
    switch(ttype){
        case 1: for(k=0;k<=15;k++){ans[k]=RL1GH(L,sa,sb,sc,k);}
            break;
        case 2: for(k=0;k<=15;k++){ans[k]=RL2GH(L,sa,sb,sc,a,k);}
            break;
        case 3: for(k=0;k<=15;k++){ans[k]=RL3GH(L,sa,sb,sc,a,b,k);}
            break;
        case 4: for(k=0;k<=15;k++){
            //cout<<"HEREB "<<k<<"\n";
            ans[k]=RL4GH(L,sa,sb,sc,a,b,c,k);}
            break;
        default: cout<<"ERROR dsg34y23ds. Abort.\n";exit(1);
    }
    //cout<<"HERE Az.\n";
    return ans;
}

void printFileAbstract(FILE* OUTFILE, SigObjectOldStyle* sigObj, vector<long long>* ghList, long long Ldenom, int targetDot, RepnumGroup* RNG){
    int i,k;int currUpto=(Ldenom+1)*targetDot;
    MyRepnumList *rpL;
    long long sa=ghList[0][0], sb=ghList[0][1], sc=ghList[0][2],
        BAa=ghList[0][3],BAb=ghList[0][4],BAc=ghList[0][5];
    rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
    rpL->printAbstractRestrictionToFile(OUTFILE,sigObj,BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,weight);
}

void freeAllQs(MySeriesTruncMod** jfqs, vector<int> used){
  int i;
  for(i=0;i<used.size();i++){
    if(used[i])if(jfqs[i]!=NULL){delete jfqs[i]; jfqs[i]=NULL;}
  }
}
void getAllQs(MySeriesTruncMod** jfqs, vector<int> used,QZSeriesWH **jfGrits,
 MyRepnumList* rpL,
 long long BAa,long long BAb,long long BAc,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,vector<int> allWeights){
  //if(VERBOSITY>4)cout<<"Freeing jfqs.\n";
  freeAllQs(jfqs,used);
  if(VERBOSITY>4)cout<<"Freed jfqs.\n";
  int i;
  for(i=0;i<used.size();i++){
    //if(VERBOSITY>4){cout<<i<<"\n";}
    if(used[i]){
     // if(VERBOSITY>4){cout<<"HERE G2: "<<used[i]<<","<<allWeights[i]<<"\n";}
      jfqs[i]=rpL->restrictGrit(jfGrits[i],BAa,BAb,BAc,Ldenom,modn,LrootPowers,
                    uptoNup,allWeights[i]);
      if(VERBOSITY>4){cout<<"allQs["<<i<<"]=";allQs[i]->printstdoutWithTruncationOrder();cout<<"\n";}
    }
  }
}
void getAllQsWithDesiredPrecision(MySeriesTruncMod** jfqs, vector<int> used,QZSeriesWH **jfGrits,
 MyRepnumList* &rpL,
 long long BAa,long long BAb,long long BAc,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,vector<int> allWeights, int desiredPrecision, RepnumGroup* RNG,
 long long sa, long long sb, long long sc){
  //if(VERBOSITY>4)cout<<"Freeing jfqs.\n";
  freeAllQs(jfqs,used);
  if(VERBOSITY>4)cout<<"Freed jfqs.\n";
  int i, ctr, rpLUpto=uptoNup, currUpto;
  for(i=0;i<used.size();i++){
    //if(VERBOSITY>4){cout<<i<<"\n";}
    if(used[i]){
     // if(VERBOSITY>4){cout<<"HERE G2: "<<used[i]<<","<<allWeights[i]<<"\n";}
      currUpto = uptoNup;
      for(ctr=0;ctr<100;ctr++){
          if(currUpto>rpLUpto){
            if(VERBOSITY>3){cout<<"Need longer rpL: rpLUpto, currUpto "<<rpLUpto<<","<<currUpto<<".\n";}
            //delete rpL;  ''DON'T DELETE!!! RNG takes care of that!
            rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
            rpLUpto=currUpto;
          }
          jfqs[i]=rpL->restrictGrit(jfGrits[i],BAa,BAb,BAc,Ldenom,modn,LrootPowers,
                    currUpto,allWeights[i]);
          jfqs[i]->normalize();
          if(jfqs[i]->len>=desiredPrecision)break;
          if(VERBOSITY>3){cout<<"Need longer jfqs[i]: currUpto, van, len = "<<currUpto<<","<<jfqs[i]->van<<","<<jfqs[i]->len<<".\n";}
          currUpto = currUpto+uptoNup;
      }
      if(VERBOSITY>4){cout<<"allQs["<<i<<"]=";allQs[i]->printstdoutWithTruncationOrder();cout<<"\n";}
    }
  }
}

MySeriesTruncMod* doOneCoset(vector<long long>* ghList, long long Ldenom, int targetDot, RepnumGroup* RNG){
    int desiredPrecision=targetDot*Ldenom;;
    int i,k;
    MyRepnumList *rpL;
    //fmpz_t L;  //Just to root out typos.
    if(VERBOSITY>3)cout<<"doOneCoset...\n";
    if(VERBOSITY>3){cout<<"theGritMonNumer=";theGritMonNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>3){cout<<"theGritMonDenom=";theGritMonDenom->printRepresentation();cout<<"\n";}
    long long sa=ghList[0][0], sb=ghList[0][1], sc=ghList[0][2],
        BAa=ghList[0][3],BAb=ghList[0][4],BAc=ghList[0][5];
    int currUpto=2*Ldenom*targetDot;
    if(VERBOSITY>1)cout<<"doOneCoset:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BAc<<"; finding mtr...\n";
    if(VERBOSITY>4){cout<<"desiredPrecision="<<desiredPrecision<<", (L,targetDot,currUpto)=("<<L<<","<<targetDot<<","<<currUpto<<").\n";}
    rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
    getAllQsWithDesiredPrecision(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,allWeights,desiredPrecision, RNG, sa,sb,sc);
    jfDenomqs=theGritMonDenom->evaluate(allQs);
    if(VERBOSITY>4){cout<<"jfDenomqs=";jfDenomqs->printstdoutWithTruncationOrder();cout<<"\n"; }
    if(VERBOSITY>4)cout<<"Putting together numerator...\n";
    numerqs=theGritMonNumer->evaluate(allQs);
    if(VERBOSITY>1){cout<<"numerqs=";numerqs->printstdoutWithTruncationOrder();cout<<"\n";
        cout<<"jfDenomqs=";jfDenomqs->printstdoutWithTruncationOrder();cout<<"\n";
        cout<<"Dividing by jfDenomqs..\n";}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"numerqs=");numerqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"jfDenomqs=");jfDenomqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    numerqs->divideBy(jfDenomqs);
    if(VERBOSITY>1){cout<<"quotient=";numerqs->printstdoutWithTruncationOrder();cout<<"\n";}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"quotientqs=");numerqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    delete jfDenomqs;
    if(VERBOSITY>4){cout<<"returning from doOneCoset function.\n";}
    return numerqs;
}

MySeriesTruncMod* doOneCosetOLD(vector<long long>* ghList, long long Ldenom, int targetDot, RepnumGroup* RNG){
    int mtr=-1;
    int i,k;
    MyRepnumList *rpL;
    //fmpz_t L;  //Just to root out typos.
    if(VERBOSITY>3)cout<<"doOneCoset...\n";
    if(VERBOSITY>3){cout<<"theGritMonNumer=";theGritMonNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>3){cout<<"theGritMonDenom=";theGritMonDenom->printRepresentation();cout<<"\n";}
    long long sa=ghList[0][0], sb=ghList[0][1], sc=ghList[0][2],
        BAa=ghList[0][3],BAb=ghList[0][4],BAc=ghList[0][5];
    int currUpto=Ldenom*targetDot;
    if(VERBOSITY>1)cout<<"doOneCoset:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BAc<<"; finding mtr...\n";
    mtr=currUpto; //TRY THIS??? 20200604
    if(0) while(mtr<0){ //SKIP THIS AND TRY IT? 20200604
        if(VERBOSITY>4)cout<<"making rpL.\n";
        rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
        if(VERBOSITY>4)cout<<"making jfDenomqs, first getAllQs.\n";
        if(VERBOSITY>4)cout.flush();
        //jfDenomqs=rpL->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,weight);
        getAllQs(allQs, usedJF, jfGrits, rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,allWeights);
        if(VERBOSITY>4)cout<<"making jfDenomqs.\n";
        jfDenomqs=theGritMonDenom->evaluate(allQs);
        if(VERBOSITY>4)cout<<"Made a putative jfDenomqs.\n";
        mtr=jfDenomqs->getVanishingOrder();
        if(mtr>=0)break;
        //if(currUpto>=maxqorder){
          //cout<<"Already doubled to the max. Must abort dgs8g8as.\n"; exit(1);
        //}
        currUpto=2*currUpto;  //Keep doubling until we get nonzero restriction.
        //if(currUpto>=maxqorder){currUpto=maxqorder;}
        if(VERBOSITY>4)cout<<"Doubling currUpto (now "<<currUpto<<") to get mtr in doOneCoset.\n";
        delete jfDenomqs;
    }
    fmpz_t uptoN;
    mtr = 12*currUpto;
    int uptoNup=targetDot*Ldenom+mtr;
    if(VERBOSITY>4){cout<<"mtr="<<mtr<<", (L,targetDot,uptoNup)=("<<L<<","<<targetDot<<","<<uptoNup<<").\n";
// NEW CODE
        rpL=RNG->getRepnumList(sa,sb,sc,uptoNup);
        getAllQs(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,allWeights);
        jfDenomqs=theGritMonDenom->evaluate(allQs);
// NEW CODE
    cout<<"jfDenomqs=";jfDenomqs->printstdoutWithTruncationOrder();cout<<"\n";
    }
    if(0) if(currUpto<uptoNup){
        if(VERBOSITY>4)cout<<"currUpto<uptoNup, making longer jfDenomqs.\n";
        delete jfDenomqs;
        rpL=RNG->getRepnumList(sa,sb,sc,uptoNup);
        getAllQs(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,allWeights);
        jfDenomqs=theGritMonDenom->evaluate(allQs);
        //jfDenomqs=RNG->getRepnumList(sa,sb,sc,uptoNup)->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,weight);
        if(VERBOSITY>4){cout<<"jfDenomqs=";jfDenomqs->printstdoutWithTruncationOrder();cout<<"\n";}
    }
    if(VERBOSITY>4)cout<<"Putting together numerator...\n";
    numerqs=theGritMonNumer->evaluate(allQs);
    if(VERBOSITY>1){cout<<"numerqs=";numerqs->printstdoutWithTruncationOrder();cout<<"\n";
        cout<<"jfDenomqs=";jfDenomqs->printstdoutWithTruncationOrder();cout<<"\n";
        cout<<"Dividing by jfDenomqs..\n";}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"numerqs=");numerqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"jfDenomqs=");jfDenomqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    numerqs->divideBy(jfDenomqs);
    if(VERBOSITY>1){cout<<"quotient=";numerqs->printstdoutWithTruncationOrder();cout<<"\n";}
    if(SMRTDIAGNOSE>1){fprintf(OUTFILE,"quotientqs=");numerqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");}
    delete jfDenomqs;
    //delete jfwt4qs;
    //delete numerqs;
    if(VERBOSITY>4){cout<<"returning from doOneCoset function.\n";}
    return numerqs;
}

