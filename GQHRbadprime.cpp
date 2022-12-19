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
#include "CNRepnumList.h"
#include "CNRepnumGroup.h"
#include "EfficientQuadratic.h"
#include "EfficientMultiPoly.h"
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
#include "CNRepnumList.cpp"
#include "CNRepnumGroup.cpp"
#include "EfficientQuadratic.cpp"
#include "EfficientMultiPoly.cpp"
#include "JacobiFormData.cpp"


using namespace std;

  const int VERBOSITY=0; ///Whether to print DIAGNOSE4 messages.
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

  EfficientMultiPoly* theGritQuoNumer, *theGritQuoDenom;
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
void readParametersFilefmpzGQ();
void myfmpz_powm(fmpz_t ans, fmpz_t g, slong e, fmpz_t m);
vector<long long> RL4GH(long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, long long c, int k);
vector<long long> RL3GH(long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, int k);
vector<long long> RL1GH(long long L, long long sa, long long sb, long long sc, int k);
vector<long long> RL2GH(long long L, long long sa, long long sb, long long sc, long long a, int k);

vector<long long>* RLT2List(int ttype, long long L, long long sa, long long sb, long long sc,
                        long long a, long long b, long long c);
MySeriesTruncMod* doOneCosetBActimeslevel(CNRepnumGroup* RNG,
    long long sa, long long sb, long long sc, 
    long long BAa, long long BAb, long long BActimeslevel, 
    long long Ldenom, fmpz_t* LrootPowers, int targetDot);
void freeAllQs(MySeriesTruncMod** jfqs, vector<int> used);
void getAllQsBActimeslevel(MySeriesTruncMod** jfqs, vector<int> used, QZSeriesWH **jfGrits,
CNRepnumList* rpL,
 long long BAa,long long BAb,long long BActimeslevel,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,int weight);
void MyExtendedGCD(int a, int b, int &c, int &d);

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
        cout<<"GritQuo Hecke Restriction for primes||level (i.e. bad primes that divide exactly)\n";
        cout << "Usage requires two or more arguments: GritQuo-file HR-file,  [uptoN], [beginCtr], [endCtr], [printEvery], [saveFileName].\n";
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
    readParametersFilefmpzGQ();


    string f;
    stringstream ffs;
    ffs<<"outputbadprime-"<<GQfileName<<"-"<<infileName;
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


    long long H, a, b, c;
    int  i, j, k;
    int maxHecke=heckeArray[0];
    for(H=1;H<numHecke;H++){
        maxHecke=max(maxHecke, heckeArray[H]);
    }

    int maxUptoDet;
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke/(sa*sc*1.0/level-sb*sb); //20190120
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke; //See notebook for theory
    SigObjectOldStyle* sigObj;
    if(SMRTDIAGNOSE){sigObj = new SigObjectOldStyle(level,maxUptoDet);}  //if diagnostics
    else{sigObj = NULL;} //new SigObjectOldStyle(level,maxUptoDet);

    maxqorder = (level*level+4*maxUptoDet)/(4*level) + uptoN; //20190104 try
    //maxqorder = (level*level+4*maxUptoDet)/(4*level);  //commented out 20211003
    maxqorder = (level*level+4*maxUptoDet)/(4*level) + 1000; //20211003 try
    maxqorder = (level*level+4*maxUptoDet)/(4*level) + 200; //20211003 try
    maxqorder = (level*level+4*maxUptoDet)/(4*level) + 2000; //new theory needed
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
    CNRepnumGroup* RNG=new CNRepnumGroup(level,1);  //mindet is 1 because all Jacobi cusp forms.
    //cout<<"HERE B4.\n";
    MySeriesTruncMod* totqs, *oneqs0, *oneqs, *partialsum;
    MySeriesTruncMod* phifqs=NULL, *origqs=NULL;
    fmpz_t detAfactor, fmpzL, tmpfmpz1, tmpfmpz2;
    fmpz_init(detAfactor);  fmpz_init(fmpzL);
    fmpz_init(tmpfmpz1);  fmpz_init(tmpfmpz2);
    int ctr, totctr, globalctr;
    //cout<<"HERE B5.\n";
    CNRepnumList* RL;
    int useShortcut, upperBound;

    int leveloverp, M, ahat, aa, cc, unused;
    //long long  Linvmodn;


    for(H=0;H<numHecke;H++){
        hecke=heckeArray[H];
        fmpz_set(modn,modnfmpzArray[H]);
        Lroot=LrootArray[H];
        L=hecke; fmpz_set_si(fmpzL, L);
        leveloverp=level/L;
        MyExtendedGCD(L, leveloverp, aa, M);
        cc = -M;
        
     //Check that L divides level simply.
        if((level%L==0) and ((level/L)%L != 0)){
           cout<<"Yes, "<<L<<" divides "<<level<<" simply.\n";
        }else{
           cout<<"Yikes, "<<L<<" does not divide "<<level<<" simply. Aborting.\n";
           exit(0);
        }
        //fmpz_invmod(tmpfmpz2, fmpzL, modn);
        //Linvmodn=fmpz_get_si(tmpfmpz2);
        //cout<<"Inverse of "<<L<<" is "<<Linvmodn<<".\n";
        //cout<<"HERE B6.\n";
        LrootPowers=new fmpz_t[L];
        for(i=0;i<L;i++){fmpz_init(LrootPowers[i]);}
        LrootPowersT2=new fmpz_t[2*L];
        for(i=0;i<2*L;i++){fmpz_init(LrootPowersT2[i]);}

        cout<<"Temporarily setting LrootPowers and LrootPowersT2 for L=1.\n";
        fmpz_set_si(LrootPowers[0],1);
        fmpz_set_si(LrootPowersT2[0],1);
        fmpz_set_si(LrootPowersT2[1],-1);

        cout<<"Making origqs... \n";

         RL=RNG->getRepnumList(sa, sb, sc,L*uptoN);
         if(RL->hasDesiredDot(uptoN)){
           cout<<"Original restriction has desired product. Continuing.\n";
         }else{
           cout<<"Original restriction does not have desired product.  Aborting.\n";exit(0);
         }
         //exit(0); //TESTING

        origqs=doOneCosetBActimeslevel(RNG, sa, sb,
                   sc, 0,0,0, 1, LrootPowers, uptoN);

        cout<<"origqs=";origqs->printstdout(); cout<<"\n";
        fprintf(OUTFILE,"origqs = (");origqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");

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
        totctr=L+L+L*L+L*L; ctr=0; globalctr=0;
        totqs=new MySeriesTruncMod(modn, uptoN*L);


        //TYPE 2
        //No shortcut unless we rewrite input to have c be c/N instead.
        //That would be a major undertaking.
        cout<<"Doing type 2 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;
        //int bega=(int)(-L/2);
        int bega=-L+1;

//cout<<"SKIPPING TYPE 2 for diagnostics purposes.\n"; if(0)
        for(a=bega;a<bega+L;a++){
         RL=RNG->getRepnumList(L*L*sa, L*(sb+a*sa), sc+level*(2*a*sb+a*a*sa),L*uptoN);
         if(!RL->hasDesiredDot(L*uptoN)){
            cout<<"Does not have desired dot using a="<<a<<", skipping.\n";
         }else

          for(c=0;c<L;c++){
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
            if(VERBOSITY>4){
              cout<<"type 2: (a,c)=("<<a<<","<<c<<")\n";
            }
            oneqs=doOneCosetBActimeslevel(RNG, L*L*sa, L*(a*sa + sb),
                   level*(a*a*sa + 2*a*sb) + sc, 
                   0,0,c, L, LrootPowers, uptoN);
            //oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
            myfmpz_powm(detAfactor, fmpzL, weight-3, modn);
            oneqs->scalarMultiplityBy(detAfactor);
                //cout<<"HERE Q1\n";
                //if(useShortcut)oneqs->scalarMultiplityBy(fmpzL);
                totqs->addBy(oneqs);
                //cout<<"HERE Q2\n";
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE2X%d=(",a);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                delete oneqs;
            }
         }
       }

        //TYPE 3
        cout<<"Doing type 3 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        if(sa%L==0){cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;}
        else{cout<<"Shortcut available.\n";useShortcut=1;upperBound=1;}
 
//cout<<"SKIPPING TYPE 3 for diagnostics purposes.\n"; if(0)
        for(b=0;b<L;b++){
         RL=RNG->getRepnumList(
            sa + b*level*M* (-2*sb + b*M*sc),
                L*(sb - b*M*sc), L*L*sc,L*uptoN);
         if(!RL->hasDesiredDot(L*uptoN)){
            cout<<"Does not have desired dot using b="<<b<<", skipping.\n";
         }else{
          //if(VERBOSITY>4){cout<<"starting loop of a from 0 to "<<upperBound<<".\n";}
          for(a=0;a<upperBound;a++){
            cout<<"TOUGH to debug when not all diagnostic messages redirect to file.\n";
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    //cout<<"\n";
                    cout<<" [ type 3 b="<<b<<"]\n";
                }
            if(VERBOSITY>4){
              cout<<"type 3: (b,c)=("<<b<<","<<c<<")\n";
            }
            oneqs=doOneCosetBActimeslevel(RNG, sa + b*level*M* (-2*sb + b*M*sc),
                L*(sb - b*M*sc),
           L*L*sc, a,0,0, L, LrootPowers, uptoN);
            //oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
            myfmpz_powm(detAfactor, fmpzL, weight-3, modn);
            oneqs->scalarMultiplityBy(detAfactor);
                //cout<<"HERE Q1\n";
                if(useShortcut)oneqs->scalarMultiplityBy(fmpzL);
                totqs->addBy(oneqs);
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 3 a="<<a<<",b="<<b<<"]\n";
                    cout.flush();
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE3X%dX%X = (",b,c);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                delete oneqs;
           }
          }}
          //delete RL;
          //if(useShortcut)delete RL;
"SOLUTION IS TO JUST ALLOW MEMORY LEAK?????\n";
"WE can track down the bug later.";
"ACTUALLY YOU MUST NOT DELETE CNRepnumList's that were gotten from CNRepnumGroup!!!";
                //cout<<"HERE F7\n";
       }

        //TYPE 4
      bool skip;
      cout<<"No shortcut available.\n";useShortcut=0;upperBound=L;
      for(a=0;a<L;a++){
        cout<<"Doing type 4 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<".\n";
        globalctr++;
        skip=0;
        if(a!=0){
          RL=RNG->getRepnumList(
              L*L*sa - 2*L*level*sb + level*sc,
                -cc*L*sa + L*sb + 2*cc*level*sb - aa*L*sc,
                cc*cc*sa*level - 2*aa*cc*L*sb*level + aa*L*sc+aa*cc*level*sc
          ,L*uptoN);
          if(!RL->hasDesiredDot(L*uptoN))skip=1;
        }else{
          RL=RNG->getRepnumList( L*sa,L*sb,L*sc ,uptoN);
          if(!RL->hasDesiredDot(uptoN))skip=1;
        }
        if(skip){
            cout<<"Does not have desired dot using a="<<a<<", skipping.\n";
        }else
        {
        if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
            ctr++;
            if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                cout<<"\n";
            }
          if(a!=0){
            MyExtendedGCD(a,L,ahat,unused);  //ahat*a = 1 mod L
            oneqs=doOneCosetBActimeslevel(RNG, 
                L*L*sa - 2*L*level*sb + level*sc,
                -cc*L*sa + L*sb + 2*cc*level*sb - aa*L*sc,
                cc*cc*sa*level - 2*aa*cc*L*sb*level + aa*L*sc+aa*cc*level*sc,
            0,-ahat,0, L, LrootPowers, uptoN);
            myfmpz_powm(detAfactor, fmpzL, weight-3, modn);
          }else{
            oneqs0=doOneCosetBActimeslevel(RNG, 
                L*sa,L*sb,L*sc,
            0,0,0, 1, LrootPowers, uptoN);
            oneqs= oneqs0->expandExpBy(L, uptoN*L); delete oneqs0;
            myfmpz_powm(detAfactor, fmpzL, 2*weight-3, modn);
          }
            //cout<<"DIAGNOSTICS A: detAfactor="<<fmpz_get_si(detAfactor)<<".\n";
            oneqs->scalarMultiplityBy(detAfactor);
            totqs->addBy(oneqs);
            if(PRINTTYPE){fprintf(OUTFILE,"TYPE1 = (");oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}

            delete oneqs;
           }}
        }


        //TYPE 1
        cout<<"Doing type 1 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        partialsum=new MySeriesTruncMod(modn, uptoN*L);
        if(sb%L==0){cout<<"No shortcut available because upper right "<<sb<<" is a multiple of "<<L<<".\n";
            useShortcut=0;upperBound=L;}
        else{cout<<"Shortcut available.\n";useShortcut=1;upperBound=1;}
        for(b=0;b<upperBound;b++)
            for(a=-((int)(L/2));a<-((int)(L/2))+L;a++)
            for(c=-((int)(L/2));c<-((int)(L/2))+L;c++){
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<". partialsum=";partialsum->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<") ("
                        <<(100*(globalctr-beginCtr+1))/(endCtr-beginCtr+1)<<"%)";}
                    cout<<" [ type 1 a="<<a<<",b="<<b<<",c="<<c<<"]\n";
                }
                //cout<<"HERE Q1\n";
            oneqs=doOneCosetBActimeslevel(RNG, sa, sb, sc,
                   a,b,c, L, LrootPowers, uptoN);
                partialsum->addBy(oneqs);
                if((beginCtr>0)&&(endCtr==beginCtr)){
                    fprintf(OUTFILE,"abc={%d,%d,%d};\n",a,b,c);
                }
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE1X%dX%dX%d = (",(a+L)%L,(b+L)%L,(c+L)%L);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                delete oneqs;
            }
        }
        if(useShortcut)partialsum->scalarMultiplityBy(fmpzL);
        myfmpz_powm(detAfactor, fmpzL, -3, modn);
        partialsum->scalarMultiplityBy(detAfactor);
        totqs->addBy(partialsum);
        if(PRINTTYPE){fprintf(OUTFILE,"TYPE1partialsum = (");partialsum->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
        delete partialsum;
//
//
// DO NOT DO A sweeping overall power of L
//        myfmpz_powm(detAfactor, fmpzL, weight-3, modn);
//        totqs->scalarMultiplityBy(detAfactor);



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
        fprintf(OUTFILE, "doSomething;\nseparator=nothing;\n\n");
        fflush(OUTFILE);



        delete totqs;
    }
return 0;

}

void readParametersFilefmpzGQ(){
  FILE* PFILE;
  int x;
  cout<<"Reading GritQuo parameters from "<<GQfileName<<" ...\n";
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
  int indSize;
  vector<int>* ind;
  fmpz_t* coeffs;

  fscanf(PFILE, "%d", &indSize); 
  if(VERBOSITY>1){cout<<"indSize="<<indSize<<"\n";}
  ind = new vector<int>[indSize];;
  coeffs=new fmpz_t[indSize];;
  for(i=0;i<indSize;i++){
    fmpz_init(coeffs[i]);
    fmpz_fread(PFILE, coeffs[i]);
    if(VERBOSITY>1){cout<<"(";fmpz_print(coeffs[i]);cout<<") ";}
    fscanf(PFILE, "%d", &len); 
    if(VERBOSITY>1){cout<<len<<":";}
    for(j=0;j<len;j++){
      fscanf(PFILE, "%d", &x); 
      if(VERBOSITY>1){cout<<x<<" ";}
      ind[i].push_back(x);
      usedJF[x]=1;
    }
    if(VERBOSITY>1){cout<<"\n";}
  }
  theGritQuoNumer = new EfficientMultiPoly(coeffs, ind, indSize);
  denomIsLinear=0;  //FUTURE: Test if denom is one Grit
    if(VERBOSITY>3){cout<<"[A] theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
  for(i=0;i<indSize;i++){
    fmpz_clear(coeffs[i]);
  }
    if(VERBOSITY>3){cout<<"[B] theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
  delete[] coeffs; delete[] ind;

  fscanf(PFILE, "%d", &indSize); 
  if(VERBOSITY>1){cout<<"indSize="<<indSize<<"\n";}
  ind = new vector<int>[indSize];;
  coeffs=new fmpz_t[indSize];;
  for(i=0;i<indSize;i++){
    fmpz_init(coeffs[i]);
    fmpz_fread(PFILE, coeffs[i]);
    if(VERBOSITY>1){cout<<"(";fmpz_print(coeffs[i]);cout<<") ";}
    fscanf(PFILE, "%d", &len); 
    if(VERBOSITY>1){cout<<len<<":";}
    for(j=0;j<len;j++){
      fscanf(PFILE, "%d", &x); 
      if(VERBOSITY>1){cout<<x<<" ";}
      ind[i].push_back(x);
      usedJF[x]=1;
    }
    if(VERBOSITY>1){cout<<"\n";}
  }
  theGritQuoDenom= new EfficientMultiPoly(coeffs, ind, indSize);
  for(i=0;i<indSize;i++){
    fmpz_clear(coeffs[i]);
  }
    if(VERBOSITY>3){cout<<"theGritQuoDenom=";theGritQuoDenom->printRepresentation();cout<<"\n";}
  delete[] coeffs; delete[] ind;
    if(VERBOSITY>0){cout<<"theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>0){cout<<"theGritQuoDenom=";theGritQuoDenom->printRepresentation();cout<<"\n";}
    //cout<<"TEST ENDS HERE\n";exit(0);
}

void readParametersFilefmpzHR(){
  FILE* PFILE;
  cout<<"Reading Hecke Restriction parameters from "<<infileName<<" ...\n";
  PFILE = fopen(infileName.c_str(),"r");
  if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
  int x, i, j, len, numer, denom;
  char s[256];
  char charvar;
  fscanf(PFILE, "%d", &weight);
  fscanf(PFILE, "%d", &level);



    fscanf(PFILE, "%d", &x);  sa=x;
    fscanf(PFILE, "%d", &x);  sb=x;
    //fscanf(PFILE, "%d", &x);  sc=x;
    fscanf(PFILE, "%s", s);
    numer=-999; denom=-999;
    sscanf(s, "%d%c%d", &numer, &charvar, &denom);
    cout<<"HERE ND: "<<numer<<","<<charvar<<","<<denom<<"\n";
    if(denom==-999){
      sscanf(s, "%d", &x); sc = x*level;
    }else{
      sc = (level*numer)/denom;
    }
    
    fscanf(PFILE, "%d", &x);  if(uptoN==0){uptoN=x;}

    //fscanf(PFILE, "%s", HRfileCode);

    fscanf(PFILE, "%d", &x);  numHecke=x;


    cout<<"weight = "<<weight<<".\n";
    cout<<"level = "<<level<<".\n";

  cout<<"restriction matrix = "<<sa<<","<<sb<<","<<sc<<"/"<<level<<".\n";
  if(sc*sa-level*sb*sb<=0){cout<<"Restriction matrix has negative determinant. Aborting.\n";exit(1);}
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
  //exit(0); //TESTING
}


void freeAllQs(MySeriesTruncMod** jfqs, vector<int> used){
  int i;
  for(i=0;i<used.size();i++){
    if(used[i])if(jfqs[i]!=NULL){delete jfqs[i]; jfqs[i]=NULL;}
  }
}
void getAllQsBActimeslevel(MySeriesTruncMod** jfqs, vector<int> used,QZSeriesWH **jfGrits,
 CNRepnumList* rpL,
 long long BAa,long long BAb,long long BActimeslevel,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,vector<int> allWeights){
  //This version has lower right corner of translate pre-multiplied by level
  if(VERBOSITY>4)cout<<"Freeing jfqs.\n";
  freeAllQs(jfqs,used);
  if(VERBOSITY>4)cout<<"Freed jfqs.\n";
  int i;
  for(i=0;i<used.size();i++){
    //cout<<i<<"\n";
    if(used[i]){
      //cout<<"HERE G2: "<<used[i]<<","<<allWeights[i]<<"\n";
      jfqs[i]=rpL->restrictGritBActimeslevel(jfGrits[i],BAa,BAb,BActimeslevel,Ldenom,modn,LrootPowers,
                    uptoNup,allWeights[i]);
      if(VERBOSITY>4){cout<<"allQs["<<i<<"]=";allQs[i]->printstdout();cout<<"\n";}
    }
  }
}

void MyExtendedGCD(int a, int b, int &c, int &d){
  if(b==0){c=1;d=0;return;}
  if(abs(b)>abs(a)){return MyExtendedGCD(b,a,d,c);}
  int r=a%b;
  if(r==0){c=0;d=1;return;}
  int x,y;
  MyExtendedGCD(b, r, x, y);
  c=y; d=x-((a-r)/b)*y;
}


MySeriesTruncMod* doOneCosetBActimeslevel(CNRepnumGroup* RNG,
    long long sa, long long sb, long long sc, 
    long long BAa, long long BAb, long long BActimeslevel, 
    long long Ldenom, fmpz_t* LrootPowers, int targetDot){
/**
 * Automatically expands to targetDot*Ldenom!!
 * This version has lower right corner of translate pre-multiplied by level
 * **/
    int mtr=-1;
    int i,k;
    CNRepnumList *rpL;
    //fmpz_t L;  //Just to root out typos.
    if(VERBOSITY>3)cout<<"doOneCosetBActimeslevel...\n";
    if(VERBOSITY>3){cout<<"theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>3){cout<<"theGritQuoDenom=";theGritQuoDenom->printRepresentation();cout<<"\n";}
    int currUpto=Ldenom*targetDot;
    if(VERBOSITY>1)cout<<"doOneCosetBActimeslevel:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BActimeslevel<<"; finding mtr...\n";
    while(mtr<0){
        if(VERBOSITY>4)cout<<"making rpL.\n";
        if(VERBOSITY>4)cout<<"doOneCosetBActimeslevel:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BActimeslevel<<","<<Ldenom<<","<<targetDot<<"\n";
        rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
        if(VERBOSITY>4)cout<<"making jfDenomqs, first getAllQsBActimeslevel.\n";
        if(VERBOSITY>4)cout.flush();
        //jfDenomqs=rpL->
            //restrictGrit(jfDenom,BAa,BAb,BActimeslevel,Ldenom,modn,LrootPowers,currUpto,weight);
        getAllQsBActimeslevel(allQs, usedJF, jfGrits, rpL,
           BAa,BAb,BActimeslevel,Ldenom,modn,LrootPowers,currUpto,allWeights);
        if(VERBOSITY>4){cout<<"making jfDenomqs.\n";cout.flush();}
        jfDenomqs=theGritQuoDenom->evaluate(allQs);
        if(VERBOSITY>4){cout<<"Made a putative jfDenomqs.\n";cout.flush();}
        mtr=jfDenomqs->getVanishingOrder();
        if(VERBOSITY>4){cout<<"mtr="<<mtr<<"\n";cout.flush();}
        if(mtr>=0)break;
        //if(currUpto>=maxqorder){
          //cout<<"Already doubled to the max. Must abort dgs8g8as.\n"; exit(1);
        //}
        currUpto=2*currUpto;  //Keep doubling until we get nonzero restriction.
        //if(currUpto>=maxqorder){currUpto=maxqorder;}
        if(VERBOSITY>4)cout<<"Doubling currUpto (now "<<currUpto<<") to get mtr in doOneCosetBActimeslevel.\n";
        delete jfDenomqs;
    }
    fmpz_t uptoN;
    int uptoNup=targetDot*Ldenom+mtr;
    if(VERBOSITY>4){cout<<"mtr="<<mtr<<", (L,targetDot,uptoNup)=("<<L<<","<<targetDot<<","<<uptoNup<<").\n";
    cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";
    cout.flush();
    }
    if(currUpto<uptoNup){
        if(VERBOSITY>4)cout<<"currUpto<uptoNup, making longer jfDenomqs.\n";
        cout.flush();
        delete jfDenomqs;
        rpL=RNG->getRepnumList(sa,sb,sc,uptoNup);
        getAllQsBActimeslevel(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BActimeslevel,Ldenom,modn,LrootPowers,uptoNup,allWeights);
        jfDenomqs=theGritQuoDenom->evaluate(allQs);

        //jfDenomqs=RNG->getRepnumList(sa,sb,sc,uptoNup)->
            //restrictGrit(jfDenom,BAa,BAb,BActimeslevel,Ldenom,modn,LrootPowers,uptoNup,weight);
        if(VERBOSITY>4){cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";}
        cout.flush();
    }


    if(VERBOSITY>4)cout<<"Putting together numerator...\n";
    numerqs=theGritQuoNumer->evaluate(allQs);
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

    if(VERBOSITY>4){cout<<"returning from doOneCosetBActimeslevel function.\n";}
    return numerqs;
}
