/*
 * 20190120 converting...
 *
 * ***/

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
//#include <mpc.h>
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

  EfficientMultiPoly* theGritQuoNumer, *theGritQuoDenom;
  int denomIsLinear;
  QZSeriesWH **jfGrits;
  QZSeriesWH **allJFs;
  MySeriesTruncMod **allQs;
  int numJFs;
  JacobiFormData** jfData;
  vector<int> usedJF, allWeights;

  int sa, sb, sc; //the restriction matrix
  vector<int> heckeArray, modnArray, LrootArray, L2rootArray;
  fmpz_t* modnfmpzArray;
  int hecke, Lroot, L2root; int numHecke;
  fmpz_t modn, tmpfmpz;
  fmpz_t* LrootPowers, *L2rootPowers;
  string outfileRoot;
  int useFile;
  //string whjfFilename;
  char whjfStr[255];
    FILE* OUTFILE;


void readParametersFilefmpzHRTp2();
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
MySeriesTruncMod* doOneCoset(vector<long long>* ghList, long long L, int targetDot, RepnumGroup* RNG);
MySeriesTruncMod* doOneCoset(RepnumGroup* RNG,
    long long sa, long long sb, long long sc, 
    long long BAa, long long BAb, long long BAc, 
    long long Ldenom, fmpz_t* LrootPowers, int targetDot);

void printFileAbstract(FILE* OUTFILE, SigObjectOldStyle* sigObj, vector<long long>* ghList,
                       long long Ldenom, int targetDot, RepnumGroup* RNG);
void freeAllQs(MySeriesTruncMod** jfqs, vector<int> used);
void getAllQs(MySeriesTruncMod** jfqs, vector<int> used, QZSeriesWH **jfGrits,
MyRepnumList* rpL,
 long long BAa,long long BAb,long long BAc,long long Ldenom,
 fmpz_t modn,
 fmpz_t* LrootPowers, int uptoNup,int weight);

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
        cout << "Usage requires two or more arguments: GritQuo-file HR-file,  [uptoN], [beginCtr], [endCtr], [printEvery], [saveFileName].\n";
        cout<<"uptoN by default is read from the HR file. Use 0 for default.\n";
        cout<<"Runs cosets numbered beginCtr to endCtr if set. Set beginCtr=0 and endCtr=0 if you want to do all cosets.\n";
        cout<<"Output progress every printEvery. Use 0 for default.\n";
        cout<<"saveFileName from file is used by default.\n";
        cout<<"Info about this machine:\n";
        printf("sizeof(slong) = %lu;\n", sizeof(slong));
        printf("sizeof(int) = %lu;\n", sizeof(int));
        printf("sizeof(long long) = %lu;\n", sizeof(long long));

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
    readParametersFilefmpzHRTp2();
    readParametersFilefmpzGQ();


    string f;
    stringstream ffs;
    ffs<<"outputTp2-"<<GQfileName<<"-"<<infileName;
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


    int jfTrunc=0;

    int H, a, b, c, i, j, k;
    int maxHecke=heckeArray[0];
    for(H=1;H<numHecke;H++){
        maxHecke=max(maxHecke, heckeArray[H]);
    }

    int maxUptoDet=123;
    //maxUptoDet = uptoN*uptoN*maxHecke*maxHecke; //See notebook for theory
    maxUptoDet = uptoN*uptoN*maxHecke*maxHecke*maxHecke*maxHecke/(sa*sc-sb*sb); //20190120
    SigObjectOldStyle* sigObj;
    if(SMRTDIAGNOSE){sigObj = new SigObjectOldStyle(level,maxUptoDet);}  //if diagnostics
    else{sigObj = NULL;} //new SigObjectOldStyle(level,maxUptoDet);

    //maxqorder = (level*level+4*maxUptoDet)/(4*level) + 10;
    maxqorder = (level*level+4*maxUptoDet)/(4*level);
    cout<<"maxUptoDet = "<<maxUptoDet<<".\n";
    cout<<"maxqorder (for expanding TBs) = "<<maxqorder<<".\n";
    //cout<<"maxqorder will be calculated exactly.\n";
    jfGrits=new QZSeriesWH*[numJFs];
    jfTrunc=maxqorder+2;    //Do theory later
    jfTrunc=maxqorder+1000;    //Do theory later
    jfTrunc=maxqorder+20;    //Do theory later
    //jfTrunc=maxqorder+200;    //Do theory later
    cout<<"Using jfTrunc="<<jfTrunc<<"\n";
    cout.flush();
    for(i=0;i<numJFs;i++){
            cout<<"Making JF #"<<i<<"\n";
            jfGrits[i]=jfData[i]->makeJF(jfTrunc);
    }


    //cout<<"HERE B1.\n";
    allQs=new MySeriesTruncMod*[numJFs];
    for(i=0;i<numJFs;i++)allQs[i]=NULL;
    //cout<<"HERE B1.\n";

    //cout<<"HERE B3.\n";
    RepnumGroup* RNG=new RepnumGroup(level,1);  //mindet is 1 because all Jacobi cusp forms.
    //cout<<"HERE B4.\n";
    MySeriesTruncMod* totqs, *oneqs, *oneqs0, *partialsum;
    MySeriesTruncMod* phifqs=NULL, *origqs=NULL;
    fmpz_t detAfactor, fmpzL;
    fmpz_init(detAfactor);  fmpz_init(fmpzL);
    int ctr, totctr, globalctr;
    vector<long long>* ghL;
    //cout<<"HERE B5.\n";
    MyRepnumList* RL;
    int useShortcut, upperBound;

    MyRepnumList *RL1, *RL3, *RL4, *RL5, *RL0;
    MyRepnumList **RL2, **RL6;

    for(H=0;H<numHecke;H++){
        hecke=heckeArray[H];
        fmpz_set(modn,modnfmpzArray[H]);
        Lroot=LrootArray[H];
        L2root=L2rootArray[H];
        L=hecke; fmpz_set_si(fmpzL, L);
    //cout<<"HERE B6.\n";
        LrootPowers=new fmpz_t[L];
        for(i=0;i<L;i++){fmpz_init(LrootPowers[i]);}
        L2rootPowers=new fmpz_t[L*L];
        for(i=0;i<L*L;i++){fmpz_init(L2rootPowers[i]);}

    //Calculate RLs.

        cout<<"Temporarily setting LrootPowers and L2rootPowers for L=1.\n";
        fmpz_set_si(LrootPowers[0],1);
        fmpz_set_si(L2rootPowers[0],1);

        cout<<"Getting ghL to make oriqqs... \n";
        ghL=RLT2List(4, 1, sa, sb, sc, 0, 0, 0);
        cout<<"Making origqs... \n";

/**THIS NEEDS TO BE UPDATED**/
        origqs = doOneCoset(ghL,1, uptoN, RNG);

        cout<<"origqs=";origqs->printstdout(); cout<<"\n";
        fprintf(OUTFILE,"origqs = (");origqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");
        if(SMRTDIAGNOSE){fprintf(OUTFILE,"origqsAbstract=");printFileAbstract(OUTFILE,sigObj,ghL,1,uptoN,RNG);}
        delete[] ghL;

//cout<<"END TEST\n"; return 0;

        fmpz_set_si(tmpfmpz,Lroot);
        for(i=0;i<L;i++){
            myfmpz_powm(LrootPowers[i], tmpfmpz, i, modn);  
        }
    //cout<<"HERE B7.\n";
        fmpz_set_si(tmpfmpz,L2root);
        for(i=0;i<L*L;i++){
            myfmpz_powm(L2rootPowers[i], tmpfmpz, i, modn);
        }
    //cout<<"HERE B8.\n";




        cout<<"**************** L="<<L<<" ***********************\n";
        totctr=1+1+L+L*L; ctr=0; globalctr=0;
        totctr=L+L*L;
	if(sa%L==0){totctr+=L*L*L;}else{totctr+=L;}
	for(i=0;i<L;i++){
	  if((level*(sc+2*i*sb+i*i*sa))%L==0){totctr+=L*L*L;}
	  else{totctr+=L;}
	}
	cout<<"NUMBER OF COSETS="<<totctr<<"\n";
        totqs=new MySeriesTruncMod(modn, uptoN*L);

        //TYPE 1 Tp2
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
            oneqs0=doOneCoset(RNG, sa, sb*L, sc*L*L, 0,0,0, 1, LrootPowers, uptoN);
            oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
            myfmpz_powm(detAfactor, fmpzL, 3*weight-6, modn);

            //cout<<"DIAGNOSTICS A: detAfactor="<<fmpz_get_si(detAfactor)<<".\n";
            oneqs->scalarMultiplityBy(detAfactor);
            totqs->addBy(oneqs);
            if(PRINTTYPE){fprintf(OUTFILE,"TYPE1 = (");oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}

            delete oneqs;
        }
	
    int skipped;
        //TYPE 2 Tp2
        cout<<"Doing type 2 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        skipped=0;
        for(a=0;a<L;a++){
         RL=RNG->getRepnumList(L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa,L*L*uptoN);
         if(RL->hasDesiredDot(uptoN)){

         if(VERBOSITY>4)cout<<"RNG->getRepnumList("<<L*L*sa<<","<< L*(sb+a*sa)<<","<< sc+2*a*sb+a*a*sa<<","<<L*uptoN<<");\n";
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 2 a="<<a<<"]\n";
                    cout.flush();
                }
		oneqs0=doOneCoset(RNG, L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa,
                   0,0,0, 1, LrootPowers, uptoN);
               oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
             myfmpz_powm(detAfactor, fmpzL, 3*weight-6, modn);

		oneqs->scalarMultiplityBy(detAfactor);
           totqs->addBy(oneqs);
                if(PRINTTYPE){fprintf(OUTFILE,"TYPE2X%dX = (",a);oneqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");}
                delete oneqs;
          }
        }else{
          skipped++;
	}
       }
	  if(VERBOSITY>2){
            cout<<"Skipped "<<skipped<<" becasue hasDesiredDot was false.\n";
          }

        //TYPE 3 Tp2
        cout<<"Doing type 3 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=1;a<L;a++){  //a!=0
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 3 a="<<a<<"]\n";
                    cout.flush();
                }
                oneqs0=doOneCoset(RNG, L*sa, L*sb, L*sc, a,0,0, L, LrootPowers, uptoN);
                oneqs= oneqs0->expandExpBy(L, uptoN*L*L); delete oneqs0;
                myfmpz_powm(detAfactor, fmpzL, 2*weight-6, modn);
                oneqs->scalarMultiplityBy(detAfactor);
                totqs->addBy(oneqs);
                delete oneqs;
            }
       }


        //TYPE 4 TP2
        cout<<"Doing type 4 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=1;a<L;a++)for(b=0;b<L;b++){  //a!=0
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 4 a="<<a<<",b="<<b<<"]\n";
                    cout.flush();
                }
                oneqs0=doOneCoset(RNG,L*sa,L*sb,L*sc,a*b*b,a*b,a, L, LrootPowers, uptoN);
                oneqs= oneqs0->expandExpBy(L, uptoN*L*L); delete oneqs0;
                myfmpz_powm(detAfactor, fmpzL, 2*weight-6, modn);
                oneqs->scalarMultiplityBy(detAfactor);
                totqs->addBy(oneqs);
                delete oneqs;
            }
       }

        
        //TYPE 5 Tp2
        cout<<"Doing type 5 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        partialsum=new MySeriesTruncMod(modn, uptoN*L*L);
        int begI, endI;
        if((sa%L)==0){
                begI=-((int)((L*L)/2));endI=begI+L*L;
            }else{
                cout<<"Can employ speed up #1!.\n";
                begI=0;endI=1;
            }
         for(a=0;a<L;a++){
           for(b=begI;b<endI;b++){
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<". partialsum=";partialsum->printstdout();
                    //if(!doAll){cout<<" (endCtr="<<endCtr<<") ("
                        //<<(100*(globalctr-beginCtr+1))/(endCtr-beginCtr+1)<<"%)";}
                    //cout<<"\n";
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 5 a="<<a<<",b="<<b<<"]\n";
                    cout.flush();
                }
                oneqs=doOneCoset(RNG, sa, L*sb, L*L*sc, b,L*a,0, L*L, L2rootPowers, uptoN);
                //No expand
                partialsum->addBy(oneqs);
                if((beginCtr>0)&&(endCtr==beginCtr)){
                    fprintf(OUTFILE,"abc={%d,%d};\n",a,b);
                }
                delete oneqs;
            }
        }
        }
        if((sa%L)==0){
            myfmpz_powm(detAfactor, fmpzL, weight-6, modn);
        }else{
            myfmpz_powm(detAfactor, fmpzL, weight-6+2, modn);
        }
        partialsum->scalarMultiplityBy(detAfactor);
        totqs->addBy(partialsum);
        delete partialsum;

        //TYPE 6 Tp2
        cout<<"Doing type 6 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        skipped=0;
        partialsum=new MySeriesTruncMod(modn, uptoN*L*L);
        for(a=0;a<L;a++){
        RL=RNG->getRepnumList(L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa,L*L*uptoN);
        if(!(RL->hasDesiredDot(L*L*uptoN))){skipped++;}
        else{
         if(((sc+2*a*sb+a*a*sa)%L)==0){
                begI=-((int)((L*L)/2));endI=begI+L*L;
            }else{
                cout<<"Can employ speed up #1!.\n";
                begI=0;endI=1;
            }
          for(b=0;b<L;b++)for(c=begI;c<endI;c++){
            globalctr++;

           if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<" [ type 6 a="<<a<<",b="<<b<<",c="<<c<<"]\n";
                    cout.flush();
                }
                oneqs=doOneCoset(RNG, L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa,
                  0,L*b, a*b+c, L*L, L2rootPowers, uptoN);
               //No expand
                if(((sc+2*a*sb+a*a*sa)%L)==0){
                }else{
                    myfmpz_powm(detAfactor, fmpzL, 2, modn);
                    oneqs->scalarMultiplityBy(detAfactor);
               }
                partialsum->addBy(oneqs);
                if((beginCtr>0)&&(endCtr==beginCtr)){
                    fprintf(OUTFILE,"abc={%d,%d};\n",a,b);
                }
                delete oneqs;
            }

          }

          }
       }

         myfmpz_powm(detAfactor, fmpzL, weight-6, modn);
        partialsum->scalarMultiplityBy(detAfactor);
        totqs->addBy(partialsum);
        delete partialsum;

        cout<<"global counter="<<globalctr<<"\n";
        cout<<"ctr counter="<<ctr<<"\n";

        cout<<"Made totqs = ";totqs->printstdout();cout<<"\n";

        fprintf(OUTFILE, "L=%d;\n", L);
        fprintf(OUTFILE, "doAll=%d;\n", doAll);
        fprintf(OUTFILE, "beginCtr=%d;\n", beginCtr);
        fprintf(OUTFILE, "endCtr=%d;\n", endCtr);
        fprintf(OUTFILE, "modn=");
        fmpz_fprint(OUTFILE,modn);
        fprintf(OUTFILE, ";\n");
        fprintf(OUTFILE, "L2root=%d;\n", L2root);
        fprintf(OUTFILE, "Lroot=%d;\n", Lroot);
        fprintf(OUTFILE, "numCosets=%d;\n", ctr);
        fprintf(OUTFILE, "sizeofslong = %lu;\n", sizeof(slong));
        fprintf(OUTFILE, "sizeofint = %lu;\n", sizeof(int));
        fprintf(OUTFILE, "sizeoflonglong = %lu;\n", sizeof(long long));

        fprintf(OUTFILE,"totqs = (");totqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");
        fprintf(OUTFILE, "doSomething;\nseparator=nothing;\n\n");
        fflush(OUTFILE);




        delete totqs;
    }
return 0;


return 0;






    fmpz_clear(modn); fmpz_clear(detAfactor); fmpz_clear(fmpzL); fmpz_clear(tmpfmpz);
    for(i=0;i<numHecke;i++){
            fmpz_clear(modnfmpzArray[i]);
    }
    delete[] modnfmpzArray;

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

void readParametersFilefmpzHRTp2(){
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
        fscanf(PFILE, "%d", &x);   L2rootArray.push_back(x);
        cout<<"L2root = "<<x<<", ";
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
  if(VERBOSITY>4)cout<<"Freeing jfqs.\n";
  freeAllQs(jfqs,used);
  if(VERBOSITY>4)cout<<"Freed jfqs.\n";
  int i;
  for(i=0;i<used.size();i++){
    //cout<<i<<"\n";
    if(used[i]){
      //cout<<"HERE G2: "<<used[i]<<","<<allWeights[i]<<"\n";
      jfqs[i]=rpL->restrictGrit(jfGrits[i],BAa,BAb,BAc,Ldenom,modn,LrootPowers,
                    uptoNup,allWeights[i]);
      if(VERBOSITY>4){cout<<"allQs["<<i<<"]=";allQs[i]->printstdout();cout<<"\n";}
    }
  }
}
MySeriesTruncMod* doOneCoset(vector<long long>* ghList, long long Ldenom, int targetDot, RepnumGroup* RNG){
/**
 * Automatically expands to targetDot*Ldenom!!
 * **/
    int mtr=-1;
    int i,k;
    MyRepnumList *rpL;
    //fmpz_t L;  //Just to root out typos.
    if(VERBOSITY>3)cout<<"doOneCoset...\n";
    if(VERBOSITY>3){cout<<"theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>3){cout<<"theGritQuoDenom=";theGritQuoDenom->printRepresentation();cout<<"\n";}
    long long sa=ghList[0][0], sb=ghList[0][1], sc=ghList[0][2],
        BAa=ghList[0][3],BAb=ghList[0][4],BAc=ghList[0][5];
    int currUpto=Ldenom*targetDot;
    if(VERBOSITY>1)cout<<"doOneCoset:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BAc<<"; finding mtr...\n";
    while(mtr<0){
        if(VERBOSITY>4)cout<<"making rpL.\n";
        rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
        if(VERBOSITY>4)cout<<"making jfDenomqs, first getAllQs.\n";
        if(VERBOSITY>4)cout.flush();
        //jfDenomqs=rpL->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,weight);
        getAllQs(allQs, usedJF, jfGrits, rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,allWeights);
        if(VERBOSITY>4)cout<<"making jfDenomqs.\n";
        jfDenomqs=theGritQuoDenom->evaluate(allQs);
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
    int uptoNup=targetDot*Ldenom+mtr;
    if(VERBOSITY>4){cout<<"mtr="<<mtr<<", (L,targetDot,uptoNup)=("<<L<<","<<targetDot<<","<<uptoNup<<").\n";
    cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";
    }
    if(currUpto<uptoNup){
        if(VERBOSITY>4)cout<<"currUpto<uptoNup, making longer jfDenomqs.\n";
        delete jfDenomqs;
        rpL=RNG->getRepnumList(sa,sb,sc,uptoNup);
        getAllQs(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,allWeights);
        jfDenomqs=theGritQuoDenom->evaluate(allQs);

        //jfDenomqs=RNG->getRepnumList(sa,sb,sc,uptoNup)->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,weight);
        if(VERBOSITY>4){cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";}
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

    if(VERBOSITY>4){cout<<"returning from doOneCoset function.\n";}
    return numerqs;
}

MySeriesTruncMod* doOneCoset(RepnumGroup* RNG,
    long long sa, long long sb, long long sc, 
    long long BAa, long long BAb, long long BAc, 
    long long Ldenom, fmpz_t* LrootPowers, int targetDot){
/**
 * Automatically expands to targetDot*Ldenom!!
 * **/
    int mtr=-1;
    int i,k;
    MyRepnumList *rpL;
    //fmpz_t L;  //Just to root out typos.
    if(VERBOSITY>3)cout<<"doOneCoset...\n";
    if(VERBOSITY>3){cout<<"theGritQuoNumer=";theGritQuoNumer->printRepresentation();cout<<"\n";}
    if(VERBOSITY>3){cout<<"theGritQuoDenom=";theGritQuoDenom->printRepresentation();cout<<"\n";}
    int currUpto=Ldenom*targetDot;
    if(VERBOSITY>1)cout<<"doOneCoset:"<<sa<<","<<sb<<","<<sc<<","<<BAa<<","<<BAb<<","<<BAc<<"; finding mtr...\n";
    while(mtr<0){
        if(VERBOSITY>4)cout<<"making rpL.\n";
        rpL=RNG->getRepnumList(sa,sb,sc,currUpto);
        if(VERBOSITY>4)cout<<"making jfDenomqs, first getAllQs.\n";
        if(VERBOSITY>4)cout.flush();
        //jfDenomqs=rpL->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,weight);
        getAllQs(allQs, usedJF, jfGrits, rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,currUpto,allWeights);
        if(VERBOSITY>4)cout<<"making jfDenomqs.\n";
        jfDenomqs=theGritQuoDenom->evaluate(allQs);
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
    int uptoNup=targetDot*Ldenom+mtr;
    if(VERBOSITY>4){cout<<"mtr="<<mtr<<", (L,targetDot,uptoNup)=("<<L<<","<<targetDot<<","<<uptoNup<<").\n";
    cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";
    }
    if(currUpto<uptoNup){
        if(VERBOSITY>4)cout<<"currUpto<uptoNup, making longer jfDenomqs.\n";
        delete jfDenomqs;
        rpL=RNG->getRepnumList(sa,sb,sc,uptoNup);
        getAllQs(allQs, usedJF,jfGrits,  rpL,
           BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,allWeights);
        jfDenomqs=theGritQuoDenom->evaluate(allQs);

        //jfDenomqs=RNG->getRepnumList(sa,sb,sc,uptoNup)->
            //restrictGrit(jfDenom,BAa,BAb,BAc,Ldenom,modn,LrootPowers,uptoNup,weight);
        if(VERBOSITY>4){cout<<"jfDenomqs=";jfDenomqs->printstdout();cout<<"\n";}
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

    if(VERBOSITY>4){cout<<"returning from doOneCoset function.\n";}
    return numerqs;
}
