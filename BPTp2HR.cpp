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
//#include <mpfr.h>
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
//#include "BPdata.h"
//#include "FindBPunified.h"
//#include "ThetaBlockList.h"
//#include "BPdataList.h"
#include "SigObjectOldStyle.h"
//#include "BPcheckCuspiness.h"
#include "SingularPartList.h"
//#include "WeightTwelveBases.h"
#include "QZWSeriesQ.h"
#include "MySeriesTruncMod.h"
#include "MyRepnumList.h"


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
//#include "ThetaBlockList.cpp"
//#include "BPdataList.cpp"
#include "SigObjectOldStyle.cpp"
//#include "BPcheckCuspiness.cpp"
//#include "SingularPartList.cpp"
//#include "WeightTwelveBases.cpp"
#include "QZWSeriesQ.cpp"
#include "MySeriesTruncMod.cpp"
#include "MyRepnumList.cpp"


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
  int level, weight, uptoN, L, mindet;
  string saveFileName; stringstream saveFileNameStream; int saveFileQ;
  int maxqorder, maxd;
  vector<int> TBphi, TBTheta; int tVal, ell;
  int sa, sb, sc; //the restriction matrix
  vector<int> heckeArray, modnArray, L2rootArray,LrootArray;
  fmpz_t* modnfmpzArray;
  int hecke, Lroot, L2root; int numHecke;
  fmpz_t modn;
  string outfileRoot;
  int useFile;
  //string whjfFilename;
  char whjfStr[255];

void readParametersFilefmpz();
void myfmpz_powm(fmpz_t ans, fmpz_t g, slong e, fmpz_t m);



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
    if(0){//testcode
        FILE* f = fopen("testdata-3.txt", "r");
        QZSeriesWH* qqq=new QZSeriesWH(f);
        fclose(f);
        cout<<qqq->getString()<<"\n";
        qqq->multiplyWith(qqq);
        FILE* g = fopen("testdata-4.txt", "w");
        qqq->saveToFile(g);
        return(0);
    }

    if(argc<2){
        cout << "Usage requires one or more arguments: filename,  [uptoN], [beginCtr], [endCtr], [printEvery], [saveFileName].\n";
        cout<<"uptoN by default is read from the file. Use 0 for default.\n";
        cout<<"Runs cosets numbered beginCtr to endCtr if set. Set beginCtr=0 and endCtr=0 if you want to do all cosets.\n";
        cout<<"Output progress every printEvery. Use 0 for default.\n";
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
    if(argc>4){
     beginCtr = atoi(argv[3]);
     endCtr = atoi(argv[4]);
     if((beginCtr>0)&&(endCtr>0)){doAll=0;}
    }
    if(argc>5){
     printEvery = atoi(argv[5]);
     if(printEvery==0){printEvery=1000;}
    }
    if(argc>6){
     saveFileNameStream<<argv[6];
     saveFileNameStream>>saveFileName;
     saveFileQ=1;
    }

    fmpz_init(modn);
    readParametersFilefmpz();

    FILE* OUTFILE;
    string f;
    stringstream ffs;
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

    int jfTrunc;    //Do theory later
    QZSeriesWH *jf, *jfBuddy, *jfpsi;
    jfpsi=NULL;
    int jfpsiTrunc=0;

    int H, a, b, c;
    MyRepnumList *RL1, *RL3, *RL4, *RL5, *RL0;
    MyRepnumList **RL2, **RL6;
    MySeriesTruncMod* phifqs=NULL;
    MySeriesTruncMod* totqs, *oneqs, *oneqs0, *partialsum;
    fmpz_t detAfactor, fmpzL;
    fmpz_init(detAfactor);  fmpz_init(fmpzL);
    int ctr, totctr, globalctr;
    int *useRL2, *useRL6;


    if(useFile){
        //FILE *ff = fopen(whjfFilename.c_str(),"r");
        FILE *ff = fopen(whjfStr,"r");
        //if(!ff){cout<<"File "<<whjfFilename<<" open failure. Abort.\n"; exit(1);}
        if(!ff){cout<<"File "<<whjfStr<<" open failure. Abort.\n"; exit(1);}
        jfpsi=new QZSeriesWH(ff);
        jfpsiTrunc=jfpsi->getMaxExponent();
    }
    for(H=0;H<numHecke;H++){
        hecke=heckeArray[H];
        //fmpz_set_si(modn,modnArray[H]);

        fmpz_set(modn,modnfmpzArray[H]);

        L2root=L2rootArray[H];
        Lroot=LrootArray[H];
        L=hecke;
        cout<<"**************** L="<<L<<" ***********************\n";
        fmpz_set_si(fmpzL,L);
        mindet=-1*level; //Can improve this later

//cout<<"HERE1:"<<sa<<","<<sb<<","<<sc<<"\n";

        RL0 = new MyRepnumList(sa, sb, sc, level, uptoN, mindet);
        maxqorder=RL0->maxq;
        RL5 = new MyRepnumList(sa, sb*L, sc*L*L, level, uptoN*L*L, mindet);
        RL1=RL5;
        maxqorder=max(maxqorder,RL1->maxq);

        RL4 = new MyRepnumList(L*sa, L*sb, L*sc, level, uptoN*L, mindet);
        RL3=RL4;
        maxqorder=max(maxqorder,RL4->maxq);
        RL2 = new MyRepnumList*[L];
        RL6 = new MyRepnumList*[L];
        useRL2 = new int[L];
        useRL6 = new int[L];

        for(a=0;a<L;a++){
            //RL2[a] = new MyRepnumList(L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa, level, uptoN, mindet);
            RL6[a] = new MyRepnumList(L*L*sa, L*(sb+a*sa), sc+2*a*sb+a*a*sa, level, uptoN*L*L, mindet);
            RL2[a]=RL6[a];
            if(RL2[a]->hasDesiredDot(uptoN)){
                useRL2[a]=1;
                maxqorder=max(maxqorder,RL2[a]->maxq);
            }else{
                useRL2[a]=0;
            }
            if(RL6[a]->hasDesiredDot(uptoN*L*L)){
                useRL6[a]=1;
                maxqorder=max(maxqorder,RL6[a]->maxq);
            }else{
                useRL6[a]=0;
            }
            //cout<<"RL3["<<a<<"]: Num Ts = "<<RL3[a]->n.size()<<"; Maxq = "<<RL3[a]->maxq<<"\n";
        }
        cout<<"Requires maxqorder = "<<maxqorder<<".\n";

        if(!useFile)
        if((jfpsi==NULL)||(jfpsiTrunc<maxqorder)){
            if(jfpsi!=NULL){delete jfpsi;}
            //TEST!!!!   maxqorder=50;  TEST if too low of a trunc is caught.  Test passed.
            jfTrunc=maxqorder+2;
            jfpsiTrunc=maxqorder;
            if(DIAGNOSE1)cout<<"Making TBphi (with target maxqorder = "<<maxqorder<<")...";
            jf = QZSeriesWH::ThetaBlock(jfTrunc*2, weight, TBphi);
            jfpsi = jf->copy();
            if(DIAGNOSE1)cout<<"Applying V2...";
            jfpsi->applyUp(level, weight, 2);
            if(DIAGNOSE1)cout<<"Making Theta Buddy...";
            jfBuddy=QZSeriesWH::ThetaBlock(jfTrunc, weight, TBTheta);
            if(DIAGNOSE1)cout<<"subtracting theta buddy...";
            jfpsi->subtractWith(jfBuddy);
            if(DIAGNOSE1)cout<<"Divide by TBphi...\n";
            jfpsi->divideBy(jf);
            if(DIAGNOSE1)cout<<"Done making weakly holomorphic weight 0 psi.\n";
            delete jf; delete jfBuddy;
        }

        RL0->prepEfficientSpecial(L, jfpsi, Lroot, modn);
        //RL1->prepEfficientSpecial(L, jfpsi, Lroot, modn);
        //RL3->prepEfficientSpecial(L, jfpsi, Lroot, modn);
        RL4->prepEfficientSpecial(L, jfpsi, Lroot, modn);
        RL5->prepEfficientSpecialGeneralPowersLen(L, jfpsi, L2root, modn,L*L);
        //cout<<"After prepEfficient:\n";
        //cout<<"RL4: Num Ts = "<<RL4->nmodL.size()<<"; Maxq = "<<RL4->maxqEff<<"\n";
        //cout<<"RL2: Num Ts = "<<RL2->nmodL.size()<<"; Maxq = "<<RL2->maxqEff<<"\n";

        if(0){
            fprintf(OUTFILE, "LrootPowers={");
            for(int i=0;i<RL4->LrootPowersLen;i++){
                if(i>0){fprintf(OUTFILE,",");}
                fmpz_fprint(OUTFILE,RL4->LrootPowers[i]);
            }
            fprintf(OUTFILE, "};\n");
            RL4->fprintfDiagnostics(OUTFILE, uptoN*L);
            return 0;
        }


        //for(a=0;a<L;a++)
        //    if(useRL2[a]){
        //    RL2[a]->prepEfficientSpecial(L, jfpsi, Lroot, modn);
            //cout<<"RL3["<<a<<"]: Num Ts = "<<RL3[a]->nmodL.size()<<"; Maxq = "<<RL3[a]->maxqEff<<"\n";
        //}
        for(a=0;a<L;a++)
            if(useRL6[a]||useRL2[a])
                {
            RL6[a]->prepEfficientSpecialGeneralPowersLen(L, jfpsi, L2root, modn, L*L);
            //cout<<"RL3["<<a<<"]: Num Ts = "<<RL3[a]->nmodL.size()<<"; Maxq = "<<RL3[a]->maxqEff<<"\n";
        }

        if(phifqs==NULL){
            cout<<"Now doing original restriction...\n";
            phifqs = RL0->oneBPCosetRestrictionPrepped(0,0,0, 1, modn, 1, uptoN+fmpz_get_si(RL0->singExp));
            cout<<"Made phifqs = ";phifqs->printstdout();cout<<"\n";
            fprintf(OUTFILE, "L=%d;\n", L);
            fprintf(OUTFILE, "modn=%ld;\n", fmpz_get_si(modn));
            fprintf(OUTFILE, "Lroot=%d;\n", Lroot);
            fprintf(OUTFILE,"phifqs = ");phifqs->printFile(OUTFILE);fprintf(OUTFILE,";\n");
            fprintf(OUTFILE, "domeSomethingphifqs;\nseparator=nothing;\n\n");

        }

        //totctr=1+L+L*L+L*L*L; ctr=0; globalctr=0;
        totctr=L+L*L+L*L*L+L*L*L*L; ctr=0; globalctr=0;
        totqs=new MySeriesTruncMod(modn, uptoN*L*L);

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
            oneqs0=RL1->oneBPCosetRestrictionPrepped(0,0,0, 1, modn, 1, uptoN);
            oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
            myfmpz_powm(detAfactor, fmpzL, 3*weight-6, modn);
            //cout<<"DIAGNOSTICS A: detAfactor="<<fmpz_get_si(detAfactor)<<".\n";
            oneqs->scalarMultiplityBy(detAfactor);
            totqs->addBy(oneqs);
            delete oneqs;
        }

        //TYPE 2
        cout<<"Doing type 2 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=0;a<L;a++)if(useRL2[a]){
         globalctr++;

            {
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<"\n";
                }
                oneqs0=RL2[a]->oneBPCosetRestrictionPrepped(0,0,0, 1, modn, 1, uptoN);
               oneqs= oneqs0->expandExpBy(L*L, uptoN*L*L); delete oneqs0;
             myfmpz_powm(detAfactor, fmpzL, 3*weight-6, modn);
                //detAfactor is 1 for Type 3
                //cout<<"DIAGNOSTICS A: type 3: ";oneqs->printstdout();cout<<".\n";
              oneqs->scalarMultiplityBy(detAfactor);
           totqs->addBy(oneqs);
          delete oneqs;
        }
       }
        }
        //TYPE 3
        cout<<"Doing type 3 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=1;a<L;a++){  //a!=0
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<"\n";
                }
                oneqs0=RL3->oneBPCosetRestrictionPrepped(a,0,0, L, modn, Lroot, uptoN*L);
                oneqs= oneqs0->expandExpBy(L, uptoN*L*L); delete oneqs0;
                myfmpz_powm(detAfactor, fmpzL, 2*weight-6, modn);
                //detAfactor is 1 for Type 2
                oneqs->scalarMultiplityBy(detAfactor);
                totqs->addBy(oneqs);
                delete oneqs;
            }
       }
        //TYPE 4
        cout<<"Doing type 4 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        for(a=1;a<L;a++)for(b=0;b<L;b++){  //a!=0
            globalctr++;
            if(doAll||((globalctr>=beginCtr)&&(globalctr<=endCtr))){
                ctr++;
                if((ctr<10)||((ctr<100)&&(ctr%10==0))||((ctr<1000)&&(ctr%100==0))||(ctr%printEvery==0)){
                    cout<<globalctr<<" of "<<totctr;
                    //cout<<<<". totqs=";totqs->printstdout();
                    if(!doAll){cout<<" (endCtr="<<endCtr<<")";}
                    cout<<"\n";
                }
                oneqs0=RL4->oneBPCosetRestrictionPrepped(a*b*b,a*b,a, L, modn, Lroot, uptoN*L);
                oneqs= oneqs0->expandExpBy(L, uptoN*L*L); delete oneqs0;
                myfmpz_powm(detAfactor, fmpzL, 2*weight-6, modn);
                //detAfactor is 1 for Type 2
                oneqs->scalarMultiplityBy(detAfactor);
                totqs->addBy(oneqs);
                delete oneqs;
            }
       }

        //TYPE 5
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
                    if(!doAll){cout<<" (endCtr="<<endCtr<<") ("
                        <<(100*(globalctr-beginCtr+1))/(endCtr-beginCtr+1)<<"%)";}
                    cout<<"\n";
                }
                oneqs=RL5->oneBPCosetRestrictionPrepped(b,L*a,0, L*L, modn, L2root, uptoN*L*L);
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

        //TYPE 6
        cout<<"Doing type 6 (globalctr="<<globalctr<<"). totqs=";totqs->printstdout();cout<<"\n";
        partialsum=new MySeriesTruncMod(modn, uptoN*L*L);
        for(a=0;a<L;a++)if(useRL6[a]){
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
                    cout<<"\n";
                }
                oneqs=RL6[a]->oneBPCosetRestrictionPrepped(0,L*b, a*b+c, L*L, modn, L2root, uptoN*L*L);
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
        fprintf(OUTFILE, "sizeofslong = %lu;\n", sizeof(slong));
        fprintf(OUTFILE, "sizeofint = %lu;\n", sizeof(int));
        fprintf(OUTFILE, "sizeoflonglong = %lu;\n", sizeof(long long));

        fprintf(OUTFILE,"totqs = (");totqs->printFile(OUTFILE);fprintf(OUTFILE,");\n");
        fprintf(OUTFILE, "domeSomething;\nseparator=nothing;\n\n");
        fflush(OUTFILE);

        for(a=0;a<L;a++){
                //cout<<"HEREA.";
            delete RL6[a];
                //cout<<"HEREB.";
        }
        delete[] RL6; //cout<<"HEREC.";
        delete[] useRL6;
        delete[] useRL2;
        delete RL0;
        delete RL4; //cout<<"HERED.";
        delete RL5; //cout<<"HEREE.";
        delete totqs; totqs=NULL;
        delete phifqs; phifqs=NULL;
    }

    //MySeriesTruncMod* phifqs = MySeriesTruncMod::doOneCosetRestriction(
    //   jfpsi, TBs, RL4, 0,0,0, 1, LthRoots, uptoN);

    if(0){
        fmpz_t testn; fmpz_init_set_si(testn,127);
        fmpz_t c; fmpz_init_set_si(c,-1);
        cout<<"Made constants.\n";
        MySeriesTruncMod* f = new MySeriesTruncMod(127, 17);
        cout<<"Made f.\n";
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing1.\n";
        fmpz_mod_poly_set_coeff_fmpz(f->p,15,c);
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing1b.\n";
        f->setMonicBinomialBPSpecial(c, 15);
        cout<<"setMonicBinomialBPSpecial(c, 2):";
        f->printstdout("y"); cout<<"\n";
        f->setMonicBinomialBPSpecial(c, -2);
        cout<<"setMonicBinomialBPSpecial(c, -2):";
        f->printstdout("y"); cout<<"\n";

        return 0;
    }
    if(0){
        int H=7;
        MyRepnumList* RL0 = new MyRepnumList(sa, sb, sc, level, 11*2*H, -1*level);
        cout<<"Num Ts = "<<RL0->n.size()<<"; ";
        cout<<"Maxq = "<<RL0->maxq<<"\n";
        //RL0->prepEfficientSpecial(H, jfpsi, 2, 127);
        //cout<<"After prepEfficientSpecial: Num Ts = "<<RL0->nmodL.size()<<"; ";
        //cout<<"Maxq = "<<RL0->maxqEff<<"\n";
    }


    if(0){
        fmpz_t testn; fmpz_init_set_si(testn,11);
        fmpz_t c; fmpz_init_set_si(c,5);
        cout<<"Made constants.\n";
        MySeriesTruncMod* f = MySeriesTruncMod::monicBinomial(testn, 4, c, 1);
        cout<<"Made f.\n";
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing1.\n";
        f->multiplyByPower(f, 2);
        cout<<"Done multiply by f^2.\n";
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing2.\n";
        f->multiplyByPower(f, -4);
        cout<<"Done multiply by f^-4.\n";
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing2.\n";
        f->multiplyByPowerOfQandChangeTrunc(2);
        cout<<"Done shifting up by 2.\n";
        f->printstdout("y"); cout<<"\n";
        cout<<"Done printing2.\n";

    }


    fmpz_clear(modn); fmpz_clear(detAfactor); fmpz_clear(fmpzL);
    int i;
    for(i=0;i<numHecke;i++){
            fmpz_clear(modnfmpzArray[i]);
    }
    delete[] modnfmpzArray;
    return 0;
}

void readParametersFilefmpz(){
  FILE* PFILE;
  cout<<"Reading parameters from "<<infileName<<" ...\n";
  PFILE = fopen(infileName.c_str(),"r");
  if(!PFILE){cout<<" no such parameters file.\n"; exit( 1); }
  int x, i, len;
  fscanf(PFILE, "%d", &weight);
  fscanf(PFILE, "%d", &level);
  fscanf(PFILE, "%d", &ell);
  TBphi.clear(); TBTheta.clear();

  //char str[255];
  if(ell<0){
    useFile=1;
    fscanf(PFILE, "%s", whjfStr);
    cout<<"Will read WHJF from file "<<whjfStr<<".\n";
   // whjfFilename=str;
    }else{
  for(i=0;i<ell;i++){
      fscanf(PFILE, "%d", &x);
        TBphi.push_back(x);
  }
  for(i=0;i<ell;i++){
    fscanf(PFILE, "%d", &x);
    TBTheta.push_back(x);
  }
    }
    fscanf(PFILE, "%d", &x);  sa=x;
    fscanf(PFILE, "%d", &x);  sb=x;
    fscanf(PFILE, "%d", &x);  sc=x;
    fscanf(PFILE, "%d", &x);  if(uptoN==0){uptoN=x;}

    char fff[255];
    fscanf(PFILE, "%s", fff);
    stringstream fffsss;
    fffsss<<fff;
    fffsss>>outfileRoot;

    fscanf(PFILE, "%d", &x);  numHecke=x;


    cout<<"weight = "<<weight<<".\n";
    cout<<"level = "<<level<<".\n";
    cout<<"ell = "<<ell<<".\n";

  cout<<"restriction matrix = "<<sa<<","<<sb<<","<<sc<<".\n";

    cout<<"Output file root = "<<outfileRoot<<"\n";

  cout<<"Desired up to dot product, aka uptoN = "<<uptoN<<".\n";

  heckeArray.clear(); modnArray.clear(); LrootArray.clear();L2rootArray.clear();
  cout<<"Will read "<<numHecke<<" sets of (hecke, modulus, Lroot):\n";
  modnfmpzArray=new fmpz_t[numHecke];
  for(i=0;i<numHecke;i++){
        fscanf(PFILE, "%d", &x);  heckeArray.push_back(x);
        cout<<"Hecke = "<<x<<", ";
        fmpz_fread(PFILE, modn);
        fmpz_init_set(modnfmpzArray[i],modn);
        cout<<"modulus = "<<fmpz_get_si(modn)<<", ";
        fscanf(PFILE, "%d", &x);   L2rootArray.push_back(x);
        cout<<"L2root = "<<x<<".\n";
        fscanf(PFILE, "%d", &x);   LrootArray.push_back(x);
        cout<<"Lroot = "<<x<<".\n";
  }
  fclose(PFILE);
}

