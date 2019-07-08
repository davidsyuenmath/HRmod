#include "MyRepnumList.h"

MyRepnumList::MyRepnumList(int a0, int b0, int c0, int lev0, int uptoDot0, int mindet0){
    initBasics(a0,b0,c0,lev0,uptoDot0,mindet0);
}
void MyRepnumList::initBasics(int a0, int b0, int c0, int lev0, int uptoDot0, int mindet0){
//This initializes to where det is at least - |mindet0|. Weird!!
    a=a0;b=b0;c=c0;lev=lev0;uptoDot=uptoDot0,mindet=mindet0; delta=abs(mindet);
    jfCoeffsSize=0; jfCoeffs=NULL; LrootPowersLen=0; LrootPowers=NULL;
    int i; long long J=uptoDot, nmax=uptoDot; fmpz_init(singExp);
    n.clear(); r.clear(); m.clear(); reducedmn.clear(); reducedr.clear();
    trList.clear();
    nmodL.clear(); rmodL.clear(); mmodL.clear(); nEff.clear(); rEff.clear();
    mEff.clear(); trEff.clear(); reducedmnEff.clear(); reducedrEff.clear();
    long long x,y,z,xmin,xmax,ymin,ymax,zmax;
    long long det=a*c-b*b, t;
    double sq; long long tr;
    long long newmn, newr, mn, zN;
    zmax=myfloor((sqrt(J*J+delta*det)+J)*a/2/det/lev);
    maxq=0;
    for(z=1;z<=zmax;z+=1){
    //cout<<"Doing z = "<<z<<"\n";
        zN=z*lev;
        t=a*a*delta+a*J*zN-det*zN*zN;
        if(t<0){cout<<"Something is wrong in theory?? df9ga9dg in MyRepnumList. Abort.\n";
                cout<<"uptoDot = "<<uptoDot<<"\n";
                cout<<"nmax = "<<nmax<<"\n";
                cout<<"zmax = "<<zmax<<"\n";
                cout<<"(a,b,c)=("<<a<<","<<b<<","<<c<<").\n";
                cout<<"z="<<z<<";det="<<det<<";zN="<<zN<<";delta="<<delta<<";\n";
                cout<<"t="<<t<<"\n";
                exit(1);}
        sq=sqrt(t);
        //cout<<"W? "<<sq<<"\n";
        ymin = myceiling(((-1.0*b*zN)-sq)/(0.5*a)-0.0001);
        ymax = myfloor(((-1.0*b*zN)+sq)/(0.5*a)+0.0001);
        //cout<<"y:"<<sq<<","<<ymin<<","<<ymax<<"\n";
        for(y=ymin;y<=ymax;y++){
            xmin=myceiling(((y*y/4.0-delta))/zN);
            if(xmin<1)xmin=1;
            xmax = (nmax - c*zN - b*y)/a;
            //cout<<"x:"<<xmin<<","<<xmax<<"\n";
            for(x=xmin;x<=xmax;x+=1){
                tr = (a* x +  b*y + c*lev*z);
    if(0)if(tr<0){
      cout<<"x,y,z,lev,a,b,c,uptoDot,mindet,sq,ymin,ymax,xmin,xmax="<<
       x<<","<<y<<","<<z<<","<<lev<<","<<a<<","<<b<<","<<c<<","<<
       uptoDot<<","<<mindet<<
       ","<<sq<<","<<ymin<<","<<ymax<<","<<xmin<<","<<xmax<<"\n";
      exit(0);
    }
                if(tr<=uptoDot){
                    n.push_back(x);
                    r.push_back(y);
                    m.push_back(z);
                    trList.push_back(tr);
                    mn=x*z;
                    myreduce(mn, y, lev, newmn, newr);
                    reducedmn.push_back(newmn);
                    reducedr.push_back(newr);
                    if(newmn>maxq){maxq=newmn;}
                    //cout<<">>>>>>Found: "<<tr<<": "<<x<<","<<y<<","<<z<<";"<<newmn<<","<<newr<<","<<(x*lev*z-y*y/4.0)<<"\n";
                    //if(tr<=0)cout<<">>>>>>Found: "<<tr<<": "<<x<<","<<y<<","<<z<<";"<<newmn<<","<<newr<<","<<(x*lev*z-y*y/4.0)<<"\n";

                }else{
                    cout<<"tr too large?!:"<<x<<","<<y<<","<<z<<"\n";
                }
            }
        }
    }
}


MyRepnumList:: ~MyRepnumList(){
    int i;//destructor
    //cout<<"HERE D1\n";
    if(jfCoeffsSize>0){
        for(i=0;i<jfCoeffsSize;i++){
            fmpz_clear(jfCoeffs[i]);
        }
        delete[] jfCoeffs;
    }
    //cout<<"HERE D20\n";
    if(LrootPowersLen>0){
        for(i=0;i<LrootPowersLen;i++){
            fmpz_clear(LrootPowers[i]);
            fmpz_clear(negPowers[i]);
        }
        delete[] LrootPowers;
        delete[] negPowers;
    }
    //cout<<"HERE D30\n";
    fmpz_clear(singExp);
    //cout<<"HERE D40\n";
}
long long MyRepnumList::myfloor(double x){
  int ans = (int)floor(x);
  if(ans>x){ans--;}
  if(ans+1<=x){ans++;}
  return ans;
}
long long MyRepnumList::myceiling(double x){
  int ans = (int)ceil(x);
  if(ans<x){ans++;}
  if(ans-1>=x){ans--;}
  return ans;
}


void MyRepnumList::myreduce(long long n, long long r, long long index, long long &nnew, long long &rnew){
  long long longp2=2L*index;
  rnew=((r%longp2)+longp2)%longp2;
  if(rnew>index)rnew-=longp2;  //rnew is between -index and index
  nnew=n-((r-rnew)/(longp2))*((r+rnew)/(2L)); //improvement to prevent overflow
  if((4L*index*nnew-rnew*rnew)!=(4L*index*n-r*r)){
    cout<<"ERROR fdhaertegert in getFC of QZSeriesWH\n";
    cout<<"index="<<index<<"\n";
    cout<<"r="<<r<<"\n";
    cout<<"n="<<n<<"\n";
    cout<<"nnew="<<nnew<<"\n";
    cout<<"rnew="<<rnew<<"\n";
    cout<<(4*index*nnew-rnew*rnew)<<","<<(4*index*n-r*r)<<"\n";
    exit(1);
  }
}
void MyRepnumList::prepEfficientSpecial(int L0, QZSeriesWH *jf, int Lroot, fmpz_t modn){
    prepEfficientSpecialGeneralPowersLen(L0, jf, Lroot, modn, L0);
}

void MyRepnumList::prepEfficientSpecialGeneralPowersLen(int L0, QZSeriesWH *jf,
                                                        int Lroot, fmpz_t modn, int powersLen){

//    if(0){"Assumes jf is at worst -1/q +...";}
    int i;
    L=L0;
    fmpz_t coeff; fmpz_init(coeff);
    if(LrootPowersLen!=powersLen){
        if(LrootPowersLen>0){
            for(i=0;i<LrootPowersLen;i++){
                fmpz_clear(LrootPowers[i]);
                fmpz_clear(negPowers[i]);
            }
            delete[] LrootPowers;
            delete[] negPowers;
        }
        LrootPowersLen=powersLen;
        LrootPowers=new fmpz_t[LrootPowersLen];
        negPowers=new fmpz_t[LrootPowersLen];
        fmpz_init_set_si(LrootPowers[0],1);
        fmpz_init_set_si(negPowers[0],-1);
        for(i=1;i<LrootPowersLen;i++){
            fmpz_init(LrootPowers[i]);
            fmpz_mul_si(LrootPowers[i],LrootPowers[i-1],Lroot);
            fmpz_mod(LrootPowers[i],LrootPowers[i],modn);
            fmpz_init(negPowers[i]);
            fmpz_neg(negPowers[i],LrootPowers[i]);
        }
    }
    if(jfCoeffsSize>0){
        for(i=0;i<jfCoeffsSize;i++){
            fmpz_clear(jfCoeffs[i]);
        }
    }
    delete[] jfCoeffs;
    int maxr0=jf->getCoeffMaxr(0);  //jf's max r.
    int maxnWhenmIs0 = (uptoDot +maxr0*b )/a;
    int maxmWhennIs0 = (uptoDot + maxr0*b)/(c*lev);
    int conceivableMore= (2*maxr0+1)*(maxnWhenmIs0+maxmWhennIs0)+maxr0+1;
    jfCoeffs=new fmpz_t[n.size()+conceivableMore];
    jfCoeffsSize=0;
    nmodL.clear();rmodL.clear();mmodL.clear();trEff.clear();reducedmnEff.clear();reducedrEff.clear();
    nEff.clear();rEff.clear();mEff.clear();
    maxqEff=0;
    long long tr; long long x,y,z,mn;
    //cout<<"uptDot="<<uptoDot<<", maxr0="<<maxr0<<".\n";
    //cout<<"maxmWhennIs0="<<maxmWhennIs0<<", maxnWhenmIs0="<<maxnWhenmIs0<<".\n";
    x=0;for(z=1;z<=maxmWhennIs0;z++){
        for(y=-maxr0;y<=maxr0;y++){
            tr = (a* x +  b*y + c*lev*z);
            if(tr<=uptoDot){
                n.push_back(x);
                r.push_back(y);
                m.push_back(z);
                mn=x*z;
                reducedmn.push_back(mn);
                reducedr.push_back(y);
            }
        }
    }
    z=0;for(x=1;x<=maxnWhenmIs0;x++){
        for(y=-maxr0;y<=maxr0;y++){
            tr = (a* x +  b*y + c*lev*z);
            if(tr<=uptoDot){
                n.push_back(x);
                r.push_back(y);
                m.push_back(z);
                mn=x*z;
                reducedmn.push_back(mn);
                reducedr.push_back(y);
            }
        }
    }
    x=0;z=0;
    for(y=-maxr0;y<=-1;y++){
        tr = (a* x +  b*y + c*lev*z);
        if(tr<=uptoDot){
            n.push_back(x);
            r.push_back(y);
            m.push_back(z);
            mn=x*z;
            reducedmn.push_back(mn);
            reducedr.push_back(y);
        }
    }
    x=-1;z=1;y=0;
    tr = (a* x +  b*y + c*lev*z);
    if(tr<=uptoDot){
        n.push_back(x);
        r.push_back(y);
        m.push_back(z);
        mn=x*z;
        reducedmn.push_back(mn);
        reducedr.push_back(y);
    }
    //Now get only those with c(mn,r)!=0
    fmpz_t tot, tmp;
    fmpz_init_set_si(tot,0);
    fmpz_init(tmp);
    //cout<<"size of slong = "<<sizeof(slong)<<".\n";
    //cout<<"size of int = "<<sizeof(int)<<".\n";
    //slong test = 2<<30;
    //cout<<test<<","<<test*test<<".\n";
    fmpz_t slongmax, slongmin;
    fmpz_init_set_si(slongmax,2);
    fmpz_pow_ui(slongmax,slongmax,sizeof(slong)*8-2);
    fmpz_init(slongmin);
    fmpz_neg(slongmin,slongmax);
    for(i=0;i<n.size();i++){
        jf->getFC_noreduce(coeff, reducedmn[i],reducedr[i]);
        if((!fmpz_is_zero(coeff))){
            //if(fmpz_cmp(coeff,slongmax)>0){cout<<"coeff ";fmpz_print(coeff);cout<<" is too big for slong! Abort.\n";exit(1);}
            //if(fmpz_cmp(coeff,slongmin)<0){cout<<"coeff ";fmpz_print(coeff);cout<<" is too negative for slong! Abort.\n";exit(1);}
            tr=a*n[i] +  b*r[i] + c*lev*m[i];
            if(tr<0){
                fmpz_mul_si(tmp,coeff,(slong)tr);
                fmpz_add(tot,tot,tmp);
            }
            nEff.push_back(n[i]);
            rEff.push_back(r[i]);
            mEff.push_back(m[i]);
            nmodL.push_back(n[i]%powersLen);
            rmodL.push_back(r[i]%powersLen);
            mmodL.push_back(m[i]%powersLen);
            trEff.push_back(tr);
            reducedmnEff.push_back(reducedmn[i]);
            reducedrEff.push_back(reducedr[i]);
            if(reducedmn[i]>maxqEff){maxqEff=reducedmn[i];}
            fmpz_init(jfCoeffs[jfCoeffsSize]);
            fmpz_set(jfCoeffs[jfCoeffsSize], coeff);
            jfCoeffsSize++;
            if(0)if(tr<=0)cout<<">>>>>>Well: "<<tr
                <<": "<<n[i]<<","<<r[i]<<","<<m[i]<<";"<<reducedmn[i]<<","<<reducedr[i]
                <<", "<<fmpz_get_si(coeff)<<"\n";
       }
    }
    fmpz_clear(coeff);
    int tBPA,tBPB,tBPC;
    jf->getABC(tBPA,tBPB,tBPC);
    BPA=tBPA;BPB=tBPB;BPC=tBPC;
    //singExp=(BPA*a+BPB*b+BPC*c)+((long long)fmpz_get_si(tot));
    fmpz_set_si(singExp,(slong)(BPA*a+BPB*b+BPC*c));
    fmpz_add(singExp,singExp,tot);
            if(fmpz_cmp(singExp,slongmax)>0){cout<<"singExp ";fmpz_print(singExp);cout<<" is too big for slong! Abort.\n";exit(1);}
            if(fmpz_cmp(singExp,slongmin)<0){cout<<"singExp ";fmpz_print(singExp);cout<<" is too negative for slong! Abort.\n";exit(1);}
    //cout<<"NEW CODE: ("<<(BPA*a+BPB*b+BPC*c)<<","<<fmpz_get_si(tot)<<","<<singExp<<").\n";

    fmpz_clear(tot);fmpz_clear(tmp);fmpz_clear(slongmax);fmpz_clear(slongmin);
    //cout<<"BP ABC = "<<BPA<<","<<BPB<<","<<BPC<<".\n";
}


MySeriesTruncMod* MyRepnumList::oneBPCosetRestrictionPrepped(long long BAa, long long BAb, long long BAc, long long Ldenom,
                                       fmpz_t modn, int oneLroot, int uptoN){
    fmpz_t tmp; fmpz_init_set_si(tmp, oneLroot);
    MySeriesTruncMod* ans=oneBPCosetRestrictionPrepped(BAa,BAb,BAc, Ldenom, modn, tmp, uptoN);
    fmpz_clear(tmp);
    return ans;
}
MySeriesTruncMod* MyRepnumList::oneBPCosetRestrictionPrepped(long long BAa, long long BAb, long long BAc,
                                                             long long Ldenom,
                                       fmpz_t modn, fmpz_t oneLroot, int uptoN0){
    if(fmpz_cmp_si(singExp,uptoN0)>0){
    //if(singExp>uptoN0){
        return new MySeriesTruncMod(modn, uptoN0);
    }
    int uptoN=uptoN0 - fmpz_get_si(singExp);
    BAa=BAa%Ldenom;    BAb=BAb%Ldenom;    BAc=BAc%Ldenom;
    int i;
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn, uptoN);
    ans->setCoeff(1, 0);
    //MySeriesTruncMod* negPart = new MySeriesTruncMod(modn, uptoN);
    //negPart->setCoeff(1, 0);
    fmpz_t coeff; fmpz_init(coeff);
    MySeriesTruncMod* tmpBinomial=new MySeriesTruncMod(modn, uptoN);
    //cout<<"HERE1. mnEff.size()="<<reducedmnEff.size()<<", jfCoeffsSize="<<jfCoeffsSize<<", ";;
    //cout<<"trEff.size()="<<trEff.size()<<"\n";;
    //cout<<"Restricting: trEff.size()="<<trEff.size()<<"\n";;
    long long e;
    //////////slong qAdjust=BPA*a+BPB*b+BPC*c;
    slong returnedAdjust;
    //long long expon;
    //cout<<"qAdjust="<<qAdjust<<".\n";
    fmpz_t fmpzexpon;
    fmpz_init(fmpzexpon);
    for(i=0;i<reducedmnEff.size();i++)if(trEff[i]<=uptoN){

        e=(((BAa*nmodL[i]+BAb*rmodL[i]+BAc*mmodL[i]*lev)%Ldenom)+Ldenom)%Ldenom;
            //cout<<"HERE. i="<<i<<", (n,r,m)="<<nEff[i]<<","<<rEff[i]<<","<<mEff[i]<<","
            //<<" tr="<<trEff[i]<<";";
        //fmpz_powm_ui(coeff,oneLroot,(ulong)e,ans->modn);
        ////////if(0)if(trEff[i]<0){
        ///////    cout<<"{tr,coeff,ppp}={"<<trEff[i]<<",";
        ///////}
        //////if(0){"Improve above to coeff be LrootArray[e]";}
               //cout<<" e="<<e<<", coeff="<<fmpz_get_si(coeff)
               //<<", c(mn,r)="<<fmpz_get_si(jfCoeffs[i])<<"\n.";
        //fmpz_neg(coeff, LrootPowers[e]);
        returnedAdjust=tmpBinomial->setMonicBinomialBPSpecial(negPowers[e],trEff[i]);
        //cout<<"returnedAdjust from tmpBinomial = "<<returnedAdjust<<","<<fmpz_get_si(jfCoeffs[i])<<",";
        ////////expon=(long long)fmpz_get_si(jfCoeffs[i]);

        ///////if(0)returnedAdjust=returnedAdjust*expon;
        //cout<<returnedAdjust<<",";
        ////////////if(0)qAdjust+=returnedAdjust;
               //cout<<"HERE5: binomial= ";tmpBinomial->printstdout();cout<<", adjustments: "<<returnedAdjust
               //<<","<<qAdjust<<".\n";
        //cout<<qAdjust<<".\n";
        //if(expon>=0){
            ans->multiplyByPower(tmpBinomial, jfCoeffs[i]);
        //}else{
        //    negPart->multiplyByPower(tmpBinomial, -expon);
        //}
               //cout<<"HERE6: ("<<fmpz_get_si(jfCoeffs[i])<<") ";
               //tmpBinomial->printstdout();cout<<", ";
               //ans->printstdout();cout<<".\n";
    }
    //cout<<"abcs: "<<a<<","<<b<<","<<c<<","<<BPA<<","<<BPB<<","<<BPC<<","<<BAa<<","<<BAb<<","<<BAc<<"\n";
    /////////if(0)if(qAdjust!=singExp){
    /////////    cout<<"ERROR daga4333 in MyRepnumList oneBPCosetRestrictionPrepped ("<<qAdjust<<","<<singExp<<").Abort.\n";
    ///////    exit(1);
    ///////}
    //ans->multiplyByPower(negPart, -1);
    ans->multiplyByPowerOfQandChangeTrunc(fmpz_get_si(singExp));
    //cout<<"initial,final qAdjust, ans->trunc="<<BPA*a+BPB*b+BPC*c<<","<<qAdjust<<","<<ans->trunc<<".\n";
    if(ans->trunc!=uptoN0+1){
        cout<<"ERROR dgasdg2t5 in MyRepnumList oneBPCosetRestrictionPrepped ("<<ans->trunc<<","<<uptoN0<<").Abort.\n";
        exit(1);
    }
    e=(((BAa*BPA+BAb*BPB+BAc*BPC)%Ldenom)+Ldenom)%Ldenom;
    //fmpz_powm_ui(coeff,oneLroot,(ulong)e,ans->modn);
    //cout<<"scalar multiply by "<<fmpz_get_si(coeff)<<".\n";
    ans->scalarMultiplityBy(LrootPowers[e]);
    //cout<<"HERE7: ";ans->printstdout();cout<<".\n";
    fmpz_clear(coeff); delete tmpBinomial;fmpz_clear(fmpzexpon);
    return ans;
}

void MyRepnumList::fprintfDiagnostics(FILE* F, int uptoN){
    fprintf(F, "BPA=%lld;\n", BPA);
    fprintf(F, "BPB=%lld;\n", BPB);
    fprintf(F, "BPC=%lld;\n", BPC);
    fprintf(F, "a=%lld;\n", a);
    fprintf(F, "b=%lld;\n", b);
    fprintf(F, "c=%lld;\n", c);
    fprintf(F, "L=%lld;\n", L);
    fprintf(F, "singExp=%ld;\n", fmpz_get_si(singExp));
    int i;
    fprintf(F, "data={\n");
    for(i=0;i<reducedmnEff.size();i++)if(trEff[i]<=uptoN){
        if(i>0){fprintf(F, ",\n");}
        //nEff, rEff, mEff, trEff, reducedmnEff, reducedrEff, jfCoeff
        fprintf(F, "{");
        fprintf(F, "%lld,", nEff[i]);
        fprintf(F, "%lld,", rEff[i]);
        fprintf(F, "%lld,", mEff[i]);
        fprintf(F, "%lld,", trEff[i]);
        fprintf(F, "%lld,", reducedmnEff[i]);
        fprintf(F, "%lld,", reducedrEff[i]);
        fmpz_fprint(F, jfCoeffs[i]);
        fprintf(F,"}");
    }
    fprintf(F, "\n};\n");
}

int MyRepnumList::hasDesiredDot(long long targetDot){
    if(targetDot>uptoDot){
        cout<<"MyRepnumList asked for targetDot ("<<targetDot<<") > uptoDot ("
            <<uptoDot<<"). Abort.\n";
        exit(1);
    }
    int i;
    for(i=0;i<trList.size();i++)if(trList[i]==targetDot){return 1;}
    return 0;
}
void MyRepnumList::printAbstractRestrictionToFile(FILE* OUTFILE, SigObjectOldStyle* sigObj,
                                                   long long BAa, long long BAb, long long BAc,
                                   long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoN, int weight){
    if(sigObj==NULL){fprintf(OUTFILE,"Indeterminate");}
    if(uptoN>uptoDot){
        initBasics(a,b,c,lev,uptoN,mindet);
        //cout<<"ERROR 4328gdgsd uptoN>uptoDot in restrictGrit ("<<uptoN<<">"<<uptoDot<<"). Abort.\n";
        //exit(1);
    } //[So make a new one!!!]
    BAa=BAa%Ldenom;    BAb=BAb%Ldenom;    BAc=BAc%Ldenom;
    int i;
    slong e;
    int ind; long long pm;
    //fmpz_t fc;fmpz_init(fc);
    cout<<"Starting printAbstractRestrictionToFile.\n";
    cout<<"(a,b,c, uptoDot)=("<<a<<","<<b<<","<<c<<","<<uptoDot<<")\n";
    fprintf(OUTFILE, "(0");
    for(i=0;i<n.size();i++)if(trList[i]<=uptoN){
        //cout<<"(n,r,m)=("<<n[i]<<","<<r[i]<<","<<m[i]<<")\n";
        e=(((BAa*n[i]+BAb*r[i]+BAc*m[i]*lev)%Ldenom)+Ldenom)%Ldenom;
        fprintf(OUTFILE,"+");
        fmpz_fprint(OUTFILE, rootPowers[e]);
        fprintf(OUTFILE,"*");
        if(weight%2==0){
            ind=sigObj->getIndexFromSig(2*m[i]*lev,r[i],2*n[i]);
            pm=1;
        }else{
            ind=sigObj->getIndexFromSig(2*m[i]*lev,r[i],2*n[i],pm);
        }
        cout<<"(n,r,m, tr, ind )=("<<n[i]<<","<<r[i]<<","<<m[i]<<","<<trList[i]<<","<<ind<<")\n";
        if(pm==1){
            fprintf(OUTFILE, "za[%d]*q^%lld\n",ind,trList[i]);
        }else{
            fprintf(OUTFILE, "(-za[%d])*q^%lld\n",ind,trList[i]);
        }
    }
    fprintf(OUTFILE,");\n");
}

MySeriesTruncMod* MyRepnumList::restrictGrit(QZSeriesWH* jf, long long BAa, long long BAb, long long BAc,
                               long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoN, int weight){
    //if((BAa==0)&&(BAb==0)&&(BAc==0))return restrictGritNoCharacter(jf,Ldenom,modn,uptoN,weight);
    //
    //cout<<"HERE 1\n";
    if(uptoN>uptoDot){
        initBasics(a,b,c,lev,uptoN,mindet);
        //cout<<"ERROR 4328gdgsd uptoN>uptoDot in restrictGrit ("<<uptoN<<">"<<uptoDot<<"). Abort.\n";
        //exit(1);
    } //[So make a new one!!!]
    BAa=BAa%Ldenom;    BAb=BAb%Ldenom;    BAc=BAc%Ldenom;
    //cout<<"HERE 2\n";
    int i;
    slong e;
    fmpz_t fc;fmpz_init(fc);
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn, uptoN,0);
    //cout<<"Starting restrictGrit.\n";
    //cout<<"HERE 3\n";
    //cout<<"n.size()="<<n.size()<<"\n";
    //cout<<"uptoN="<<uptoN<<"\n";
    for(i=0;i<n.size();i++)if(trList[i]<=uptoN){
        e=(((BAa*n[i]+BAb*r[i]+BAc*m[i]*lev)%Ldenom)+Ldenom)%Ldenom;
        //add oneLroot^e * GritFC(n,r,m) *q^trList
        //cout<<"e="<<e<<"\n";
        //cout<<i<<" Getting getFCGrit.\n";
        //cout<<"nrm size="<<n.size()<<".\n";
        //cout<<"(n,r,m)=("<<n[i]<<","<<r[i]<<","<<m[i]<<")\n";
        //cout<<"(a,b,c,uptoDot)=("<<a<<","<<b<<","<<c<<","<<uptoDot<<")\n";
        //cout<<"(n,r,m,weight,lev)=("<<n[i]<<","<<r[i]<<","<<m[i]<<","<<weight<<","<<lev<<")\n";
        //cout<<jf->getString();cout<<"\n";
        jf->getFCGrit(fc,n[i],r[i],m[i],weight,lev);
        //cout<<"multiply.\n";
        //cout<<"fc=";fmpz_print(fc);cout<<", rootPowers[e]=";fmpz_print(rootPowers[e]);cout<<"\n";
        fmpz_mul(fc,fc,rootPowers[e]);
        //cout<<"fc=";fmpz_print(fc);cout<<"\n";
        //cout<<"updateCoeffAdd.\n";
        //cout<<"a,b,c,lev="<<a<<","<<b<<","<<c<<","<<lev<<"\n";
        //cout<<"n[i]="<<n[i]<<"\n";
        //cout<<"r[i]="<<r[i]<<"\n";
        //cout<<"m[i]="<<m[i]<<"\n";
        //cout<<"trList[i]="<<trList[i]<<"\n";
        //cout<<"rootPowers["<<jj<<"]=";fmpz_print(rootPowers[e]);cout<<"\n";
        //cout.flush();
        if(!fmpz_is_zero(fc)){ans->updateCoeffAdd(fc, trList[i]);}
        //cout<<"HERE 3.5\n";
    }
    //cout<<"Made restrictGrit.\n";
    //cout<<"HERE 4\n";

    fmpz_clear(fc);
    return ans;
}

MySeriesTruncMod* MyRepnumList::restrictGritNoCharacter(QZSeriesWH* jf,
                                    long long Ldenom, fmpz_t modn, int uptoN, int weight){
    //
    if(uptoN>uptoDot){cout<<"ERROR 4328gdgsd uptoN>uptoDot in restrictGrit ("<<uptoN<<">"<<uptoDot<<"). Abort.\n";
        exit(1);initBasics(a,b,c,lev,uptoN,mindet);
    }
    int i;
    fmpz_t fc;fmpz_init(fc);
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn, uptoN,0);
    for(i=0;i<n.size();i++)if(trList[i]<=uptoN){
        jf->getFCGrit(fc,n[i],r[i],m[i],weight,lev);
        ans->updateCoeffAdd(fc, trList[i]);
    }

    fmpz_clear(fc);
    return ans;
}

MySeriesTruncMod* MyRepnumList::restrictGritTweak(QZSeriesWH* jf, long long BAa, long long BAb, long long BAc,
                                    long long Ldenom, fmpz_t modn, fmpz_t* rootPowers, int uptoNx2, int weight){
    /**

    **/
    if(uptoNx2>uptoDot){cout<<"ERROR fg3whsdg uptoNx2>uptoDot in restrictGritTweak ("<<uptoNx2<<">"<<uptoDot<<"). Abort.\n";
        exit(1);initBasics(a,b,c,lev,uptoNx2,mindet);
    }
    BAa=BAa%Ldenom;    BAb=BAb%Ldenom;    BAc=BAc%Ldenom;
    int i;
    slong e;
    fmpz_t fc;fmpz_init(fc);
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn, uptoNx2,1);
    //cout<<"HERE Tweak A0: n.size()="<<n.size()<<"\n";
    //cout<<"HERE Tweak A0: (a,b,c, uptoDot)=("<<a<<","<<b<<","<<c<<","<<uptoDot<<")\n";
    for(i=0;i<n.size();i++)if(trList[i]<=uptoNx2){
        //cout<<"HERE Tweak A1: (n,r,m)=("<<n[i]<<","<<r[i]<<","<<m[i]<<")\n";
      if((n[i]%2!=0)&&(r[i]%2!=0)&&(m[i]%2!=0)){
        //cout<<"HERE Tweak A2: (n,r,m)=("<<n[i]<<","<<r[i]<<","<<m[i]<<")";

        e=(((BAa*n[i]+BAb*r[i]+BAc*m[i]*lev)%Ldenom)+Ldenom)%Ldenom;
        //add oneLroot^e * GritFC(n,r,m) *q^trList
        jf->getFCGritTweak(fc,n[i],r[i],m[i],weight,lev);
        fmpz_mul(fc,fc,rootPowers[e]);
        //cout<<": ";fmpz_print(fc);cout<<"\n";
        ans->updateCoeffAdd(fc, trList[i]);
      }
    }

    fmpz_clear(fc);
    return ans;
}

int MyRepnumList::getMinDotGrit(QZSeriesWH* jf, int weight){
    //Find min dot where Grit(jf) is not zero
    //IF all zero, then redo this MyRepnumList with twice uptoDot. Refactor initialization routine.
    cout<<"getMinDotGrit to be programmed.\n";exit(1);
    int i, ans=-1;
    fmpz_t fc;fmpz_init(fc);
    while(ans<0){
        for(i=0;i<n.size();i++){
            //add oneLroot^e * GritFC(n,r,m) *q^trList
            jf->getFCGrit(fc,n[i],r[i],m[i],weight,lev);
            if(!fmpz_is_zero(fc)){
                if((ans<0)||(trList[i]<ans)){ans=trList[i];}
            }
        }
        if(ans<0){
            initBasics(a,b,c,lev,2*uptoDot,mindet);
        }
    }
    fmpz_clear(fc);
    return ans;
}
