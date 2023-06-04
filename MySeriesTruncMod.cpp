#include "MySeriesTruncMod.h"

int MySeriesTruncMod::checkSameTrunc=1;

MySeriesTruncMod::MySeriesTruncMod(fmpz_t modn0, int upto){
    fmpz_init(modn);
    fmpz_set(modn, modn0);
    fmpz_mod_ctx_init(modnctx,modn);
    fmpz_mod_poly_init(p, modnctx);
    fmpz_mod_poly_zero(p, modnctx);
    len=upto+1;
    van=0; //Could be len also, but we choose van=0.  Do not change
    normalized=0;
    inSqrtQ=0;
}
MySeriesTruncMod::MySeriesTruncMod(int modn0, int upto){
    fmpz_init(modn);
    fmpz_set_si(modn, modn0);
    fmpz_mod_ctx_init(modnctx,modn);
    fmpz_mod_poly_init(p, modnctx);
    fmpz_mod_poly_zero(p, modnctx);
    len=upto+1;
    van=0; 
    normalized=0;
    inSqrtQ=0;
}
MySeriesTruncMod::MySeriesTruncMod(int modn0, int upto, int inSqrtQ0){
    fmpz_init(modn);
    fmpz_set_si(modn, modn0);
    fmpz_mod_ctx_init(modnctx,modn);
    fmpz_mod_poly_init(p, modnctx);
    fmpz_mod_poly_zero(p, modnctx);
    len=upto+1;
    van=0; 
    normalized=0;
    inSqrtQ=inSqrtQ0;
}
MySeriesTruncMod::MySeriesTruncMod(fmpz_t modn0, int upto, int inSqrtQ0){
    fmpz_init(modn);
    fmpz_set(modn, modn0);
    fmpz_mod_ctx_init(modnctx,modn);
    fmpz_mod_poly_init(p, modnctx);
    fmpz_mod_poly_zero(p, modnctx);
    len=upto+1;
    van=0; 
    normalized=0;
    inSqrtQ=inSqrtQ0;
}

MySeriesTruncMod::MySeriesTruncMod(MySeriesTruncMod*f){
    fmpz_init(modn);
    fmpz_set(modn, f->modn);
    fmpz_mod_ctx_init(modnctx,modn);
    fmpz_mod_poly_init(p, modnctx);
    fmpz_mod_poly_set(p, f->p, modnctx);
    len=f->len;
    van=f->van;
    normalized=f->normalized;
    inSqrtQ=f->inSqrtQ;
}

MySeriesTruncMod:: ~MySeriesTruncMod(){//destructor
    fmpz_clear(modn);
    fmpz_mod_ctx_clear(modnctx);
    fmpz_mod_poly_clear(p, modnctx);
}
slong MySeriesTruncMod::getTrunc(){
    return van + len;
}


void MySeriesTruncMod::checkSameinSqrtQ(MySeriesTruncMod *f, MySeriesTruncMod *g){
    if(f->inSqrtQ!=g->inSqrtQ){
        cout<<"ERROR rh4572253 in MySeriesTruncMod, attempt to operate with two different inSqrtQ:"
                <<f->inSqrtQ<<","<<g->inSqrtQ<<". Abort.\n";
        exit(1);
    }
}

int MySeriesTruncMod::getVanishingOrder(){
    //SENTINNEL VALUE OF -1 for all zero coeffs.
    normalize();
    if(len==0)return -1;
    return van;
}
int MySeriesTruncMod::getVanishingOrderNoSentinnel(){
    //NO MORE SENTINNEL VALUE OF -1 for all zero coeffs.  Just return the known order of vanishin.
    normalize();
    return van;
}
int MySeriesTruncMod::couldBeZero(){
    normalize();
    return (len==0);
}
void  MySeriesTruncMod::normalize(){
    if(normalized)return;
    normalized=1;
    fmpz_t tmp;fmpz_init(tmp);
    int newvan = van + len;
    //for(i=1;i<trunc;i++){  This was original.  Shouldn't this be i=0; ?? Fixed 20200602
    int i;
    for(i=0;i<len;i++){
        fmpz_mod_poly_get_coeff_fmpz(tmp,p,i, modnctx);
        if(!fmpz_is_zero(tmp)){
            newvan =  i+van;
            break;
        }
    }
    fmpz_clear(tmp);
    if(newvan==van)return;
    if(newvan==van+len){//zero poly
        van = van + len;
        len = 0;
        fmpz_mod_poly_zero(p, modnctx); //maybe not necessary?
        return;
    }
    // p should be adjusted down by v-van and van=v;
    fmpz_mod_poly_t tmpp;
    fmpz_mod_poly_init(tmpp, modnctx);
    fmpz_mod_poly_set(tmpp, p, modnctx);
    fmpz_mod_poly_shift_right(p, tmpp, newvan-van, modnctx);
    fmpz_mod_poly_clear(tmpp, modnctx);
    len = len -(newvan-van);
    van = newvan;
}

MySeriesTruncMod* MySeriesTruncMod::monicBinomial(fmpz_t modn0, int upto, fmpz_t oneCoeff, int e){
    //initializes to 1 + c * x^e
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn0, upto);
    fmpz_mod_poly_set_coeff_ui(ans->p, 0, 1, modnctx);
    if(e<=upto){fmpz_mod_poly_set_coeff_fmpz(ans->p, e, oneCoeff, modnctx);}
    if(e<1){cout<<"ERROR czvzv in MySeriesTruncMod. Abort.\n";exit(1);}
    return ans;
}  //Initializes to 1 + c*x^e
int MySeriesTruncMod::setMonicBinomialBPSpecial(fmpz_t oneCoeff, int e){
    //initializes to 1 + c * x^e if e>=0   return 0
    //else  x^(-e) + c   if e<0            return e
    //cout<<"setMonicBinomialBPSpecial (oneCoeff, exponent, trunc, modn)=("<<fmpz_get_si(oneCoeff)
    //    <<","<<e<<","<<trunc<<","<<fmpz_get_si(modn)<<") returning...\n ";
    if(inSqrtQ){cout<<"Error 4t8asf0 in setMonicBinomialBPSpecial. Abort.\n";exit(1);}
    len=getTrunc();
    van=0;
    fmpz_mod_poly_zero(p, modnctx);
    //cout<<"[STEP zeroing out] ";printstdout();cout<<"\n";
    if(e>0){
        fmpz_mod_poly_set_coeff_ui(p, 0, 1, modnctx);
        //cout<<"[STEP set coeff 0 to 1] ";printstdout();cout<<"\n";
        if(e<len){
            fmpz_mod_poly_set_coeff_fmpz(p, e, oneCoeff, modnctx);
            //cout<<"[STEP set coeff e to oneCoeff] ";printstdout();cout<<"\n";
        }
        //cout<<"[expon>0] ";printstdout();cout<<"\n";
        return 0;
    }else if(e<0){
        fmpz_mod_poly_set_coeff_fmpz(p, 0, oneCoeff, modnctx);
        //cout<<"oneCoeff="<<fmpz_get_si(oneCoeff)<<"\n";
        //cout<<"[STEP Z2 set coeff 0 to oneCoeff] ";printstdout();cout<<"\n";
        if(-e<len){
            fmpz_mod_poly_set_coeff_ui(p, -e, 1, modnctx);
            //cout<<"[STEP set coeff -e to 1] ";printstdout();cout<<"\n";
        }
        //cout<<"[expon<0] ";printstdout();cout<<"\n";
        return e;
    }else{
        fmpz_t temp;
        fmpz_init_set_si(temp, 1);
        fmpz_add(temp, temp, oneCoeff);
        fmpz_mod_poly_set_coeff_fmpz(p, 0, temp, modnctx);
        //cout<<"[STEP set coeff 0 to temp=1+oneCoeff] ";printstdout();cout<<"\n";
        //cout<<"[expon o/w] ";printstdout();cout<<"\n";
        fmpz_clear(temp);
        return 0;
    }
}
MySeriesTruncMod* MySeriesTruncMod::monomial(fmpz_t modn0, int upto, fmpz_t oneCoeff, int e){
    //initializes to c * x^e
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn0, upto);
    if(e<0){cout<<"ERROR cvrgaaa in MySeriesTruncMod. Abort.\n";exit(1);}
    if(e<=upto){fmpz_mod_poly_set_coeff_fmpz(ans->p, e, oneCoeff, modnctx);}
    ans->normalize();
    return ans;
}  //Initializes to c*x^e

void MySeriesTruncMod::getFullPoly(fmpz_mod_poly_t ans){
    if(van==0){
        fmpz_mod_poly_set(ans, p, modnctx);
    }else{
        fmpz_mod_poly_shift_left(ans, p, van, modnctx);
    }
}

MySeriesTruncMod* MySeriesTruncMod::add(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkSameinSqrtQ(f,g);
    fmpz_mod_poly_t fp, gp, ansp;
    fmpz_mod_poly_init(fp, f->modnctx);
    fmpz_mod_poly_init(gp, g->modnctx);
    f->getFullPoly(fp);
    g->getFullPoly(gp);
    slong e = f->getTrunc(); if(g->getTrunc()<e){e=g->getTrunc();}
    fmpz_mod_poly_init(ansp, f->modnctx);
    fmpz_mod_poly_add(ansp, fp, gp, f->modnctx);
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, e-1, f->inSqrtQ);
    fmpz_mod_poly_set(ans->p, ansp, f->modnctx);
    ans->normalize();
    fmpz_mod_poly_clear(fp, f->modnctx);
    fmpz_mod_poly_clear(gp, f->modnctx);
    fmpz_mod_poly_clear(ansp, f->modnctx);
    return ans;
}
void MySeriesTruncMod::addAndAdjustTrunc(MySeriesTruncMod* g){
    addBy(g);
}
void MySeriesTruncMod::addBy(MySeriesTruncMod* g){
    normalized=0;
    checkSameinSqrtQ(this,g);
    fmpz_mod_poly_t fp, gp;
    fmpz_mod_poly_init(fp, modnctx);
    fmpz_mod_poly_init(gp, g->modnctx);
    getFullPoly(fp);
    g->getFullPoly(gp);
    slong e = getTrunc(); if(g->getTrunc()<e){e=g->getTrunc();}
    fmpz_mod_poly_add(p, fp, gp, modnctx);
    van=0;
    len=e;
    normalize();
    fmpz_mod_poly_clear(fp, modnctx);
    fmpz_mod_poly_clear(gp, modnctx);
}
MySeriesTruncMod* MySeriesTruncMod::scalarMultiplity(MySeriesTruncMod* g, fmpz_t a){
    MySeriesTruncMod* ans = new MySeriesTruncMod(g->modn, g->getTrunc()-1, g->inSqrtQ);
    fmpz_mod_poly_scalar_mul_fmpz(ans->p, g->p, a, g->modnctx);
    return ans;
}
void MySeriesTruncMod::scalarMultiplyBy(fmpz_t a){
    normalized=0;
    fmpz_mod_poly_scalar_mul_fmpz(p, p, a, modnctx);
}
void MySeriesTruncMod::scalarMultiplityBy(fmpz_t a){
    normalized=0;
    fmpz_mod_poly_scalar_mul_fmpz(p, p, a, modnctx);
}
MySeriesTruncMod* MySeriesTruncMod::multiply(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkSameinSqrtQ(f,g);
    f->normalize();
    g->normalize();
    slong e = f->len; if(g->len<e){e=g->len;}
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, f->van+g->van+e, f->inSqrtQ);
    ans->van=f->van+g->van;
    ans->len = e;
    ans->normalized=1;
    // NOTE: result of multiplication is automatically normalized.
    if(e>0)fmpz_mod_poly_mullow(ans->p, f->p, g->p, e, f->modnctx);
    return ans;
}
void MySeriesTruncMod::multiplyBy(MySeriesTruncMod* g){
    checkSameinSqrtQ(this,g);
    normalize();
    g->normalize();
    // NOTE: result of multiplication is automatically normalized.
    slong e = len; if(g->len<e){e=g->len;}
    van=van+g->van;
    len = e;
    if(e==0){fmpz_mod_poly_zero(p, modnctx);}
    else{fmpz_mod_poly_mullow(p, p, g->p, e, g->modnctx);}
}
MySeriesTruncMod* MySeriesTruncMod::pow(MySeriesTruncMod* f, slong e){
    f->normalize();
    if((e<0)&&(f->van>0)){
      cout<<"Cannot take inverse of series of vanishing order > 0\n"; exit(1);
    }
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, e*f->van+f->len-1, f->inSqrtQ);
    // DO WE NEED SPECIAL CASES WHEN f is possibly zero? or when e is 0?
    if(e>=0){
        fmpz_mod_poly_pow_trunc_binexp(ans->p, f->p, e,f->len, f->modnctx);
        ans->van = e*f->van;
        ans->len = f->len;
    }else{ // e<0 and f->van==0 here.
        ans->len=f->len; ans->van=0;
        fmpz_mod_poly_pow_trunc_binexp(ans->p, f->p, -e,ans->len, f->modnctx);
        fmpz_mod_poly_inv_series(ans->p, ans->p, ans->len, f->modnctx);
    } 
    ans->normalized=1;
    // NOTE: result of multiplication is automatically normalized.
    return ans;
}
void MySeriesTruncMod::multiplyByPower(MySeriesTruncMod* f, slong e){
    checkSameinSqrtQ(this,f);
    normalize();
    f->normalize();
    // NOTE: result of multiplication is automatically normalized.
    if(e==0)return;
    if(e==1){
        multiplyBy(f);return;
    }
    MySeriesTruncMod* tmp = MySeriesTruncMod::pow(f,e);
    multiplyBy(tmp);
    delete tmp;
}

void MySeriesTruncMod::getCoeff(fmpz_t c, int e){
    if(e<van){
        fmpz_set_si(c, 0);
        return;
    }
    if(e>=getTrunc()){cout<<"ERROR 32tfgdf: Request for coeff is beyond valid trunc ("<<e<<","<<getTrunc()<<").  Abort.\n";exit(1);}
    fmpz_mod_poly_get_coeff_fmpz(c, p, e-van, modnctx);
}
void MySeriesTruncMod::setCoeff(fmpz_t c, int e){
    normalized=0;
    if(e>=getTrunc()){cout<<"ERROR adfshye45y: Attempt to set coeff beyond valid trunc ("<<e<<","<<getTrunc()<<").  Abort.\n";exit(1);}
    if(e<van){
        fmpz_mod_poly_shift_left(p, p, van-e, modnctx);
        len = len +(van-e);
        van = e;
        fmpz_mod_poly_set_coeff_fmpz(p, e-van, c, modnctx);
    }else{
        fmpz_mod_poly_set_coeff_fmpz(p, e-van, c, modnctx);
    }
    // SAVE TIME BY not always normalizing.
    //normalize();
}
void MySeriesTruncMod::setCoeff(ulong c, int e){
    normalized=0;
    if(e>=getTrunc()){cout<<"ERROR jgfdgr8lll: Attempt to set coeff beyond valid trunc ("<<e<<","<<getTrunc()<<").  Abort.\n";exit(1);}
    if(e<van){
        fmpz_mod_poly_shift_left(p, p, van-e, modnctx);
        len = len +(van-e);
        van = e;
        fmpz_mod_poly_set_coeff_ui(p, e-van, c, modnctx);
    }else{
        fmpz_mod_poly_set_coeff_ui(p, e-van, c, modnctx);
    }
    // SAVE TIME BY not always normalizing.
    //normalize();
}
void MySeriesTruncMod::updateCoeffAdd(fmpz_t c, int e){
    if(e>=getTrunc()){cout<<"ERROR dawesdgag: Attempt to update coeff beyond valid trunc. Abort.\n";exit(1);}
    fmpz_t tmp; fmpz_init(tmp);
    getCoeff(tmp, e);
    fmpz_add(tmp,tmp,c);
    setCoeff(tmp, e);
    fmpz_clear(tmp);
    // SAVE TIME BY not always normalizing.
}
void MySeriesTruncMod::printstdout(const char *x){
    if(van==0){
        fmpz_mod_poly_print_pretty(p, x, modnctx);
    }else{
        fmpz_mod_poly_t temp;
        fmpz_mod_poly_init(temp, modnctx);
        getFullPoly(temp);
        fmpz_mod_poly_print_pretty(temp, x, modnctx);
        fmpz_mod_poly_clear(temp, modnctx);
    }
}
void MySeriesTruncMod::printFile(FILE *outf, const char *x){
    if(van==0){
        fmpz_mod_poly_fprint_pretty(outf, p, x, modnctx);
    }else{
        fmpz_mod_poly_t temp;
        fmpz_mod_poly_init(temp, modnctx);
        getFullPoly(temp);
        fmpz_mod_poly_fprint_pretty(outf, temp, x, modnctx);
        fmpz_mod_poly_clear(temp, modnctx);
    }
}
void MySeriesTruncMod::printstdout(){
    if(inSqrtQ){
        printstdout("(q^(1/2))");
    }else{
        printstdout("q");
   }
}
void MySeriesTruncMod::printstdoutWithTruncationOrder(){
    if(inSqrtQ){
        printstdout("(q^(1/2))");
    }else{
        printstdout("q");
   }
   cout<<"+O(q^"<<getTrunc()<<")";
}
void MySeriesTruncMod::printFile(FILE *outf){
    if(inSqrtQ){
        printFile(outf, "(q^(1/2))");
    }else{
        printFile(outf, "q");
   }
}

void MySeriesTruncMod::multiplyByPowerOfQandChangeTrunc(int e){
    if(e<0){
        cout<<"ERROR 834502 in MySeriesTruncMod. exponent e = "<<e<<" should be nonnegative. Abort.";
        exit(1);
    }
    if(e==0)return;
    van = van +e;
}

MySeriesTruncMod*  MySeriesTruncMod::expandExpBy(int factor, int newupto){
    MySeriesTruncMod* ans=new MySeriesTruncMod(modn, newupto,inSqrtQ);
    fmpz_t tmp;fmpz_init(tmp);
    int i, newi;
    for(i=0;i<getTrunc();i++){
        newi=i*factor;
        if(newi<=newupto){
            getCoeff(tmp,i);
            fmpz_mod_poly_set_coeff_fmpz(ans->p,newi,tmp, modnctx);
        }
    }
    fmpz_clear(tmp);
    ans->normalize();
    return ans;
}
MySeriesTruncMod*  MySeriesTruncMod::contractExpBy(int factor){
    MySeriesTruncMod* ans=new MySeriesTruncMod(modn, (getTrunc()-1)/factor,inSqrtQ);
    fmpz_t tmp;fmpz_init(tmp);
    int i, newi;
    for(i=0;i<getTrunc();i++){
        if((i%factor)==0){
            newi=i/factor;
            getCoeff(tmp,i);
            fmpz_mod_poly_set_coeff_fmpz(ans->p,newi,tmp, modnctx);
        }else{ //Check that coefficient is zero!
            getCoeff(tmp,i);
            if(!fmpz_is_zero(tmp)){
                cout<<"ERROR in contractExpBy (factor="<<factor<<"): ";
                printstdout();
                cout<<"\n";
                exit(1);
            }
        }
    }
    ans->normalize();
    fmpz_clear(tmp);
    return ans;
}

void MySeriesTruncMod::setUp(){
    fmpz_init_set_ui(ulongmax,2);
    fmpz_pow_ui(ulongmax,ulongmax,sizeof(ulong)*8-2);
}
MySeriesTruncMod* MySeriesTruncMod::pow(MySeriesTruncMod* f, fmpz_t e){
    f->normalize();
    if(!ulongMaxSetUp){
        ulongMaxSetUp=1;
        setUp();
    }
    if(fmpz_cmp_si(e,0)<0){
        if(f->van!=0){cout<<"Can't raise to negative power when van("<<f->van<<") > 0. Abort.\n";exit(1);}
        fmpz_t nege;
        fmpz_init(nege);
        fmpz_neg(nege, e);
        MySeriesTruncMod* ans = MySeriesTruncMod::pow(f, nege);
        fmpz_mod_poly_inv_series(ans->p, ans->p, ans->len, f->modnctx);
        fmpz_clear(nege);
        return ans;
    }
    //Now assume e>0
    if(fmpz_cmp(e, ulongmax)<0){
        return MySeriesTruncMod::pow(f, fmpz_get_si(e));
    }
    fmpz_t q, r;
    fmpz_init(q);fmpz_init(r);
    fmpz_fdiv_qr(q, r, e, ulongmax);
    MySeriesTruncMod *ans;
    ans=pow(f, q);
    fmpz_mod_poly_pow_trunc_binexp(ans->p, ans->p, fmpz_get_ui(ulongmax),f->len, f->modnctx);
    ans->van = fmpz_get_ui(ulongmax)*ans->van;

    fmpz_mod_poly_t tmp;
    fmpz_mod_poly_init(tmp, f->modnctx);
    fmpz_mod_poly_pow_trunc_binexp(tmp, f->p, fmpz_get_ui(r),f->len, f->modnctx);
    fmpz_mod_poly_mullow(ans->p, ans->p, tmp, f->len, f->modnctx);
    ans->van = ans->van + fmpz_get_ui(r) * f->van;

    fmpz_clear(q); fmpz_clear(r);fmpz_mod_poly_clear(tmp, f->modnctx);
    return ans;
}


void MySeriesTruncMod::multiplyByPower(MySeriesTruncMod* f, fmpz_t e){
    checkSameinSqrtQ(f,this);
    MySeriesTruncMod* tmp = MySeriesTruncMod::pow(f,e);
    multiplyBy(tmp);
    delete tmp;
}

void MySeriesTruncMod::clearSqrts(){ //convert inSqrtQ=1 to 0.
    if(!inSqrtQ){return;}
    int i;
    fmpz_t tmp;fmpz_init(tmp);
    for(i=1;i<getTrunc();i+=2){
        getCoeff(tmp,i);
        if(!fmpz_is_zero(tmp)){
            cout<<"ERROR 125sd23 in clearSqrts: cannot clearSqrts on polynomial:\n";
            printstdout();
            cout<<". Abort.\n";
            exit(1);
        }
    }
    for(i=2;i<getTrunc();i+=2){
        getCoeff(tmp,i);
        fmpz_mod_poly_set_coeff_fmpz(p,i/2,tmp, modnctx);
    }
    len=(len+1)/2;
    van=0;
    fmpz_mod_poly_truncate(p, len, modnctx);
    inSqrtQ=0;
    normalize();
    fmpz_clear(tmp);
}

void MySeriesTruncMod::shift(int up){
    normalize();
    if(up+van>=0){
        van=van+up;
        return;
    }
    cout<<"ERROR dg8isasugdkgk in shift. (van,len,up)=("<<van<<","<<len<<","<<up<<").  No precision left. \n"; exit(1);
    // DON'T SHIFT OFF PRECISION. fmpz_mod_poly_shift_right(p,p,-(up+origvan), modnctx);
}

MySeriesTruncMod* MySeriesTruncMod::makeCopy(){
    MySeriesTruncMod* newg = new MySeriesTruncMod(modn, getTrunc()-1, inSqrtQ);
    fmpz_mod_poly_set(newg->p, p, modnctx);
    newg->van=van;
    newg->len=len;
    newg->normalized=normalized;
    return newg;
}

MySeriesTruncMod* MySeriesTruncMod::divide(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkSameinSqrtQ(f,g);
    f->normalize();
    g->normalize();
    if(g->couldBeZero()){
        cout<<"ERROR 2883fafdas: attempted division by possibly zero polynomial (van,len)=("<<g->van<<","<<g->len<<"). Abort.\n"; exit(1);
    }
    int fVanishingOrder=f->getVanishingOrderNoSentinnel();
    int gVanishingOrder=g->getVanishingOrderNoSentinnel();
    if(gVanishingOrder>fVanishingOrder){
        cout<<"ERROR dgas2t23: attempted division by polynomial with greater vanishing order. Abort.\n";
        cout<<"(fvan,gvan)=("<<fVanishingOrder<<","<<gVanishingOrder<<").\n";
        exit(1);
    }
    if(0){
      cout<<"divide: numer=";
      f->printstdout();
      cout<<"\n denom=";
      g->printstdout();
      cout<<"\n";
      cout<<"(fvan,gvan)=("<<fVanishingOrder<<","<<gVanishingOrder<<").\n";
      cout<<"f(van,len)=("<<f->van<<","<<f->len<<")\n";
      cout<<"g(van,len)=("<<g->van<<","<<g->len<<")\n";
    }
    MySeriesTruncMod* newg = g->makeCopy();
    newg->shift(-gVanishingOrder);
    MySeriesTruncMod* ans = f->makeCopy();
    ans->shift(-gVanishingOrder);
    ans->multiplyByPower(newg,-1);
    delete newg;
    return ans;
}
void MySeriesTruncMod::divideBy(MySeriesTruncMod* g){
    checkSameinSqrtQ(this,g);
    normalize();
    g->normalize();
    if(g->couldBeZero()){
        cout<<"ERROR 4t8gasgkkdk: attempted division by possibly zero polynomial (van,len)=("<<g->van<<","<<g->len<<"). Abort.\n";
        cout<<"numer=";
        printstdoutWithTruncationOrder();
        cout<<"\n";
        cout<<"denom=";
        g->printstdoutWithTruncationOrder();
        cout<<"\n";
        exit(1);
    }
    int fVanishingOrder=getVanishingOrderNoSentinnel();
    int gVanishingOrder=g->getVanishingOrderNoSentinnel();
    if(0){
      cout<<"divideBy: numer=";
      printstdout();
      cout<<"\n denom=";
      g->printstdout();
      cout<<"\n";
      cout<<"(fvan,gvan)=("<<fVanishingOrder<<","<<gVanishingOrder<<").\n";
      cout<<"f(van,len)=("<<van<<","<<len<<")\n";
      cout<<"g(van,len)=("<<g->van<<","<<g->len<<")\n";
    }

    if(gVanishingOrder>fVanishingOrder){
        cout<<"ERROR hhfdhfdhsgfgf: attempted division by polynomial with greater vanishing order. Abort.\n";
        cout<<"(fvan,gvan)=("<<fVanishingOrder<<","<<gVanishingOrder<<").\n";
        exit(1);
    }
    MySeriesTruncMod* newg = g->makeCopy();
    newg->shift(-gVanishingOrder);
    shift(-gVanishingOrder);
    multiplyByPower(newg,-1);
    delete newg;
}



