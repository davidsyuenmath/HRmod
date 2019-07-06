#include "MySeriesTruncMod.h"

int MySeriesTruncMod::checkSameTrunc=1;

MySeriesTruncMod::MySeriesTruncMod(fmpz_t modn0, int upto){
    fmpz_init(modn);
    fmpz_set(modn, modn0);
    fmpz_mod_poly_init(p, modn);
    fmpz_mod_poly_zero(p);
    trunc=upto+1;
    inSqrtQ=0;
}
MySeriesTruncMod::MySeriesTruncMod(int modn0, int upto){
    fmpz_init(modn);
    fmpz_set_si(modn, modn0);
    fmpz_mod_poly_init(p, modn);
    fmpz_mod_poly_zero(p);
    trunc=upto+1;
    inSqrtQ=0;
}
MySeriesTruncMod::MySeriesTruncMod(int modn0, int upto, int inSqrtQ0){
    fmpz_init(modn);
    fmpz_set_si(modn, modn0);
    fmpz_mod_poly_init(p, modn);
    fmpz_mod_poly_zero(p);
    trunc=upto+1;
    inSqrtQ=inSqrtQ0;
}
MySeriesTruncMod::MySeriesTruncMod(fmpz_t modn0, int upto, int inSqrtQ0){
    fmpz_init(modn);
    fmpz_set(modn, modn0);
    fmpz_mod_poly_init(p, modn);
    fmpz_mod_poly_zero(p);
    trunc=upto+1;
    inSqrtQ=inSqrtQ0;
}

MySeriesTruncMod::MySeriesTruncMod(MySeriesTruncMod*f){
    fmpz_init(modn);
    fmpz_set(modn, f->modn);
    fmpz_mod_poly_init(p, modn);
    fmpz_mod_poly_set(p, f->p);
    trunc=f->trunc;
    inSqrtQ=f->inSqrtQ;
}

MySeriesTruncMod:: ~MySeriesTruncMod(){//destructor
    fmpz_clear(modn);
    fmpz_mod_poly_clear(p);

}




void MySeriesTruncMod::checkOperands(MySeriesTruncMod *f, MySeriesTruncMod *g){
    if(f->trunc!=g->trunc){
        if(checkSameTrunc){
            cout<<"ERROR 2411A in MySeriesTruncMod, attempt to operate with two different truncs:"
                <<f->trunc<<","<<g->trunc<<". Abort.\n";
            exit(1);
        }else{
            if(f->trunc<g->trunc){g->trunc=f->trunc;}else{f->trunc=g->trunc;}
        }
    }
    checkSameinSqrtQ(f,g);
}
void MySeriesTruncMod::checkSameinSqrtQ(MySeriesTruncMod *f, MySeriesTruncMod *g){
    if(f->inSqrtQ!=g->inSqrtQ){
        cout<<"ERROR rh4572253 in MySeriesTruncMod, attempt to operate with two different inSqrtQ:"
                <<f->inSqrtQ<<","<<g->inSqrtQ<<". Abort.\n";
        exit(1);
    }
}


MySeriesTruncMod* MySeriesTruncMod::monicBinomial(fmpz_t modn0, int upto, fmpz_t oneCoeff, int e){
    //initializes to 1 + c * x^e
    MySeriesTruncMod* ans = new MySeriesTruncMod(modn0, upto);
    fmpz_mod_poly_set_coeff_ui(ans->p, 0, 1);
    if(e<=upto){fmpz_mod_poly_set_coeff_fmpz(ans->p, e, oneCoeff);}
    if(e<1){cout<<"ERROR czvzv in MySeriesTruncMod. Abort.\n";exit(1);}
    return ans;
}  //Initializes to 1 + c*x^e
int MySeriesTruncMod::setMonicBinomialBPSpecial(fmpz_t oneCoeff, int e){
    //initializes to 1 + c * x^e if e>=0   return 0
    //else  x^(-e) + c   if e<0            return e
    //cout<<"setMonicBinomialBPSpecial (oneCoeff, exponent, trunc, modn)=("<<fmpz_get_si(oneCoeff)
    //    <<","<<e<<","<<trunc<<","<<fmpz_get_si(modn)<<") returning...\n ";
    if(inSqrtQ){cout<<"Error 4t8asf0 in setMonicBinomialBPSpecial. Abort.\n";exit(1);}
    fmpz_mod_poly_zero(p);
    //cout<<"[STEP zeroing out] ";printstdout();cout<<"\n";
    if(e>0){
        fmpz_mod_poly_set_coeff_ui(p, 0, 1);
        //cout<<"[STEP set coeff 0 to 1] ";printstdout();cout<<"\n";
        if(e<trunc){
            fmpz_mod_poly_set_coeff_fmpz(p, e, oneCoeff);
            //cout<<"[STEP set coeff e to oneCoeff] ";printstdout();cout<<"\n";
        }
        //cout<<"[expon>0] ";printstdout();cout<<"\n";
        return 0;
    }else if(e<0){
        fmpz_mod_poly_set_coeff_fmpz(p, 0, oneCoeff);
        //cout<<"oneCoeff="<<fmpz_get_si(oneCoeff)<<"\n";
        //cout<<"[STEP Z2 set coeff 0 to oneCoeff] ";printstdout();cout<<"\n";
        if(-e<trunc){
            fmpz_mod_poly_set_coeff_ui(p, -e, 1);
            //cout<<"[STEP set coeff -e to 1] ";printstdout();cout<<"\n";
        }
        //cout<<"[expon<0] ";printstdout();cout<<"\n";
        return e;
    }else{
        fmpz_t temp;
        fmpz_init_set_si(temp, 1);
        fmpz_add(temp, temp, oneCoeff);
        fmpz_mod_poly_set_coeff_fmpz(p, 0, temp);
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
    if(e<=upto){fmpz_mod_poly_set_coeff_fmpz(ans->p, e, oneCoeff);}
    return ans;
}  //Initializes to c*x^e

MySeriesTruncMod* MySeriesTruncMod::add(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkOperands(f,g);
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, f->trunc-1, f->inSqrtQ);
    fmpz_mod_poly_add(ans->p,f->p,g->p);
    return ans;
}
void MySeriesTruncMod::addAndAdjustTrunc(MySeriesTruncMod* g){
    fmpz_mod_poly_add(p,p,g->p);
    if(trunc!=g->trunc){
        int e=min(trunc,g->trunc);
        fmpz_mod_poly_truncate(p, e);
        trunc=e;
    }

}
void MySeriesTruncMod::addBy(MySeriesTruncMod* g){
    checkOperands(this, g);
    fmpz_mod_poly_add(p,p,g->p);
}
MySeriesTruncMod* MySeriesTruncMod::scalarMultiplity(MySeriesTruncMod* g, fmpz_t a){
    MySeriesTruncMod* ans = new MySeriesTruncMod(g->modn, g->trunc-1, g->inSqrtQ);
    fmpz_mod_poly_scalar_mul_fmpz(ans->p, g->p, a);
    return ans;
}
void MySeriesTruncMod::scalarMultiplityBy(fmpz_t a){
    fmpz_mod_poly_scalar_mul_fmpz(p, p, a);
}
MySeriesTruncMod* MySeriesTruncMod::multiply(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkOperands(f,g);
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, f->trunc-1, f->inSqrtQ);
    fmpz_mod_poly_mullow(ans->p, f->p, g->p, f->trunc);
    return ans;
}
void MySeriesTruncMod::multiplyBy(MySeriesTruncMod* g){
    checkOperands(this, g);
    fmpz_mod_poly_mullow(p, p, g->p,trunc);
}
MySeriesTruncMod* MySeriesTruncMod::pow(MySeriesTruncMod* f, slong e){
    MySeriesTruncMod* ans = new MySeriesTruncMod(f->modn, f->trunc-1, f->inSqrtQ);
    if(e>=0){
        fmpz_mod_poly_pow_trunc_binexp(ans->p, f->p, e,f->trunc);
    }else{
        fmpz_mod_poly_pow_trunc_binexp(ans->p, f->p, -e,f->trunc);
        fmpz_mod_poly_inv_series(ans->p, ans->p, ans->trunc);
    }
    return ans;
}
void MySeriesTruncMod::multiplyByPower(MySeriesTruncMod* f, slong e){
    checkOperands(this, f);

    if(e==0)return;
    if(e==1){
        multiplyBy(f);return;
    }
    MySeriesTruncMod* tmp = MySeriesTruncMod::pow(f,e);
    multiplyBy(tmp);
    delete tmp;
}

void MySeriesTruncMod::getCoeff(fmpz_t c, int e){
    if(e>=trunc){cout<<"ERROR 32tfgdf: Request for coeff is beyond valid trunc. Abort.\n";exit(1);}
    fmpz_mod_poly_get_coeff_fmpz(c, p, e);
}
void MySeriesTruncMod::setCoeff(fmpz_t c, int e){
    if(e>=trunc){cout<<"ERROR adfshye45y: Attempt to set coeff beyond valid trunc. Abort.\n";exit(1);}
    fmpz_mod_poly_set_coeff_fmpz(p, e, c);
}
void MySeriesTruncMod::setCoeff(ulong c, int e){
    if(e>=trunc){cout<<"ERROR jtrj4532: Attempt to set coeff beyond valid trunc. Abort.\n";exit(1);}
    fmpz_mod_poly_set_coeff_ui(p, e, c);
}
void MySeriesTruncMod::updateCoeffAdd(fmpz_t c, int e){
    if(e>=trunc){cout<<"ERROR dawesdgag: Attempt to update coeff beyond valid trunc. Abort.\n";exit(1);}
    fmpz_t tmp; fmpz_init(tmp);
    fmpz_mod_poly_get_coeff_fmpz(tmp, p, e);
    fmpz_add(tmp,tmp,c);
    fmpz_mod_poly_set_coeff_fmpz(p, e, tmp);
    fmpz_clear(tmp);
}
void MySeriesTruncMod::printstdout(const char *x){
    fmpz_mod_poly_print_pretty(p, x);
}
void MySeriesTruncMod::printFile(FILE *outf, const char *x){
    fmpz_mod_poly_fprint_pretty(outf, p, x);
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
   cout<<"+O(q^"<<trunc<<")";
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
    trunc = trunc+e;
    fmpz_mod_poly_shift_left(p, p, e);
}

MySeriesTruncMod*  MySeriesTruncMod::expandExpBy(int factor, int newupto){
    MySeriesTruncMod* ans=new MySeriesTruncMod(modn, newupto,inSqrtQ);
    fmpz_t tmp;fmpz_init(tmp);
    int i, newi;
    for(i=0;i<trunc;i++){
        newi=i*factor;
        if(newi<=newupto){
            fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
            fmpz_mod_poly_set_coeff_fmpz(ans->p,newi,tmp);
        }
    }
    fmpz_clear(tmp);
    return ans;
}
MySeriesTruncMod*  MySeriesTruncMod::contractExpBy(int factor){
    MySeriesTruncMod* ans=new MySeriesTruncMod(modn, (trunc-1)/factor,inSqrtQ);
    fmpz_t tmp;fmpz_init(tmp);
    int i, newi;
    for(i=0;i<trunc;i++){
        if((i%factor)==0){
            newi=i/factor;
            fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
            fmpz_mod_poly_set_coeff_fmpz(ans->p,newi,tmp);
        }else{ //Check that coefficient is zero!
            fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
            if(!fmpz_is_zero(tmp)){
                cout<<"ERROR in contractExpBy (factor="<<factor<<"): ";
                printstdout();
                cout<<"\n";
                exit(1);
            }
        }
    }
    fmpz_clear(tmp);
    return ans;
}

void MySeriesTruncMod::setUp(){
    fmpz_init_set_ui(ulongmax,2);
    fmpz_pow_ui(ulongmax,ulongmax,sizeof(ulong)*8-2);
}
MySeriesTruncMod* MySeriesTruncMod::pow(MySeriesTruncMod* f, fmpz_t e){
    if(!ulongMaxSetUp){
        ulongMaxSetUp=1;
        setUp();
    }
    if(fmpz_cmp_si(e,0)<0){
        fmpz_t nege;
        fmpz_init(nege);
        fmpz_neg(nege, e);
        MySeriesTruncMod* ans = MySeriesTruncMod::pow(f, nege);
        fmpz_mod_poly_inv_series(ans->p, ans->p, ans->trunc);
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
    fmpz_mod_poly_pow_trunc_binexp(ans->p, ans->p, fmpz_get_ui(ulongmax),f->trunc);

    fmpz_mod_poly_t tmp;
    fmpz_mod_poly_init(tmp, f->modn);
    fmpz_mod_poly_pow_trunc_binexp(tmp, f->p, fmpz_get_ui(r),f->trunc);
    fmpz_mod_poly_mullow(ans->p, ans->p, tmp, f->trunc);

    fmpz_clear(q); fmpz_clear(r);fmpz_mod_poly_clear(tmp);
    return ans;
}


void MySeriesTruncMod::multiplyByPower(MySeriesTruncMod* f, fmpz_t e){
    checkOperands(this, f);
    MySeriesTruncMod* tmp = MySeriesTruncMod::pow(f,e);
    multiplyBy(tmp);
    delete tmp;
}

void MySeriesTruncMod::clearSqrts(){ //convert inSqrtQ=1 to 0.
    if(!inSqrtQ){return;}
    int i;
    fmpz_t tmp;fmpz_init(tmp);
    for(i=1;i<trunc;i+=2){
        fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
        if(!fmpz_is_zero(tmp)){
            cout<<"ERROR 125sd23 in clearSqrts: cannot clearSqrts on polynomial:\n";
            printstdout();
            cout<<". Abort.\n";
            exit(1);
        }
    }
    for(i=2;i<trunc;i+=2){
        fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
        fmpz_mod_poly_set_coeff_fmpz(p,i/2,tmp);
    }
    trunc=(trunc+1)/2;
    fmpz_mod_poly_truncate(p, trunc);
    inSqrtQ=0;
    fmpz_clear(tmp);
}

int MySeriesTruncMod::getVanishingOrder(){
    int i;
    fmpz_t tmp;fmpz_init(tmp);
    for(i=1;i<trunc;i++){
        fmpz_mod_poly_get_coeff_fmpz(tmp,p,i);
        if(!fmpz_is_zero(tmp)){
            fmpz_clear(tmp);
            return i;
        }
    }
    fmpz_clear(tmp);
    return -1; //Sentinnel for zero polynomial.
}

void MySeriesTruncMod::shift(int up){
    trunc=trunc+up;
    if(up>=0){
        fmpz_mod_poly_shift_left(p, p, up);
        return;
    }
    if(trunc<=0){
        trunc=0;
    }
    fmpz_mod_poly_shift_right(p,p,-up);
}

MySeriesTruncMod* MySeriesTruncMod::makeCopy(){
    MySeriesTruncMod* newg = new MySeriesTruncMod(modn, trunc-1, inSqrtQ);
    fmpz_mod_poly_set(newg->p, p);
    return newg;
}

MySeriesTruncMod* MySeriesTruncMod::divide(MySeriesTruncMod* f, MySeriesTruncMod* g){
    checkSameinSqrtQ(f,g);
    int fVanishingOrder=f->getVanishingOrder();
    int gVanishingOrder=g->getVanishingOrder();
    if(gVanishingOrder==-1){
        cout<<"ERROR 32t88fsd: attempted division by zero polynomial. Abort.\n"; exit(1);
    }
    if(gVanishingOrder>fVanishingOrder){
        cout<<"ERROR dgas2t23: attempted division by polynomial with greater vanishing order. Abort.\n"; exit(1);
    }
    MySeriesTruncMod* newg = g->makeCopy();
    newg->shift(-gVanishingOrder);
    MySeriesTruncMod* ans = f->makeCopy();
    ans->shift(-gVanishingOrder);
    ans->multiplyByPower(newg,-1);
    return ans;
}
void MySeriesTruncMod::divideBy(MySeriesTruncMod* g){
    checkSameinSqrtQ(this,g);
    int fVanishingOrder=getVanishingOrder();
    int gVanishingOrder=g->getVanishingOrder();
    if(gVanishingOrder==-1){
        cout<<"ERROR dfaewtwe: attempted division by zero polynomial. Abort.\n"; exit(1);
    }
    if(gVanishingOrder>fVanishingOrder)if(fVanishingOrder!=-1){
        cout<<"ERROR 4hdgdfga: attempted division by polynomial with greater vanishing order. Abort.\n"; exit(1);
    }
    MySeriesTruncMod* newg = g->makeCopy();
    newg->shift(-gVanishingOrder);
    shift(-gVanishingOrder);
    multiplyByPower(newg,-1);
}



