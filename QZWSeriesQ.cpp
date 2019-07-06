#include "QZWSeriesQ.h"

QZWSeriesQ::QZWSeriesQ(){ //constructor does nothing
  coeffArraySize=0;
  trunc=-1;
  defaultQTrunc=0;
  fmpz_init(denom);
  fmpz_set_ui(denom,1);
}
QZWSeriesQ::QZWSeriesQ(int wTrunc, int qTrunc){
  coeffArraySize=wTrunc+1;
  trunc = wTrunc;
  defaultQTrunc=qTrunc;
  coeff=new QZSeriesWH*[coeffArraySize];
  int i;
  for(i=0;i<coeffArraySize;i++){
    coeff[i]=new QZSeriesWH(defaultQTrunc);
  fmpz_init(denom);
  fmpz_set_ui(denom,1);
  }
}

QZWSeriesQ:: ~QZWSeriesQ(){//destructor
  destroy();
}

void QZWSeriesQ::destroy(){
    //cout<<"QZWSeriesQ destroy activated. (trunc,coeffArraySize)=("<<trunc<<
    //     ","<<coeffArraySize<<")\n";
  int i;
  if(coeffArraySize>0){
    for(i=0;i<coeffArraySize;i++){
        delete coeff[i];
    }
    delete[] coeff;
  }
  coeffArraySize=0;
  trunc=-1;
  fmpz_clear(denom);
}
void QZWSeriesQ::resizeArray(int newsize){
  if(newsize<coeffArraySize)return;
  QZSeriesWH** newArray=new QZSeriesWH*[newsize];
  int i;
  for(i=0;i<coeffArraySize;i++)newArray[i]=coeff[i];
  for(i=coeffArraySize;i<newsize;i++)newArray[i]=new QZSeriesWH(defaultQTrunc);
  delete[] coeff;
  coeff=newArray;
  coeffArraySize=newsize;
}



void QZWSeriesQ::shiftExponents(int shift){

  if(shift>0){//multiply by w^shift
    resizeArray(trunc+shift+1);
    int i;
    for(i=trunc+1;i<trunc+shift;i++)delete coeff[i];
    for(i=trunc;i>=0;i--)coeff[i+shift]=coeff[i];
    for(i=0;i<shift;i++)coeff[i]=new QZSeriesWH(defaultQTrunc);
    trunc=trunc+shift;
  }else if(shift<0){
    int i;
    for(i=0;i<-shift;i++)delete coeff[i];
    for(i=-shift;i<=trunc;i++)coeff[i+shift]=coeff[i];
    for(i=trunc+shift+1;i<=trunc;i++)coeff[i]=new QZSeriesWH(defaultQTrunc);
    trunc=trunc+shift;
  }
}

string QZWSeriesQ::getString(){
    return getString("q", "z", "w");
}

string QZWSeriesQ::getString(string qvar, string zvar, string wvar){
  string ans, str;
  stringstream ss;
  slong  n;
  ans="1/";
  ans=ans+fmpz_get_str(NULL,10,denom);
  ans=ans+"(";
  int isEmpty=1;
  for(n=0;n<=trunc;n++){
   if(!(coeff[n]->isZero())){
    if(isEmpty){isEmpty=0;}
    else{ans = ans + " + ";}
    ans = ans + "("+(coeff[n]->getString(qvar, zvar)) + ")";
    if(n==0){

    }else{
      ans = ans + "*" + wvar;
      if(n!=1){
        ans = ans + "^";
        ss.clear();
        ss<<n;
        ss>>str;
        //to_string(n) works on some compilers.
        ans = ans + str;
      }
    }
   }
  }
  ans=ans+")";
  if(isEmpty){ans="0";}
  return ans;
}


QZWSeriesQ* QZWSeriesQ::copy(){
  QZWSeriesQ* ans=new QZWSeriesQ();
  ans->coeffArraySize=trunc+1;
  ans->trunc=trunc;
  ans->defaultQTrunc=defaultQTrunc;
  fmpz_set(ans->denom,denom);
  ans->coeff=new QZSeriesWH*[coeffArraySize];
  int i;
  for(i=0;i<=trunc;i++){
    ans->coeff[i]=coeff[i]->copy();
  }
  return ans;
}


QZWSeriesQ* QZWSeriesQ::add(QZWSeriesQ* a, QZWSeriesQ* b){
  QZWSeriesQ *ans = a->copy();
  ans->addWith(b);
  return ans;
}
void QZWSeriesQ::addWith(QZWSeriesQ* b){
  //decide if we need to shorten this series
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  if(fmpz_equal(denom,b->denom)){
    for(i=0;i<=trunc;i++)coeff[i]->addWith(b->coeff[i]);
  }else{
    fmpz_t lcm, m, mb;
    fmpz_init(lcm);fmpz_init(m);fmpz_init(mb);
    fmpz_lcm(lcm, denom, b->denom);
    fmpz_divexact(m, lcm, denom);
    fmpz_divexact(mb, lcm, b->denom);
    multiplyByScalar(m);
    for(i=0;i<=trunc;i++)coeff[i]->addScalarMultipleWith(mb, b->coeff[i]);
    fmpz_set(denom, lcm);
    fmpz_clear(lcm);fmpz_clear(m);fmpz_clear(mb);
  }
}

void QZWSeriesQ::multiplyByScalar(fmpz_t m){
    int i;
    for(i=0;i<=trunc;i++)coeff[i]->multiplyByScalar(m);
}

void QZWSeriesQ::subtractWith(QZWSeriesQ* b){
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  if(fmpz_equal(denom,b->denom)){
    for(i=0;i<=trunc;i++)coeff[i]->subtractWith(b->coeff[i]);
  }else{
    fmpz_t lcm, m, mb;
    fmpz_init(lcm);fmpz_init(m);fmpz_init(mb);
    fmpz_lcm(lcm, denom, b->denom);
    fmpz_divexact(m, lcm, denom);
    fmpz_divexact(mb, lcm, b->denom);
    fmpz_neg(mb,mb);
    multiplyByScalar(m);
    for(i=0;i<=trunc;i++)coeff[i]->addScalarMultipleWith(mb, b->coeff[i]);
    fmpz_set(denom, lcm);
    fmpz_clear(lcm);fmpz_clear(m);fmpz_clear(mb);
  }
}


QZWSeriesQ* QZWSeriesQ::multiply(QZWSeriesQ* a, QZWSeriesQ *b){
  int newtrunc = a->trunc;
  if(a->trunc > b->trunc){newtrunc=b->trunc;}
  QZWSeriesQ *ans = new QZWSeriesQ(newtrunc, a->defaultQTrunc);
  fmpz_mul(ans->denom,a->denom,b->denom);
  QZSeriesWH* tmp;
  int n, i,j;
  for(n=0;n<=newtrunc;n++){
    for(i=0;i<=n;i++){
      j=n-i;
      tmp = QZSeriesWH::multiply(a->coeff[i], b->coeff[j]);
      ans->coeff[n]->addWith(tmp);
      delete tmp;
    }
  }
  return ans;
}


void QZWSeriesQ::multiplyWith(QZWSeriesQ* b){
  QZWSeriesQ *ans = multiply(this, b);
  //Memory leak fixed 11-22-2016 Need the following destroy()
  destroy();
  coeffArraySize = ans->coeffArraySize;
  coeff = ans->coeff;
  fmpz_init(denom);
  fmpz_set(denom,ans->denom);
  trunc = ans->trunc;
  defaultQTrunc=ans->defaultQTrunc;
  ans->coeffArraySize=0;
  delete ans;
}


QZWSeriesQ* QZWSeriesQ::addScalarMultiple(QZWSeriesQ *a, fmpz_t m, QZWSeriesQ* b){
  QZWSeriesQ *ans = a->copy();
  ans->addScalarMultipleWith(m, b);
  return ans;
}
void QZWSeriesQ::addScalarMultipleWith(fmpz_t c, QZWSeriesQ* b){
  //decide if we need to shorten this series
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  if(fmpz_equal(denom,b->denom)){
    for(i=0;i<=trunc;i++)coeff[i]->addScalarMultipleWith(c,b->coeff[i]);
  }else{
    fmpz_t lcm, m, mb;
    fmpz_init(lcm);fmpz_init(m);fmpz_init(mb);
    fmpz_lcm(lcm, denom, b->denom);
    fmpz_divexact(m, lcm, denom);
    fmpz_divexact(mb, lcm, b->denom);
    fmpz_mul(mb,mb,c);
    multiplyByScalar(m);
    for(i=0;i<=trunc;i++)coeff[i]->addScalarMultipleWith(mb, b->coeff[i]);
    fmpz_set(denom, lcm);
    fmpz_clear(lcm);fmpz_clear(m);fmpz_clear(mb);
  }
}


void QZWSeriesQ::multiplyByQZSeriesWH(QZSeriesWH* b){
  int i;
  for(i=0;i<=trunc;i++){
    coeff[i]->multiplyWith(b);
  }
}


bool QZWSeriesQ::divideByDenominator(){ //returns 1 if divide exactly
    bool success=true;
    int i;
    for(i=0;i<=trunc;i++){
        success = success && (coeff[i]->divideByIntegerWithCheck(denom));
    }
    fmpz_set_ui(denom,1);
    return success;
}
void QZWSeriesQ::divideByInteger(fmpz_t c){ //returns 1 if divide exactly
    fmpz_mul(denom,denom,c);
}


void QZWSeriesQ::negate(){
  int i;
  for(i=0;i<=trunc;i++){
    coeff[i]->negate();
  }
}

QZWSeriesQ* QZWSeriesQ::makeBorchLift(QZSeriesWH* psi, QZSeriesWH* tb, int level, int leadingJC, int numJC){
    //cout<<"MBL 1:\n";
    if(Diagnostics::saveDiagnosticFile)psi->saveToFile(Diagnostics::f, "JFpsi");
    QZWSeriesQ* grit = makeGritLift(psi, level, 0, numJC-1);
    if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLift");

    //cout<<"MBL 2:\n";
    grit->negate();
    //if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLiftNeg");

    //cout<<"MBL 3:\n";
    QZWSeriesQ* bp =grit->exp(numJC-1);
    //if(Diagnostics::saveDiagnosticFile)bp->saveToFile(Diagnostics::f, "JFbpExpansionStep1");
    //cout<<"MBL 4:\n";
    bp->multiplyByQZSeriesWH(tb);

    //cout<<"MBL 5:\n";
    cout<<"bp->denom="<<fmpz_get_str(NULL,10,bp->denom)<<"\n";
    //if(Diagnostics::saveDiagnosticFile)bp->saveToFile(Diagnostics::f, "JFbpExpansionStep2");

    if(!(bp->divideByDenominator())){
        cout<<"Yikes, BorchLift not integral!! Abort.\n";
        //cout<<bp->getString()<<"\n";
        exit(1);
    };
        //Make sure the first term is tb, because of possible BS with qTrunc's.
    delete bp->coeff[0]; bp->coeff[0] = tb->copy();

    bp->shiftExponents(leadingJC);
    if(Diagnostics::saveDiagnosticFile)bp->saveToFile(Diagnostics::f, "JFbpExpansion");

    delete grit;
    return bp;
}
QZWSeriesQ* QZWSeriesQ::makeGritLift(QZSeriesWH* phi, int level, int weight, int numJC){
    QZWSeriesQ* grit = new QZWSeriesQ(numJC, phi->getMaxExponent());
    //if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLiftStep1");

    int i;
    for(i=1;i<=numJC;i++){
        delete grit->coeff[i];
        grit->coeff[i] = phi->copy();
        grit->coeff[i]->applyUp(level, weight, i);
    }
    //if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLiftStep2");

    if(weight==0){
        fmpz_t lcm; fmpz_init(lcm);
        fmpz_t c; fmpz_init(c);
        fmpz_set_ui(lcm,1);
        for(i=1;i<=numJC;i++){fmpz_set_si(c,i);fmpz_lcm(lcm,lcm,c);}
        for(i=1;i<=numJC;i++){
            fmpz_set_si(c,i);
            fmpz_divexact(c, lcm, c);
            //cout<<"(i, lcm, c)=("<<i<<","<<fmpz_get_str(NULL,10,lcm)<<","<<fmpz_get_str(NULL,10,c)<<")\n";
            grit->coeff[i]->multiplyByScalar(c);
        }
        //if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLiftStep3");

        grit->divideByInteger(lcm);
        fmpz_clear(lcm); fmpz_clear(c);
    }
    //if(Diagnostics::saveDiagnosticFile)grit->saveToFile(Diagnostics::f, "JFgritLiftStep4");

    return grit;
}
QZWSeriesQ* QZWSeriesQ::exp(int terms){
    QZWSeriesQ* ans = copy();
    QZWSeriesQ* temp, *power;
    temp =new QZWSeriesQ(trunc, defaultQTrunc);
    delete temp->coeff[0];
    temp->coeff[0] = new QZSeriesWH(defaultQTrunc, 1);
    ans->addWith(temp);
    delete temp;

    //cout<<"EXP 1:\n";

    fmpz_t nfactorial; fmpz_init(nfactorial);

    fmpz_set_ui(nfactorial,2);
    power = multiply(this, this);
    temp =power->copy();
    temp->divideByInteger(nfactorial);
    ans->addWith(temp);
    delete temp;
    //cout<<"EXP 2:\n";

    int n;
    for(n=3;n<=terms;n++){

        fmpz_mul_si(nfactorial, nfactorial,n);
        //cout<<"EXP n, n!: "<<n<<","<<fmpz_get_str(NULL,10,nfactorial)<<"\n";

        power->multiplyWith(this);
        temp=power->copy();
        temp->divideByInteger(nfactorial);
        ans->addWith(temp);
        delete temp;
    }
    //cout<<"EXP 4:\n";

    delete power;     cout<<"EXP 6:\n";

    fmpz_clear(nfactorial);
    //cout<<"EXP 7:\n";

    return ans;
}
int QZWSeriesQ::getFC(fmpz_t fc, int n, int r, int m, int level){
    if(m>trunc){
        cout<<"QZWSeries getFC: (n,r,m)=("<<n<<","<<r<<","<<m<<") exceeds wtrunc of "<<trunc<<"\n";

            return 0;
    }
    if(n>coeff[m]->getMaxExponent()){
        cout<<"QZWSeries getFC: (n,r,m)=("<<n<<","<<r<<","<<m<<
            ") exceeds q-maxExponent of "<<coeff[m]->getMaxExponent()<<"\n";

            return 0;
    }
    coeff[m]->getFC(fc, n, r, m*level);
        cout<<"QZWSeries getFC: (n,r,m)=("<<n<<","<<r<<","<<m<<") = "<<fmpz_get_str(NULL,10,fc)<<"\n";

    return 1;
}
void QZWSeriesQ::saveToFile(ofstream &f, string varname){
    f<<varname<<"=(";
    f<< getString();
    f<<");\n";
    f.flush();
}
