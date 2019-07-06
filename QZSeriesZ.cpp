#include "QZSeriesZ.h"

QZSeriesZ::QZSeriesZ(){ //constructor does nothing
  coeffArraySize=0;
  trunc=-1;
}
QZSeriesZ::QZSeriesZ(int trunc0){
  coeffArraySize=trunc0+1;
  trunc = trunc0;
  coeff=new LaurentZ*[coeffArraySize];
  int i;
  for(i=0;i<coeffArraySize;i++){
    coeff[i]=new LaurentZ(0);
  }
}

QZSeriesZ::QZSeriesZ(int trunc0, int topOne, int topCoeff){
  coeffArraySize=trunc0+1;
  trunc = trunc0;
  coeff=new LaurentZ*[coeffArraySize];
  int i;
  for(i=0;i<coeffArraySize;i++){
    if(i==topOne){
      coeff[i]=new LaurentZ(topCoeff); //Note topCoeff has precedence over constant if topOne==0
    }else if(i==0){
      coeff[i] = new LaurentZ(1);
    }else{
      coeff[i]=new LaurentZ(0);
    }
  }
}

QZSeriesZ* QZSeriesZ::OnePlusSkip(int trunc0, int skip){
  vector<int> c;
  int i;
  for(i=0;i<=trunc0;i++){
    if(i%skip==0){
      c.push_back(1);
    }else{
      c.push_back(0);
    }
  }
  return new QZSeriesZ(trunc0, c);
}

QZSeriesZ:: ~QZSeriesZ(){//destructor
  destroy();
}
void QZSeriesZ::destroy(){
    //cout<<"QZSeriesZ destroy activated.\n";
  int i;
  if(coeffArraySize>0){
    for(i=0;i<coeffArraySize;i++){
        delete coeff[i];
    }
    delete[] coeff;
  }
  coeffArraySize=0;
  trunc=-1;
}

QZSeriesZ::QZSeriesZ(int trunc0, vector<int> c){
  coeffArraySize=trunc0+1;
  trunc=trunc0;
  coeff=new LaurentZ*[coeffArraySize];
  int i;
  for(i=0;i<c.size();i++){
    if(i<coeffArraySize)coeff[i]=new LaurentZ(c[i]);
  }
  for(i=c.size();i<coeffArraySize;i++){
    coeff[i]=new LaurentZ(0);
  }
}
QZSeriesZ* QZSeriesZ::EtaFunctionWithoutFractionalExponent(int trunc0){
  //The efficiency of this can be improved greatly by multiplying vector<slong> first
  //and then using the construction (trunc0, vector<slong>)
  //but let us not worry about that now.
  QZSeriesZ* ans=new QZSeriesZ(trunc0,1,-1);
  int i;
  QZSeriesZ* tmp;
  for(i=2;i<=trunc0;i++){
    tmp=new QZSeriesZ(trunc0,i, -1);
    ans->multiplyWith(tmp);
    delete tmp;
  }
  return ans;
}
QZSeriesZ* QZSeriesZ::InverseEtaFunctionWithoutFractionalExponent(int trunc0){
  //The efficiency of this can be improved greatly by multiplying vector<slong> first
  //and then using the construction (trunc0, vector<slong>)
  //but let us not worry about that now.
  QZSeriesZ* ans=OnePlusSkip(trunc0,1);
  int i;
  QZSeriesZ* tmp;
  for(i=2;i<=trunc0;i++){
    tmp=OnePlusSkip(trunc0,i);
    ans->multiplyWith(tmp);
    delete tmp;
  }
  return ans;
}
//QZSeriesZ* QZSeriesZ::ThetaFunctionWithoutFractionalExponentWithoutBTB(int trunc0, int d){
//    cout<<"Not yet implemented.\n";exit(1);
//}

void QZSeriesZ::resizeArray(int newsize){
  if(newsize<coeffArraySize)return;
  LaurentZ** newArray=new LaurentZ*[newsize];
  int i;
  for(i=0;i<coeffArraySize;i++)newArray[i]=coeff[i];
  for(i=coeffArraySize;i<newsize;i++)newArray[i]=new LaurentZ(0);
  delete[] coeff;
  coeff=newArray;
  coeffArraySize=newsize;
}

void QZSeriesZ::shiftExponents(int shift){

  if(shift>0){//multiply by q^shift
    resizeArray(trunc+shift+1);
    int i;
    for(i=trunc+1;i<trunc+shift;i++)delete coeff[i];
    for(i=trunc;i>=0;i--)coeff[i+shift]=coeff[i];
    for(i=0;i<shift;i++)coeff[i]=new LaurentZ(0);
    trunc=trunc+shift;
  }else if(shift<0){
    int i;
    for(i=0;i<-shift;i++)delete coeff[i];
    for(i=-shift;i<=trunc;i++)coeff[i+shift]=coeff[i];
    for(i=trunc+shift+1;i<=trunc;i++)coeff[i]=new LaurentZ(0);
    trunc=trunc+shift;
  }
}


bool QZSeriesZ::isZero(){

  bool ans=true;
  int i;
  for(i=0;i<=trunc;i++){
    if(! coeff[i]->isZero()){
        ans=false;
        break;
    }
  }
  return ans;
}

string QZSeriesZ::getString(){
    return getString("q", "z");
}

string QZSeriesZ::getString(string qvar, string zvar){
  string ans, str;
  stringstream ss;
  slong  n;
  ans="";
  int isEmpty=1;
  for(n=0;n<=trunc;n++){
   if(!(coeff[n]->isZero())){
    if(isEmpty){isEmpty=0;}
    else{ans = ans + " + ";}
    ans = ans + "("+(coeff[n]->getString(zvar)) + ")";
    if(n==0){

    }else{
      ans = ans + "*" + qvar;
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
  if(isEmpty){ans="0";}
  return ans;
}

QZSeriesZ* QZSeriesZ::copy(){
  QZSeriesZ* ans=new QZSeriesZ();
  ans->coeffArraySize=trunc+1;
  ans->trunc=trunc;
  ans->coeff=new LaurentZ*[coeffArraySize];
  int i;
  for(i=0;i<=trunc;i++){
    ans->coeff[i]=coeff[i]->copy();
  }
  return ans;
}

QZSeriesZ* QZSeriesZ::add(QZSeriesZ* a, QZSeriesZ* b){
  QZSeriesZ *ans = a->copy();
  ans->addWith(b);
  return ans;
}
void QZSeriesZ::addWith(QZSeriesZ* b){
  //decide if we need to shorten this series
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  for(i=0;i<=trunc;i++)coeff[i]->addWith(b->coeff[i]);
}
void QZSeriesZ::subtractWith(QZSeriesZ* b){
  //decide if we need to shorten this series
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  for(i=0;i<=trunc;i++)coeff[i]->subtractWith(b->coeff[i]);
}

QZSeriesZ* QZSeriesZ::multiply(QZSeriesZ* a, QZSeriesZ *b){
  int newtrunc = a->trunc;
  if(a->trunc > b->trunc){newtrunc=b->trunc;}
  QZSeriesZ *ans = new QZSeriesZ(newtrunc);
  LaurentZ* tmp = new LaurentZ();
  int n, i,j;
  for(n=0;n<=newtrunc;n++){
    for(i=0;i<=n;i++){
      j=n-i;
      LaurentZ::multiply(tmp, a->coeff[i], b->coeff[j]);
      ans->coeff[n]->addWith(tmp);
    }
  }
  delete tmp;
  return ans;
}

void QZSeriesZ::multiplyWith(QZSeriesZ* b){
  QZSeriesZ *ans = multiply(this, b);
  //Memory leak fixed 11-22-2016 Need the following destroy()
  destroy();
  coeffArraySize = ans->coeffArraySize;
  coeff = ans->coeff;
  trunc = ans->trunc;
  ans->coeffArraySize=0;
  delete ans;
}


QZSeriesZ* QZSeriesZ::addScalarMultiple(QZSeriesZ *a, fmpz_t m, QZSeriesZ* b){
  QZSeriesZ *ans = a->copy();
  ans->addScalarMultipleWith(m, b);
  return ans;
}
void QZSeriesZ::addScalarMultipleWith(fmpz_t m, QZSeriesZ* b){
  //decide if we need to shorten this series
  if(trunc > b->trunc){trunc=b->trunc;}
  int i;
  for(i=0;i<=trunc;i++)coeff[i]->addScalarMultipleWith(m,b->coeff[i]);
}
QZSeriesZ* QZSeriesZ::divideByMonic(QZSeriesZ* a, QZSeriesZ *b){
  if(b->trunc<0){cout<<"ERROR 38tdfgdg in QZSeriesZ.\n";exit(1);}
  int leadingCoeffIsUnit = b->coeff[0]->isUnit();
  if(leadingCoeffIsUnit!=1){   //This tests for monic (==-1 would be -1).
    cout<<"ERROR in QZSeriesZ::divideByMonic. Divisor not monic.\n";
    cout<<"denom = "<<b->getString()<<"\n";
    exit(1);
  }
  int newtrunc = a->trunc;
  if(a->trunc > b->trunc){newtrunc=b->trunc;}
  QZSeriesZ *ans = new QZSeriesZ();
  ans->coeffArraySize=newtrunc+1;
  ans->trunc=newtrunc;
  ans->coeff=new LaurentZ*[ans->coeffArraySize];
  LaurentZ* newcoeff, *tmp;
  tmp = new LaurentZ(0);
  int n, i,j;
  ans->coeff[0] = a->coeff[0]->copy();

  for(n=1;n<=newtrunc;n++){
    newcoeff = a->coeff[n]->copy();
    for(j=0;j<=n-1;j++){
      LaurentZ::multiply(tmp, ans->coeff[j], b->coeff[n-j]);
      newcoeff->subtractWith(tmp);
    }
    ans->coeff[n]=newcoeff;
  }
  delete tmp;
  return ans;
}

int  QZSeriesZ::divideByLaurentZ(LaurentZ* b){
  //returns 1 if no remainders, returns 0 otherwise.
  int noRemainder=1;
  int i, res;
  for(i=0;i<=trunc;i++){
    res = coeff[i]->divideWith(b);
    if(res==0)noRemainder=0;
  }
  return noRemainder;
}

QZSeriesZ* QZSeriesZ::multiplyByLaurentZ(QZSeriesZ *a, LaurentZ* b){
  QZSeriesZ *ans = a->copy();
  ans->multiplyByLaurentZ(b);
  return ans;

}
void QZSeriesZ::multiplyByLaurentZ(LaurentZ* b){
  int i;
  for(i=0;i<=trunc;i++){
    coeff[i]->multiplyWith(b);
  }
}

bool QZSeriesZ::divideByIntegerWithCheck(fmpz_t c){ //returns 1 if divide exactly
    bool success=true;
    int i;
    for(i=0;i<=trunc;i++){
        success = success && (coeff[i]->divideByIntegerWithCheck(c));
    }
    return success;
}
void QZSeriesZ::divideByIntegerWithoutCheck(fmpz_t c){ //assumes division will be exact
    int i;
    for(i=0;i<=trunc;i++){
        coeff[i]->divideByIntegerWithoutCheck(c);
    }
}
void QZSeriesZ::negate(){
  int i;
  for(i=0;i<=trunc;i++){
    coeff[i]->negate();
  }
}

void QZSeriesZ::multiplyByScalar(fmpz_t m){
  int i;
  for(i=0;i<=trunc;i++){
    coeff[i]->scalarMultiply(m);
  }
}

QZSeriesZ* QZSeriesZ::E4(int trunc0){
  QZSeriesZ* ans=new QZSeriesZ(trunc0,0,1);
  int i;
  fmpz_t mult, tmp;
  fmpz_init(mult); fmpz_init(tmp);
  fmpz_set_si(mult,240); //-504 for #6
  for(i=1;i<=trunc0;i++){
    sigma(tmp, 3, i);
    fmpz_mul(tmp, tmp, mult);
    ans->coeff[i] = new LaurentZ(tmp);
  }
  fmpz_clear(mult); fmpz_clear(tmp);
  return ans;
}
QZSeriesZ* QZSeriesZ::E6(int trunc0){
  QZSeriesZ* ans=new QZSeriesZ(trunc0,0,1);
  int i;
  fmpz_t mult, tmp;
  fmpz_init(mult); fmpz_init(tmp);
  fmpz_set_si(mult,-504); //-504 for #6
  for(i=1;i<=trunc0;i++){
    sigma(tmp, 5, i);
    fmpz_mul(tmp, tmp, mult);
    ans->coeff[i] = new LaurentZ(tmp);
  }
  fmpz_clear(mult); fmpz_clear(tmp);
  return ans;
}

void QZSeriesZ::sigma(fmpz_t result, int k, int n){
    fmpz_set_si(result, 0);
    fmpz_t tmp;
    fmpz_init(tmp);
    int d;
    for(d=1;d<=n;d++){
        if((n%d)==0){
            fmpz_set_si(tmp, d);
            fmpz_pow_ui(tmp, tmp, k);
            fmpz_add(result, result, tmp);
        }
    }
    fmpz_clear(tmp);
}
void QZSeriesZ::saveToFile(FILE *f){
    fprintf(f, "%d\n", trunc);
    int i;
    for(i=0;i<coeffArraySize;i++){
        coeff[i]->saveToFile(f);
        fprintf(f, "\n");
    }
}
QZSeriesZ::QZSeriesZ(FILE *f){
    fscanf(f, "%d", &trunc);
    coeffArraySize=trunc+1;
    coeff=new LaurentZ*[coeffArraySize];
    int i;
    for(i=0;i<coeffArraySize;i++){
        coeff[i]=new LaurentZ(f);
    }
}

