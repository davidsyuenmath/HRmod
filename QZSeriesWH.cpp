/***
  encodes a Weakly Holomorphic series in q,z.
  20181220: Fixed ThetaBlockWithOps.
***/

#include "QZSeriesWH.h"

fmpz_mat_t QZSeriesWH::diagMat;
int QZSeriesWH::printGG=0;
int QZSeriesWH::etaPowerArraySize=0;
int QZSeriesWH::E4PowerArraySize=0;
int QZSeriesWH::E6PowerArraySize=0;
int QZSeriesWH::thetaArraySize=0;
int* QZSeriesWH::thetaPowerArraySize=NULL;
QZSeriesWH** QZSeriesWH::etaPowers=NULL;
QZSeriesWH** QZSeriesWH::E4Powers=NULL;
QZSeriesWH** QZSeriesWH::E6Powers=NULL;
int* QZSeriesWH::etaTruncations=NULL;
int* QZSeriesWH::E4Truncations=NULL;
int* QZSeriesWH::E6Truncations=NULL;
QZSeriesWH*** QZSeriesWH::thetaPowers=NULL;
int** QZSeriesWH::thetaTruncations=NULL;


QZSeriesWH::QZSeriesWH(){ //does nothing!
}

QZSeriesWH::QZSeriesWH(int trunc0){
  qzs=new QZSeriesZ(trunc0);
  qMinExp=0;
  qExpAdjustOver24=0;
  zExpAdjustOver2=0;
}
QZSeriesWH::QZSeriesWH(int trunc0, int initalCoeff){
  qzs=new QZSeriesZ(trunc0,0,initalCoeff);
  qMinExp=0;
  qExpAdjustOver24=0;
  zExpAdjustOver2=0;
}

QZSeriesWH:: ~QZSeriesWH(){//destructor
  destroy();
}
void QZSeriesWH::destroy(){
    //cout<<"QZSeriesWH destroy activated.\n";
    if(qzs!=NULL){    delete qzs;}
    qzs=NULL;
}


QZSeriesWH::QZSeriesWH(int trunc0, vector<int> c){
  qzs=new QZSeriesZ(trunc0,c);
  qMinExp=0;
  qExpAdjustOver24=0;
  zExpAdjustOver2=0;
}


QZSeriesWH* QZSeriesWH::EtaFunction(int trunc0){
    QZSeriesWH* ans = new QZSeriesWH();
    ans->qzs=QZSeriesZ::EtaFunctionWithoutFractionalExponent(trunc0);
    ans->qMinExp=0;
    ans->qExpAdjustOver24=1;
    ans->zExpAdjustOver2=0;
    return ans;
}

QZSeriesWH* QZSeriesWH::InverseEtaFunction(int trunc0){
    QZSeriesWH* ans = new QZSeriesWH();
    ans->qzs=QZSeriesZ::InverseEtaFunctionWithoutFractionalExponent(trunc0);
    ans->qMinExp=0;
    ans->qExpAdjustOver24=-1;
    ans->zExpAdjustOver2=0;
    return ans;
}
QZSeriesWH* QZSeriesWH::EtaFunctionPower(int trunc0, int pow){
    //This could be made more efficient by using binary expansion of pw
    //But don't care right now.

    if(pow>0){
        if(pow<etaPowerArraySize){
            if(etaPowers[pow]!=NULL){
                if(etaTruncations[pow]>=trunc0){
                    return etaPowers[pow]->copy();
                }else{
                    delete etaPowers[pow]; etaPowers[pow]=NULL;etaTruncations[pow]=0;
                }
            }
        }
    }

    if(pow==0){return new QZSeriesWH(trunc0,1);}
    QZSeriesWH* ans;
    if(pow>0){ans=EtaFunction(trunc0);}else{ans=InverseEtaFunction(trunc0);}
    if(abs(pow)==1)return ans;
    int i;
    QZSeriesWH* tmp = ans->copy();
    for(i=2;i<=abs(pow);i++)ans->multiplyWith(tmp);
    ans->normalize();
    delete tmp;
    if(pow>0){
        if(pow<etaPowerArraySize){
            if(etaPowers[pow]==NULL){
                etaPowers[pow]=ans->copy();
                etaTruncations[pow]=etaPowers[pow]->getMaxExponent();
            }
        }
    }

    return ans;
}
QZSeriesWH* QZSeriesWH::ThetaFunctionD(int trunc0, int d){
    if(d<0)d=-d;
    if(d==0){cout<<"ERROR 349gsdg ThetaFunctionD called with d=0.Abort.\n";exit(1);}
    QZSeriesWH* ans = new QZSeriesWH();
    ans->qzs=new QZSeriesZ(trunc0);
    ans->qMinExp=0;
    ans->qExpAdjustOver24=3;
    ans->zExpAdjustOver2=d;
    int i,n,r;
    LaurentZ* tmpL;
    for(i=0;i<=trunc0+2;i++){
        n=(i*(i+1))/2;
        if(n>trunc0)break;
        r=i*d;
        if(i%2==0){
            tmpL=new LaurentZ(1,r);
        }else{
            tmpL=new LaurentZ(-1,r);
        }
        ans->qzs->coeff[n]->addWith(tmpL);
        delete tmpL;
    }
    for(i=-1;i>=-trunc0-2;i--){
        n=(i*(i+1))/2;
        if(n>trunc0)break;
        r=i*d;
        if(i%2==0){
            tmpL=new LaurentZ(1,r);
        }else{
            tmpL=new LaurentZ(-1,r);
        }
        ans->qzs->coeff[n]->addWith(tmpL);
        delete tmpL;
    }
    return ans;
}

void QZSeriesWH::normalize(){
  int deg = qzs->trunc;
  if(deg<0){ //series is zero
    qMinExp=0; qExpAdjustOver24=0; zExpAdjustOver2=0;
    return;
  }
  if(qExpAdjustOver24<0){
    while(qExpAdjustOver24<0){
      qMinExp-=1;
      qExpAdjustOver24+=24;
    }
  }
  if(qExpAdjustOver24>=24){
    qMinExp+=(qExpAdjustOver24/24); qExpAdjustOver24=qExpAdjustOver24%24;
  }
  int zAdj=0;
  if(zExpAdjustOver2<0){
    while(zExpAdjustOver2<0){
      zExpAdjustOver2+=2;
      zAdj-=1;
    }
  }
  if(zExpAdjustOver2>=2){
    while(zExpAdjustOver2>=2){
      zExpAdjustOver2-=2;
      zAdj+=1;
    }
  }
  if(zAdj!=0){
    vector<int> c; c.push_back(1);
    LaurentZ* tmpL=new LaurentZ(c, zAdj);
    qzs->multiplyByLaurentZ(tmpL);
  }
  normalizePartially();
}
void QZSeriesWH::normalizePartially(){
  int deg = qzs->trunc;
  int i, ord=0;
  for(i=0;i<=deg;i++){
    if(qzs->coeff[i]->isZero()){
      //if(deg==0){//series is zero
        //qMinExp=0; qExpAdjustOver24=0; zExpAdjustOver2=0;
        //Do not change because trunc is important
       // return;
      //}
    }else{
      ord=i;
      qzs->shiftExponents(-ord);
      qMinExp+=ord;
      return;
    }
  }
  //Really no error. This means the series is zero.
  return;

   cout<<"ERROR 325yt75 in QZSeriesWH normalize() ord not found. deg = "<<deg
   <<". Abort.\n";
   cout<<getStringNoNormalize()<<"\n";
   exit(1);
}
void QZSeriesWH::shiftExponents(int shift){
  if(shift<0){cout<<"ERROR 4t306w4309: cannot call with shift<0. Abort.\n";exit(1);}
  if(shift>0){
    qzs->shiftExponents(shift);
    qMinExp-=shift;
  }
}



bool QZSeriesWH::isZero(){
  return qzs->isZero();
}


string QZSeriesWH::getString(){
    return getString("q", "z");
}
string QZSeriesWH::getStringNoNormalize(){
    return getStringNoNormalize("q", "z");
}

string QZSeriesWH::getString(string qvar, string zvar){
    normalize();
    return getStringNoNormalize(qvar,zvar);
}
string QZSeriesWH::getStringNoNormalize(string qvar, string zvar){
  string ans, str;
  stringstream ss;
  slong  n;
  ans="";
  int needParentheses=0;
  if(qExpAdjustOver24!=0){
    ss.clear();
    ss<<qExpAdjustOver24;
    ss>>str;
    ans = ans + "q^(" + str + "/24)*";
    needParentheses=1;
  }
  if(zExpAdjustOver2!=0){
    ss.clear();
    ss<<zExpAdjustOver2;
    ss>>str;
    ans = ans + "z^(" + str + "/2)*";
    needParentheses=1;
  }
  if(qMinExp!=0){
    ss.clear();
    ss<<qMinExp;
    ss>>str;
    ans = ans + "q^(" + str + ")*";
    needParentheses=1;
  }
  if(needParentheses){ans = ans + "(";}
  ans = ans + qzs->getString(qvar, zvar);
  if(needParentheses){ans = ans + ")";}
  return ans;
}

QZSeriesWH* QZSeriesWH::copy(){
  return new QZSeriesWH(this);
}

QZSeriesWH::QZSeriesWH(QZSeriesWH* a){
  qExpAdjustOver24=a->qExpAdjustOver24;
  zExpAdjustOver2=a->zExpAdjustOver2;
  qMinExp=a->qMinExp;
  qzs=a->qzs->copy();
}

QZSeriesWH* QZSeriesWH::multiply(QZSeriesWH* a, QZSeriesWH *b){
  QZSeriesWH *ans = new QZSeriesWH();
  a->normalizePartially();b->normalizePartially();
  ans->qzs = QZSeriesZ::multiply(a->qzs, b->qzs);
  ans->qMinExp =a->qMinExp + b->qMinExp;
  ans->qExpAdjustOver24 =a->qExpAdjustOver24 + b->qExpAdjustOver24;
  ans->zExpAdjustOver2 =a->zExpAdjustOver2 + b->zExpAdjustOver2;
  return ans;
}

void QZSeriesWH::multiplyWith(QZSeriesWH* b){
  normalizePartially();b->normalizePartially();
  qzs->multiplyWith(b->qzs);
  qMinExp =qMinExp + b->qMinExp;
  qExpAdjustOver24 =qExpAdjustOver24 + b->qExpAdjustOver24;
  zExpAdjustOver2 =zExpAdjustOver2 + b->zExpAdjustOver2;
}


QZSeriesWH* QZSeriesWH::add(QZSeriesWH* a, QZSeriesWH* b){
  QZSeriesWH *ans = a->copy();
  ans->addWith(b);
  return ans;
}

void QZSeriesWH::addWith(QZSeriesWH* b){
  //Decide on common denominator, so to speak.
  normalize(); b->normalize();
  if(qExpAdjustOver24!=b->qExpAdjustOver24){
    cout<<"ERROR 325fdsfs: cannot addWith two series with different qExpAdjustOver24: "
      <<qExpAdjustOver24<<","<<b->qExpAdjustOver24<<". Abort.\n";
    exit(1);
  }
  if(zExpAdjustOver2!=b->zExpAdjustOver2){
    cout<<"ERROR hsfdh323: cannot addWith two series with different zExpAdjustOver2: "
      <<zExpAdjustOver2<<","<<b->zExpAdjustOver2<<". Abort.\n";
    exit(1);
  }
  if(qMinExp > b->qMinExp){shiftExponents(qMinExp - b->qMinExp);}
  if(qMinExp < b->qMinExp){b->shiftExponents(b->qMinExp - qMinExp);}
  if(qMinExp != b->qMinExp){
    cout<<"ERROR jgfjghg in QZSeriesWH::addWith. Abort.\n"; exit(1);
  }
  qzs->addWith(b->qzs);
  normalize();
}
void QZSeriesWH::subtractWith(QZSeriesWH* b){
  //Decide on common denominator, so to speak.
  normalize(); b->normalize();
  if(qExpAdjustOver24!=b->qExpAdjustOver24){
    cout<<"ERROR 34yergw: cannot subtractWith two series with different qExpAdjustOver24: "
      <<qExpAdjustOver24<<","<<b->qExpAdjustOver24<<". Abort.\n";
    exit(1);
  }
  if(zExpAdjustOver2!=b->zExpAdjustOver2){
    cout<<"ERROR hsadh523td: cannot subtractWith two series with different zExpAdjustOver2: "
      <<zExpAdjustOver2<<","<<b->zExpAdjustOver2<<". Abort.\n";
    exit(1);
  }
  if(qMinExp > b->qMinExp){shiftExponents(qMinExp - b->qMinExp);}
  if(qMinExp < b->qMinExp){b->shiftExponents(b->qMinExp - qMinExp);}
  if(qMinExp != b->qMinExp){
    cout<<"ERROR jdfjdjaaa in QZSeriesWH::subtractWith. Abort.\n"; exit(1);
  }
  qzs->subtractWith(b->qzs);
  normalize();
}

QZSeriesWH* QZSeriesWH::addScalarMultiple(QZSeriesWH *a, fmpz_t m, QZSeriesWH* b){
  QZSeriesWH *ans = a->copy();
  ans->addScalarMultipleWith(m, b);
  return ans;
}
void QZSeriesWH::addScalarMultipleWith(fmpz_t m, QZSeriesWH* b){
  //Decide on common denominator, so to speak.
  if(fmpz_is_zero(m))return;

  normalize(); b->normalize();
  if(qExpAdjustOver24!=b->qExpAdjustOver24){
    cout<<"ERROR jlululhl: cannot addScalarMultipleWith two series with different qExpAdjustOver24: "
      <<qExpAdjustOver24<<","<<b->qExpAdjustOver24<<". Abort.\n";
    exit(1);
  }
  if(zExpAdjustOver2!=b->zExpAdjustOver2){
    cout<<"ERROR 67kygk7p;: cannot addScalarMultipleWith two series with different zExpAdjustOver2: "
      <<zExpAdjustOver2<<","<<b->zExpAdjustOver2<<". Abort.\n";
    exit(1);
  }
  if(qMinExp > b->qMinExp){shiftExponents(qMinExp - b->qMinExp);}
  if(qMinExp < b->qMinExp){b->shiftExponents(b->qMinExp - qMinExp);}
  if(qMinExp != b->qMinExp){
    cout<<"ERROR nnygntu in QZSeriesWH::addScalarMultipleWith. Abort.\n"; exit(1);
  }
  qzs->addScalarMultipleWith(m,b->qzs);
  normalize();
  b->normalize();
}

QZSeriesWH* QZSeriesWH::divideByMonic(QZSeriesWH* a, QZSeriesWH *b){
  QZSeriesWH *ans = new QZSeriesWH();
  a->normalizePartially();b->normalizePartially();
  ans->qzs = QZSeriesZ::divideByMonic(a->qzs, b->qzs);
  ans->qMinExp =a->qMinExp - b->qMinExp;
  ans->qExpAdjustOver24 =a->qExpAdjustOver24 - b->qExpAdjustOver24;
  ans->zExpAdjustOver2 =a->zExpAdjustOver2 - b->zExpAdjustOver2;
  return ans;
}

void QZSeriesWH::divideByMonic(QZSeriesWH* b){
  normalizePartially();b->normalizePartially();
  QZSeriesZ* quot = QZSeriesZ::divideByMonic(qzs, b->qzs);
  delete qzs;
  qzs=quot;
  qMinExp =qMinExp - b->qMinExp;
  qExpAdjustOver24 =qExpAdjustOver24 - b->qExpAdjustOver24;
  zExpAdjustOver2 =zExpAdjustOver2 - b->zExpAdjustOver2;
}


int  QZSeriesWH::divideByLaurentZ(LaurentZ* b){
  //returns 1 if no remainders, returns 0 otherwise.
  return qzs->divideByLaurentZ(b);
}

QZSeriesWH* QZSeriesWH::multiplyByLaurentZ(QZSeriesWH *a, LaurentZ* b){
  QZSeriesWH *ans = a->copy();
  ans->multiplyByLaurentZ(b);
  return ans;

}


void QZSeriesWH::multiplyByLaurentZ(LaurentZ* b){
  qzs->multiplyByLaurentZ(b);
}
void QZSeriesWH::negate(){
  qzs->negate();
}
int QZSeriesWH::getMaxExponent(){
    normalize();
    return qMinExp+qzs->trunc;
}
int QZSeriesWH::getMinExponent(){
    normalize();
    return qMinExp;
}
int QZSeriesWH::divideByIntegerWithCheck(fmpz_t c){ //returns 1 if divide exactly
    return qzs->divideByIntegerWithCheck(c);
}
void QZSeriesWH::divideByIntegerWithoutCheck(fmpz_t c){ //assumes division will be exact
    qzs->divideByIntegerWithCheck(c);
}
int QZSeriesWH::divideBy(QZSeriesWH* b){ //returns 1 if successful
    QZSeriesWH* tmp = b->copy();
    int success=1;
    tmp->normalize();
    if(tmp->isZero()){cout<<"ERROR fdkgjsdflg division by zero in divideBy. Abort.\n";exit(1);}
    LaurentZ* lead = tmp->qzs->coeff[0]->copy();
    success = tmp->divideByLaurentZ(lead);
    if(!success){cout<<"WARNING: denominator not evenly divisible by leading term.\n";}
    int success2 = (divideByLaurentZ(lead));  //Was missing = sign until Dec 29, 2016.
    if(!success2){cout<<"WARNING: numerator not evenly divisible by leading term.\n";
        cout<<"lead="<<lead->getString()<<";\n";
        //cout<<"newnumer="<<getString()<<";\n";
    }
    success=success&&success2;
    delete lead;
    divideByMonic(tmp);
    delete tmp;
    return success;
}
QZSeriesWH* QZSeriesWH::linearCombo(QZSeriesWH** qq, fmpz_t* c, int len){
    if(len<=0){cout<<"ERROR fh8hf84 linearCombo called with zero len. Abort.\n";exit(1);}
    int trunc, i;
    trunc=qq[0]->getMaxExponent();
    for(i=1;i<len;i++){
        trunc = min(trunc, qq[i]->getMaxExponent());
    }
    QZSeriesWH* ans = new QZSeriesWH(trunc);
    for(i=0;i<len;i++){
        ans->addScalarMultipleWith(c[i],qq[i]);
    }
    return ans;
}

QZSeriesWH* QZSeriesWH::rationalLinearCombo(QZSeriesWH** qq, fmpq_t* c, int len){
    cout<<"rationalLinearCombo To be implemented.\n";exit(1);
}
int QZSeriesWH::signum(int x){
    if(x==0)return 0;
    if(x>0)return 1;
    return -1;
}
QZSeriesWH* QZSeriesWH::ThetaBlock(int trunc0, int weight, vector<int> d){
    return ThetaBlockEfficient(trunc0, weight, d);
}
QZSeriesWH* QZSeriesWH::ThetaBlockTweak(int trunc0, int weight, vector<int> d){
    return ThetaBlockEfficient(trunc0, weight, d, 1);
}
QZSeriesWH* QZSeriesWH::ThetaBlockInefficient(int trunc0, int weight, vector<int> d){
    int ell, t,i;
    ell=0;
    for(i=0;i<d.size();i++){ell+=signum(d[i]);}
    t=2*weight-ell;
    if((3*ell+t)%24!=0){
        cout<<"Theta block does not have 24|3L+t (weight, L,t)=("<<weight<<","<<ell<<","<<t<<"). Abort.\n";exit(1);
    }
    QZSeriesWH* ans = QZSeriesWH::EtaFunctionPower(trunc0, t);
    QZSeriesWH* th;
    for(i=0;i<d.size();i++)if(d[i]>0){
        th=ThetaFunctionD(trunc0, d[i]);
        ans->multiplyWith(th);
        delete th;
    }
    int success;
    for(i=0;i<d.size();i++)if(d[i]<0){
        th=ThetaFunctionD(trunc0, -d[i]);
        success=ans->divideBy(th);
        if(!success){
            cout<<"ERROR 94rjdfgsd in Theta Block with denominator";
            for(int j=0;j<d.size();j++)cout<<","<<d[j];
            cout<<". Abort.\n"; exit(1);
        }
        delete th;
    }
    ans->normalize();
    if(ans->qExpAdjustOver24!=0){
        cout<<"ERROR fdkhd0hdgk qExpAdjustOver24 is not zero in making a ThetaBlock. Abort.\n"; exit(1);
    }
    if(ans->zExpAdjustOver2!=0){
        cout<<"ERROR dgasdgadgadg zExpAdjustOver2 is not zero in making a ThetaBlock. Abort.\n"; exit(1);
    }
    //ans->truncate(trunc0);
    return ans;
}
QZSeriesWH* QZSeriesWH::ThetaBlockWithOps(int trunc0, int weight, vector<int> b){
    //Using encoding of b:  numOps type num ... type num d1 d2 ... dL
    vector<int>d, ops;
    int numOps=b[0];
    int i;
    for(i=0;i<numOps;i++){
        ops.push_back(b[2*i+1]);
        ops.push_back(b[2*i+2]);
    }
    for(i=2*numOps+1;i<b.size();i++){
        d.push_back(b[i]);
    }
    if(numOps==0){return ThetaBlock(trunc0, weight, d);}
    return ThetaBlockWithOps(trunc0, weight, d, ops);
}
QZSeriesWH* QZSeriesWH::ThetaBlockWithOps(int trunc0, int weight, vector<int> d, vector<int> heckeOps){
    int trunc; int i, type, L, index=0;
    for(i=0;i<d.size();i++){
        index = index + d[i]*d[i]*signum(d[i]);
    }
    if(index%2!=0){
        cout<<"ERROR ds8g8fd in ThetaBlockWithOps: sum of squares not even. Abort.\n";exit(1);
    }
    index=index/2;
    int m=index;
    for(i=0;i<heckeOps.size();i+=2){
        type=heckeOps[i]; L=heckeOps[i+1];
        if(type==1){//hecke up
            m = m*L;
        }else if(type==2){//hecke down
            m = m/L;
        }else if(type==3){//double up
            //trunc = trunc;  //do nothing
            m = m*(L*L);
        }else{
            cout<<"ERROR gdjtyie in ThetaBlockWithOps, illegal value of type: "<<type<<". Abort.\n";exit(1);
        }
    }
    int resultIndex = m; //This is the rssulting index
    trunc=trunc0;
    for(i=heckeOps.size()-2;i>=0;i-=2){ //Run through Hecke ops backwards
        type=heckeOps[i]; L=heckeOps[i+1];
        if(type==1){//hecke up
            trunc=trunc*L;
            m=m/L;
        }else if(type==2){//hecke down
            trunc=myHeckeDownMaxq(m,trunc,L);
            m = m*L;
        }else if(type==3){//double up
            //trunc = trunc;  //do nothing
            m = m/(L*L);
        }else{
            cout<<"ERROR gdjtyie in ThetaBlockWithOps, illegal value of type: "<<type<<". Abort.\n";exit(1);
        }
    }
    if(m!=index){
        cout<<"ERRO dfkasdf; in ThetaBlockWithOps. A double check of m==index ("<<m<<","<<index<<") failed. Abort.\n";exit(1);
    }
    QZSeriesWH* ans = ThetaBlock(trunc, weight, d); //need current index!!
    for(i=0;i<heckeOps.size();i+=2){
        type=heckeOps[i]; L=heckeOps[i+1];
        if(type==1){//hecke up
            ans->applyUp(index,weight,L);
            index = index*L;
        }else if(type==2){//hecke down
            ans->applyDown(index,weight,L);
            index=index/L;
        }else if(type==3){//double up
            ans->applyDoubleUp(L);
            index = index*L*L;
        }else{
            cout<<"ERROR ufgiifd in ThetaBlockWithOps, illegal value of type: "<<type<<". Abort.\n";exit(1);
        }
    }
    if(ans->getMaxExponent()<trunc0){
        cout<<"ERROR t439wts in ThetaBlockWithOps. Did not get desired max q exponent ("
        <<ans->getMaxExponent()<<","<<trunc0<<"). Abort\n"; exit(1);
    }
    //ans->truncate(trunc0);
    return ans;
}
void QZSeriesWH::getFC_noreduce(fmpz_t fc, slong n, slong r){
    //cout<<"What???\n\n";
    //cout<<"n="<<n<<"\n";
    //cout<<"r="<<r<<"\n";
    //cout.flush();//sleep(1);
    //cout<<"qMinExp="<<qMinExp<<"\n";
    //cout<<"What??????\n\n";
    //cout<<"getFC_noreduce: (n,r, qMinExp) = ("<<n<<","<<r<<","<<qMinExp<<")\n";
    slong qexp = n - qMinExp;
    if(qexp<0){
        fmpz_set_si(fc, 0); return;
    }
                if(printGG==-24){cout<<"GGd1 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

    if(qexp>qzs->trunc){
        cout<<"getFC_noreduce: (n,r, qMinExp) = ("<<n<<","<<r<<","<<qMinExp<<")\n";
        cout<<"ERROR 515445u7 in getFC_noreduce: insufficient length.\n"
        <<"(n, qMinExp, trunc) = ("<<n<<","<<qMinExp<<","<<qzs->trunc<<"). Abort.\n";
        exit(1);
    }
                if(printGG==-24){
                        cout<<"(qexp, r, fc)<<"<<qexp<<","<<r<<","<<fmpz_get_str(NULL,10,fc)<<"\n";
                        cout<<"GGd2 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}
                if(printGG==-24){
                        cout<<"GGe1 (coeff[qexp])<<"<<qexp<<","<<qzs->coeff[qexp]->getString()<<"\n";
               }

    qzs->coeff[qexp]->getCoeff(fc, r);

               if(printGG==-24){
                        cout<<"(qexp, r, fc)<<"<<qexp<<","<<r<<","<<fmpz_get_str(NULL,10,fc)<<"\n";
                        cout<<"GGd3 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";
               }
               //if(printGG==-24){
               //         fmpz_set_ui(fc, 12000+100*qexp+r);
               //         cout<<"(fc)<<"<<fmpz_get_str(NULL,10,fc)<<"\n";
               //         cout<<"GGd4 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

}
void QZSeriesWH::getFC(fmpz_t fc, slong n, slong r, slong index){
  //cout<<"HERE Y1 (test code).\n";
  //getFC_noreduce(fc,0,0);
  //fmpz_set_si(fc, 3);
  //cout<<"HERE Y1 end.\n";

  //cout<<"getFC: (n,r,index) = ("<<n<<","<<r<<","<<index<<")\n";
  slong longp2=2L*index;
  slong rnew=((r%longp2)+longp2)%longp2;
  if(rnew>index)rnew-=longp2;  //rnew is between -index and index
  slong nnew=n-((r-rnew)/(longp2))*((r+rnew)/(2L)); //improvement to prevent overflow
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
            //if(printGG==-24){cout<<"GGc diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}
  //cout<<"HERE 1 getFC: (nnew,rnew,index) = ("<<nnew<<","<<rnew<<","<<index<<")\n";
  //cout<<"HERE 1 getFC: (nnew,rnew,index,qminExp) = ("<<nnew<<","<<rnew<<","<<index<<","<<qMinExp<<")\n";

  //fmpz_set_si(fc,7);
  //return;

  //getFC_noreduce(fc,0,0);
  //cout<<"HERE 2 getFC: (nnew,rnew,index) = ("<<nnew<<","<<rnew<<","<<index<<")\n";
  getFC_noreduce(fc, nnew, rnew);

            //if(printGG==-24){cout<<"GGc diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

}
/***
void QZSeriesWH::getFC_Dr(fmpz_t fc, slong D, slong r, slong index){
  //D = 4 * index * n - r^2
  //so  n = (D+r*r)/(4*index)
  slong n = (D + r*r)/(4*index);
  slong longp2 = 2*index;
  r=((r%longp2)+longp2)%longp2;
  if(r>index)r-=longp2;  //r is between -index and index
  if(D!= 4*index*n-r*r){
    cout<<"ERROR dfegs in getFC_Dr of QZSeriesWH\n";
    cout<<"index="<<index<<"\n";
    cout<<"D="<<D<<"\n";
    cout<<"r="<<r<<"\n";
    cout<<"n="<<n<<"\n";
    cout<<(D)<<","<<(4*index*n-r*r)<<"\n";
    exit(1);
  }
}
***/

void QZSeriesWH::getFC_Dr(fmpz_t fc, slong D, slong r, slong index){
  //D = 4 * index * n - r^2
  //so  n = (D+r*r)/(4*index)
  //cout<<"HEREAA1: (D, r)=("<<D<<","<<r<<")...";
  slong longp2 = 2*index;
  r=((r%longp2)+longp2)%longp2;
  if(r>index)r-=longp2;  //r is between -index and index
  slong n = (D + r*r)/(4*index);
  if(D!= 4*index*n-r*r){
    cout<<"ERROR dfegs in getFC_Dr of QZSeriesWH\n";
    cout<<"index="<<index<<"\n";
    cout<<"D="<<D<<"\n";
    cout<<"r="<<r<<"\n";
    cout<<"n="<<n<<"\n";
    cout<<(D)<<","<<(4*index*n-r*r)<<"\n";
    exit(1);
  }
  getFC(fc,n,r,index);
  //cout<<"(n, r,index)=("<<n<<","<<r<<","<<index<<")=> "<<fmpz_get_str(NULL,10,fc)<<"\n";
}

void QZSeriesWH::applyDoubleUp(int L){
    //Replace z by L*z
    int n;
    QZSeriesZ* tmp = qzs;
    qzs=new QZSeriesZ(tmp->trunc);
    for(n=0;n<=qzs->trunc;n++){
        delete qzs->coeff[n];
        qzs->coeff[n]=LaurentZ::expand(tmp->coeff[n],L);
    }
}
LaurentZ* QZSeriesWH::getCoeffCopy(int n){
    int i = n-qMinExp;
    if(i<0)return new LaurentZ(0);  //These are all copies!!
    if(i>qzs->trunc){
        cout<<"ERROR h9fdh9ds Requested getCoeff is beyond truncation.\n"
        <<"that is, request n ("<<n<<") is > maxExponent ("<<getMaxExponent()<<"). Abort.\n";
        exit(1);
    }
    return qzs->coeff[i]->copy();
}
int QZSeriesWH::getCoeffMaxr(int n){
    int i = n-qMinExp;
    if(i<0)return 0; //These are all copies!!
    if(i>qzs->trunc){
        cout<<"ERROR 7ifghdsfjd Requested getCoeffMaxr is beyond truncation.\n"
        <<"that is, request n ("<<n<<") is > maxExponent ("<<getMaxExponent()<<"). Abort.\n";
        exit(1);
    }
    return qzs->coeff[i]->getMaxExponent();
}

void QZSeriesWH::applyUp(int index, int weight, int L){
    fmpz_t factor;
    LaurentZ* tmpL;
    if(weight<0){
        cout<<"ERROR applyUp only avail for weight>=1 ("<<weight<<"). Abort.\n";exit(1);
    }
    /**
    IF WEIGHT ZERO, THEN RETURNS L*answer
    It is up to the user to remember to divide the answer by L!!!!
    **/
    fmpz_init(factor); //fmpz_pow_ui(factor,factor,weight-2);
    int d,n,newMinExp, newMaxExp,i;
    newMaxExp=getMaxExponent()/L;
    newMinExp=min(0,L*qMinExp);
    QZSeriesZ* ans = new QZSeriesZ(newMaxExp - newMinExp);
    for(n=newMinExp;n<=newMaxExp;n++){
        i=n-newMinExp;
        for(d=1;d<=L;d++)if(L%d==0){
            if(d==1){
                delete ans->coeff[i];
                if(weight>0){
                    ans->coeff[i]=getCoeffCopy(L*n);
                }else{
                    ans->coeff[i]=getCoeffCopy(L*n);
                    fmpz_set_si(factor,L/1);
                    ans->coeff[i]->scalarMultiply(factor);
                }
                //cout<<"HERE SS 1:"<<ans->coeff[i]->checkSymmetry()<<"\n";
            //}else if((L*n)%(d*d)==0){ THIS WAS WRONG, CAUSE OF A BUG
            }else if((n%d)==0){
                if(weight>0){
                    fmpz_set_si(factor,d);
                    fmpz_pow_ui(factor, factor, weight-1);
                }else{
                    fmpz_set_si(factor,L/d); //For weight zero. Remember to divide answer by L
                }
                tmpL = getCoeffCopy((L*n)/(d*d));
                tmpL->expand(d);//Expand, not shrink
                //cout<<"HERE SS 2:"<<tmpL->checkSymmetry()<<"\n";
                tmpL->scalarMultiply(factor);
                //cout<<"HERE SS 2:"<<tmpL->checkSymmetry()<<"\n";
                //cout<<"coeff[i]="<<ans->coeff[i]->getString()<<"\n";
                ans->coeff[i]->addWith(tmpL);
                //cout<<"HERE SS 4:"<<ans->coeff[i]->checkSymmetry()<<"\n";
                //cout<<"tmpL="<<tmpL->getString()<"\n";
                //cout<<"coeff[i]="<<ans->coeff[i]->getString()<<"\n";

                delete tmpL;
            }
        }
    }
    fmpz_clear(factor);
    delete qzs;
    qzs=ans;
    qMinExp=newMinExp;
    normalize();
}
long long QZSeriesWH:: mysqrtfloor(long long x){
  if(x<0){
    cout<<"ERROR in mysqrtfloor 3r232r2 of Nums. x="<<x<<"\n";exit(1);
  }
  long long y = (long long)sqrt((double)x);
  if(((y*y)>x)||((y+1)*(y+1)<=x)){
    cout<<"ERROR in mysqrtfloor 248tdf of Nums. x="<<x<<",y="<<y<<"\n";exit(1);
  }
  return y;
}
int QZSeriesWH::myHeckeDownMaxq(int targetm, int n,int h){
  return max(h*n, n/h + (h*targetm)/4 + 1);
}
int QZSeriesWH::myHeckeDownMaxqInverse(int targetm, int n,int h){
  return max(0,min(n/h, h*(n - ((h*targetm)/4 + 1))));
}
void QZSeriesWH::applyDown(int index, int weight, int L){ //prime only, assumes holomorphic
    //cout<<"HERE T00 \n";
    fmpz_t tmpL, factor;
    //cout<<"HERE T0\n";
    fmpz_init(tmpL); fmpz_set_si(tmpL, L);
    if(weight<2){
        cout<<"ERROR applyDown only avail for weight>=2 ("<<weight<<"). Abort.\n";exit(1);
    }
    fmpz_init(factor); fmpz_pow_ui(factor,tmpL,weight-2);
    //cout<<"HERE TA\n";
    if(!fmpz_is_prime(tmpL)){
        fmpz_clear(tmpL);
        cout<<"applyDown Hecke is only implemented for prime L. Abort.\n";exit(1);
    }
    //cout<<"HERE TB\n";
    normalize();
    //if((qExpAdjustOver24!=0)||(zExpAdjustOver2!=0)||(qMinExp<0)
    slong n, r, rmax,a,t ,m ;
    m=index/L; //m is new index in formula from JacobiForm.cpp
    //cout<<"HERE TC\n";
    int newmaxq=myHeckeDownMaxqInverse(m, getMaxExponent(),L);
    if(newmaxq<=0){
        cout<<"ERROR f8gsdgd in applyDown: newmaxq="<<newmaxq<<". Abort.\n"; exit(1);
    }
    //cout<<"HERE TD newmaxq="<<newmaxq<<"\n";
    QZSeriesZ* ans = new QZSeriesZ(newmaxq);
    //cout<<"HERE TE\n";
    //ans->coeff[1];   ////What is this for????????????? DXXX
    //cout<<"HERE TE2\n";
    for(n=0;n<=ans->trunc;n++){
        rmax=mysqrtfloor(4*m*n);
        for(r=-rmax;r<=rmax;r++){
            //cout<<"HERE TF0: (n,r)=("<<n<<","<<r<<").";
            //cout<<"HERE TF: getFC(tmpL,L*n,L*r,index):("<<L*n<<","<<L*r<<","<<index<<")\n";
            getFC(tmpL,L*n,L*r,index);
            //if(n==1&&r==10){cout<<">>> tmpL="<<fmpz_get_str(NULL,10,tmpL)<<".\n";}
            //cout<<"TF gotten.\n";
            ans->coeff[n]->addToTerm(tmpL, r);
            //cout<<"TF coeff set.\n";
            for(a=0;a<L;a++){
                t=n+a*r+a*a*m;
                if(t%L==0){
                    //cout<<"HERE TG: getFC(tmpL,t/L,r+2*a*m,index):("<<t/L<<","<<r+2*a*m<<","<<index<<")\n";
                    getFC(tmpL,t/L,r+2*a*m,index);
                    //cout<<"TG gotten.\n";
                    //if(n==1&&r==10){cout<<">>> tmpL="<<fmpz_get_str(NULL,10,tmpL)<<".\n";}
                    //if(n==1&&r==10){cout<<">>> "<<getCoeff(14)->getString()<<"<<<\n";}
                    fmpz_mul(tmpL,tmpL,factor);
                    //if(n==1&&r==10){cout<<">>> tmpL="<<fmpz_get_str(NULL,10,tmpL)<<".\n";}
                   ans->coeff[n]->addToTerm(tmpL, r);
                    //cout<<"TG coeff set.\n";
                }
            }
            //if(n==1&&r==10){exit(1);}
        }
    }


    fmpz_clear(tmpL);
    fmpz_clear(factor);
    delete qzs;
    qzs=ans;
    qMinExp=0;qExpAdjustOver24=0;zExpAdjustOver2=0;
    normalize();

}

void QZSeriesWH::saveToFile(ofstream &f, string varname){
    f<<varname<<"=";
    f<< getString();
    f<<";\n";
    f.flush();
}

int QZSeriesWH::checkSymmetry(){
    //Only symmetry, not antisymmetry
    int i, ans;
    ans=1;
    for(i=0;i<=qzs->trunc;i++){
        if(!(qzs->coeff[i]->checkSymmetry())){
            ans=0;break;
        }
    }
    return ans;
}

void QZSeriesWH::putIntoDRList(DRList* dr, int index, int uptoN){
    int n, r, d;
    if(uptoN>getMaxExponent()){
        cout<<"ERROR hhsfh99a: putIntoDRList has uptoN>maxExp ("<<uptoN<<","<<getMaxExponent()<<").\n";
        exit(1);
    }
    fmpz_t c;  fmpz_init(c); //Can't believe I forgot to initialize this c 2/5/2017
    for(n=getMinExponent();n<=getMaxExponent();n++){
        for(r=0;r<=getCoeffMaxr(n);r++){ //Can be more efficient?
            d=4*index*n-r*r;
            if(d<0){
                getFC(c, n, r, index);
                if(!fmpz_is_zero(c))
                {
                    dr->insert(d, r%(2*index));
                }
            }
        }
    }
    fmpz_clear(c);
}

int QZSeriesWH::hasNonzeroNewDR(DRList* dr, int index, int uptoN){
    int n, r, d;
    if(uptoN>getMaxExponent()){
        cout<<"ERROR hhsfh99a: putIntoDRList has uptoN>maxExp ("<<uptoN<<","<<getMaxExponent()<<").\n";
        exit(1);
    }
    fmpz_t c;  fmpz_init(c); //Can't believe I forgot to initialize this c 2/5/2017
    for(n=getMinExponent();n<=getMaxExponent();n++){
        for(r=0;r<=getCoeffMaxr(n);r++){ //Can be more efficient?
            d=4*index*n-r*r;
            if(d<0){
                getFC(c, n, r, index);
                if(!fmpz_is_zero(c))
                {
                    if(!dr->alreadyOnList(d,r)){

                            fmpz_clear(c);
                            return 1;
                                           //dr->insert(d, r%(2*index));
                    }
                }
            }
        }
    }
    fmpz_clear(c);
    return 0;
}

void QZSeriesWH::getHumbertMultiplicity(fmpz_t h, int minD, int d, int r, int index){
    fmpz_set_si(h,0);
    fmpz_t c; fmpz_init(c);
    int n;
    //cout<<"getHumbertMultiplicity:\n";
    //cout<<"HEREAA2: (minD,d,r,index)=("<<minD<<","<<d<<","<<r<<","<<index<<").\n";
    if(printGG==-24){cout<<"GGa diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

    for(n=1;n<=mysqrtfloor((-minD)/(-d));n++){
            if(printGG==-24){cout<<"GGb diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

        getFC_Dr(c, n*n*d, n*r, index);

            if(printGG==-24){cout<<"GGb diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

        fmpz_add(h,h,c);

            if(printGG==-24){cout<<"GGb diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,0))<<"\n";}

        //cout<<"n,h,c="<<n<<","<<fmpz_get_str(NULL,10,h)<<","<<fmpz_get_str(NULL,10,c)<<".\n";
    }
    //cout<<"done TT\n";
    fmpz_clear(c);
}
void QZSeriesWH::makeHumbertVector(fmpz_mat_t M, DRList* dr, int index, int colIndex){
    int i, d, r;
    //fmpz_t c; fmpz_t cc; fmpz_init(cc);
    //fmpz_init(c);
    //cout<<"HEREAA3: makeHumbertVector: colIndex = "<<colIndex<<"\n";
    //printGG=-24;
    //for(i=0;i<dr->getLength();i++){

    for(i=0;i<dr->getLength();i++){  //DEBUG MODE!!! starting at 1 instead of 0.

        //fmpz_init(c);
        d=dr->DList[i]; r=dr->RList[i];
        //cout<<"HHz(i,colIndex,entry): "<<i<<","<<colIndex<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,i,colIndex))<<"\n";
        //cout<<"HHz(0,0,entry): "<<0<<","<<0<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,0,colIndex))<<"\n";
        //if(d==-24){
                getHumbertMultiplicity(fmpz_mat_entry(M,i,colIndex), dr->getMinD(), d, r, index);

                //if(i==0)if(colIndex==0){fmpz_set(Diagnostics::HumVecEntry00, fmpz_mat_entry(M,i,colIndex));
                //cout<<"HumVecEntry00 set to "<<fmpz_get_str(NULL,10,Diagnostics::HumVecEntry00)<<"\n";}

                //getHumbertMultiplicity(cc, dr->getMinD(), d, r, index);
                //fmpz_set_si(c, 100*d+r);
        //}else{getHumbertMultiplicity(cc, dr->getMinD(), d, r, index);}
        //cout<<"HH:d,r,c="<<d<<","<<r<<","<<fmpz_get_str(NULL,10,cc)<<".\n";
        //cout<<"HHa(i,colIndex,entry): "<<i<<","<<colIndex<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,i,colIndex))<<"\n";
        //cout<<"HHa(0,0,entry): "<<0<<","<<0<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,0,colIndex))<<"\n";
        //if(d==-24){
        if(printGG==-24){cout<<"FFa diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,colIndex))<<"\n";}
        //fmpz_set(fmpz_mat_entry(M,i,colIndex),cc);
        //}else{fmpz_set(fmpz_mat_entry(M,i,colIndex),cc);}
        //cout<<"HHb(i,colIndex,entry): "<<i<<","<<colIndex<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,i,colIndex))<<"\n";
        if(printGG==-24){cout<<"FFb diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,colIndex))<<"\n";}
        //cout<<"HHb(0,0,entry): "<<0<<","<<0<<","<<fmpz_get_str(NULL,10,fmpz_mat_entry(M,0,colIndex))<<"\n";
         if(printGG==-24){cout<<"FFc diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,colIndex))<<"\n";}
       //     fmpz_clear(c);
         if(printGG==-24){cout<<"FFd diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,colIndex))<<"\n";}
       //if(i==0){fmpz_set(fmpz_mat_entry(diagMat,i,colIndex), fmpz_mat_entry(M,i,colIndex));}
         if(printGG==-24){cout<<"FFe diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(diagMat,0,colIndex))<<"\n";}
   }
    //fmpz_clear(c);
    //fmpz_clear(cc);
}
void QZSeriesWH::makeHumbertVector(fmpz_mat_t M, DRList* dr, int index){
    makeHumbertVector(M, dr, index, 0);
}
void QZSeriesWH::truncate(int newtrunc){
    int currDegree=getMaxExponent();
    if(newtrunc>currDegree){cout<<"ERROR in truncate request 348dsgf in QZSeriesWH. Abort.\n";exit(1);}
    qzs->trunc = newtrunc - qMinExp;
    if(getMaxExponent()!=newtrunc){
        cout<<"Logical error in truncate dg8afdhgdafg. Abort.\n";exit(1);
    }
}
void QZSeriesWH::makeSingularVector(fmpz_mat_t M, DRList* dr, int index, int colIndex){
    int i, d, r;
    fmpz_t c;  fmpz_init(c); //Can't believe I forgot to initialize this c 2/5/2017
    //cout<<"dr->getMinD()="<<dr->getMinD()<<"\n";
    for(i=0;i<dr->getLength();i++){
        d=dr->DList[i]; r=dr->RList[i];
        getFC_Dr(c, d, r, index);
        fmpz_set(fmpz_mat_entry(M,i,colIndex),c);
        //cout<<"HERE HH: "<<fmpz_get_str(NULL,10,c)<<"\n";
    }
    fmpz_clear(c);
}
void QZSeriesWH::makeSingularVector(fmpz_mat_t M, DRList* dr, int index){
    makeSingularVector(M, dr, index, 0);
}
double QZSeriesWH::getOrd(int index){
    int n, r, d;
    double ord = 0;
    double tempord;
    int uptoN = index/4;
    if(uptoN>getMaxExponent()){
        cout<<"ERROR gjkrdfig34jadgasdg: getOrd requires uptoN=index/4 has uptoN>maxExp ("<<uptoN<<","<<getMaxExponent()<<").\n";
        exit(1);
    }
    fmpz_t c;  fmpz_init(c);
    for(n=getMinExponent();n<=getMaxExponent();n++){
        for(r=0;r<=getCoeffMaxr(n);r++){ //Can be more efficient?
            d=4*index*n-r*r;
            if(d<0){
                getFC(c, n, r, index);
                if(!fmpz_is_zero(c))
                {
                   tempord = d / (4.0*index);
                   if(tempord<ord){ord=tempord;}
                }
            }
        }
    }
    fmpz_clear(c);
    return ord;
}
int QZSeriesWH::getSymmetry(int index){
    int i, ibegin, fs;
    ibegin=getMinExponent();
    if(ibegin>=0)return 1;
    ibegin=-ibegin;
    ibegin=((int)sqrt(ibegin));
    if((ibegin+1)*(ibegin+1)<=-getMinExponent()){ibegin=ibegin+1;}//EXTRA CAUTIOUS because of the (int)sqrt()
    fmpz_t tot, c;
    fmpz_init(tot);fmpz_init(c);
    fmpz_set_ui(tot,0);
    for(i=ibegin;i>=1;i--){
        getFC(c, -(i*i),0,index);
        fmpz_add(tot, tot, c);
    }
    if(fmpz_is_even(tot)){fs=1;}else{fs=-1;}
    fmpz_clear(tot);
    fmpz_clear(c);
    return fs;
}

void QZSeriesWH::multiplyByScalar(fmpz_t m){
  qzs->multiplyByScalar(m);

}

void QZSeriesWH::initPowers(){

    //cout<<"QZSeriesWH::initPowers() started.\n";
    etaPowerArraySize=1000;
    E4PowerArraySize=100;
    E6PowerArraySize=100;
    etaPowers = new QZSeriesWH*[etaPowerArraySize];
    E4Powers = new QZSeriesWH*[E4PowerArraySize];
    E6Powers = new QZSeriesWH*[E6PowerArraySize];
    etaTruncations = new int[etaPowerArraySize];
    E4Truncations = new int[E4PowerArraySize];
    E6Truncations = new int[E6PowerArraySize];
    int i,j;
    for(i=0;i<etaPowerArraySize;i++){etaPowers[i]=NULL;etaTruncations[i]=0;}
    for(i=0;i<E4PowerArraySize;i++){E4Powers[i]=NULL;E4Truncations[i]=0;}
    for(i=0;i<E6PowerArraySize;i++){E6Powers[i]=NULL;E6Truncations[i]=0;}
    thetaArraySize=50;
    thetaPowers = new QZSeriesWH**[thetaArraySize];
    thetaTruncations= new int*[thetaArraySize];
    thetaPowerArraySize = new int[thetaArraySize];
    for(i=0;i<thetaArraySize;i++){
        thetaPowerArraySize[i]=100;
        thetaPowers[i]=new QZSeriesWH*[thetaPowerArraySize[i]];
        thetaTruncations[i]=new int[thetaPowerArraySize[i]];
        for(j=0;j<thetaPowerArraySize[i];j++){
            thetaPowers[i][j]=NULL;
            thetaTruncations[i][j]=0;
        }
    }
    //cout<<"QZSeriesWH::initPowers() done.\n";

}
QZSeriesWH* QZSeriesWH::ThetaFunctionDPower(int trunc0, int d, int pow){

    if(pow<0){
        cout<<"ThetaFunctionDPower only implemented for pow>=0. Abort.\n";exit(1);
    }
    if(pow==0){return new QZSeriesWH(trunc0,1);}

    if(d<thetaArraySize){
        if(pow<thetaPowerArraySize[d]){
            if(thetaPowers[d][pow]!=NULL){
                if(thetaTruncations[d][pow]>=trunc0){
                    //cout<<"Yay, efficient reuse of thetaPowers["<<d<<"]["<<pow<<"].\n";
                    return thetaPowers[d][pow]->copy();
                }else{
                    //cout<<"Aw, need to increase thetaTruncation ("<<thetaTruncations[d][pow]
                    //<<","<<trunc0<<") at thetaPowers["<<d<<"]["<<pow<<"].\n";
                    thetaTruncations[d][pow]=0;
                    delete thetaPowers[d][pow]; thetaPowers[d][pow]=NULL;
                }
            }
        }
    }

    //THIS CAN BE WAY MORE EFFICIENT!
    QZSeriesWH* ans = ThetaFunctionD(trunc0,d);
    if(pow>1){
        int i;
        QZSeriesWH* tmp = ans->copy();
        for(i=2;i<=pow;i++){
                ans->multiplyWith(tmp);
                if(d<thetaArraySize){
                    if(i<thetaPowerArraySize[d]){
                        if(thetaPowers[d][i]==NULL){
                            thetaPowers[d][i]=ans->copy();
                            thetaTruncations[d][i]=trunc0;
                        }else{
                            if(thetaTruncations[d][i]<trunc0){
                                delete thetaPowers[d][i];
                                thetaPowers[d][i]=ans->copy();
                                thetaTruncations[d][i]=trunc0;
                            }
                        }
                    }
                }
        }
        ans->normalize();
        delete tmp;
    }
    if(d<thetaArraySize){
        if(pow<thetaPowerArraySize[d]){
            if(thetaPowers[d][pow]==NULL){
                thetaPowers[d][pow]=ans->copy();
                thetaTruncations[d][pow]=trunc0;

            }
        }
    }

    return ans;

}

QZSeriesWH* QZSeriesWH::ThetaBlockEfficient(int trunc0, int weight, vector<int> d){
    return ThetaBlockEfficient(trunc0, weight, d, 0);
}
QZSeriesWH* QZSeriesWH::ThetaBlockEfficient(int trunc0, int weight, vector<int> d, int allowTweak){
    int ell, t,i;
    ell=0;
    for(i=0;i<d.size();i++){ell+=signum(d[i]);}
    t=2*weight-ell;
    if(!allowTweak)if((3*ell+t)%24!=0){
        cout<<"Theta block does not have 24|3L+t (weight, L,t)=("<<weight<<","<<ell<<","<<t<<"). Abort.\n";exit(1);
    }
    QZSeriesWH* ans = QZSeriesWH::EtaFunctionPower(trunc0, t);
    int maxD=0, absD, oneD, pow;
    for(i=0;i<d.size();i++){
        absD=abs(d[i]);
        if(absD>maxD){maxD=absD;}
    }
    if(maxD<=0){
        cout<<"ERROR fdhsafwf in ThetaBlockEfficient maxD<=0 ("<<maxD<<"). Abort.\n"; exit(1);
    }
    QZSeriesWH* th;
    for(oneD=1;oneD<=maxD;oneD++){
        pow=0;
        for(i=0;i<d.size();i++)if(d[i]==oneD){pow++;}
        if(pow>0){
            th=ThetaFunctionDPower(trunc0, oneD, pow);;
            ans->multiplyWith(th);
            delete th;
        }
    }
    int success;
    for(oneD=1;oneD<=maxD;oneD++){
        pow=0;
        for(i=0;i<d.size();i++)if(d[i]==-oneD){pow++;}
        if(pow>0){
            th=ThetaFunctionDPower(trunc0, oneD, pow);;
            success=ans->divideBy(th);
            if(!success){
                cout<<"ERROR dagaeyeru4333 in ThetaBlockEfficient: ";
                for(int j=0;j<d.size();j++)cout<<","<<d[j];
                cout<<". Abort.\n"; exit(1);
            }
            delete th;
        }
    }

    ans->normalize();
    if(!allowTweak)if(ans->qExpAdjustOver24!=0){
        cout<<"ERROR fdkhd0hdgk qExpAdjustOver24 is not zero in making a ThetaBlock. Abort.\n"; exit(1);
    }
    if(!allowTweak)if(ans->zExpAdjustOver2!=0){
        cout<<"ERROR dgasdgadgadg zExpAdjustOver2 is not zero in making a ThetaBlock. Abort.\n"; exit(1);
    }
    if(allowTweak)if(ans->qExpAdjustOver24!=12){
        cout<<"ERROR fdaf qExpAdjustOver24 is not 12 in making a ThetaBlock with Tweak. Abort.\n"; exit(1);
    }
    if(allowTweak)if(ans->zExpAdjustOver2!=1){
        cout<<"ERROR dsafasdf zExpAdjustOver2 is not 1 in making a ThetaBlock with Tweak. Abort.\n"; exit(1);
    }
    return ans;
}

QZSeriesWH* QZSeriesWH::ThetaBlockWithOpsjfFormat(int targetIndex, int trunc0, int weight, vector<int> d, vector<int> heckeOps){
    int trunc; int i, m, type, L, index;

    int resultIndex = targetIndex; //This is the rssulting index
    m=targetIndex;
    trunc=trunc0;
    for(i=heckeOps.size()-2;i>=0;i-=2){ //Run through Hecke ops backwards
        type=heckeOps[i]; L=heckeOps[i+1];
        if(type==1){//hecke up
            trunc=trunc*L;
            m=m/L;
        }else if(type==2){//hecke down
            trunc=myHeckeDownMaxq(m,trunc,L);
            m = m*L;
        }else if(type==3){//double up
            //trunc = trunc;  //do nothing
            m = m/(L*L);
        }else{
            cout<<"ERROR gdjtyie in ThetaBlockWithOps, illegal value of type: "<<type<<". Abort.\n";exit(1);
        }
    }
    index=m;

    QZSeriesWH* ans = ThetaBlockProductjfFormat(trunc, d); //need current index!!

    for(i=0;i<heckeOps.size();i+=2){
        type=heckeOps[i]; L=heckeOps[i+1];
        if(type==1){//hecke up
            ans->applyUp(index,weight,L);
            index = index*L;
        }else if(type==2){//hecke down
            ans->applyDown(index,weight,L);
            index=index/L;
        }else if(type==3){//double up
            ans->applyDoubleUp(L);
            index = index*L*L;
        }else{
            cout<<"ERROR ufgiifd in ThetaBlockWithOps, illegal value of type: "<<type<<". Abort.\n";exit(1);
        }
    }
    if(ans->getMaxExponent()<trunc0){
        cout<<"ERROR t439wts in ThetaBlockWithOps. Did not get desired max q exponent ("
        <<ans->getMaxExponent()<<","<<trunc0<<"). Abort\n"; exit(1);
    }
    //ans->truncate(trunc0);
    return ans;
}

QZSeriesWH* QZSeriesWH::ThetaBlockProductjfFormat(int trunc, vector<int> d){
    int n, type, ell, tVal, power, i, j, ctr, numProd;
    vector<int>tb;
    ctr=0;
    numProd=d[ctr];ctr++;
    QZSeriesWH *ans, *tmp;
    for(n=0;n<numProd;n++){
        type=d[ctr];ctr++;
        if(type==4){
            power=d[ctr];ctr++;
            if(n==0){ans=E4Power(trunc, power);}
            else{tmp=E4Power(trunc, power);ans->multiplyWith(tmp);delete tmp;}
        }else if(type==6){
            power=d[ctr];ctr++;
            if(n==0){ans=E6Power(trunc, power);}
            else{tmp=E6Power(trunc, power);ans->multiplyWith(tmp);delete tmp;}
        }else if(type==123){ //theta block
            ell=d[ctr];ctr++;
            tVal=d[ctr];ctr++;
            tb.clear();
            for(j=0;j<ell;j++){
                tb.push_back(d[ctr]);ctr++;
            }
            tmp=ThetaBlock(trunc, (ell+tVal)/2, tb);
            if(n==0){ans=tmp;}
            else{ans->multiplyWith(tmp);delete tmp;}

        }else{
            cout<<"NO SUCH TYPE. ABORT.\n";exit(1);
        }
    }
    return ans;
}
QZSeriesWH* QZSeriesWH::E4Power(int trunc0, int pow){
    //This could be made more efficient by using binary expansion of pw
    //But don't care right now.
    if(pow<0){cout<<"E4Power not implemented for negative powers. Abort.\n";exit(1);}
    if(pow>0){
        if(pow<E4PowerArraySize){
            if(E4Powers[pow]!=NULL){
                if(E4Truncations[pow]>=trunc0){
                    return E4Powers[pow]->copy();
                }else{
                    delete E4Powers[pow]; E4Powers[pow]=NULL;E4Truncations[pow]=0;
                }
            }
        }
    }

    if(pow==0){return new QZSeriesWH(trunc0,1);}
    QZSeriesWH* ans;
    ans=E4(trunc0);
    if(abs(pow)==1)return ans;
    int i;
    QZSeriesWH* tmp = ans->copy();
    for(i=2;i<=abs(pow);i++)ans->multiplyWith(tmp);
    ans->normalize();
    delete tmp;
    if(pow>0){
        if(pow<E4PowerArraySize){
            if(E4Powers[pow]==NULL){
                E4Powers[pow]=ans->copy();
                E4Truncations[pow]=E4Powers[pow]->getMaxExponent();
            }
        }
    }

    return ans;
}
QZSeriesWH* QZSeriesWH::E6Power(int trunc0, int pow){
    //This could be made more efficient by using binary expansion of pw
    //But don't care right now.
    if(pow<0){cout<<"E6Power not implemented for negative powers. Abort.\n";exit(1);}
    if(pow>0){
        if(pow<E6PowerArraySize){
            if(E6Powers[pow]!=NULL){
                if(E6Truncations[pow]>=trunc0){
                    return E6Powers[pow]->copy();
                }else{
                    delete E6Powers[pow]; E6Powers[pow]=NULL;E6Truncations[pow]=0;
                }
            }
        }
    }

    if(pow==0){return new QZSeriesWH(trunc0,1);}
    QZSeriesWH* ans;
    ans=E6(trunc0);
    if(abs(pow)==1)return ans;
    int i;
    QZSeriesWH* tmp = ans->copy();
    for(i=2;i<=abs(pow);i++)ans->multiplyWith(tmp);
    ans->normalize();
    delete tmp;
    if(pow>0){
        if(pow<E6PowerArraySize){
            if(E6Powers[pow]==NULL){
                E6Powers[pow]=ans->copy();
                E6Truncations[pow]=E6Powers[pow]->getMaxExponent();
            }
        }
    }

    return ans;
}
QZSeriesWH* QZSeriesWH::E4(int trunc0){
    QZSeriesWH* ans = new QZSeriesWH();
    ans->qzs=QZSeriesZ::E4(trunc0);
    ans->qMinExp=0;
    ans->qExpAdjustOver24=0;
    ans->zExpAdjustOver2=0;
    return ans;
}
QZSeriesWH* QZSeriesWH::E6(int trunc0){
    QZSeriesWH* ans = new QZSeriesWH();
    ans->qzs=QZSeriesZ::E6(trunc0);
    ans->qMinExp=0;
    ans->qExpAdjustOver24=0;
    ans->zExpAdjustOver2=0;
    return ans;
}

void QZSeriesWH::getABC(int&A, int&B, int&C){
    fmpz_t tmp;
    fmpz_init_set_si(tmp, 0);
    int r, c; A=0; B=0; C=0;
    for(r=0;r<=getCoeffMaxr(0);r++){
        getFC_noreduce(tmp, 0, r);
        c=fmpz_get_si(tmp);
        if(r==0){
            A+=c;
        }else{
            A+=2*c;
            B+=r*c;
            C+=r*r*c;
        }
    }
    if(A%24==0){A=A/24;}else{cout<<"WARNING: dfat43y A is not integral ("<<A<<"/24), aborting to be safe.";exit(1);}
    if(B%2==0){B=B/2;}else{cout<<"WARNING: fhehsdg B is not integral ("<<B<<"/2), aborting to be safe.";exit(1);}
    if(C%2==0){C=C/2;}else{cout<<"WARNING: dg5u4gda A is not integral ("<<C<<"/2), aborting to be safe.";exit(1);}
    fmpz_clear(tmp);
}

void QZSeriesWH::saveToFile(FILE *f){
    fprintf(f, "%d\n", qMinExp);
    fprintf(f, "%d\n", qExpAdjustOver24);
    fprintf(f, "%d\n", zExpAdjustOver2);
    qzs->saveToFile(f);
}
QZSeriesWH::QZSeriesWH(FILE *f){
    fscanf(f, "%d", &qMinExp);
    fscanf(f, "%d", &qExpAdjustOver24);
    fscanf(f, "%d", &zExpAdjustOver2);
    qzs=new QZSeriesZ(f);
}

QZSeriesWH* QZSeriesWH::makeTBEcombo(ThetaBlockList *tbList,
                                     fmpz_t* combo, fmpz_t denom, int maxdegree){
    int tbi, i, success;
    QZSeriesWH* jftb, *jftbV2, *ans;
    ans=new QZSeriesWH(maxdegree);
    vector<int> d;
    cout<<"Using theta blocks of weight "<<tbList->weight<<", index "<<tbList->index<<".\n";
    cout<<"Reading "<<tbList->num<<":";
    for(tbi=0;tbi<tbList->num;tbi++){
        cout<<" "<<tbi;
        if(fmpz_is_zero(combo[tbi])){
            //cout<<"Skipping because of zero combo coefficient.\n";
        }else{
            d=tbList->d[tbi];
            jftb = QZSeriesWH::ThetaBlock(maxdegree*2,tbList->weight,d);
            jftbV2 = jftb->copy();
            jftbV2->applyUp(tbList->index,tbList->weight,2);
            success=jftbV2->divideBy(jftb);  //Later check if this was possible at all.
            if(!success){
                cout<<"ERROR32twe  jf|V2 not divisible by jf (must be a TBWD?). Abort. \n";
                exit(1);
            }
            ans->addScalarMultipleWith(combo[tbi],jftbV2);
            delete jftb; delete jftbV2;
        }
    }
    cout<<"\n";
    success=ans->divideByIntegerWithCheck(denom);
    if(!success){
        cout<<"j6l;iret ans not divisible by denom ("<<fmpz_get_str(NULL, 10,denom)
            <<"). Abort.\n";
        exit(1);
    }
    return ans;
}
vector<slong> QZSeriesWH ::getDivisors(slong a){
  vector<slong>v;
  slong i;
  for(i=1;i<=a;i++){
    //very inefficient!!
    if((a%i)==0)v.push_back(i);
  }
  return v;
}
slong QZSeriesWH ::GCDL(slong a, slong b){
  slong tmp;
  if(a<0)a=-a;
  if(b<0)b=-b;
  if(a>b){tmp=a;a=b;b=tmp;};
  while(a!=0) {
    tmp=b%a;
    b=a;
    a=tmp;
  }
  return b;
}

void QZSeriesWH::getFCGrit(fmpz_t fc, slong n, slong r, slong m, slong weight, slong level){
    vector<slong> d=getDivisors(GCDL(GCDL(n, r), m));
    fmpz_zero(fc);
    int i,j;
    fmpz_t tmp; fmpz_init(tmp);
    //cout<<"div.size="<<d.size()<<"\n";
    if(weight<1){cout<<"getFCGrit not implemented for weight ("<<weight<<") < 1. Abort.\n";exit(1);}
    for(i=0;i<d.size();i++){
        //ans+=mypower((int)d[i],weight-1)*((long long)getcnr((m/d[i])*(n/d[i]),r/d[i], (long long)p));
        //cout<<"HERE Z1.\n";
        getFC(tmp,(m/d[i])*(n/d[i]),r/d[i], level);
        //cout<<"HERE Z2.\n";
        if(d[i]!=1)for(j=0;j<weight-1;j++){
            fmpz_mul_si(tmp,tmp,d[i]);
        }
        fmpz_add(fc,fc,tmp);
    }
    fmpz_clear(tmp);
}
void QZSeriesWH::getFCGritTweak(fmpz_t fc, slong nx2, slong rx2, slong mx2, slong weight, slong level){
    /***
    ONLY FOR ODD SQUARE FREE LEVEL
    ***/
    vector<slong> d=getDivisors(GCDL(GCDL(nx2, rx2), mx2));
    fmpz_zero(fc);
    int i,j;
    fmpz_t tmp; fmpz_init(tmp);
    for(i=0;i<d.size();i++){
        //ans+=mypower((int)d[i],weight-1)*((long long)getcnrx2((mx2/d[i])*(nx2/d[i]),rx2/d[i], (long long)p));
    getFCx2(tmp,(mx2/d[i])*(nx2/d[i]),rx2/d[i], level);
        if(d[i]!=1)for(j=0;j<weight-1;j++){
            fmpz_mul_si(tmp,tmp,d[i]);
        }
        fmpz_add(fc,fc,tmp);
    }
    fmpz_clear(tmp);
}
void QZSeriesWH::getFCx2(fmpz_t fc, slong n, slong r, slong index){
    //Assumes index is odd
    //cout<<"getFCx2: (n,r,m) = ("<<n<<","<<r<<","<<index<<")\n";
    slong longp2=2L*index;
    slong rnew=((r%longp2)+longp2)%longp2;
    if(rnew>index)rnew-=longp2;  //rnew is between -index and index
    slong nnew=n-((r-rnew)/(longp2))*((r+rnew)/(2L)); //improvement to prevent overflow
    if((4L*index*nnew-rnew*rnew)!=(4L*index*n-r*r)){
        cout<<"ERROR dsafasdgew in getFCx2 of QZSeriesWH\n";
        cout<<"index="<<index<<"\n";
        cout<<"r="<<r<<"\n";
        cout<<"n="<<n<<"\n";
        cout<<"nnew="<<nnew<<"\n";
        cout<<"rnew="<<rnew<<"\n";
        cout<<(4*index*nnew-rnew*rnew)<<","<<(4*index*n-r*r)<<"\n";
        exit(1);
    }
    if((qExpAdjustOver24==12)&&(zExpAdjustOver2==1)){
        if(nnew%2!=0)if(rnew%2!=0){
            slong nlookup = (nnew-1)/2, rlookup = (rnew-1)/2;
            getFC_noreduce(fc, nlookup, rlookup);
            //long long a = (rnewx2 - rx2)/longp2;
            //long long sign=1;
            //if((a%2)!=0)sign=-1;
            //return sign*c[nlookup][rlookup+mr[nlookup]];
            if(((rnew-r)/(2*index)%2)!=0){
                fmpz_neg(fc,fc);
            }
            return;
        }
        cout<<"index="<<index<<"\n";
        cout<<"r="<<r<<"\n";
        cout<<"n="<<n<<"\n";
        cout<<"nnew="<<nnew<<"\n";
        cout<<"rnew="<<rnew<<"\n";
        cout<<(4*index*nnew-rnew*rnew)<<","<<(4*index*n-r*r)<<"\n";
        cout<<"CAUTION: getFCx2 called when nnew,rnew not both odd. Abort to be safe.\n";exit(1);
    }
    cout<<"getFCx2 called for series that has qExpAdjustOver24,zExpAdjustOver2 ="<<
        qExpAdjustOver24<<","<<zExpAdjustOver2<<". Abort to be safe.\n"; exit(1);

}


/******

*****/
