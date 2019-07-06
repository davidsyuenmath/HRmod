#include "LaurentZ.h"

LaurentZ::LaurentZ(){
  vector<int> coeffs;
  coeffs.push_back(0);
  init(coeffs, 0);
  normalize();
}
LaurentZ::LaurentZ(int c){
  vector<int> coeffs;
  coeffs.push_back(c);
  init(coeffs, 0);
  normalize();
}
LaurentZ::LaurentZ(fmpz_t c){
  minexp=0;
  fmpz_poly_init(p);
  fmpz_poly_set_fmpz(p, c);
  fmpz_init(tmp);
}
LaurentZ::LaurentZ(int c, int minexp0){
  vector<int> coeffs;
  coeffs.push_back(c);
  init(coeffs, minexp0);
  normalize();
}
LaurentZ::LaurentZ(vector<int> coeffs){
  init(coeffs, 0);
  normalize();
}

LaurentZ::LaurentZ(vector<int> coeffs, int minexp0){
  init(coeffs, minexp0);
  normalize();
}

void LaurentZ::init(vector<int> coeffs, int minexp0){
  minexp=minexp0;
  fmpz_poly_init(p);
  int i;
  for(i=0;i<coeffs.size();i++){
    fmpz_poly_set_coeff_si(p, i, coeffs[i]);
  }
  fmpz_init(tmp);
}

LaurentZ:: ~LaurentZ(){//destructor
  fmpz_poly_clear(p);
  fmpz_clear(tmp);
}

void LaurentZ::normalize(){
  //Adjust so constant term of p is nonzero (unless p is zero)
  int deg = fmpz_poly_degree(p);
  int i,ord=0;
  if(deg<0){ //poly is zero
    minexp=0;
    return;
  }
  for(i=0;i<=deg;i++){
    fmpz_poly_get_coeff_fmpz(tmp, p, i);
    if(fmpz_is_zero(tmp)){
      if(deg==0){
        //poly is zero!
        minexp=0;
        return;
      }//else keep going to find first nonzero coeff
    }else{
      ord=i;
      shiftExponents(-ord);
      return;
    }
  }
  cout<<"ERROR in LaurentZ normalize() ord not found. deg = "<<deg<<". Abort.\n"; exit(1);
}

void LaurentZ::shiftExponents(int shift){ //SHOULD Use internal shifting routines!

  if(shift>0){ //multiply by z^shift
    fmpz_poly_shift_left(p, p, shift);
    minexp=minexp-shift;
  }else if(shift<0){
    fmpz_poly_shift_right(p, p, -shift);
    minexp=minexp-shift;
  }
}

bool LaurentZ::isZero(){
  return fmpz_poly_is_zero(p);
}

string LaurentZ::getString(){
    return getString("z");
}
string LaurentZ::getString(string var){
  if(fmpz_poly_is_zero(p))return "0";
  string ans, str;
  stringstream ss;
  slong i, n;
  ans="";
  int isEmpty=1;
  for(i=fmpz_poly_degree(p);i>=0;i--){
    n=i+minexp;
    fmpz_poly_get_coeff_fmpz(tmp, p, i);
    if(!fmpz_is_zero(tmp)){
      if(isEmpty){isEmpty=0;}else{ans = ans + " + ";}
      ans = ans + fmpz_get_str(NULL,10,tmp);
      if(n==0){

      }else{
        ans = ans + "*" + var;
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

LaurentZ* LaurentZ::copy(){
  LaurentZ* ans=new LaurentZ();
  fmpz_poly_set(ans->p, p);
  ans->minexp = minexp;
  return ans;
}
LaurentZ* LaurentZ::multiply(LaurentZ* a, LaurentZ *b){
  LaurentZ *ans = a->copy();
  ans->multiplyWith(b);
  return ans;
}
void LaurentZ::multiplyWith(LaurentZ* b){
  fmpz_poly_mul(p, p, b->p);
  minexp = minexp + b->minexp;
}
LaurentZ* LaurentZ::add(LaurentZ* a, LaurentZ* b){
  LaurentZ *ans = a->copy();
  ans->addWith(b);
  return ans;
}
void LaurentZ::addWith(LaurentZ* b){
  //Decide on common denominator, so to speak.
  if(minexp > b->minexp){shiftExponents(minexp - b->minexp);}
  if(minexp < b->minexp){b->shiftExponents(b->minexp - minexp);}
  if(minexp != b->minexp){
    cout<<"ERROR 438t534t in LaurentZ::addWith. Abort.\n"; exit(1);
  }
  fmpz_poly_add(p, p, b->p);
  normalize();
  b->normalize();
}
void LaurentZ::subtractWith(LaurentZ* b){
  //Decide on common denominator, so to speak.
  if(minexp > b->minexp){shiftExponents(minexp - b->minexp);}
  if(minexp < b->minexp){b->shiftExponents(b->minexp - minexp);}
  if(minexp != b->minexp){
    cout<<"ERROR 438t534t in LaurentZ::addWith. Abort.\n"; exit(1);
  }
  fmpz_poly_sub(p, p, b->p);
  normalize();
  b->normalize();
}

LaurentZ* LaurentZ::addScalarMultiple(LaurentZ *a, fmpz_t m, LaurentZ* b){
  LaurentZ *ans = a->copy();
  ans->addScalarMultipleWith(m, b);
  return ans;
}
void LaurentZ::addScalarMultipleWith(fmpz_t m, LaurentZ* b){
  //Decide on common denominator, so to speak.
  if(minexp > b->minexp){shiftExponents(minexp - b->minexp);}
  if(minexp < b->minexp){b->shiftExponents(b->minexp - minexp);}
  if(minexp != b->minexp){
    cout<<"ERROR 438t534t in LaurentZ::addWith. Abort.\n"; exit(1);
  }
  fmpz_poly_scalar_addmul_fmpz(p, b->p, m);
  normalize();
  b->normalize();
}
int LaurentZ::divideWith(LaurentZ* b){
  //returns 1 if success and no remainder
  //CAUTION: divrem does NOT allow for aliasing!! must use fresh variables!
  //cout<<"p=";fmpz_poly_print(p); cout<<endl;
  //cout<<"b->p=";fmpz_poly_print(b->p); cout<<endl;

  fmpz_poly_t rem;
  fmpz_poly_init(rem);
  fmpz_poly_t div;
  fmpz_poly_init(div);
  fmpz_poly_set(div, p);
  fmpz_poly_divrem(p, rem, div, b->p);
  //cout<<"ZZZ p="<<fmpz_poly_print(p); cout<<endl;

  minexp = minexp - b->minexp;
  int remIsZero = fmpz_poly_is_zero(rem);
  fmpz_poly_clear(rem);
  fmpz_poly_clear(div);
  return remIsZero;
}
void LaurentZ::setRemainder(fmpz_poly_t r, LaurentZ* b){
  //b must be doubly monic
  //Assumes you can use aliasing with rem. I.e. rem(r,r,b->p) is okay. Let's hope so.
  b->normalize();

  if(!(fmpz_is_pm1((fmpz_poly_lead(b->p))))){
    cout<<"In LaurentZ::setRemainder: Attempted remainder division by non-monic polynomial. Abort.\n";
    exit(1);
  }
  //cout<<"OK1 ";
  if(minexp>=0){
    fmpz_poly_shift_left(r, p, minexp);
    fmpz_poly_rem(r, r, b->p);
    return;
  }
  //cout<<"OK2 ";
  int constCoeffIsOne=1;
  fmpz_poly_get_coeff_fmpz(tmp, b->p, 0);
  if(!(fmpz_is_pm1(tmp))){
    cout<<"In LaurentZ::setRemainder: Attempted remainder division by non-constant-monic polynomial. Abort.\n";
    exit(1);
  }
  if(!(fmpz_is_one(tmp))){
    constCoeffIsOne=0;
  }
  LaurentZ* ans=new LaurentZ();
  fmpz_poly_rem(ans->p, p, b->p);
  ans->minexp = minexp;
  ans->normalize();
  //ans now holds the remainder, but it is Laurent.
  slong currexp=ans->minexp; //As a means to diagnose an infinite loop.
  while(ans->minexp<0){
    fmpz_poly_get_coeff_fmpz(tmp, ans->p, 0);
    if(constCoeffIsOne){fmpz_poly_scalar_submul_fmpz(ans->p,b->p,tmp);}
    else{fmpz_poly_scalar_addmul_fmpz(ans->p,b->p,tmp);}
    ans->normalize();
    if(ans->minexp<=currexp){
      cout<<"ERROR in LaurentZ::setRemainder, minexp should have increased.\n"<<
        "(currexp, minexp) = ("<<currexp<<","<<ans->minexp<<"). Abort.\n";
      exit(1);
    }
    currexp=ans->minexp;
  }
  fmpz_poly_shift_left(r, ans->p, ans->minexp);
  fmpz_poly_rem(r, r, b->p);

  //CHECK
  //LaurentZ* remainder=new LaurentZ(0);
  //remainder->minexp=0;
  //fmpz_poly_set(remainder->p,r);
  //LaurentZ* orig = copy();
  //orig->subtractWith(remainder); //This should be a multiple of b
  //int result=orig->divideWith(b);
  //cout<<"RESULT OF DIVISION = "<<result<<".\n";
  //delete remainder; delete orig;
  //
  delete ans;
}
LaurentZ* LaurentZ::scalarMultiply(LaurentZ* b, fmpz_t m){
  LaurentZ *ans = b->copy();
  ans->scalarMultiply(m);
  ans->normalize();
  return ans;
}
void LaurentZ::scalarMultiply(fmpz_t m){
  fmpz_poly_scalar_mul_fmpz(p, p, m);
}
void LaurentZ::multiply(LaurentZ* res, LaurentZ* a, LaurentZ *b){
  res->minexp = a->minexp + b->minexp;
  fmpz_poly_mul(res->p,a->p,b->p);
}

int LaurentZ::isUnit(){
  if(minexp==0){
    if(fmpz_poly_is_unit(p)){
      if(fmpz_poly_is_one(p)){return 1;}
      else{return -1;}
    }else{
    return 0;}
  }
  return 0;
}
LaurentZ* LaurentZ::BinomialZnMinusOne(int d){
    if(d<0){d=-d;}
    if(d==0){cout<<"ERROR 238dsf00 d=0 in BinomialZnMinusOne. Abort\n";exit(1);}
    vector<int> c; int i;
    c.push_back(1);
    for(i=1;i<d;i++){c.push_back(0);}
    c.push_back(-1);
    return new LaurentZ(c);
}
LaurentZ* LaurentZ::BinomialZnPlusZMinusn(int d){
    if(d<0){d=-d;}
    if(d==0){cout<<"ERROR 3523yhghf d=0 in BinomialZnPlusZMinusn. Abort\n";exit(1);}
    vector<int> c; int i;
    c.push_back(1);
    for(i=1;i<2*d;i++){c.push_back(0);}
    c.push_back(1);
    return new LaurentZ(c, -d);
}
LaurentZ* LaurentZ::BabyThetaBlock(vector<int> d){
    return BTB(d);
}
LaurentZ* LaurentZ::BTB(vector<int> d){
    vector<int> p, n; int i, sum=0;
    for(i=0;i<d.size();i++){
      sum+=d[i];
      if(d[i]==0){cout<<"ERROR kkht4332y6 d[i]=0 in BTB. Abort\n";exit(1);}
      if(d[i]>0){p.push_back(d[i]);}else{n.push_back(d[i]);}
    }
    if(sum%2!=0){
        cout<<"ERROR 32sd8aagd in BTB, sum of d[] not even. Abort\n";exit(1);
    }
    LaurentZ* ans = new LaurentZ(1, -sum/2);
    LaurentZ* tmpL;
    for(i=0;i<p.size();i++){
      tmpL=LaurentZ::BinomialZnMinusOne(p[i]);
      ans->multiplyWith(tmpL);
      delete tmpL;
    }
    int success;
    for(i=0;i<n.size();i++){
      tmpL=LaurentZ::BinomialZnMinusOne(-n[i]);
      success = ans->divideWith(tmpL);
      delete tmpL;
      if(success==0){
        cout<<"INVALID d-vector for BTB:\n";
        for(int j=0;j<d.size();j++){cout<<d[i]<<",";}
        cout<<" Abort.\n";
        exit(1);
      }
    }
    return ans;
}
LaurentZ* LaurentZ::Germ(int weight, vector<int> d){
    LaurentZ* ans = new LaurentZ(2*weight);
    int i;
    LaurentZ* tmpL;
    for(i=0;i<d.size();i++){
      tmpL=LaurentZ::BinomialZnPlusZMinusn(abs(d[i]));
      if(d[i]>0){
        ans->addWith(tmpL);
      }else{
        ans->subtractWith(tmpL);
      }
      delete tmpL;
    }
    return ans;
}
bool LaurentZ::divideByIntegerWithCheck(fmpz_t c){ //returns 1 if divide exactly
    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);
    fmpz_poly_scalar_mod_fmpz(tmp, p, c);
    bool success = fmpz_poly_is_zero(tmp);
    fmpz_poly_scalar_divexact_fmpz(p, p, c);
    fmpz_poly_clear(tmp);
    return success;
}
void LaurentZ::divideByIntegerWithoutCheck(fmpz_t c){ //assumes division will be exact
    fmpz_poly_scalar_divexact_fmpz(p, p, c);
}
void LaurentZ::getCoeff(fmpz_t coeff, int n){
    fmpz_t desperateTry;
    fmpz_init(desperateTry);
    if(QZSeriesWH::printGG==-24){

                        cout<<"In LaurentZ::getCoeff: p="<<fmpz_poly_get_str_pretty(p,"x")<<"\n";
                        cout<<"(minexp, n, coeff)<<"<<minexp<<","<<n<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf1 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
    if(n-minexp<0){
                if(QZSeriesWH::printGG==-24){
                        cout<<"(minexp, n, coeff)<<"<<minexp<<","<<n<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf2 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
        fmpz_set_si(coeff,0);
        if(QZSeriesWH::printGG==-24){
                        cout<<"(minexp, n, coeff)<<"<<minexp<<","<<n<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf3 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
    }else{
            if(QZSeriesWH::printGG==-24){
                        cout<<"(minexp, n, n-minexp, coeff)<<"<<minexp<<","<<n<<","<<n-minexp<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf4 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
        if(QZSeriesWH::printGG!=-24){}
        if(1){

                //fmpz_poly_get_coeff_fmpz(coeff,p,n-minexp);
                if(QZSeriesWH::printGG==-24){cout<<"GGz1 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";}
                fmpz_poly_get_coeff_fmpz(desperateTry,p,n-minexp);
                if(QZSeriesWH::printGG==-24){cout<<"GGz2 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";}
                fmpz_set(coeff,desperateTry);
                if(QZSeriesWH::printGG==-24){
                    cout<<"GGz3 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
                    fmpz_set_si(coeff,314159);
                    cout<<"GGz4 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
                    fmpz_poly_get_coeff_fmpz(coeff,p,n-minexp+1);
                    cout<<"GGz5 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
                    fmpz_set(coeff,desperateTry);
                    cout<<"GGz6 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
                    fmpz_poly_get_coeff_fmpz(fmpz_mat_entry(QZSeriesWH::diagMat,3,1),p,n-minexp);
                    cout<<"GGz7 diagMat(3,1): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,3,1))<<"\n";
                    cout<<"GGz7 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
                    fmpz_poly_get_coeff_fmpz(fmpz_mat_entry(QZSeriesWH::diagMat,4,1),p,n-minexp+1);
                    cout<<"GGz8 diagMat(4,1): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,4,1))<<"\n";
                    cout<<"GGz7 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
       }else{
            cout<<"GGff setting coeff manually.\n";
            fmpz_set_si(coeff, 1000*n+minexp);
        }

            if(QZSeriesWH::printGG==-24){
                        cout<<"(minexp, n, coeff)<<"<<minexp<<","<<n<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf5 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
    }
        if(QZSeriesWH::printGG==-24){
                        cout<<"(minexp, n, coeff)<<"<<minexp<<","<<n<<","<<fmpz_get_str(NULL,10,coeff)<<"\n";
                        cout<<"GGf6 diagMat(0,0): "<<fmpz_get_str(NULL,10,fmpz_mat_entry(QZSeriesWH::diagMat,0,0))<<"\n";
               }
    fmpz_clear(desperateTry);
}
LaurentZ* LaurentZ::expand(LaurentZ* b, int L){
    LaurentZ* ans = new LaurentZ(0);
    int i;
    ans->minexp = L* b->minexp;
    for(i=0;i<=fmpz_poly_degree(b->p);i++){
        fmpz_poly_get_coeff_fmpz(ans->tmp,b->p,i);
        fmpz_poly_set_coeff_fmpz(ans->p,L*i,ans->tmp);
    }
    return ans;
}
void LaurentZ::expand(int L){
    LaurentZ* ans = expand(this,L);
    minexp=ans->minexp;
    fmpz_poly_set(p, ans->p);
    delete ans;
}
LaurentZ* LaurentZ::shrink(LaurentZ* b, int L){
    LaurentZ* ans = new LaurentZ(0);
    int i, shift;
    shift = ((b->minexp%L)+L)%L;
    //cout<<"b->minexp="<<b->minexp<<"\n";
    //cout<<"shift="<<shift<<"\n";
    b->shiftExponents(shift);
    //cout<<"b->minexp="<<b->minexp<<"\n";
    if(b->minexp%L !=0){
        cout<<"ERROR 4dgjfdaa in shrink. Abort.\n";exit(1);
    }
    ans->minexp = b->minexp/L;
    for(i=0;i<=fmpz_poly_degree(b->p);i+=L){
        fmpz_poly_get_coeff_fmpz(ans->tmp,b->p,i);
        fmpz_poly_set_coeff_fmpz(ans->p,i/L,ans->tmp);
    }
    return ans;
}
void LaurentZ::shrink(int L){
    LaurentZ* ans = shrink(this,L);
    minexp=ans->minexp;
    fmpz_poly_set(p, ans->p);
    delete ans;
}

void LaurentZ::addToTerm(fmpz_t c, int n){
    //D'uh - you have to adjust minexp if need be!!
    //cout<<"UA (n,minexp)=("<<n<<","<<minexp<<")\n";
    int t=n-minexp; //The exponent of p to change.
    if(t<0){
        shiftExponents(-t);
        fmpz_poly_set_coeff_fmpz(p,0,c);

    }else{
        fmpz_poly_get_coeff_fmpz(tmp,p,t);
        //cout<<"UB\n";
        fmpz_add(tmp,tmp,c);
        //cout<<"UC\n";
        fmpz_poly_set_coeff_fmpz(p,t,tmp);
        //cout<<"UD\n";
    }
}

int LaurentZ::getMaxExponent(){
    return fmpz_poly_degree(p)+minexp;
}

int LaurentZ::equals(LaurentZ* b){
    normalize();
    b->normalize();
    if(minexp!=b->minexp){return 0;}
    if(fmpz_poly_equal(p, b->p)){return 1;}else{return 0;}
}

int LaurentZ::checkSymmetry(){
    normalize();
    if(fmpz_poly_is_zero(p)){return 1;}
    if(fmpz_poly_degree(p)!=-2*minexp){
        cout<<"minexp="<<minexp<<"\n";
        cout<<getString()<<"\n";
        return 0;
    }
    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);
    fmpz_poly_reverse(tmp, p, fmpz_poly_length(p));
    int ans = fmpz_poly_equal(p, tmp);
    fmpz_poly_clear(tmp);
    return ans;
}
void LaurentZ::negate(){
    fmpz_poly_scalar_mul_si(p,p,-1);
}

void LaurentZ::saveToFile(FILE *f){
    fprintf(f, "%d ", minexp);
    fprintf(f, "%ld", fmpz_poly_length(p));
    slong i;
    for(i=0;i<fmpz_poly_length(p);i++){
        fprintf(f, " ");
        fmpz_poly_get_coeff_fmpz(tmp, p, i);
        fmpz_fprint(f, tmp);
    }
    fprintf(f, "\n");
}

LaurentZ::LaurentZ(FILE *f){
    fmpz_poly_init(p);
    fmpz_init(tmp);
    fscanf(f, "%d", &minexp);
    int len;
    fscanf(f, "%d", &len);
    slong i;
    for(i=0;i<len;i++){
        fmpz_fread(f, tmp);
        fmpz_poly_set_coeff_fmpz(p, i, tmp);
    }
}
