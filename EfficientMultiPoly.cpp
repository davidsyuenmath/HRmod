#include "EfficientMultiPoly.h"

EfficientMultiPoly::EfficientMultiPoly(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0
    ){
  init(coeffs0, ind0, indSize0, -1);
}
EfficientMultiPoly::EfficientMultiPoly(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0, int oneExisting0
    ){
  init(coeffs0, ind0, indSize0, oneExisting0);
}

void EfficientMultiPoly::init(
        fmpz_t* coeffs0, vector<int>* ind0, int indSize0, int oneExisting0
    ){
  oneExisting=oneExisting0;
  doRecurse=0; polyA=NULL; polyB=NULL;
  maxTBnum=0;
  indSize=indSize0;
  coeffs=new fmpz_t[indSize];
  ind=new vector<int>[indSize];
  int i, j;
  for(i=0;i<indSize;i++){
    fmpz_init_set(coeffs[i],coeffs0[i]);
    ind[i]=ind0[i];
    for(j=0;j<ind[i].size();j++){
      if(ind[i][j]>maxTBnum)maxTBnum=ind[i][j];
    }
    if(oneExisting<0){
      if(ind[i].size()>0)oneExisting=ind[i][0];
    }
  }
  makeEfficient();
  //printRepresentation();
  if(oneExisting<0){cout<<"ERROR in instantiating EfficientMultiPoly dfadfs. Could not find an existent used index. Abort.\n";exit(1);}
}

EfficientMultiPoly::~EfficientMultiPoly(){
  int i;
  for(i=0;i<indSize;i++){
    fmpz_clear(coeffs[i]);
  }
  delete[] coeffs;
  delete[] ind;
  if(doRecurse){
    delete polyA; delete polyB; doRecurse=0;
  }
}

void EfficientMultiPoly::makeEfficient(){
  //find candidate splittingIndex
  int best, i, c, count;
  best=0; splittingIndex=-1; doRecurse=0;
  for(c=0;c<maxTBnum;c++){
    count=0;
    for(i=0;i<indSize;i++){if(contains(ind[i],c))count++;}
    if(count>best){best=count;splittingIndex=c;}
  }
  if(best>1){
    vector<int>* indA = new vector<int>[indSize];
    vector<int>* indB = new vector<int>[indSize];
    fmpz_t* cA = new fmpz_t[indSize];
    fmpz_t* cB = new fmpz_t[indSize];
    int indASize=0, indBSize=0;
    for(i=0;i<indSize;i++){
      if(contains(ind[i],splittingIndex)){
        indA[indASize]=copyWithDeleteOne(ind[i],splittingIndex);
        fmpz_init_set(cA[indASize],coeffs[i]);
        indASize++;
      }else{
        indB[indBSize]=ind[i];
        fmpz_init_set(cB[indBSize],coeffs[i]);
        indBSize++;
      }
    }
    doRecurse=1;
    polyA=new EfficientMultiPoly(cA, indA, indASize, oneExisting);
    polyB=new EfficientMultiPoly(cB, indB, indBSize, oneExisting);
    for(i=0;i<indASize;i++){fmpz_clear(cA[i]);}
    for(i=0;i<indASize;i++){fmpz_clear(cB[i]);}
    delete[] indA; delete[] cA;
    delete[] indB; delete[] cB;
  }
}
int  EfficientMultiPoly::contains(vector<int>v, int a){
  int i;
  for(i=0;i<v.size();i++){if(v[i]==a)return 1;}
  return 0;
}
vector<int>  EfficientMultiPoly::copyWithDeleteOne(vector<int>v, int a){
  vector<int> ans;
  int deleted=0;
  for(int i=0;i<v.size();i++){
    if(v[i]!=a){ans.push_back(v[i]);}
    else{
      if(deleted){ans.push_back(v[i]);}
      else deleted=1;
    }
  }
  return ans;
}

MySeriesTruncMod* EfficientMultiPoly::evaluate(MySeriesTruncMod**qs){
  //cout<<"HERE Eff 2: doRecurse="<<doRecurse<<"\n";
  //printRepresentation();
  //cout<<"\n";
  if(!doRecurse)return evaluateNoFrills(qs);
  //cout<<"HERE Eff 3: Doing recursion.\n";
  MySeriesTruncMod* ans, *tmpqsB;
  ans=polyA->evaluate(qs);
  tmpqsB=polyB->evaluate(qs);
  ans->multiplyBy(qs[splittingIndex]);
  ans->addBy(tmpqsB);
  delete tmpqsB;
  return ans;
}

MySeriesTruncMod* EfficientMultiPoly::evaluateNoFrills(MySeriesTruncMod**qs){
  int COMMENT=0;
  if(COMMENT){cout<<"HERE Eff 0: oneExisting="<<oneExisting<<"\n";}
  MySeriesTruncMod*ans = new MySeriesTruncMod(qs[oneExisting]->modn, qs[oneExisting]->trunc-1);
  int i,j;
  MySeriesTruncMod* tmpqs;
  if(COMMENT){cout<<"HERE Eff 1: indSize="<<indSize<<"\n";}
  if(COMMENT){cout<<"So far [A] ans=";ans->printstdoutWithTruncationOrder();cout<<"\n";}
  for(i=0;i<indSize;i++){
    if(COMMENT){cout<<"i, ind[i].size()="<<i<<","<<ind[i].size()<<"\n";}
    if(ind[i].size()==0){ //Constant
      ans->updateCoeffAdd(coeffs[i],0);
      if(COMMENT){cout<<"So far [B] ans=";ans->printstdoutWithTruncationOrder();cout<<"\n";}
    }else{
      if(COMMENT){cout<<"ind["<<i<<"]["<<0<<"]="<<ind[i][0]<<"\n";}
      tmpqs=new MySeriesTruncMod(qs[ind[i][0]]);
      for(j=1;j<ind[i].size();j++){
        if(COMMENT){cout<<"ind["<<i<<"]["<<j<<"]="<<ind[i][j]<<"\n";}
        tmpqs->multiplyBy(qs[ind[i][j]]);
      }
      tmpqs->scalarMultiplyBy(coeffs[i]);
      ans->addBy(tmpqs);
      if(COMMENT){cout<<"So far [C] ans=";ans->printstdoutWithTruncationOrder();cout<<"\n";}
      delete tmpqs;
    }
  }
  return ans;
}

void EfficientMultiPoly::printRepresentation(){
  int i,j;
  int COMMENT=0;
  if(COMMENT){cout<<"printRepresentation\n";}
  if(COMMENT){cout<<"indSize,doRecurse="<<indSize<<","<<doRecurse<<"\n";}
  if(!doRecurse){
    for(i=0;i<indSize;i++){
      if(COMMENT){cout<<"i="<<i<<"\n";}
      if(COMMENT){cout<<"ind[i].size()="<<ind[i].size()<<"\n";}
      if(i>0)if(fmpz_cmp_si(coeffs[i],0)>=0)cout<<"+";
      fmpz_print(coeffs[i]);
      cout<<"[";
      for(j=0;j<ind[i].size();j++){
        if(j>0)cout<<" ";
        cout<<ind[i][j];
      }
      cout<<"]";
    }
    return;
  }
  cout<<"["<<splittingIndex<<"]*(";
  polyA->printRepresentation();
  cout<<") + ";
  polyB->printRepresentation();
}

