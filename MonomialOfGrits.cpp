#include "MonomialOfGrits.h"

MonomialOfGrits::MonomialOfGrits( vector<int> ind0, vector<int>exponent0){
  ind=ind0;
  exponent=exponent0;
}

MonomialOfGrits::~MonomialOfGrits(){
}


MySeriesTruncMod* MonomialOfGrits::evaluate(MySeriesTruncMod**qs){
  return evaluate(qs, qs[0]->getTrunc()-1);
}
MySeriesTruncMod* MonomialOfGrits::evaluate(MySeriesTruncMod**qs, int upto){
  int COMMENT=0;
  MySeriesTruncMod*ans;
  ans = new MySeriesTruncMod(qs[0]->modn, qs[0]->getTrunc()-1);
  ans->setCoeff((ulong)1, 0);
  if(COMMENT){cout<<"HERE D: ans->van, ans->len="<<ans->van<<","<<ans->len<<"\n";
     cout<<"  and  qs[0]->getTrunc() = "<< qs[0]->getTrunc()<<"\n";
     ans->printstdoutWithTruncationOrder();cout<<"\n";
  }
  if(ind.size()==0){ //RETURN 1
    return ans;
  }
  int i,j;
  if(COMMENT){cout<<"HERE x 1: num inds="<<ind.size()<<"\n";}
  if(COMMENT){cout<<"So far [A] ans=";ans->printstdoutWithTruncationOrder();cout<<"\n";}
  for(i=0;i<ind.size();i++){
    if(exponent[i]<0){
      cout<<"Exponent "<<exponent[i]<<" not permitted in MonomialOfGrits::evaluate. Abort\n";
    }
    if(exponent[i]==0){//Do nothing
    }else{
      if(COMMENT){cout<<"ind["<<i<<"], exponent="<<ind[i]<<","<<exponent[i]<<"\n";}
      ans->multiplyByPower(qs[ind[i]], exponent[i]);
      if(COMMENT){cout<<"So far [B] ans=";ans->printstdoutWithTruncationOrder();cout<<"\n";}
    }
  }
  return ans;
}

void MonomialOfGrits::printRepresentation(){
  int i,j;
  int COMMENT=0;
  if(COMMENT){cout<<"printRepresentation\n";}
  if(COMMENT){cout<<"indSize="<<ind.size()<<"\n";}
  for(i=0;i<ind.size();i++){
      if(i>0)cout<<"*";
      cout<<"["<<ind[i]<<"]^"<<exponent[i];
  }
  return;
}

