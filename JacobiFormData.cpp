#include "JacobiFormData.h"

JacobiFormData::JacobiFormData( int tee0, vector<int> d0 ){
  tee=tee0;
  d=d0;
  jfType=0;
  int i;
  weight=tee;
  index=0;
  for(i=0;i<d.size();i++){
    if(d[i]==0){
       cout<<"LIKELY ERROR df8asd8: d=0 in JacobiFormData. Aborting to be safe.\n";
       exit(1);
    }else if(d[i]>0){
       weight++;
       index+=d[i]*d[i];
    }else{
       weight--;
       index-=d[i]*d[i];
    }
  } 
  if(weight%2!=0){
       cout<<"LIKELY ERROR 48fgik: nonintegral weight in JacobiFormData. Aborting to be safe.\n";
       exit(1);
  }else weight=weight/2;
  if(index%2!=0){
       cout<<"LIKELY ERROR 49yhdrgsd: nonintegral index in JacobiFormData. Aborting to be safe.\n";
       exit(1);
  }else index=index/2;
}
JacobiFormData::JacobiFormData(
        int jfType0, int weight0, int index0, int tee0, vector<int> d0 
    ){
  tee=tee0;
  d=d0;
  jfType=jfType0;
  weight=weight0;
  index=index0;
}
JacobiFormData::~JacobiFormData(){
}

QZSeriesWH* JacobiFormData::makeJF(int trunc){
  return QZSeriesWH::ThetaBlockWithOps(trunc, weight, d);
}

