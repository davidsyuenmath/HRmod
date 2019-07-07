#include "EfficientQuadratic.h"

EfficientQuadratic::EfficientQuadratic(vector<int>coeff, vector<int>TBi, vector<int>TBi2){
    int i,j, m,a,b;
    m=0;
    for(i=0;i<coeff.size();i++){
        if(TBi[i]>m){m=TBi[i];}
        if(TBi2[i]>m){m=TBi2[i];}
    }
    numTB=m+1;
    qc=new int*[numTB];
    for(i=0;i<numTB;i++){
        qc[i]=new int[numTB];
        for(j=0;j<numTB;j++){
            qc[i][j]=0;
        }
    }
    for(i=0;i<coeff.size();i++){
        if(TBi[i]<=TBi2[i]){a=TBi[i];b=TBi2[i];}
        else{b=TBi[i];a=TBi2[i];}
        qc[a][b]+=coeff[i];
    }
    fmpz_init(tmp);
}


EfficientQuadratic:: ~EfficientQuadratic(){
    int i;
    if(numTB>0){
        for(i=0;i<numTB;i++){
            delete[] qc[i];
        }
        delete[] qc;
    }
    fmpz_clear(tmp);
}

MySeriesTruncMod* EfficientQuadratic::evaluate(MySeriesTruncMod**qs){
    MySeriesTruncMod* ans, *t, *s;
    int i,j; int notYetStarted=1;  int isZero;
    for(i=0;i<numTB;i++){
        isZero=1;
        for(j=0;j<numTB;j++){
            if(qc[i][j]!=0){
                if(isZero){
                    s = new MySeriesTruncMod(qs[j]);
                    if(qc[i][j]!=1){
                        fmpz_set_si(tmp,qc[i][j]);
                        s->scalarMultiplityBy(tmp);
                    }
                    isZero=0;
                }else{
                    if(qc[i][j]==1){
                        s->addBy(qs[j]);
                    }else{
                        t=new MySeriesTruncMod(qs[j]);
                        fmpz_set_si(tmp, qc[i][j]);
                        t->scalarMultiplityBy(tmp);
                        s->addBy(t);
                        delete t;
                    }
                }
            }
        }
        if(!isZero){
            s->multiplyBy(qs[i]);
            if(notYetStarted){
                ans=s;
                notYetStarted=0;
            }else{
                ans->addBy(s);
                delete s;
            }
        }
    }
    if(notYetStarted){
        ans=new MySeriesTruncMod(qs[0]->modn,qs[0]->trunc-1);
    }
    return ans;

}

