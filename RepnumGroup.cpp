#include "RepnumGroup.h"

RepnumGroup::RepnumGroup(int level0, int mindet0){
    level=level0;
    mindet=mindet0;
    a.clear();b.clear();c.clear();rpn.clear();
}
RepnumGroup::RepnumGroup(int level0){
    level=level0;
    mindet=1;
    a.clear();b.clear();c.clear();rpn.clear();
}

RepnumGroup:: ~RepnumGroup(){
    int i;//destructor
    for(i=0;i<rpn.size();i++){
        delete rpn[i];
    }
    a.clear();b.clear();c.clear();rpn.clear();
}

int RepnumGroup::findIndex(slong sa, slong sb, slong sc){
    int i;
    for(i=0;i<a.size();i++){
        if(sa==a[i])if(sb==b[i])if(sc==c[i]){
            return i;
        }
    }
    return -1; //NOT FOUND
}

void RepnumGroup::makeRepnumList(slong sa, slong sb, slong sc, int uptoDot){
    //Check if a,b,c already exists
    int ii=findIndex(sa,sb,sc);
    if(ii>=0)if(rpn[ii]->uptoDot>=uptoDot){return;}
    MyRepnumList* rrr=new MyRepnumList(sa,sb,sc,level,uptoDot,mindet);
    if(ii>=0){
        delete rpn[ii];
        rpn[ii]=rrr;
        return;
    }else{
        a.push_back(sa);
        b.push_back(sb);
        c.push_back(sc);
        rpn.push_back(rrr);
        return;
    }
}

MyRepnumList* RepnumGroup::getRepnumList(slong sa, slong sb, slong sc, int uptoDot){
    int ii=findIndex(sa,sb,sc);
    if(ii>=0)if(rpn[ii]->uptoDot>=uptoDot){return rpn[ii];}
    makeRepnumList(sa,sb,sc,uptoDot);
    ii=findIndex(sa,sb,sc);
    if(ii>=0)if(rpn[ii]->uptoDot>=uptoDot){return rpn[ii];}
    cout<<"ERROR g8ads8gasg0 in getRepnumList. This should not have happened. Abort.\n";exit(1);
}
