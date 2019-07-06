#include "DRList.h"

DRList::DRList(){
    minDcalculated=0;
}
DRList::~DRList(){
}
void DRList::clear(){
    DList.clear(); RList.clear(); minDcalculated=0;
}
int DRList::alreadyOnList(int d, int r){
    for(int i=0;i<DList.size();i++){
        if((d==DList[i])&&(r==RList[i])){
            return 1;
        }
    }
    return 0;
}

void DRList::insert(int d, int r){
    if(!alreadyOnList(d,r)){
        DList.push_back(d);
        RList.push_back(r);
        if(minDcalculated){
            minD=min(minD, d);
        }else{
            minD=d;
            minDcalculated=1;
        }
    }
}
void DRList::insert(vector<int> d, vector<int> r)
{
    for(int i=0;i<d.size();i++){
        insert(d[i],r[i]);
    }

}
vector<int> DRList::getDList(){
    return DList;
}
vector<int> DRList::getRList(){
    return RList;
}
int DRList::getLength(){
    return DList.size();
}

int DRList::getMinD(){
    return minD;
}
void DRList::saveToFile(ofstream &f, string varname){
    f<<varname<<"={";
    int i;
    for(i=0;i<DList.size();i++){
        if(i>0){f<<",\n";}
        f<<"{"<<DList[i]<<","<<RList[i]<<"}";
    }
    f<<"\n};\n";
    f.flush();
}
void DRList::print(){
    cout<<"DRList={";
    int i;
    for(i=0;i<DList.size();i++){
        if(i>0){cout<<",\n";}
        cout<<"{"<<DList[i]<<","<<RList[i]<<"}";
    }
    cout<<"\n};\n";
}
