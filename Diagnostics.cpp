#include "Diagnostics.h"

ofstream Diagnostics::f;
int Diagnostics::checkSymmetry;
int Diagnostics::saveDiagnosticFile;
int Diagnostics::debugMode=0;
int Diagnostics::skipMathematica=0;
int Diagnostics::useNegativeSingularOnly=0;
fmpz_t Diagnostics::HumVecEntry00;

void Diagnostics::init(){
    fmpz_init(HumVecEntry00);
    init("diagnostic-file.ma");
}

void Diagnostics::init(string fname){
    f.open(fname.c_str());
    if(f.fail()){cout<<"File write open failure: "<<fname<<"\n";exit(1);}
    checkSymmetry=1;
    saveDiagnosticFile=1;

}
void Diagnostics::printVector(vector<int> c){
    cout<<"{";
    for(int i=0;i<c.size();i++){
        if(i>0)cout<<",";
        cout<<c[i];
    }
    cout<<"}";
}

string Diagnostics::vectorToString(vector<int> c){
    string ans; stringstream ss;
    ss<<"{";
    for(int i=0;i<c.size();i++){
        if(i>0)ss<<",";
        ss<<c[i];
    }
    ss<<"}";
    ss>>ans;
    return ans;
}
