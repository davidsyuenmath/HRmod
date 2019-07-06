#include "ThetaBlockList.h"




ThetaBlockList::ThetaBlockList(int weight, int index, int ell){
    init("",weight,index,ell);
}
ThetaBlockList::ThetaBlockList(string directory, int weight, int index, int ell){
    init(directory,weight,index,ell);
}
ThetaBlockList::ThetaBlockList(string directory, int weight, int index, int ell, int withDenom){
    if(withDenom==0){
        init(directory,weight,index,ell);
    }else if(withDenom>0){
        initWD(directory,weight,index,ell,withDenom);
    }else{
        cout<<"Custom TBs with negative withDenom not implemented. Abort.\n"; exit(1);
        //initWDcustom(directory,weight,index,ell,withDenom);
    }
}

void ThetaBlockList::init(string directory, int weight0, int index0, int ell){
    weight=weight0; index=index0;
  ifstream TBFILE;
  string f;
  stringstream ffs;
  if(directory!=""){
    ffs<<directory<<"/";
  }
  ffs<<"allTB-"<<weight<<"-"<<index<<"-"<<ell<<".ma";
  ffs>>f;
  string command="";
  command = "perl grabComSkipCommentsAppendMinusOne.pl " + f;
  //if(Diagnostics::debugMode){command = "perl ..\\..\\..\\grabComSkipCommentsAppendMinusOne.pl " + f;}
  cout<<"B: Executing command "<<command<<".\n";
  system(command.c_str());
  f=f+".txt";
  cout<<"Trying to read spanning blocks from "<<f<<" ...";
  TBFILE.open(f.c_str());
  if(TBFILE.fail()){cout<<" no such AllTB file.\n";
     system("dir");
     exit( 1); }
  int x;
  TBFILE>>x;  //should be weight
  if(x!=weight){
    cout<<"First number in file does not match weight "<<weight<<"\n";exit(1);
  }
  TBFILE>>x;  //should be index
  if(x!=index){
    cout<<"Second number in file does not match index "<<index<<"\n";exit(1);
  }
  TBFILE>>x;  //should be ell
  if(x!=ell){
    cout<<"Third number in file does not match ell "<<ell<<"\n";exit(1);
  }
  arraysize=100;
  d=new vector<int>[arraysize];
  num=0;
  cout<<"Reading theta block d-values...\n";
  //vector<int>* spanningBlocks = new vector<int>[numBlocks];
  vector<int>c; int i;
  while(1){
    TBFILE>>x;
    if(x==-1){break;}
    c.clear();
    c.push_back(x);
    for(i=1;i<ell;i++){
        TBFILE>>x;
        c.push_back(x);
    }
    insertThetaBlock(c);
  }
  cout<<"done. "<<num<<" blocks read.\n";
}
ThetaBlockList::~ThetaBlockList(){//destructor
    destroy();
}
void ThetaBlockList::destroy(){
    delete[] d;
}

void ThetaBlockList::insertThetaBlock(vector<int> c){
    if(num>=arraysize){
        increaseSize();
    }
    d[num]=c;
    num=num+1;
}

void ThetaBlockList::increaseSize(){
    int newarraysize=2*arraysize+1;
    vector<int>* newd = new vector<int>[newarraysize];
    for(int i=0;i<num;i++){
        newd[i]=d[i];
    }
    delete[] d;
    d=newd;
    arraysize=newarraysize;
}

void ThetaBlockList::initWD(string directory, int weight0, int index0, int ell, int withDenom){
        weight=weight0; index=index0;
  ifstream TBFILE;
  string f;
  stringstream ffs;
  if(directory!=""){
    ffs<<directory<<"/";
  }
  if(withDenom==1){
    ffs<<"allTBWD-reformatted-"<<weight<<"-"<<index<<"-"<<ell<<".ma";
    ffs>>f;
  }else if(withDenom==2){
    ffs<<"allTBWD-only-"<<weight<<"-"<<index<<"-"<<ell<<".ma";
    ffs>>f;
  }else{
    ffs<<"allTBWD-custom"<<withDenom<<"-"<<weight<<"-"<<index<<"-"<<ell<<".ma";
    ffs>>f;
  }
  string command="";
  command = "perl grabComSkipCommentsAppendMinusOne.pl " + f;
  //if(Diagnostics::debugMode){command = "perl ..\\..\\..\\grabComSkipCommentsAppendMinusOne.pl " + f;}
  cout<<"A: Executing command "<<command<<".\n";
  system(command.c_str());
  f=f+".txt";
  cout<<"Trying to read spanning blocks from "<<f<<" ...";
  TBFILE.open(f.c_str());
  if(TBFILE.fail()){cout<<" no such JBasis file.\n";;
     system("dir");
      exit( 1); }
  int x;
  TBFILE>>x;  //should be weight
  if(x!=weight){
    cout<<"First number in file does not match weight "<<weight<<"\n";exit(1);
  }
  TBFILE>>x;  //should be index
  if(x!=index){
    cout<<"Second number in file does not match index "<<index<<"\n";exit(1);
  }
  TBFILE>>x;  //should be ell
  if(x!=ell){
    cout<<"Third number in file does not match ell "<<ell<<"\n";exit(1);
  }
  arraysize=100;
  d=new vector<int>[arraysize];
  num=0;
  cout<<"Reading theta block d-values...\n";
  //vector<int>* spanningBlocks = new vector<int>[numBlocks];
  vector<int>c; int i;
  int len;
  while(1){
    TBFILE>>len;
    if(len==-1){break;}
    c.clear();
    for(i=0;i<len;i++){
        TBFILE>>x;
        c.push_back(x);
    }
    insertThetaBlock(c);
  }
  cout<<"done. "<<num<<" blocks read.\n";
}
