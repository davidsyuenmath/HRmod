/****
  Encodes a set of (should be linearly independent) QZSeriesWH Jacobi form expansions
  9/7/2017: add capability to read a JBasis-wt-index.jf file
  1/7/2018: less verbose :-)
  1/9/2018: finally programmed general myLinearSolveOverIntegers
****/


#include "JacobiFormBasis.h"

JacobiFormBasis::JacobiFormBasis(){ //does nothing other than to make num=0
    num=0;
    jfArray=NULL;
    tbArray=NULL;
    maxrMade=0; DrMade=0;
}
JacobiFormBasis::JacobiFormBasis(int index0, int trunc){ //weight 0 nothing else
    num=0; weight=0; index=index0; maxDegree=trunc;
    jfArray=NULL;
    tbArray=NULL;
    maxrMade=0; DrMade=0;
}
JacobiFormBasis::JacobiFormBasis(int weight0, int index0, int trunc){
    weight=weight0; index=index0; maxDegree=trunc;
    maxrMade=0; DrMade=0;
    tbArray=NULL;
    jfArray=readJBasisFile(weight, index, trunc, num,NULL);
}
JacobiFormBasis::JacobiFormBasis(int weight0, int index0, int trunc,
                                 fmpz_t* combo){
    weight=weight0; index=index0; maxDegree=trunc;
    maxrMade=0; DrMade=0;
    tbArray=NULL;
    jfArray=readJBasisFile(weight, index, trunc, num,combo);
}
JacobiFormBasis::JacobiFormBasis(string directory, int weight0, int index0, int trunc){
    weight=weight0; index=index0; maxDegree=trunc;
    maxrMade=0; DrMade=0;tbArray=NULL;
    jfArray=readJBasisFile(directory,weight, index, trunc, num,NULL);}
JacobiFormBasis::JacobiFormBasis(string directory, int weight0, int index0, int trunc,
                                 fmpz_t* combo){
    weight=weight0; index=index0; maxDegree=trunc;
    maxrMade=0; DrMade=0;tbArray=NULL;
    cout<<"Hah?"<<directory<<"?\n";
    jfArray=readJBasisFile(directory,weight, index, trunc, num, combo);
}


QZSeriesWH** JacobiFormBasis::readJBasisFile(int weight, int index, int trunc, int &num){
    return readJBasisFile(".",weight,index,trunc,num,NULL);
}
QZSeriesWH** JacobiFormBasis::readJBasisFile(int weight, int index, int trunc, int &num,
                                 fmpz_t* combo){
    return readJBasisFile(".",weight,index,trunc,num,combo);
}

QZSeriesWH** JacobiFormBasis::readJBasisFile(string directory, int weight, int index, int trunc, int &numBlocks){
    return readJBasisFile(directory,weight,index,trunc,numBlocks,NULL);
}

QZSeriesWH** JacobiFormBasis::readJBasisFile(string directory, int weight, int index, int trunc, int &numBlocks,
                                 fmpz_t* combo){
  ifstream TBFILE;
  string f;
  stringstream ffs;
  if(directory!=""){
    ffs<<directory<<"/";
  }
  ffs<<"JBasis-"<<index<<"-"<<weight<<".ma.txt";
  ffs>>f;
  cout<<"Trying to read spanning blocks from "<<f<<" ...";
  TBFILE.open(f.c_str());
  if(TBFILE.fail()){
          cout<<" no such JBasis file.  Will try jf file.\n";
        return readjfJBasisFile(directory, weight, index, trunc, numBlocks);
  }
  int x, tVal, ell;
  TBFILE>>x;  //should be index
  if(x!=index){
    cout<<"First number in file does not match index "<<index<<"\n";exit(1);
  }
  TBFILE>>x;  //should be weight
  if(x!=weight){
    cout<<"Second number in file does not match weight "<<weight<<"\n";exit(1);
  }
  TBFILE>>x;  //should be t
  tVal=x;
  cout<<"tVal = "<<tVal<<"\n";
  //if(x!=tee){
  //  cout<<"Second number in file does not match t "<<tee<<"\n";exit(1);
  //}
  TBFILE>>x;  //should be ell
  ell=x;
  cout<<"ell = "<<ell<<"\n";
  //if(x!=ell){
  //  cout<<"Second number in file does not match ell "<<ell<<"\n";exit(1);
  // }
  TBFILE>>x;  //number of blocks to calculate
  numBlocks=x;
  //  cout<<"Number of blocks to calculate: "<<numBlocks<<"\n";

  cout<<"Reading theta block d-values...\n";
  //vector<int>* spanningBlocks = new vector<int>[numBlocks];
  vector<int> tb;
  int i, j, numOps;
  vector<int> h;
  QZSeriesWH** jfb =new QZSeriesWH*[numBlocks];

  for(i=0;i<numBlocks;i++){
    TBFILE>>numOps;
    h.clear(); tb.clear();
    cout<<"TB#"<<i<<"(of "<<numBlocks<<"): ";
    for(j=0;j<2*numOps;j++){
      TBFILE>>x;
      h.push_back(x);
      cout<<x<<",";
    }
    cout<<"//";
    for(j=0;j<ell;j++){
      TBFILE>>x;
      tb.push_back(x);
      cout<<x<<",";
    }
    cout<<"\n";

    if((combo==NULL)||(!fmpz_is_zero(combo[i]))){
        jfb[i]=QZSeriesWH::ThetaBlockWithOps(trunc,weight,tb,h);
    }else{
        cout<<"Skipping because combo coefficient is known to be zero.\n";
        jfb[i]=NULL;
    }
  }
  TBFILE.close();
  cout<<"done.\n";

  return jfb;
}

QZSeriesWH** JacobiFormBasis::readjfJBasisFile(string directory, int weight, int index, int trunc, int &numBlocks){
    return readjfJBasisFile(directory,weight,index,trunc,numBlocks,NULL);
}
QZSeriesWH** JacobiFormBasis::readjfJBasisFile(string directory, int weight, int index, int trunc, int &numBlocks,
                                               fmpz_t* combo){
  ifstream TBFILE;
  string f;
  stringstream ffs;
  cout<<"La la la#"<<directory<<"#.\n";
  if(directory!=""){
    ffs<<directory<<"/";
  }
  ffs<<"JBasis-"<<index<<"-"<<weight<<".jf";
  ffs>>f;
  cout<<"Trying to read spanning blocks from "<<f<<" ...";
  TBFILE.open(f.c_str());
  if(TBFILE.fail()){cout<<" no such JBasis file.\n"; exit( 1); }
  int x, tVal, ell;
  TBFILE>>x;  //should be index
  if(x!=index){
    cout<<"First number in file does not match index "<<index<<"\n";exit(1);
  }
  TBFILE>>x;  //should be weight
  if(x!=weight){
    cout<<"Second number in file does not match weight "<<weight<<"\n";exit(1);
  }

  TBFILE>>x;  //number of blocks to calculate
  numBlocks=x;
  //  cout<<"Number of blocks to calculate: "<<numBlocks<<"\n";

  cout<<"Reading theta block d-values...\n";
  //vector<int>* spanningBlocks = new vector<int>[numBlocks];
  vector<int> tb;
  int i, j, n, numOps, numProd, type;
  vector<int> h;
  QZSeriesWH** jfb =new QZSeriesWH*[numBlocks];

  for(i=0;i<numBlocks;i++){
    TBFILE>>numOps;
    h.clear(); tb.clear();
    cout<<"TB#"<<i<<"(of "<<numBlocks<<"): ";
    for(j=0;j<2*numOps;j++){
      TBFILE>>x;
      h.push_back(x);
      cout<<x<<",";
    }
    cout<<"//";

    TBFILE>>numProd;
    tb.push_back(numProd);
    cout<<numProd<<":";
    for(n=0;n<numProd;n++){
        TBFILE>>type;
        tb.push_back(type);
        cout<<"("<<type<<",";
        if((type==4)||(type==6)){
            TBFILE>>x;  //Exponent of E4 or E6
            tb.push_back(x);
            cout<<x<<")";
        }else if(type==123){ //theta block
            TBFILE>>ell;
            tb.push_back(ell);
            TBFILE>>tVal;
            tb.push_back(tVal);
            cout<<ell<<","<<tVal<<":";
            for(j=0;j<ell;j++){
                TBFILE>>x;
                tb.push_back(x);
                cout<<x;
                if(j<ell-1){cout<<",";}
            }
            cout<<")";
        }else{
            cout<<"NO SUCH TYPE. ABORT.\n";exit(1);
        }
        if(n<numProd-1){cout<<",";}
    }

    cout<<"\n";
    if((combo==NULL)||(!fmpz_is_zero(combo[i]))){
        jfb[i]=QZSeriesWH::ThetaBlockWithOpsjfFormat(index,trunc,weight,tb,h);
    }else{
        cout<<"Skipping because combo coefficient is known to be zero.\n";
        jfb[i]=NULL;
    }

   }
  TBFILE.close();
  cout<<"done.\n";

  //exit(0);

  return jfb;
}

JacobiFormBasis::~JacobiFormBasis(){ //destructor
    destroy();
}
void JacobiFormBasis::destroy(){
    int i;
    if((num>0)&&(jfArray!=NULL)){
        for(i=0;i<num;i++){
            if(jfArray[i]!=NULL){delete jfArray[i];}
        }
        delete[] jfArray;
    }
    if((num>0)&&(tbArray!=NULL)){delete[] tbArray;}
    num=0;
}
void JacobiFormBasis::divideBy(QZSeriesWH* jf, int weight0, int index0){
    int i;
    cout<<"divideBy...";
    for(i=0;i<num;i++){
        cout<<i<<",";
        if(Diagnostics::saveDiagnosticFile){
            if(i==7){
                Diagnostics::f<<"badIndex7="<<i<<";\n";
                jfArray[i]->saveToFile(Diagnostics::f,"badjfArrayJF7");
                jf->saveToFile(Diagnostics::f,"attemptedJF7");
            }
        }

        if(!(jfArray[i]->divideBy(jf))){
            cout<<"ERROR 5r9gdfgs JacobiFormBasis divideBy did not divide evenly at ["<<i<<"]. Aborting to be safe.\n";
            if(Diagnostics::saveDiagnosticFile){
                Diagnostics::f<<"badIndex="<<i<<";\n";
                jfArray[i]->saveToFile(Diagnostics::f,"badjfArrayJF");
                jf->saveToFile(Diagnostics::f,"attemptedJF");
            }
            exit(1);
        };
        if(i==0){maxDegree=jfArray[i]->getMaxExponent();}
        else{maxDegree=min(maxDegree, jfArray[i]->getMaxExponent());}
    }
    weight=weight-weight0;
    index=index-index0;
    maxrMade=0; DrMade=0;
    cout<<"\n";
}
void JacobiFormBasis::myInvert(fmpz_mat_t A){
    fmpz_t den; fmpz_init(den);
    fmpz_mat_inv(A, den, A);
    if(!fmpz_is_pm1(den)){cout<<"ERROR gjdtykb in myInvert: input is not unimodular. Abort.\n";exit(1);}
    fmpz_clear(den);
    if(!fmpz_is_one(den)){
        fmpz_mat_scalar_mul_si(A,A,-1);
    }
}
void JacobiFormBasis::mySetIdentityMatrix(fmpz_mat_t A){
    int n=fmpz_mat_nrows(A);
    int i,j;
    for(i=0;i<n;i++)for(j=0;j<n;j++){
        if(i==j){fmpz_set_si(fmpz_mat_entry(A,i,j),1);}
        else{fmpz_set_si(fmpz_mat_entry(A,i,j),0);}
    }
}
int JacobiFormBasis::myDiagonalForm(fmpz_mat_t A){
    int m=fmpz_mat_nrows(A), n=fmpz_mat_ncols(A);
    int i,j;
    for(i=0;i<m;i++)for(j=0;j<n;j++){
        if(i!=j)if(!fmpz_is_zero(fmpz_mat_entry(A,i,j))){
            return 0;
        }
    }
    return 1;
}

void JacobiFormBasis::mySmithNormalForm(fmpz_mat_t U, fmpz_mat_t D, fmpz_mat_t V, fmpz_mat_t M){
    //returns UM=DV
    //Assumes U,A,V are allocated correctly.
    int m=fmpz_mat_nrows(M);
    int n=fmpz_mat_ncols(M);
    if((m!=fmpz_mat_nrows(U))||(m!=fmpz_mat_ncols(U))||
       (m!=fmpz_mat_nrows(D))||(n!=fmpz_mat_ncols(D))||
       (n!=fmpz_mat_nrows(V))||(n!=fmpz_mat_ncols(V))){
        cout<<"ERROR 4t8sgdg Incorrectly allocated matrix inputs to mySmithNormalForm. Abort.\n";exit(1);
    }
    mySetIdentityMatrix(U);
    mySetIdentityMatrix(V);
    fmpz_mat_t tmpU, tmpV, A, B, At, Bt;
    fmpz_mat_init(tmpU,m,m);
    fmpz_mat_init(tmpV,n,n);
    fmpz_mat_init(A,m,n);
    fmpz_mat_init(B,n,m);
    fmpz_mat_init(At,n,m);
    fmpz_mat_init(Bt,m,n);
    int i;
    fmpz_mat_set(Bt,M);
    cout<<"Do Smith normal form on "<<m<<" x "<<n<<" matrix. This could take a while.\n";
    cout<<"Using inefficient algorithm of doing Hermite normal form repeatedly.\n";
    for(i=0;i<m+n;i++){//The theory is that it should take max(m,n) steps at most.
        cout<<"HNF #"<<i<<".\n";
        fmpz_mat_hnf_transform(A,tmpU,Bt);
        fmpz_mat_mul(U,tmpU,U);
        fmpz_mat_set(D,A);
        if(myDiagonalForm(D))break;
        fmpz_mat_transpose(At,A);
        cout<<"HNF #"<<i<<" on transpose.\n";
        fmpz_mat_hnf_transform(B,tmpV,At);
        fmpz_mat_transpose(Bt,B);
        fmpz_mat_transpose(tmpV,tmpV);
        myInvert(tmpV);
        fmpz_mat_mul(V,tmpV,V);
        fmpz_mat_set(D,Bt);
        if(myDiagonalForm(D))break;
    }
    if(!myDiagonalForm(D)){
        cout<<"ERROR kg;6vsgs mySmithNormForm. D is not SNF. Abort.\n";exit(1);
    }
    //double check
    cout<<"SNF found. Now double checking...";
    fmpz_mat_mul(A, D, V);
    fmpz_mat_mul(Bt,U,M);
    if(!fmpz_mat_equal(A,Bt)){
        cout<<"ERROR f8h=j=fd mySmithNormForm. UM=DV failed. Abort.\n";exit(1);
    }
    cout<<"checked.\n";
    fmpz_mat_clear(A);fmpz_mat_clear(B);fmpz_mat_clear(tmpU);fmpz_mat_clear(tmpV);
    fmpz_mat_clear(At);fmpz_mat_clear(Bt);
}
int JacobiFormBasis::getMinExponent(){
    if(num==0)return 0;
    int i, minqMinExp;
    minqMinExp=jfArray[0]->getMinExponent();
    for(i=1;i<num;i++){
        minqMinExp=min(minqMinExp,jfArray[i]->getMinExponent());
    }
    return minqMinExp;
}
void JacobiFormBasis::makeMaxrList(){
    nList.clear(); maxrList.clear(); nsingList.clear(); rsingList.clear();
    if(num==0)return;
    int n, maxr, i, minqMinExp,r,beginr;
    minqMinExp=getMinExponent();
    //cout<<"HERE BB: maxDegree="<<maxDegree<<"\n";
    /*** Maybe for singular part go up to index/4 only???? ***/
    int singUptoN=index/4;
    //if(maxDegree<singUptoN){cout<<"Program is set to abort if maxDegree<floor(N/4) ("<<maxDegree<<","<<singUptoN<<"). Abort.\n";
    //    exit(1);}
    for(n=minqMinExp;n<=maxDegree;n++){
        nList.push_back(n);
        maxr=jfArray[0]->getCoeffMaxr(n);
        for(i=1;i<num;i++){
            maxr=max(maxr,jfArray[i]->getCoeffMaxr(n));
        }
        //cout<<"HERE BB 2: (n,maxr,index,4n*index-maxr^2)=("<<n<<","<<maxr<<","<<index<<","<<4*n*index-maxr*maxr<<")\n";
        maxrList.push_back(maxr);
        //if(n<=singUptoN)   //Note needed any more?
        {
            if(n<0){beginr=0;}
            //else{beginr=(QZSeriesWH::mysqrtfloor(4*index*n))+1;}
            else{beginr=(QZSeriesWH::mysqrtfloor(4*index*n));}
            for(r=beginr;r<=maxr;r++){
                if(((!Diagnostics::useNegativeSingularOnly)&&
                   (4*index*n-r*r<0))||
                   (4*index*n-r*r<=0)){
                        nsingList.push_back(n); rsingList.push_back(r);
                   }
            }
        }
    }
    maxrMade=1;
}

JacobiFormBasis* JacobiFormBasis::getWeightZeroIntegerBasis(){
    int i, j,n,r;
    int nsingNum, nsingCtr;

    if(maxDegree<(index/4)){
        cout<<"Insufficient maxDegree to guarantee weight 0 integer basis. (index, index/4, maxDegree)=("
            <<index<<","<<index/4<<","<<maxDegree<<").\n";
        exit(1);
    }
    JacobiFormBasis* ans = new JacobiFormBasis();
    if(num==0){
        ans->num=0; ans->weight=weight; ans->index=index; ans->jfArray=NULL;
        return ans;
    }
    //cout<<"HERE U: "<<maxrMade<<"\n";
    cout<<"maxDegree="<<maxDegree<<"\n";
    //for(i=0;i<num;i++){
    //    cout<<"b"<<i<<"="<<jfArray[i]->getString()<<";\n";
    //}
    if(!maxrMade){makeMaxrList();}
    //if(1){makeMaxrList();}
    nsingNum=0;
    for(j=0;j<nsingList.size();j++){
        n=nsingList[j]; r=rsingList[j];
        //NEW CODE HERE March 4, 2017 to only use those with n<=index/4
        if(n<=(index/4)){\
            nsingNum++;
        }
    }
    //cout<<"HERE U2 nsingList size: "<<nsingList.size()<<", use only "<<nsingNum<<"\n";
    fmpz_mat_t M, U, D, V;
    fmpz_mat_init(M, num, nsingNum);
    fmpz_mat_init(D, num, nsingNum);
    fmpz_mat_init(U, num, num);
    fmpz_mat_init(V, nsingNum, nsingNum);
    fmpz_t c; fmpz_init(c);
    nsingCtr=0;
    for(j=0;j<nsingList.size();j++){
        n=nsingList[j]; r=rsingList[j];
        //NEW CODE HERE DEC 29, 2016 to only use those with n<=index/4

        if(n<=(index/4)){
            //cout<<"Used (n,r)=("<<n<<","<<r<<")\n";
            for(i=0;i<num;i++){
                jfArray[i]->getFC(c, n, r, index);
                fmpz_set(fmpz_mat_entry(M,i,nsingCtr),c);
            }
            nsingCtr++;
        }else{
            //cout<<"Not used (n,r)=("<<n<<",*)\n";
        }
    }
    //cout<<"M = ";fmpz_mat_print_pretty(M);cout<<endl;
    //saveMatToFile(Diagnostics::f, M, "MM");
//printMat(M,"M");
    mySmithNormalForm(U,D,V,M);
//printMat(U,"U");
//printMat(D,"D");
//printMat(V,"V");


    //cout<<"U = ";fmpz_mat_print_pretty(U);cout<<endl;
    //cout<<"D = ";fmpz_mat_print_pretty(D);cout<<endl;
    //cout<<"V = ";fmpz_mat_print_pretty(V);cout<<endl;
    //saveMatToFile(Diagnostics::f, U, "UU");
    //saveMatToFile(Diagnostics::f, D, "DD");
    //saveMatToFile(Diagnostics::f, V, "VV");

    if(fmpz_mat_rank(D)!=num){
        cout<<"CAUTION: alleged basis is not full rank ("<<
        fmpz_mat_rank(D)<<","<<num<<")in getWeightZeroIntegerBasis. Aborting to be safe.\n";exit(1);
    }
    ans->num=num;
    ans->maxDegree=maxDegree;
    ans->weight=weight;
    ans->index=index;
    ans->jfArray = new QZSeriesWH*[num];
    cout<<"Making linear combos (of "<<num<<"): ";
    for(i=0;i<num;i++){
        cout<<i<<",";
        ans->jfArray[i] = makeLinearCombo(U, i, fmpz_mat_entry(D,i,i));

        //stringstream ss; string str;
        //ss<<"allegedB["<<i<<"]";
        //ss>>str;
        //ans->jfArray[i]->saveToFile(Diagnostics::f, str);

       // jfArray[i]->getFC(c,4,17,index);  //Yhese were only used for diagnostics?
       // ans->jfArray[i]->getFC(c,4,17,index);
    }
    cout<<"done\n";
    fmpz_clear(c);
    fmpz_mat_clear(M); fmpz_mat_clear(U); fmpz_mat_clear(D); fmpz_mat_clear(V);
    return ans;
}
QZSeriesWH* JacobiFormBasis::makeLinearCombo(fmpz_mat_t M, int row,  fmpz_t denom){
    QZSeriesWH* ans = new QZSeriesWH(maxDegree);
    int j;
    for(j=0;j<num;j++){
        //cout<<"    HERE W: "<<j<<"\n";
        //cout<<"ans="<<ans->getString()<<";\n";

        ans->addScalarMultipleWith(fmpz_mat_entry(M, row, j), jfArray[j]);
    }
    //cout<<"ans="<<ans->getString()<<"\n";
    if(ans->divideByIntegerWithCheck(denom)){
        return ans;
    }else{
        cout<<"ans="<<ans->getString()<<";\n";
        cout<<"denom="<<fmpz_get_str(NULL,10,denom)<<"\n";
        //cout<<"ans="<<ans->getString()<<"\n";
        for(j=0;j<nsingList.size();j++){cout<<"("<<nsingList[j]<<","<<rsingList[j]<<")\n";}
        cout<<"num="<<num<<"\n";
        cout<<"row="<<row<<"\n";
        cout<<"M_rows="<<fmpz_mat_nrows(M)<<"\n";
        cout<<"M_cols="<<fmpz_mat_ncols(M)<<"\n";
        fmpz_mat_print(M);cout<<"\n";
        cout<<"ERROR xxxgf88glgf in makeLinearCombo: not evenly divisible by integer.\n";exit(1);
    }
}
QZSeriesWH* JacobiFormBasis::makeLinearCombo(fmpz_mat_t M, int row){
    QZSeriesWH* ans = new QZSeriesWH(maxDegree);
    int j;
    for(j=0;j<num;j++){
        //cout<<"    HERE W: "<<j<<"\n";

        ans->addScalarMultipleWith(fmpz_mat_entry(M, row, j), jfArray[j]);
    }
    return ans;
}
QZSeriesWH* JacobiFormBasis::makeLinearComboColumn(fmpz_mat_t M, int col){
    QZSeriesWH* ans = new QZSeriesWH(maxDegree);
    int j;
    for(j=0;j<num;j++){
        //cout<<"    HERE W: "<<j<<"\n";

        ans->addScalarMultipleWith(fmpz_mat_entry(M, j, col), jfArray[j]);
    }
    return ans;
}
QZSeriesWH* JacobiFormBasis::makeLinearCombo(fmpz_t* combo, fmpz_t denom){
    QZSeriesWH* ans = new QZSeriesWH(maxDegree);
    int j;
    for(j=0;j<num;j++){
        //cout<<"    HERE W: "<<j<<"\n";

        ans->addScalarMultipleWith(combo[j], jfArray[j]);
    }
    //cout<<"ans="<<ans->getString()<<"\n";
    if(ans->divideByIntegerWithCheck(denom)){
        return ans;
    }else{
        cout<<"denom="<<fmpz_get_str(NULL,10,denom)<<"\n";
        for(j=0;j<num;j++){cout<<fmpz_get_str(NULL,10,combo[j])<<",";}
        cout<<"\nERROR zzzgf88glgf in makeLinearCombo: not evenly divisible by integer.\n";exit(1);
    }
}

QZSeriesWH* JacobiFormBasis::makeLinearCombo(fmpz_t* combo){
    QZSeriesWH* ans = new QZSeriesWH(maxDegree);
    int j;
    for(j=0;j<num;j++){
        //cout<<"    HERE W: "<<j<<"\n";

        ans->addScalarMultipleWith(combo[j], jfArray[j]);
    }
    return ans;
}

void JacobiFormBasis::saveToFile(ofstream &f, string varname){
    f<<varname<<"={";
    int i;
    for(i=0;i<num;i++){
        if(i>0){f<<",\n";}
        f<< jfArray[i]->getString();
    }
    f<<"\n};\n";
    f.flush();
}

void JacobiFormBasis::saveMatToFile(ofstream &f, fmpz_mat_t m, string varname){
    f<<varname<<"={";
    int i,j;
    for(i=0;i<fmpz_mat_nrows(m);i++){
        if(i>0){f<<",\n";}
        f<<"{";
        for(j=0;j<fmpz_mat_ncols(m);j++){
            if(j>0){f<<",\n";}
            f<< fmpz_get_str(NULL,10,fmpz_mat_entry(m,i,j));
        }
        f<<"}";
    }
    f<<"\n};\n";
    f.flush();
}
void JacobiFormBasis::printMat(fmpz_mat_t m, string varname){
    cout<<varname<<"={";
    int i,j;
    for(i=0;i<fmpz_mat_nrows(m);i++){
        if(i>0){cout<<",\n";}
        cout<<"{";
        for(j=0;j<fmpz_mat_ncols(m);j++){
            if(j>0){cout<<",\n";}
            cout<< fmpz_get_str(NULL,10,fmpz_mat_entry(m,i,j));
        }
        cout<<"}";
    }
    cout<<"\n};\n";
}
JacobiFormBasis* JacobiFormBasis::getNullSpace(int targetqExponent){
    return getNullSpace(targetqExponent, targetqExponent);
}
JacobiFormBasis* JacobiFormBasis::getNullSpaceUpto(int targetqExponent){
    return getNullSpace(getMinExponent(), targetqExponent);
}
JacobiFormBasis* JacobiFormBasis::getNullSpaceAll(){
    return getNullSpace(getMinExponent(), maxDegree);
}
int JacobiFormBasis::getMaxr(int n){
    int i;
    if(n<getMinExponent()){
        return -1;
    }
    for(i=0;i<nList.size();i++){
        if(nList[i]==n)return maxrList[i];
    }
    cout<<"ERROR fdh5dfjg getMaxr failed for requested n="<<n<<", (min,max)=("<<
    getMinExponent()<<","<<maxDegree<<"). Abort.\n";
    exit(1);
}
JacobiFormBasis* JacobiFormBasis::getNullSpace(int targetqExponentFrom, int targetqExponentTo){
    if(!maxrMade){makeMaxrList();}
    int numCoeffs=0;
    int i,j,n,r,ctr;
    for(n=targetqExponentFrom;n<=targetqExponentTo;n++){
        numCoeffs+=(getMaxr(n)+1);
    }
    //cout<<"YY numCoeffs="<<numCoeffs<<".\n";
    fmpz_mat_t M; fmpz_t c;
    fmpz_init(c);
    fmpz_mat_init(M, numCoeffs, num);
    ctr=0;
    for(n=targetqExponentFrom;n<=targetqExponentTo;n++){
        for(r=0;r<=getMaxr(n);r++){
            for(i=0;i<num;i++){
                jfArray[i]->getFC(c,n,r,index);
                if(ctr>=numCoeffs){cout<<"ERROR ggasd9999 ctr>=numCoeffs: "<<ctr<<">="<<numCoeffs<<". Abort.\n";exit(1);}
                fmpz_set(fmpz_mat_entry(M,ctr,i),c);
            }
            ctr++;
        }
    }
    //cout<<"YY ctr="<<ctr<<".\n";
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f,M,"nspMat");
    JacobiFormBasis* ans;
    fmpz_mat_t nsp;
    fmpz_mat_init(nsp, num,num);
    int rk;
    rk=fmpz_mat_nullspace(nsp, M);
    //cout<<"YY rk="<<rk<<".\n";
    ans = new JacobiFormBasis();
    ans->num=rk;
    ans->maxDegree=maxDegree;
    ans->weight=weight;
    ans->index=index;
    if(rk>0){
        ans->jfArray = new QZSeriesWH*[rk];
        for(i=0;i<rk;i++){
            ans->jfArray[i] = makeLinearComboColumn(nsp, i);
        }
    }
    fmpz_clear(c);
    fmpz_mat_clear(M); fmpz_mat_clear(nsp);
    return ans;
}

fmpz_t* JacobiFormBasis::myLinearSolve(fmpz_mat_t M, fmpz_t* b, fmpz_t denom){
    int nrows, ncols, i,j, rk;
    nrows=fmpz_mat_nrows(M); ncols=fmpz_mat_ncols(M);
    cout<<"(nrows,ncols)=("<<nrows<<","<<ncols<<")\n";
    fmpz_mat_t augMat; //Augmented matrix
    fmpz_mat_init(augMat, nrows, ncols+1);
    //make augMat and row reduce
    for(i=0;i<nrows;i++)for(j=0;j<ncols;j++){
        fmpz_set(fmpz_mat_entry(augMat,i,j),fmpz_mat_entry(M,i,j));
    }
    for(i=0;i<nrows;i++){
        fmpz_set(fmpz_mat_entry(augMat,i,ncols),b[i]);
    }
    fmpz_t tmpdenom;
    fmpz_init(tmpdenom);
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f,augMat,"augMat");
    rk=fmpz_mat_rref(augMat,tmpdenom,augMat);  //tmpdenom is irrelevant for solving a solution
    cout<<"Rank="<<rk<<".\n";
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f,augMat,"augMatrref");
    if(Diagnostics::saveDiagnosticFile){Diagnostics::f<<"tmpdenom="<<fmpz_get_str(NULL,10,tmpdenom)<<";\n";}
    //check consistency
    int isZeroRow=1;
    for(j=0;j<ncols;j++){
        if(!fmpz_is_zero(fmpz_mat_entry(augMat,rk-1,j))){
            isZeroRow=0; break;
        }
    }
    if(isZeroRow){ //Then system is inconsistent!
        if(fmpz_is_zero(fmpz_mat_entry(augMat,rk-1,ncols))){
            cout<<"ERROR hgjjgjkqq: row rk should not be all zero in myLinearSolve. Abort.\n";exit(1);
        }
        fmpz_mat_clear(augMat); fmpz_clear(tmpdenom);
        return NULL;
    }
    //There is a solution!!
    fmpz_t* ans = new fmpz_t[ncols];
    fmpz_t* dens = new fmpz_t[ncols];
    fmpz_t gcd; fmpz_init(gcd);
    for(i=0;i<ncols;i++){fmpz_init(ans[i]); fmpz_init(dens[i]);}
    i=0;j=0;
    while(j<ncols){
        //cout<<"(i,j)=("<<i<<","<<j<<")\n";
        if(i>=rk){ //only free parameters left
            fmpz_set_si(ans[j],0);
            fmpz_set_si(dens[j],1);
            j=j+1;
        }else if(!fmpz_is_zero(fmpz_mat_entry(augMat,i,j))){
            //jth column is a pivot at ith row
            fmpz_set(ans[j], fmpz_mat_entry(augMat,i,ncols));
            fmpz_set(dens[j], fmpz_mat_entry(augMat,i,j));
            fmpz_gcd(gcd, ans[j], dens[j]);
            fmpz_divexact(ans[j],ans[j],gcd); fmpz_divexact(dens[j],dens[j],gcd);
            i=i+1; j=j+1;
        }else{ //zero entry found.
            fmpz_set_si(ans[j],0);
            fmpz_set_si(dens[j],1);
            j=j+1;
        }
    }
    if(i!=rk){
        cout<<"ERROR gf4t434: Something wrong with theory in algorithm of myLinearSolve. Abort.\n";exit(1);
    }
    //cout<<"AAC 6\n";

    //Get common denominator of dens.
    fmpz_set_si(denom,1);
    for(j=0;j<ncols;j++){fmpz_lcm(denom,denom,dens[j]);}
    for(j=0;j<ncols;j++){
        fmpz_mul(ans[j],ans[j],denom);
        fmpz_divexact(ans[j],ans[j],dens[j]);
    }
    fmpz_mat_clear(augMat); fmpz_clear(tmpdenom);fmpz_clear(gcd);
    for(j=0;j<ncols;j++){fmpz_clear(dens[j]);} delete[] dens;
    return ans;
}

fmpz_t* JacobiFormBasis::myLinearSolveOverIntegersOLD(fmpz_mat_t M, fmpz_t* b){
    int nrows, ncols, i,j, rk;
    nrows=fmpz_mat_nrows(M); ncols=fmpz_mat_ncols(M);
    //cout<<"(nrows,ncols)=("<<nrows<<","<<ncols<<")\n";
    fmpz_t* ans;
    fmpz_mat_t U, Mtr, H, Utr, Htr, y, x;
    fmpz_mat_init(Mtr, ncols, nrows);
    fmpz_mat_init(H, ncols, nrows);
    fmpz_mat_init(Htr, nrows, ncols);
    fmpz_mat_init(U, ncols, ncols);
    fmpz_mat_init(Utr, ncols, ncols);
    fmpz_mat_init(y, ncols, 1);
    fmpz_mat_init(x, ncols, 1);

    fmpz_mat_transpose(Mtr, M);
    fmpz_mat_hnf_transform(H, U, Mtr);
    fmpz_mat_transpose(Htr, H);

    fmpz_t denom;    fmpz_init(denom);
    fmpz_t* sol = myLinearSolve(Htr, b, denom);
    if(sol==NULL){
        ans=NULL;
    }else if(!fmpz_is_one(denom)){
        ////We need to do Smith Normal Form in this case to be sure there are no integer solutions!
        ////Realized this 4-1-2017 !!
        cout<<"ERROR 523842395: Row reduction yielded a noninteger rational solution (denom="<<fmpz_get_si(denom)
            <<"). Need to program Smith normal form! Abort.\n";
        exit(1);
        ans=NULL;
    }else{
        for(j=0;j<ncols;j++){
            fmpz_set(fmpz_mat_entry(y,j,0),sol[j]);
            //cout<<"sol[j]="<<fmpz_get_str(NULL,10,sol[j])<<"\n";
        }
        fmpz_mat_transpose(Utr,U);
        fmpz_mat_mul(x,Utr,y);
        ans=new fmpz_t[ncols];
        for(j=0;j<ncols;j++){
            fmpz_init(ans[j]);
            fmpz_set(ans[j],fmpz_mat_entry(x,j,0));
        }
        if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, x, "xsol");
    }

    //Clean up
    if(sol!=NULL){
        for(j=0;j<ncols;j++){fmpz_clear(sol[j]);}
        delete[] sol;
    }
    fmpz_clear(denom);
    fmpz_mat_clear(Mtr);
    fmpz_mat_clear(H);
    fmpz_mat_clear(Htr);
    fmpz_mat_clear(U);
    fmpz_mat_clear(Utr);
    fmpz_mat_clear(y);
    fmpz_mat_clear(x);
    return ans;
}

QZSeriesWH* JacobiFormBasis::getOneSolution(LaurentZ* targetLaurentPoly, int targetqExponent){
    //Returns NULL if no solution

    if(!maxrMade){makeMaxrList();}
    int maxr = max(getMaxr(targetqExponent), targetLaurentPoly->getMaxExponent());
    //cout<<"AAB maxr="<<maxr<<"\n";

    int numCoeffs=maxr+1;
    int i,j,n,r,ctr;
    fmpz_mat_t M; fmpz_t c;
    fmpz_init(c);
    fmpz_mat_init(M, numCoeffs, num);
    fmpz_t* b = new fmpz_t[numCoeffs];
    for(i=0;i<numCoeffs;i++){fmpz_init(b[i]);}
    ctr=0;
    n=targetqExponent;
    for(r=0;r<=maxr;r++){
            targetLaurentPoly->getCoeff(c, r);
            fmpz_set(b[r],c);
            for(i=0;i<num;i++){
                jfArray[i]->getFC(c,n,r,index);
                fmpz_set(fmpz_mat_entry(M,ctr,i),c);
            }
            ctr++;
    }
    //cout<<"YY ctr="<<ctr<<".\n";
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f,M,"getOneSolutionM");

    //cout<<"AAB 3\n";
    fmpz_t* combo;

    //combo = myLinearSolveOverIntegersOLD(M, b);
    combo = myLinearSolveOverIntegers(M, b);

    QZSeriesWH* ans;
    if(combo==NULL){
        ans=NULL;
    }else{
        ans =makeLinearCombo(combo);
        //check answer?
        LaurentZ* tCopy=ans->getCoeffCopy(targetqExponent);
        if(targetLaurentPoly->equals(tCopy)){delete tCopy;}
        else{cout<<"ERROR hghddaaa: Solution check at getOneSolution failed:\n"<<
            targetLaurentPoly->getString()<<"\nversus\n"<<
            tCopy->getString()<<".\nAbort.\n";
            cout<<"whole ans = "<<ans->getString()<<"\n";
            exit(1);
        }
        for(i=0;i<num;i++){fmpz_clear(combo[i]);} delete[] combo;
    }
    fmpz_clear(c);
    fmpz_mat_clear(M);
    for(i=0;i<numCoeffs;i++){fmpz_clear(b[i]);} delete[] b;
    return ans;
}
int JacobiFormBasis::checkSymmetry(){
    int i;
    for(i=0;i<num;i++){
        if(!(jfArray[i]->checkSymmetry())){
            return 0;
        }
    }
    return 1;
}

void JacobiFormBasis::putIntoDRList(DRList* dr, int uptoN){
    int i;
    for(i=0;i<num;i++){
        jfArray[i]->putIntoDRList(dr,index,uptoN);
    }
}

void JacobiFormBasis::makeHumbertMatrix(fmpz_mat_t M, DRList* dr){
    int j;
    for(j=0;j<num;j++){
        jfArray[j]->makeHumbertVector(M, dr, index,j);
    }
}
void JacobiFormBasis::makeSingularMatrix(fmpz_mat_t M, DRList* dr){
    int j;
    for(j=0;j<num;j++){
        jfArray[j]->makeSingularVector(M, dr, index,j);
    }
}


int JacobiFormBasis::allNonnegative(fmpz_mat_t M){
    int i, j;
    for(i=0;i<fmpz_mat_nrows(M);i++)for(j=0;j<fmpz_mat_ncols(M);j++){
        if(fmpz_cmp_si(fmpz_mat_entry(M,i,j),0)<0){return 0;}
    }
    return 1;
}

JacobiFormBasis* JacobiFormBasis::getGoodJFs(fmpz_mat_t HumMatrix,fmpz_mat_t HumVec, QZSeriesWH* combinedSol){
    JacobiFormBasis* ans = new JacobiFormBasis();
    ans->weight=weight; ans->index=index; ans->maxDegree=maxDegree;
    if(num==0){
        if(allNonnegative(HumVec)){
            ans->num=1;
            ans->jfArray = new QZSeriesWH*[1];
            ans->jfArray[0] = combinedSol->copy();
            ans->maxDegree = combinedSol->getMaxExponent();
        }
        return ans;
    }
    if(Diagnostics::skipMathematica){
        ans->num=0;
        return ans;
    }

    AXBsolver* axb = new AXBsolver(HumMatrix,HumVec, fmpz_mat_nrows(HumMatrix), num);
    axb->solve();
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f,axb->solM,"solM");
    cout<<"solM saved.\n";
    ans->num = axb->numSolutions;
    if(ans->num>0){
        ans->jfArray = new QZSeriesWH*[ans->num];
        int i;
        for(i=0;i<ans->num;i++){
            ans->jfArray[i] = makeLinearCombo(axb->solM,i);
            ans->jfArray[i]->addWith(combinedSol);
            //double check nonnegativity of Humbert multiplicities? Nah.
        }
    }
    delete axb;
    return ans;
}

JacobiFormBasis* JacobiFormBasis::subspaceMultipleOf(LaurentZ* pz){
    int numExponents = maxDegree-getMinExponent()+1;
    pz->normalize();
    int deg = fmpz_poly_degree(pz->p);
    fmpz_mat_t M;
    fmpz_mat_init(M, numExponents*deg, num);
    fmpz_poly_t rem; fmpz_poly_init(rem);
    LaurentZ* coeff;
    fmpz_t c; fmpz_init(c);
    int i, n,r, ctr;
    for(i=0;i<num;i++){
        jfArray[i]->truncate(maxDegree);
        ctr=0;
        for(n=getMinExponent();n<=maxDegree;n++){
            coeff=jfArray[i]->getCoeffCopy(n);
            coeff->setRemainder(rem, pz);
            delete coeff;
            for(r=0;r<deg;r++){
                fmpz_poly_get_coeff_fmpz(c, rem, r);
                if(ctr>=fmpz_mat_nrows(M)){cout<<"ERROR 845dfagaa ctr>=numCoeffs: "<<ctr<<">="<<fmpz_mat_nrows(M)<<". Abort.\n";exit(1);}
                fmpz_set(fmpz_mat_entry(M,ctr,i),c);
                ctr++;
            }
        }
    }
    cout<<"(nfrom, nto, deg, ctr, nrows, ncols)=("<<getMinExponent()<<","<<maxDegree<<","
        <<deg<<","<<ctr<<","<<fmpz_mat_nrows(M)<<","<<fmpz_mat_ncols(M)<<")\n";
    //saveMatToFile(Diagnostics::f,M,"nspMat");
    JacobiFormBasis* ans;
    fmpz_mat_t nsp;
    fmpz_mat_init(nsp, num,num);
    int rk;
    rk=fmpz_mat_nullspace(nsp, M);
    fmpz_mat_t chMat;
    fmpz_mat_init(chMat, fmpz_mat_nrows(M), num);
    fmpz_mat_mul(chMat, M, nsp);
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, chMat, "chMat");
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, M, "M");
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, nsp, "nsp");
    if(Diagnostics::saveDiagnosticFile)saveToFile(Diagnostics::f, "bbb");
    //cout<<"pz="<<pz->getString()<<"\n";

    ans = new JacobiFormBasis();
    ans->num=rk;
    ans->maxDegree=maxDegree;
    ans->weight=weight;
    ans->index=index;
    int succ; QZSeriesWH* ch;
    if(rk>0){
        ans->jfArray = new QZSeriesWH*[rk];
        for(i=0;i<rk;i++){
        //for(i=rk-1;i>=0;i--){
            //cout<<"ZZZ0 : "<<i<<"\n";
            ans->jfArray[i] = makeLinearComboColumn(nsp, i);
            //CHECK
            ch = ans->jfArray[i]->copy();
            succ=ch->divideByLaurentZ(pz);
            //cout<<"ZZZZZ ["<<i<<"] succ="<<succ<<", and maxExponent="<<ans->jfArray[i]->getMaxExponent()<<".\n";
            delete ch;
            //cout<<"ch deleted.\n";
        }

    }
    fmpz_clear(c);
    fmpz_mat_clear(M); fmpz_mat_clear(nsp); fmpz_poly_clear(rem); fmpz_mat_clear(chMat);
    return ans;
}

fmpz_t* JacobiFormBasis::matchSingularPart(QZSeriesWH* jfweak, fmpz_t denom){
    DRList *dr = new DRList();
    if(num==0)return 0;
    int uptoN;
    if(weight==0){
        uptoN=index/4;
        cout<<"In matchSingularPart, weight=0, so using uptoN=index/4="<<uptoN<<"\n";
    }else{
        uptoN=jfweak->getMaxExponent();
        cout<<"In matchSingularPart, weight="<<weight<<", so using uptoN=maxExponent="<<uptoN<<"\n";
    }
    putIntoDRList(dr, uptoN);
    jfweak->putIntoDRList(dr, index, uptoN);
    fmpz_mat_t singMatrix;
    fmpz_mat_t singVec;
    fmpz_mat_init(singVec, dr->getLength(),1);
    cout<<"Make singular vector and matrix.\n";
    jfweak->makeSingularVector(singVec, dr, index);
    fmpz_mat_init(singMatrix, dr->getLength(), num);
    makeSingularMatrix(singMatrix, dr);
    fmpz_t*ans;
    ans = myLinearSolve(singMatrix, singVec, denom);
    fmpz_mat_clear(singMatrix);
    fmpz_mat_clear(singVec);
    delete dr;
    return ans;
}
fmpz_t* JacobiFormBasis::myLinearSolve(fmpz_mat_t M, fmpz_mat_t b, fmpz_t denom){
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, M, "LSM");
    if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, b, "LSb");

    int n=fmpz_mat_nrows(b);
    fmpz_t* tmp = new fmpz_t[n];
    int i;
    for(i=0;i<n;i++){
        fmpz_init(tmp[i]);
        fmpz_set(tmp[i],fmpz_mat_entry(b,i,0));
    }
    fmpz_t* ans=myLinearSolve(M, tmp, denom);
    for(i=0;i<n;i++){
        fmpz_clear(tmp[i]);
    }
    delete[] tmp;
    return ans;
}

void JacobiFormBasis::weight0AppendMorejfv2overjf(ThetaBlockList* tbList){
    int tbi, i, success;
    QZSeriesWH** newjfArray = new QZSeriesWH*[num+tbList->num];
    vector<int>* newtbArray = new vector<int>[num+tbList->num];
    if(num>0){
        for(i=0;i<num;i++){
            newjfArray[i]=jfArray[i];
            newtbArray[i]=tbArray[i];
        }
    }
    QZSeriesWH* jftb, *jftbV2, *jftbV2overjftb;
    vector<int> d;
    cout<<"Using theta blocks of weight "<<tbList->weight<<", index "<<tbList->index<<".\n";
    cout<<"Reading "<<tbList->num<<":";
    for(tbi=0;tbi<tbList->num;tbi++){
        d=tbList->d[tbi];
        cout<<" "<<tbi;
        jftb = QZSeriesWH::ThetaBlock(maxDegree*2,tbList->weight,d);
        jftbV2 = jftb->copy();
        jftbV2->applyUp(tbList->index,tbList->weight,2);
        jftbV2overjftb = jftbV2->copy();
        success=jftbV2overjftb->divideBy(jftb);  //Later check if this was possible at all.
        if(success){
            newjfArray[num]=jftbV2overjftb;
            newtbArray[num]=d;
            num++;
        }else{
            cout<<"jftbV2overjftb not added because jf|V2 not divisible by jf (must be a TBWD).\n";
        }
        delete jftb; delete jftbV2;
    }
    cout<<"\n";
    delete[] jfArray; jfArray=newjfArray;
    delete[] tbArray; tbArray=newtbArray;
}

fmpz_t* JacobiFormBasis::myLinearSolveOverIntegers(fmpz_mat_t M, fmpz_t* b){
    int nrows, ncols, i,j, rk, hasSol;
    nrows=fmpz_mat_nrows(M); ncols=fmpz_mat_ncols(M);
    //cout<<"(nrows,ncols)=("<<nrows<<","<<ncols<<")\n";
    fmpz_t* ans; ans=NULL;
    fmpz_mat_t U, D, V, bMat, c, x, y;
    fmpz_mat_init(U, nrows, nrows);
    fmpz_mat_init(D, nrows, ncols);
    fmpz_mat_init(V, ncols, ncols);
    fmpz_mat_init(bMat, nrows, 1);
    fmpz_mat_init(c, nrows, 1);
    fmpz_mat_init(y, ncols, 1);
    fmpz_mat_init(x, ncols, 1);

    fmpz_t denom;    fmpz_init(denom);
        //cout<<"Dimensions: ";
        //cout<<"M ("<<fmpz_mat_nrows(M)<<","<<fmpz_mat_ncols(M)<<"), ";
        //cout<<"U ("<<fmpz_mat_nrows(U)<<","<<fmpz_mat_ncols(U)<<"), ";
        //cout<<"D ("<<fmpz_mat_nrows(D)<<","<<fmpz_mat_ncols(D)<<"), ";
        //cout<<"V ("<<fmpz_mat_nrows(V)<<","<<fmpz_mat_ncols(V)<<"), ";
        //cout<<"x ("<<fmpz_mat_nrows(x)<<","<<fmpz_mat_ncols(x)<<"), ";
        //cout<<"y ("<<fmpz_mat_nrows(y)<<","<<fmpz_mat_ncols(y)<<"), ";
        //cout<<"c ("<<fmpz_mat_nrows(c)<<","<<fmpz_mat_ncols(c)<<"), ";
        //cout<<"bMat ("<<fmpz_mat_nrows(bMat)<<","<<fmpz_mat_ncols(bMat)<<").\n";

    mySmithNormalForm(U, D, V, M); //UM=DV
    for(i=0;i<nrows;i++){
        fmpz_set(fmpz_mat_entry(bMat,i,0),b[i]);
    }
    fmpz_mat_mul(c, U,bMat);
    hasSol=1;
    if(nrows<=ncols){
        for(i=0;i<nrows;i++){
            //cout<<"Testing divisibility ("<<fmpz_get_si(fmpz_mat_entry(c,i,0))
            //    <<","<<fmpz_get_si(fmpz_mat_entry(D,i,i))<<").\n";
            if(!myAttemptDivide(fmpz_mat_entry(y,i,0),fmpz_mat_entry(c,i,0),fmpz_mat_entry(D,i,i))){
                hasSol=0;break;
            }
        }
        if(hasSol){
            for(i=nrows;i<ncols;i++)fmpz_set_si(fmpz_mat_entry(y,i,0),0);
        }
    }else{
        for(i=0;i<ncols;i++){
            if(!myAttemptDivide(fmpz_mat_entry(y,i,0),fmpz_mat_entry(c,i,0),fmpz_mat_entry(D,i,i))){
                hasSol=0;break;
            }
        }
        if(hasSol){
            for(i=ncols;i<nrows;i++)if(!fmpz_is_zero(fmpz_mat_entry(c,i,0))){
                hasSol=0;break;
            }
        }
    }
    if(hasSol){
        fmpz_mat_solve(x, denom, V, y);
        if(!fmpz_is_pm1(denom)){
            cout<<"ERROR hadfga268889 in JacobiFormBasis.cpp: Something is wrong denom ("<<fmpz_get_si(denom)<<") should be pm1. Abort.\n";
            exit(1);
            ans=NULL;
        }else if(fmpz_equal_si(denom,-1)){
            fmpz_mat_neg(x, x);
        }
        ans=new fmpz_t[ncols];
        for(j=0;j<ncols;j++){
            fmpz_init(ans[j]);
            fmpz_set(ans[j],fmpz_mat_entry(x,j,0));
        }
        if(Diagnostics::saveDiagnosticFile)saveMatToFile(Diagnostics::f, x, "xsol");
        //Double check Mx=b?
        //cout<<"Double checking myLinearSolveOverIntegers.\n";
        fmpz_mat_mul(c, M, x);
        //cout<<"Done multiplying, now checking equality.\n";
        if(!fmpz_mat_equal(c,bMat)){
            cout<<"ERROR ahadfg34 in JacobiFormBasis.cpp myLinearSolveOverIntegers: double check solution failed!!!. Abort.\n";
            exit(1);
        }
    }
    //cout<<"Cleaning up.\n";
    fmpz_clear(denom);
    fmpz_mat_clear(U);
    fmpz_mat_clear(D);
    fmpz_mat_clear(V);
    fmpz_mat_clear(bMat);
    fmpz_mat_clear(c);
    fmpz_mat_clear(y);
    fmpz_mat_clear(x);
    return ans;
}

int JacobiFormBasis::myAttemptDivide(fmpz_t x, fmpz_t a, fmpz_t b){  //Solves bx=a
    if(fmpz_is_zero(b)){
        if(fmpz_is_zero(a)){
            fmpz_set_si(x,0);
            return 1;
        }else{
            return 0;
        }
    }
    int ans;
    if(fmpz_cmp_si(b,0)>0){
        ans=fmpz_divisible(a,b);
    }else{
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_neg(tmp, b);
        int ans=fmpz_divisible(a,tmp);
        fmpz_clear(tmp);
    }
    if(ans){fmpz_divexact(x, a, b);}
    return ans;
}
