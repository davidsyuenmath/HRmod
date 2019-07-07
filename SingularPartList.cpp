#include "SingularPartList.h"

SingularPartList::SingularPartList(int index0){
    index=index0;
    rowLen.clear();
    maxRowLen=0;
    matSize=100;
    mat=new fmpz_t*[matSize];
    tbArray=new vector<int>[matSize];
    dr=new DRList();
    fmatMade=0;
}
SingularPartList::~SingularPartList(){
    if(matSize>0){
        int i,j;
        for(i=0;i<rowLen.size();i++){
            for(j=0;j<rowLen[i];j++){
                fmpz_clear(mat[i][j]);
            }
            delete[] mat[i];
        }
        delete[] mat;
        delete[] tbArray;
        matSize=0;
        tbArray=NULL;
        mat=NULL;
        if(fmatMade){
            fmatMade=0;
            fmpz_mat_clear(fmat);
        }
    }
} //destructor

int SingularPartList::getNumRows(){
    return rowLen.size();
}
void SingularPartList::     increaseMatSize(){
    int newmatSize=matSize+100;
    fmpz_t** newmat = new fmpz_t*[newmatSize];
    vector<int>* newtbArray=new vector<int>[newmatSize];
    int i;
    for(i=0;i<rowLen.size();i++){
        newmat[i]=mat[i];
        newtbArray[i]=tbArray[i];
    }
    delete[] mat;
    mat=newmat; matSize=newmatSize; tbArray=newtbArray;
}

void SingularPartList::     appendRow(QZSeriesWH* jf, int index, int uptoN){
    vector<int> tb; tb.clear();
    appendRow(jf, index, uptoN, tb);
}

void SingularPartList::     appendRow(QZSeriesWH* jf, int index, int uptoN, vector<int> tb){
    //cout<<"HERE A: appending row.\n";
    jf->putIntoDRList(dr, index, uptoN);
    int drLen=dr->getLength();
    maxRowLen=drLen;
    int i, d, r;
    fmpz_t c; fmpz_init(c);
    fmpz_t *v = new fmpz_t[drLen];
    for(i=0;i<dr->getLength();i++){
        d=dr->DList[i]; r=dr->RList[i];
        //jf->getHumbertMultiplicity(c, dr->getMinD(), d, r, index);
        jf->getFC_Dr(c, d, r, index);
        fmpz_init(v[i]);
        fmpz_set(v[i],c);
        //cout<<"HERE HH: "<<fmpz_get_str(NULL,10,c)<<"\n";
    }
    fmpz_clear(c);
    //Now put v into mat
    i=rowLen.size();
    if(i>=matSize){increaseMatSize();}
    mat[i]=v;
    tbArray[i]=tb;
    rowLen.push_back(drLen);
    clearfmat();
}
void SingularPartList::checkfmat(){
    if(!fmatMade){
        makeMat(fmat, 0, getNumRows());
        fmatMade=1;
    }
}
void SingularPartList::clearfmat(){
    if(fmatMade){
        fmpz_mat_clear(fmat);
        fmatMade=0;
    }
}

void SingularPartList::     makeMat(fmpz_mat_t &m, int beginRowNumber, int numRows){
    if(beginRowNumber+numRows>getNumRows()){
        cout<<"ERROR adsfhg8e in SingularPartList: Requested makeMat size exceeds number of rows. Abort.\n";exit(1);
    }
    fmpz_mat_init(m,maxRowLen, numRows);
    int i, j,k;
    for(i=0;i<numRows;i++){
        k=i+beginRowNumber;
        for(j=0;j<maxRowLen;j++){
            if(j<rowLen[k]){
                fmpz_set(fmpz_mat_entry(m,j,i),mat[k][j]);
            }else{
                fmpz_set_si(fmpz_mat_entry(m,j,i),0);
            }
        }
    }
}
fmpz_t* SingularPartList::     getVec(int rowNumber){
    fmpz_t* v = new fmpz_t[maxRowLen];
    int i;
    for(i=0;i<maxRowLen;i++){
        fmpz_init(v[i]);
        if(i<rowLen[rowNumber]){
            fmpz_set(v[i],mat[rowNumber][i]);
        }else{
            fmpz_set_si(v[i],0);
        }
    }
    return v;
}




void SingularPartList:: saveToFile(ofstream &f, string varname){
    f<<varname<<"={";
    int i, j;
    for(i=0;i<rowLen.size();i++){
        if(i>0){f<<",\n";}
        f<<"{";
        for(j=0;j<rowLen[i];j++){
            if(j>0){f<<",";}
            f<<fmpz_get_str(NULL,10,mat[i][j]);
        }
        f<<"}";
    }
    f<<"\n};\n";
    f.flush();
}
void SingularPartList::appendRows(ThetaBlockList* tbList, int index, int uptoN){
    int tbi, i, success;
    vector<int>* tb;
    QZSeriesWH* jftb, *jftbV2, *jftbV2overjftb;
    vector<int> d;
    cout<<"Using theta blocks of weight "<<tbList->weight<<", index "<<tbList->index<<".\n";
    cout<<"Reading :-)  "<<tbList->num<<":";
    for(tbi=0;tbi<tbList->num;tbi++){
        d=tbList->d[tbi];
        cout<<" "<<tbi;
        jftb = QZSeriesWH::ThetaBlock(uptoN*2,tbList->weight,d);
        jftbV2 = jftb->copy();
        jftbV2->applyUp(tbList->index,tbList->weight,2);
        jftbV2overjftb = jftbV2->copy();
        success=jftbV2overjftb->divideBy(jftb);  //Later check if this was possible at all.
        if(success){
            appendRow(jftbV2overjftb,tbList->index, index/4 ,d); //NOTE index/4!
        }else{
            cout<<"\n jftbV2overjftb not added because jf|V2 not divisible by jf (must have been a TBWD).\n";
        }
        delete jftb; delete jftbV2; delete jftbV2overjftb;
    }
    cout<<"\n";
}

fmpz_t* SingularPartList::matchSingularPart(QZSeriesWH *jf, fmpz_t denom){
    cout<<"Attempting to matchSingularPart with combo of SingularPartList.\n";
    /** Check if jf adds anything nonzero entry to a new entry of dr?**/
    if(jf->hasNonzeroNewDR(dr, index, index/4 ))
    {
        return NULL;
    }
    /**Get sing vec from jf and attempt to solve**/
    checkfmat();
    fmpz_mat_t singVec;
    fmpz_mat_init(singVec, dr->getLength(),1);
    cout<<"Make singular vector and matrix in SingularPartList.\n";
    jf->makeSingularVector(singVec, dr, index);
    fmpz_t*ans;
    ans = JacobiFormBasis::myLinearSolve(fmat, singVec, denom);
    fmpz_mat_clear(singVec);
    return ans;
}

void SingularPartList::appendRows(ThetaBlockList* tbList, ThetaBlockList *tbList2, int index, int uptoN){
    /** TO BE PROGRAMMED**/
    int tbi, i, success;
    vector<int>* tb;
    QZSeriesWH* jftb, *jftbV2, *jftbV2overjftb;
    vector<int> d;
    cout<<"Using theta blocks of weight "<<tbList->weight<<", index "<<tbList->index<<".\n";
    cout<<"Reading :-)  "<<tbList->num<<":";
    for(tbi=0;tbi<tbList->num;tbi++){
        d=tbList->d[tbi];
        cout<<" "<<tbi;
        jftb = QZSeriesWH::ThetaBlock(uptoN*2,tbList->weight,d);
        jftbV2 = jftb->copy();
        jftbV2->applyUp(tbList->index,tbList->weight,2);
        jftbV2overjftb = jftbV2->copy();
        success=jftbV2overjftb->divideBy(jftb);  //Later check if this was possible at all.
        if(success){
            appendRow(jftbV2overjftb,tbList->index, index/4 ,d); //NOTE index/4!
        }else{
            cout<<"\n jftbV2overjftb not added because jf|V2 not divisible by jf (must have been a TBWD).\n";
        }
        delete jftb; delete jftbV2; delete jftbV2overjftb;
    }
    cout<<"\n";
}
