/****
  Finds integer solutions to A x + b >= 0
****/


#include "AXBsolver.h"
#include "time.h"


AXBsolver::AXBsolver(fmpz_mat_t A0, fmpz_t* b0, int nrows0, int ncols0){
  numSolutions=0;
  nrows=nrows0;ncols=ncols0;
  fmpz_mat_init(A, nrows, ncols);
  fmpz_mat_set(A, A0);
  b=new fmpz_t[nrows];
  for(int i=0;i<nrows;i++){
    fmpz_init(b[i]);
    fmpz_set(b[i],b0[i]);
  }
}
AXBsolver::AXBsolver(fmpz_mat_t A0, fmpz_mat_t b0, int nrows0, int ncols0){
  numSolutions=0;
  nrows=nrows0;ncols=ncols0;
  fmpz_mat_init(A, nrows, ncols);
  fmpz_mat_set(A, A0);
  b=new fmpz_t[nrows];
  for(int i=0;i<nrows;i++){
    fmpz_init(b[i]);
    fmpz_set(b[i], fmpz_mat_entry(b0,i,0));
  }
}

AXBsolver::AXBsolver(int** A0, int* b0, int nrows0, int ncols0){
  numSolutions=0;
  int i,j;
  nrows=nrows0;ncols=ncols0;
  fmpz_mat_init(A, nrows, ncols);
  b=new fmpz_t[nrows];
  for(i=0;i<nrows;i++){
    fmpz_init(b[i]);
    fmpz_set_si(b[i],b0[i]);
    for(j=0;j<ncols;j++){
      fmpz_set_si(fmpz_mat_entry(A,i,j),A0[i][j]);
    }
  }
}

AXBsolver:: ~AXBsolver(){//destructor
    destroy();
}
void AXBsolver::destroy(){
  fmpz_mat_clear(A);
  fmpz_mat_clear(solM);
  int i;
  for(i=0;i<nrows;i++){
    fmpz_clear(b[i]);
  }
}

void AXBsolver::solve(){
    //Stores solutions in long long ** sols with int numSolutions solutions.
    //Find unused temporary filename
    srand(time(NULL));
    int randInt = rand()%10000;
    stringstream ss;
    string fname;

    ss.clear();
    ss<<"temp-input-"<<randInt<<".ma";
    ss>>fname;
    while(std::ifstream(fname.c_str())){
        randInt++;
        ss.clear();
        ss<<"temp-input-"<<randInt<<".ma";
        ss>>fname;
    }
    ofstream PFILE;
    PFILE.open(fname.c_str());
    if(PFILE.fail()){cout<<"File "<<fname<<" failed to open for writing.\n"; exit( 1); }
    PFILE<<"A = {";
    int i,j;
    for(i=0;i<nrows;i++){
        if(i>0){PFILE<<",\n";}
        PFILE<<"{";
        for(j=0;j<ncols;j++){
            if(j>0){PFILE<<",";}
            PFILE<<fmpz_get_str(NULL,10,fmpz_mat_entry(A,i,j));
        }
        PFILE<<"}";
    }
    PFILE<<"\n};\n";
    PFILE<<"b = {";
    for(i=0;i<nrows;i++){
        if(i>0){PFILE<<",\n";}
        PFILE<<fmpz_get_str(NULL,10,b[i]);
    }
    PFILE<<"};\n";
    PFILE.close();
    string command="";
    command = "perl SolveAXB.pl " + fname;
    system(command.c_str());
    ss.clear();
    ss<<fname<<"-output.ma.txt";
    string infileName;
    ss>>infileName;
    FILE* INFILE;
    INFILE = fopen(infileName.c_str(), "r");
    if(INFILE==NULL){cout<<"File "<<infileName<<" failed to open for reading.\n"; exit(1); }
    int success;
    fscanf(INFILE,"%i", &success);
    if(success==0){
        cout<<"YIKES success=0 in file "<<infileName<<". Aborting to be safe.\n";exit(1);
    }
    fscanf(INFILE," %i", &numSolutions);
    int x;
    fscanf(INFILE," %i", &x);
    if(x!=ncols){
        cout<<"ERROR sfdhhh ncols does not match that in file "<<infileName<<
        " ("<<ncols<<","<<x<<"). Abort.\n";exit(1);
    }
    fmpz_mat_init(solM, numSolutions, ncols);
    fmpz_t c;
    for(i=0;i<numSolutions;i++){
        for(j=0;j<ncols;j++){
            fmpz_fread(INFILE, c);
            fmpz_set(fmpz_mat_entry(solM,i,j),c);
        }
    }
    fmpz_clear(c);
    fclose(INFILE);

    string randIntString;
    ss.clear();
    ss<<randInt;
    ss>>randIntString;

    command = "perl deleteTempFiles.pl " + randIntString;
    system(command.c_str());

}
/*****************

*****/
