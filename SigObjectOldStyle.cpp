/***********
 * NOTE: sigs are stored as a,b,c where 2p|a.
 * July 3, 2010 - put in odd weight capable Reduce with sign
 * Should put in a getUnreducedSig feature later
 * ALL odd plus or minus signs are long long.
 * Need long long PX2Reduce version
 * indices go from 1 to SIGLEN!!
 ************/

#include "SigObjectOldStyle.h"

void SigObjectOldStyle ::destroy(){
  delete[] u;
  delete[] w;
  delete[] z;
  delete[] Ninv;
  delete[] Ngcd;
  int i;
  for(i=1;i<SIGLEN+1;i++){
    delete sigs[i];
  }
  delete sigs;
  for(i=0;i<=MAXDET;i++){
    delete startend[i];
  }
  delete startend;
  for(i=0;i<vset.size();i++){
    delete[] vset[i];
  }
  delete[] divisorToIndex;
  for(i=0;i<vsetv2.size();i++){
    delete[] vsetv2[i];
  }

}

SigObjectOldStyle ::SigObjectOldStyle(int p0, int d, string filename){
  uptoDet = d;
  p = p0;
  p2 = 2*p;
  verbose=1;

  initFirst();

  filename = ""; //Just to get rid of unused variable warning.
  int i,j;
  string sigf;
  ifstream SIGFILE;
  stringstream ss;
  ss<<"siginfo-"<<p<<"-"<<uptoDet<<".ma.txt";
  ss>>sigf;
  int x;
  cout<<"Reading file: "<<sigf<<"...\n";
  SIGFILE.open(sigf.c_str());
  if(SIGFILE.fail()){cout<<"File read open failure huh ("<<sigf<<")\n";exit(1);}
  SIGFILE>>x;  //should be p
  if(x!=p){
    cout<<"First number in file does not match "<<p<<"\n";exit(1);
  }
  SIGFILE>>x;
  if(x!=uptoDet){
    cout<<"uptoDet number in file does not match "<<x<<"\n";exit(1);
  }
  SIGFILE>>SIGLEN;
  cout<<"Number of sigs to be read: "<<SIGLEN<<"\n";

  cout<<"Reading sigs from file...\n";

  sigs=new int*[SIGLEN+1];
  for(i=1;i<=SIGLEN;i++){
    sigs[i]=new int[5];
    for(j=0;j<5;j++){
      SIGFILE>>sigs[i][j];
    }
  }
  SIGFILE.close();
  cout<<"File read.  Now creating startend array.\n";
  initRest();
}

void SigObjectOldStyle ::initRest(){
  int i;

  MAXDET=quickDet3(sigs[SIGLEN]);
  startend=new int*[MAXDET+1];
  for(i=0;i<=MAXDET;i++){
    startend[i]=new int[2];
    startend[i][0]=-1;
  }
  int currdet=quickDet3(sigs[1]);
  int nextdet;
  int currStart=1;
  int currI=1;
  while(true){
    currI++;
    if(currI>SIGLEN){
      startend[currdet][0]=currStart;
      startend[currdet][1]=currI-1;
      break;
    }
    nextdet=quickDet3(sigs[currI]);
    if(nextdet!=currdet){
       startend[currdet][0]=currStart;
       startend[currdet][1]=currI-1;
       currStart=currI;
       currdet=nextdet;
    }
  }
}

void SigObjectOldStyle ::initFirst(){
  int i, j;
  u=new int[2];
  w=new int[2];
  z=new int[2];

  if(verbose)cout<<"Now making gcd and inv array.\n";

  Ninv = new int[p];  // c*Ninv[c] == Ngcd[c]
  Ngcd = new int[p];

  for(i=0;i<p;i++){
    Ngcd[i] = MyGCD(i,p)%p;
  }
  for(i=0;i<p;i++){
    Ninv[i]=0;
    if(Ngcd[i]!=0){
     for(j=1;j<p;j++)if(Ngcd[j]==1){
      if((i*j%p)==Ngcd[i]){Ninv[i]=j;break;}
     }
     if(Ninv[i]==0){cout<<"Ninv error at i="<<i<<"\n";exit(1);}
    }
  }

  if(verbose)cout<<"Now making vset.\n";
  //Make reduced vset and at the same time,
  //Make quick references
  //Use new arrays
  //int* divisorToIndex;  //Set to -1 if not divisor
  //vector<int> vsetv1;
  //vector<int*> vsetv2;
  //The above are only for divisors not 0,1.
  int *k;
  k=new int[2];
  k[0]=0;
  k[1]=1;
  vset.push_back(k);
  divisorToIndex = new int[p];
  for(i=0;i<p;i++){
    k=new int[2];
    k[0]=1;
    k[1]=i;
    vset.push_back(k);
    divisorToIndex[i]=-1; //initialize to default
  }
  int ctr=0;
  int *vList; int r, a, b;
  for(i=2;i<p;i++){
    if((p%i)==0){
      divisorToIndex[i]=ctr;
      vsetv1.push_back(i);
      vList = new int[p];
      vsetv2.push_back(vList);
      for(j=0;j<p;j++)if(MyGCD(i,j)==1){

      //cout<<"HERE B";

        r=j; //best so far.
        for(a=p/i+1; a<=p; a+=p/i)if(Ngcd[a]==1){
          if(a*i%p !=  i){cout<<"ERROR 34252000das.\n";exit(1);}
          b = a*j%p;
          if(b<r)r=b;
        }
        //r is now lowest equivalent
          vsetv2[ctr][j]=r;
          if(r==j){  //Then new!!
            k=new int[2];
            k[0]=i;
            k[1]=r;
            vset.push_back(k);
          }

      }
      ctr++;
    }
  }
  if(verbose)cout<<"Length of vset = "<<vset.size()<<"\n";


}

void SigObjectOldStyle ::printOUTarray(int* v, int n){
  int i;
  for(i=0;i<n;i++)cout<<v[i]<<",";
}

int SigObjectOldStyle ::quickDet3(int *a){
  return a[0]*a[2]-a[1]*a[1];
}


void SigObjectOldStyle ::P1Reduce(int* v){
//cout<<"P1ReduceIN: ";printOUTarray(v,2);cout<<"\n";
  v[0]=(v[0]%p+p)%p; v[1]=(v[1]%p+p)%p;
  int t, b,c;
  if(v[0]==0){
    if(v[1]==0){cout<<"ERROR 435u23rkse in SigObjectOldStyle.\n";exit(1);}
    //t=v[1]/Ngcd[v[1]];
    //v[0]=0;
    //v[1]=(v[1]*Ninv[t])%p;
    //I don't know what I was thinking above here??
    v[1]=1;
  }else{
    b=v[0]; c=v[1]; //for debugging
    t=Ninv[v[0]];
    v[0] = (v[0]*t)%p;
    v[1] = (v[1]*t)%p;
    if(v[0]!=1){
      if((divisorToIndex[v[0]]==-1)
      ||(MyGCD(v[0],v[1])!=1)
      ){cout<<"ERROR 3e2t23253 in SigObjectOldStyle ("
        <<b<<","<<c<<","<<v[0]<<","<<v[1]<<") t="<<t<<" ["<<Ngcd[b]<<","<<Ninv[b]<<"].\n";exit(1);}
      //Find equiv class...
      v[1]=vsetv2[divisorToIndex[v[0]]][v[1]];
    }
  }
  //if(Ninv[v[0]]!=v[0]){cout<<"ERROR atwtwfzeer43334 in SigObjectOldStyle.\n";exit(1);}
  //What the heck is this above check?
  //When is that ever true??
  //cout<<"P1ReduceOUT: ";printOUTarray(v,2);cout<<"\n";
}

void SigObjectOldStyle ::ReduceUPart(int a, int b, int c, int* v){
  int x,y;
//cout<<"v: ";printOUTarray(v,2);cout<<"\n";
  P1Reduce(v);
//cout<<"v: ";printOUTarray(v,2);cout<<"\n";
  if(v[0]==0)return;  //That's minimal!
  w[0]=v[0];w[1]=v[1];
  x=v[0];y=v[1];
  if(b==0){
    u[0]=x;u[1]=(p-y);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if((a==c)&&(b==0)){
    u[0]=y;u[1]=(p-x);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    // if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if(a==c){
    u[0]=y;u[1]=x;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if(a==2*b){
    u[0]=(x+y);u[1]=(p-y);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if((a==c)&&(c==2*b)){
    u[0]=(p-y);u[1]=x+y;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if((a==c)&&(c==2*b)){
    u[0]=(2*p-x-y);u[1]=x;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  if(c==2*b){
    u[0]=p-x;u[1]=x+y;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];}
    //if(u[0]==0){v[0]=0;v[1]=1;return;}
    //else if(u[1]<w[1]){w[1]=u[1];}
  }
  v[0]=w[0];v[1]=w[1];
}

void SigObjectOldStyle ::PX2Reduce(int a, int b, int c, int x, int y,  int* v5){
//input 3-vector v, output 5-vector x
  //int a=v[0],b=v[1],c=v[2];
 // int x=v2[0]; int y=v2[1];
 // cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  int t;
  bool goon=true;
  while(goon){
    if(a>c){
      t=a;a=c;c=t;
      t=x;x=y;y=t;
    }else if(b<0){
      b=-b;y=(p-y)%p;
    }else if(2*b>a){
      t=(2*b+a)/(2*a);
      c=c-2*t*b+t*t*a; b=b-t*a;
      x = (x+t*y)%p;
    }else{goon=false;}
  }
  z[0]=x;z[1]=y;
//  cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  ReduceUPart(a, b, c, z);
//  printOUTarray(z,2);cout<<"\n";
  v5[0]=a; v5[1]=b; v5[2]=c; v5[3]=z[0]; v5[4]=z[1];
}

void SigObjectOldStyle ::PX2ReduceLong(long long a, long long b, long long c, int x, int y,  int* v5){
//input 3-vector v, output 5-vector x
  //int a=v[0],b=v[1],c=v[2];
 // int x=v2[0]; int y=v2[1];
 // cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  long long t;  int u;
  bool goon=true;
  while(goon){
    if(a>c){
      t=a;a=c;c=t;
      u=x;x=y;y=u;
    }else if(b<0){
      b=-b;y=(p-y)%p;
    }else if(2L*b>a){
      t=(2L*b+a)/(2L*a);
      c=c-2L*t*b+t*t*a; b=b-t*a;
      x = (x+t*y)%p;
    }else{goon=false;}
  }
  z[0]=x;z[1]=y;
//  cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  ReduceUPart(a, b, c, z);
//  printOUTarray(z,2);cout<<"\n";
  v5[0]=a; v5[1]=b; v5[2]=c; v5[3]=z[0]; v5[4]=z[1];
}

int SigObjectOldStyle ::getIndexFromSig(int *v){
  long long d=((long long)v[0])*((long long)v[2])-((long long)v[1])*((long long)v[1]);
  if((d<0)||(d>uptoDet)){
      cout<<"ERROR 43523i4i3254 in GetIndexFromSig ("<<d<<")\n";
      printOUTarray(v,5);  cout<<"---------\n";
      exit(1);
  }
  int* m = startend[d];
  if(m[0]==-1){cout<<"ERROR 4261113 in GetIndexFromSig ("<<d<<")\n";
      printOUTarray(v,5);  cout<<"---------\n";
      exit(1);
  }
  int i;
  for(i=m[0];i<=m[1];i++){
    if(v[0]==sigs[i][0])
    if(v[1]==sigs[i][1])
    if(v[2]==sigs[i][2])
    if(v[3]==sigs[i][3])
    if(v[4]==sigs[i][4])
    return i;
  }
  cout<<"ERROR 8756321 in GetIndexFromSig - sig not found ("<<d<<")\n";
  printOUTarray(v,5);
  cout<<"---------\n";
  printOUTarray(m,2);
  cout<<"---------\n";
  printOUTarray(sigs[m[0]],5);
  cout<<"---------\n";
  exit(1);
  return -1;
}

int SigObjectOldStyle ::MyGCD(int a, int b){
  if(a<0)a=-a;  if(b<0)b=-b;
  if(a<b)return  MyGCD0(a,b);
  else if(a==b) return a;
  else return MyGCD0(b,a);
}
int SigObjectOldStyle ::MyGCD0(int a, int b){
  //a<b guaranteed
  int t;
  while(a!=0){
    t=a;
    a = b %a;
    b = t;
  }
  return b;
}

int SigObjectOldStyle ::getIndexFromSig(long long a, long long b, long long c){
  int *v = new int[5];
  PX2ReduceLong(a,b,c,1,0,v);
  //cout<<"abc="<<a<<","<<b<<","<<c<<";";  printOUTarray(v,5); cout<<"\n";
  int ans = getIndexFromSig(v);
  delete[] v;
  return ans;

}

int SigObjectOldStyle ::getIndexFromDet(int d){
  if(d>MAXDET){
    cout<<"ERROR 352952ersr33 d exceeds MAXDET ("<<d<<","<<MAXDET<<")\n";
    exit(1);
  }
  int i;
  for(i=d;i>=1;i--){
    if(startend[i][0]>=0)return startend[i][1];
  }
  return 1;
}

int SigObjectOldStyle ::getDetFromIndex(int in){
  if(in>SIGLEN){
    cout<<"ERROR 2532634ewt4 in exceeds SIGLEN\n";
    exit(1);
  }
  return quickDet3(sigs[in]);
}

SigObjectOldStyle ::SigObjectOldStyle(int p0, int d0, int verbose0){
  verbose=verbose0;
  initSigObjectOldStyle(p0, d0);
}

SigObjectOldStyle ::SigObjectOldStyle(int p0, int d0){
  verbose=0;
  initSigObjectOldStyle(p0, d0);
}

void SigObjectOldStyle ::initSigObjectOldStyle(int p0, int d0){
  uptoDet = d0;
  p = p0;
  p2=2*p;
  initFirst();

  cout<<"Making a SigObjectOldStyle with uptoDet="<<uptoDet<<"\n";

  int a,b,c, maxc,maxa,maxb,i,j;
  int *t;

  vector<int*> *tempsigs;
  int maxn=1000000;
  int currn=0; int d;
  tempsigs=new vector<int*>[uptoDet+1];

  int **zset; int znum;
  zset=new int*[vset.size()];
  for(i=0;i<vset.size();i++)zset[i]=new int[5];
  int * s=new int[5];

  maxa=mysqrtfloor((4*uptoDet)/3);
  for(a=2;a<=maxa;a+=2){
    if(verbose)cout<<a<<" of "<<maxa<<"\n";
    maxb=a/2;
    for(b=0;b<=maxb;b++){
      maxc=(b*b+uptoDet)/a;
      for(c=a;c<=maxc;c+=2){
        znum=0;
        for(i=0;i<vset.size();i++){


          if((QFormEval(a,b,c,vset[i][0],vset[i][1])%p2) == 0){
            PX2Reduce(a,b,c,vset[i][0],vset[i][1],s);


            if((QFormEval(s[0],s[1],s[2],s[3],s[4])%p2) != 0){
               cout<<a<<","<<b<<","<<c<<","<<vset[i][0]<<","<<vset[i][1]<<
               ","<<QFormEval(a,b,c,vset[i][0],vset[i][1])<<"\n";
               cout<<s[0]<<","<<s[1]<<","<<s[2]<<","<<s[3]<<","<<s[4]<<
               ","<<QFormEval(s[0],s[1],s[2],s[3],s[4])<<"\n";
               cout<<"ERROR 43rhf2342 in SigObjectOldStyle constructor.  Ouch!\n";  exit(1);

            }

            if(!alreadyInList(zset,znum,s,5)){
              for(j=0;j<5;j++)zset[znum][j]=s[j];
              znum++;
            }
          }
        }
        for(i=0;i<znum;i++){
          d=quickDet3(zset[i]);
          t=new int[5];
          for(j=0;j<5;j++)t[j]=zset[i][j];
          tempsigs[d].push_back(t);
          currn++;
        }
        /*
        for(i=0;i<znum;i++){
          if(currn>=maxn){cout<<maxn<<" maxn exceeded for sigs ERROR 252352. Aborting.\n";exit(1);}
          sigs[currn]=new int[6];
          for(j=0;j<5;j++)sigs[currn][j]=zset[i][j];
          sigs[currn][5]=quickDet3(zset[i]);
          currn++;
        }
        */
      }
    }
  }
    cout<<"SIGLEN="<<currn<<"\n";
  //cout<<"Sorting...\n";


  /*
  Save p
  Save uptoDet
  Save SIGLEN
  Save the whole orderedGoodSigs
  */


  SIGLEN=currn;
  sigs=new int*[SIGLEN+1];
  for(i=1;i<=SIGLEN;i++){
    sigs[i]=new int[5];
  }


  int ctr=1; int*k;  int ell;
  for(i=0;i<=uptoDet;i++){
    for(j=0;j<tempsigs[i].size();j++){
      k=tempsigs[i][j];
      for(ell=0;ell<5;ell++){
         sigs[ctr][ell]=k[ell];
      }
      ctr++;
      delete[] tempsigs[i][j];
    }
  }
  delete[] tempsigs;



  initRest();

}

long long SigObjectOldStyle ::QFormEval(long long a, long long b, long long c,
 long long x, long long y){
  return a*x*x+2*b*x*y+c*y*y;
}

int SigObjectOldStyle ::myfloor(double x){
  if(x>0)return((int)x);
  else{
    int x1=(int)(2-x);
    return ((int)(x+x1))-x1;
  }
}
long long SigObjectOldStyle ::myDet3(int a, int b, int c){
  return ((long long)a)*((long long)c)-((long long)b)*((long long)b);
}

int SigObjectOldStyle ::myceiling(double x){
  return -myfloor(-x);
}

void SigObjectOldStyle ::mycopy(int a[][2], int b[][2]){
  int i,j;
  for(i=0;i<2;i++)for(j=0;j<2;j++)b[i][j]=a[i][j];
}
void SigObjectOldStyle ::swap(int &x, int &y){
  int t=x;
  x=y;
  y=t;
}
long long SigObjectOldStyle ::mysqrtfloor(long long  x){
  long long ans;
  ans = (long long)sqrt((double)x);
  if(ans*ans>x)ans--;
  if(ans*ans>x){cout<<"ERROR 126011351 in mysqrtfloor("<<x<<")\n";exit(1);}
  if((ans+1)*(ans+1)<=x){cout<<"ERROR 14513455 in mysqrtfloor("<<x<<")\n";exit(1);}
  return ans;
}

int SigObjectOldStyle ::eqvecQ(int*x, int*y, int n){
  int i;
  for(i=0;i<n;i++){
    if(x[i]!=y[i])return 0;
  }
  return 1;
}

int SigObjectOldStyle ::mycomp(int*x, int*y, int n){
  int i;
  if(x[n]<y[n])return -1;
  if(x[n]>y[n])return 1;
  for(i=0;i<n;i++){
    if(x[i]<y[i])return -1;
    if(x[i]>y[i])return 1;
  }
  return 0;
}

int SigObjectOldStyle ::alreadyInList(int** x, int n, int *y, int len){
  int i;
  for(i=0;i<n;i++){
    if(eqvecQ(x[i],y,len))return 1;
  }
  return 0;
}

int SigObjectOldStyle ::twinIndex(int in){
  long long* unr = new long long[5];
  int ans;
  if((in<1)||(in>SIGLEN)){
    cout<<"ERROR 435twef3rr22340 in SigObjectOldStyle twinIndex.\n";exit(1);
  }
  //cout<<"HERE ZA: "<<sigs[in][0]<<","<<sigs[in][1]<<","<<sigs[in][2]<<","<<sigs[in][3]<<","<<
  //   sigs[in][4]<<"\n";

  doStandardUNRnoverify(sigs[in][0],sigs[in][1],sigs[in][2],sigs[in][3],
     sigs[in][4],unr);
  //cout<<"HERE ZB: "<<unr[0]<<","<<unr[1]<<","<<unr[2]<<"\n";
  //
  ans = getIndexFromSig(unr[2]*p, -unr[1], unr[0]/p);
  delete[] unr;
  return ans;
}

void SigObjectOldStyle ::doStandardUNRnoverify(int A, int B, int C, int v1, int v2,
         long long* unr){
  int w1,w2;
  long long  mat[2][2], U[2][2], T[2][2], best[2][2], temp[2][2];

  /*********
  long long a=(long long)A,b=(long long)B,c=(long long)C,d=(long long)v1,e=(long long)v2;
  if(a*d*d+2L*b*d*e+c*e*e>(MAXINT)){cout<<"OUCH 42461061 in doStandardUNR\n";exit(1);}
  **********/

  if(v1==1){w1=1;w2=0;}
  else if((v1==0)&&(v2==1)){w1=0;w2=1;}
  else MyExtendedGCD(v1,v2,w1,w2);
  //else{cout<<"ERROR e34116372 in doStandardUNR. ("<<v1<<","<<v2<<")\n";exit(1);}
  setMat22L(U, v1, v2, -w2, w1);
  setTranspose22L(U, T);
  setMat22L(mat, A, B, B, C);
  mult22L(U, mat, temp);
  mult22L(temp, T, best);

  unr[0]=best[0][0]; unr[1]=best[0][1]; unr[2]=best[1][1];
  //verify by reducing... DON'T... Would need long long version. Do this later.
  /******************
  PX2ReduceL(unr[0],unr[1],unr[2], 1,0, v5);
  if((v5[0]==A)&&(v5[1]==B)&&(v5[2]==C)&&(v5[3]==v1)&&(v5[4]==v2)){
    return;
  }
  else{
    cout<<"ERROR 432811106 in doStandardUNR.\n";exit(1);
  }
  ********************/


}
void SigObjectOldStyle ::mult22L(long long a[][2], long long b[][2], long long c[][2]){
  c[0][0]=a[0][0]* b[0][0] + a[0][1]* b[1][0];
  c[0][1]=a[0][0]* b[0][1] + a[0][1]* b[1][1];
  c[1][0]=a[1][0]* b[0][0] + a[1][1]* b[1][0];
  c[1][1]=a[1][0]* b[0][1] + a[1][1]* b[1][1];
}
void SigObjectOldStyle ::setMat22L(long long a[][2], int x, int y, int z, int w){
  a[0][0]=(long long)x;
  a[0][1]=(long long)y;
  a[1][0]=(long long)z;
  a[1][1]=(long long)w;
}
void SigObjectOldStyle ::setTranspose22L(long long b[][2], long long a[][2]){
  a[0][0]=b[0][0];
  a[0][1]=b[1][0];
  a[1][0]=b[0][1];
  a[1][1]=b[1][1];
}
void SigObjectOldStyle ::MyExtendedGCD(int a, int b, int &c, int &d){
  if(b==0){c=1;d=0;return;}
  if(abs(b)>abs(a)){return MyExtendedGCD(b,a,d,c);}
  int r=a%b;
  if(r==0){c=0;d=1;return;}
  int x,y;
  MyExtendedGCD(b, r, x, y);
  c=y; d=x-((a-r)/b)*y;
}


//********* July 3, 2010 *********************

void SigObjectOldStyle ::ReduceUPart(int a, int b, int c, int* v, long long &pm){
  int x,y;
  long long bestpm=pm;
  P1Reduce(v);
  if(v[0]==0)return;  //That's minimal!
  w[0]=v[0];w[1]=v[1];
  x=v[0];y=v[1];
  if(b==0){
    u[0]=x;u[1]=(p-y);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=-pm;}
  }
  if((a==c)&&(b==0)){
    u[0]=y;u[1]=(p-x);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=pm;}
  }
  if(a==c){
    u[0]=y;u[1]=x;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=-pm;}
  }
  if(a==2*b){
    u[0]=(x+y);u[1]=(p-y);
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=-pm;}
  }
  if((a==c)&&(c==2*b)){
    u[0]=(p-y);u[1]=x+y;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=pm;}
  }
  if((a==c)&&(c==2*b)){
    u[0]=(2*p-x-y);u[1]=x;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=pm;}
  }
  if(c==2*b){
    u[0]=p-x;u[1]=x+y;
    P1Reduce(u);
    if((u[0]<w[0])||((u[0]==w[0])&&(u[1]<w[1]))){w[0]=u[0];w[1]=u[1];bestpm=-pm;}
  }
  v[0]=w[0];v[1]=w[1];
  pm=bestpm;
}

void SigObjectOldStyle ::PX2ReduceLong(long long a, long long b, long long c, int x, int y,  int* v5,
    long long &pm){
//input 3-vector v, output 5-vector x
  //int a=v[0],b=v[1],c=v[2];
 // int x=v2[0]; int y=v2[1];
 // cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  long long t;  int u;
  bool goon=true;
  pm=1;
  while(goon){
    if(a>c){
      t=a;a=c;c=t;
      u=x;x=y;y=u;
      pm=-pm;
    }else if(b<0){
      b=-b;y=(p-y)%p;
      pm=-pm;
    }else if(2L*b>a){
      t=(2L*b+a)/(2L*a);
      c=c-2L*t*b+t*t*a; b=b-t*a;
      x = (x+t*y)%p;
    }else{goon=false;}
  }
  z[0]=x;z[1]=y;
//  cout<<a<<","<<b<<","<<c<<","<<x<<","<<y<<",\n";
  ReduceUPart(a, b, c, z, pm);
//  printOUTarray(z,2);cout<<"\n";
  v5[0]=a; v5[1]=b; v5[2]=c; v5[3]=z[0]; v5[4]=z[1];
}

int SigObjectOldStyle ::getIndexFromSig(long long a, long long b, long long c, long long &pm){
  int *v = new int[5];
  PX2ReduceLong(a,b,c,1,0,v, pm);
  //cout<<"abc="<<a<<","<<b<<","<<c<<";";  printOUTarray(v,5); cout<<"\n";
  int ans = getIndexFromSig(v);
  delete[] v;
  return ans;

}

int SigObjectOldStyle ::twinIndex(int in, long long &pm){
  long long* unr = new long long[5];
  int ans;
  if((in<1)||(in>SIGLEN)){
    cout<<"ERROR 435twef3rr22340 in SigObjectOldStyle twinIndex.\n";exit(1);
  }
  //cout<<"HERE ZA: "<<sigs[in][0]<<","<<sigs[in][1]<<","<<sigs[in][2]<<","<<sigs[in][3]<<","<<
  //   sigs[in][4]<<"\n";

  doStandardUNRnoverify(sigs[in][0],sigs[in][1],sigs[in][2],sigs[in][3],
     sigs[in][4],unr);
  //cout<<"HERE ZB: "<<unr[0]<<","<<unr[1]<<","<<unr[2]<<"\n";
  //
  ans = getIndexFromSig(unr[2]*p, -unr[1], unr[0]/p, pm);
  delete[] unr;
  return ans;
}

