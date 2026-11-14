// S328: reproducible finite certificate. Not an exhaustive order-nine census.
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <algorithm>
using namespace std;
constexpr int LO=69, HI=79115;
struct W { int n=0,j=0,r=0,h=0,l=0; array<int,9> a{}; };
W WIT[HI+1];
unsigned COVER[HI+1]{};
void require(bool good){if(!good)throw runtime_error("Fermat range replay failed");}
bool reach_all(const vector<int>&A,int r,bool converse){
 int n=A.size(),v=(1<<n)-1,x=1<<r;
 while(true){
  int y=x;
  for(int s=x;s;s&=s-1){int i=__builtin_ctz(s);y|=(converse?(v^A[i]^(1<<i)):A[i]);}
  if(y==x)return x==v;
  x=y;
 }
}
void scan(const vector<int>&A){
 const int n=A.size(),v=(1<<n)-1;
 if(!reach_all(A,0,false)||!reach_all(A,0,true))return;
 int E[512][9]{},H[512]{};H[0]=1;
 for(int i=0;i<n;i++)E[1<<i][i]=1;
 for(int M=1;M<=v;M++)
  for(int i=0;i<n;i++)if(E[M][i]){
   int c=E[M][i];H[M]+=c;
   for(int r=A[i]&~M;r;r&=r-1)E[M|(r&-r)][__builtin_ctz(r)]+=c;
  }
 memset(E,0,sizeof(E));
 for(int i=0;i<n;i++)E[1<<i][i]=1;
 for(int M=1;M<=v;M++)
  for(int i=0;i<n;i++)if(E[M][i]){
   int c=E[M][i];
   for(int r=(v^A[i]^(1<<i))&~M;r;r&=r-1)E[M|(r&-r)][__builtin_ctz(r)]+=c;
  }
 for(int r=0;r<n;r++){
  int M=v^(1<<r),L=H[M];
  for(int Q=M;Q;Q=(Q-1)&M){
   int S=0;for(int b=Q&A[r];b;b&=b-1)S+=E[Q][__builtin_ctz(b)];
   L+=H[M^Q]*S;
  }
  for(int j=1,bit=0;j<=4;j*=2,bit++){
   int64_t N=H[v]+(int64_t(1)<<j)*L;
   if(N>HI)break;
   if(N<LO)continue;
   COVER[N]|=1u<<bit;
   if(!WIT[N].n||WIT[N].n>n||(WIT[N].n==n&&WIT[N].j>j)){
    auto&w=WIT[N];w.n=n;w.j=j;w.r=r;w.h=H[v];w.l=L;
    copy(A.begin(),A.end(),w.a.begin());
   }
  }
 }
}
vector<int> expand(const W&w){
 const int n=w.n,j=w.j,M=(1<<(n+j+1))-1,b=M^((1<<n)-1);
 vector<int>T(n+j+1);
 for(int i=0;i<n;i++)T[i]=w.a[i]|(i==w.r?(b^(1<<(n+j))):b);
 for(int i=n;i<n+j+1;i++)T[i]=(b&~((1<<(i+1))-1))|(i==n+j?(1<<w.r):0);
 return T;
}
uint64_t direct(const vector<int>&T){
 int n=T.size(),N=1<<n;vector<uint64_t> E(size_t(N)*n);for(int i=0;i<n;i++)E[(1<<i)*n+i]=1;
 for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(E[s*n+i])
  for(int b=T[i]&~s;b;b&=b-1)E[(s|(b&-b))*n+__builtin_ctz(b)]+=E[s*n+i];
 uint64_t h=0;for(int i=0;i<n;i++)h+=E[(N-1)*n+i];return h;
}
int main(int argc,char**argv){
 for(int n=1;n<=9;n++){
  uint64_t seed=325;
  const int count=n<=6?(1<<(n*(n-1)/2)):200000;
  for(int k=0;k<count;k++){
   seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;
   uint64_t code=n<=6?uint64_t(k):seed;
   vector<int>A(n);
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,code>>=1)
    if(code&1)A[i]|=1<<j;else A[j]|=1<<i;
   scan(A);
  }
 }
 int counts[10]{},blocks[17]{},sample[5]{};
 uint64_t digest=14695981039346656037ull;
 for(int N=LO;N<=HI;N+=2){
  auto&w=WIT[N];require(w.n&&w.n<=9&&(w.j&(w.j-1))==0&&w.j<=4);
  require(w.h+(1LL<<w.j)*w.l==N);
  counts[w.n]++;blocks[w.j]++;
  const int b=__builtin_ctz(w.j);
  if(!sample[b])sample[b]=N;
  for(int x:{N,w.n,w.j,w.r,w.h,w.l}){digest^=x;digest*=1099511628211ULL;}
  for(int i=0;i<w.n;i++){digest^=w.a[i];digest*=1099511628211ULL;}
 }
 for(int N:{69,2887,30013,30029,30047,59999,79115}){
  auto T=expand(WIT[N]);require(reach_all(T,0,false)&&reach_all(T,0,true)&&direct(T)==uint64_t(N));
 }
 for(int N:sample)if(N){auto T=expand(WIT[N]);require(direct(T)==uint64_t(N));}
 cout<<"strong Fermat coverage "<<LO<<".."<<HI<<" (all "<<(HI-LO)/2+1<<" odds); digest "<<digest<<'\n';
 cout<<"witnessed preferred source orders:";for(int n=1;n<=9;n++)cout<<' '<<counts[n];cout<<'\n';
 cout<<"witnessed preferred Fermat exponents:";for(int j=1;j<=4;j*=2)cout<<' '<<blocks[j];cout<<'\n';
 cout<<"first unwitnessed targets at progressively larger cutoffs:";for(int b=0;b<3;b++){
  int N=LO;while(N<=HI&&(COVER[N]&((1u<<(b+1))-1)))N+=2;cout<<' '<<N;
 }cout<<'\n';
 for(int i=1;i<argc;i++){
  int t=stoi(argv[i]);require(t>=LO&&t<=HI&&t%2);
  auto&w=WIT[t];cout<<"target "<<t<<": source";
  for(int k=0;k<w.n;k++)cout<<' '<<w.a[k];
  cout<<"; root "<<w.r<<"; j "<<w.j<<"; H "<<w.h<<"; L "<<w.l<<'\n';
 }

}
