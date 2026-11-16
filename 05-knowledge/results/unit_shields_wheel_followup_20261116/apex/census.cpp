#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
using std::array; using std::vector; using U=uint32_t;

struct V {
  U s[2][2][3]={}; // first ends bad, last active, middle activation/exits
  V() = default;
};
V step(V const& old, bool a, bool b) {
  V v;
  for(int f=0;f<2;f++) for(int l=0;l<2;l++) for(int m=0;m<3;m++){
    U x=old.s[f][l][m];
    v.s[a][l][m]+=x;
    if(a||m) v.s[f][l][b?2:1]+=x;
    if(b||l) v.s[f][1][m]+=x;
  }
  return v;
}
struct Pair{U ab,ba;};
Pair read(V const& ab,V const& ba) {
  return {
    ab.s[0][0][0]+ab.s[0][1][0]+ab.s[0][0][1]+ab.s[0][1][1]+ba.s[0][0][1]+ba.s[0][1][1],
    ba.s[0][0][0]+ba.s[0][1][0]+ba.s[0][0][1]+ba.s[0][1][1]+ab.s[0][0][1]+ab.s[0][1][1]
  };
}
Pair eval(int n,U a,U b){
  V x,y; x.s[0][0][0]=y.s[0][0][0]=1;
  for(int i=0;i<n;i++){
    bool ai=a>>i&1,bi=b>>i&1;
    x=step(x,ai,bi); y=step(y,bi,ai);
  }
  return read(x,y);
}
U literal(int n,U a,U b,bool ab){
  int m=n+2,N=1<<m;
  vector<U> A(m);
  for(int i=0;i<n;i++){
    for(int j=i+1;j<n;j++) A[i]|=1<<j;
    if(!(a>>i&1)) A[i]|=1<<n;
    if(!(b>>i&1)) A[i]|=1<<(n+1);
  }
  A[n]=a+(ab?(1<<(n+1)):0);
  A[n+1]=b+(ab?0:(1<<n));
  vector<uint64_t> F(N*m,0);
  for(int i=0;i<m;i++) F[(1<<i)*m+i]=1;
  for(int mask=1;mask<N;mask++){
    for(int bits=mask;bits;bits&=bits-1){
      int i=__builtin_ctz(bits),rem=mask^(1<<i);
      for(int t=rem&A[i];t;t&=t-1) F[mask*m+i]+=F[rem*m+__builtin_ctz(t)];
    }
  }
  U h=0;
  for(int i=0;i<m;i++) h+=F[(N-1)*m+i];
  return h;
}
void verify_small(){
 for(int n=0;n<=5;n++)for(U a=0;a<(1u<<n);a++)for(U b=0;b<(1u<<n);b++){
  Pair t=eval(n,a,b);
  if(t.ab!=literal(n,a,b,true)||t.ba!=literal(n,a,b,false)) {std::cerr<<n<<" "<<a<<" "<<b<<" "<<t.ab<<" "<<t.ba<<" "<<literal(n,a,b,true)<<" "<<literal(n,a,b,false)<<"\n";std::abort();}
 }
 std::cout<<"Verified every tape, both apex arcs, through spine length 5 by literal subset DP.\n";
}
struct Census {
 int n,k,cap;
 U b;
 vector<uint8_t> ab,ba;
 vector<array<U,3>> watch;
 uint64_t words=0, transitions=0;
 U maxab=0,maxba=0;
 std::map<U,array<U,3>> hits;
 Census(int n0,int k0):n(n0),k(k0),cap(1) {
   for(int i=0;i<n;i++) cap*=3;
   ab.assign(cap,255);ba.assign(cap,255);
 }
 void push(int i,U a,V const& x,V const& y){
  if(i<n){
    bool bi=b>>i&1;
    push(i+1,a,step(x,0,bi),step(y,bi,0));
    push(i+1,a|(1u<<i),step(x,1,bi),step(y,bi,1));
    return;
  }
  Pair t=read(x,y); words++;
  if(!(t.ab&1) || !(t.ba&1) ||t.ab>2*cap||t.ba>2*cap) std::abort();
  int w=__builtin_popcount(b);
  ab[t.ab/2]=std::min((int)ab[t.ab/2],w);
  ba[t.ba/2]=std::min((int)ba[t.ba/2],w);
  maxab=std::max(t.ab,maxab);maxba=std::max(t.ba,maxba);
  for (auto target: {256521,524287,524289,260909,521307,6425}) {
   if(t.ab==target && !hits.count(target)) hits[target]={a,b,1};
   if(t.ba==target && !hits.count(target)) hits[target]={a,b,0};
  }
 }
 void combined(){
  for(int i=0;i<cap;i++) ab[i]=std::min(ab[i],ba[i]);
 }
 void summary(vector<uint8_t> const& ax, char const* tag){
   std::cout<<"spine "<<n<<", "<<tag<<"\n";
   for(int m=0;m<=k;m++) {
    int nfound=0,lowmiss=-1;
    for(int i=0;i<cap;i++){
      if(ax[i]<=m) nfound++;
      else if(lowmiss==-1&&i!=3&&i!=10) lowmiss=2*i+1;
    }
    std::cout<<"  |S_B|<="<<m<<" distinct="<<nfound<<", first miss other than 7,21="<<lowmiss<<"\n";
   }
   for(int exp: {16,17,18,19,20}) {
    int count[4]={},samples=0;
    std::cout<<"  below 2^"<<exp<<" visits by minimum support ";
    for(int i=0;i<(1<<(exp-1))&&i<cap;i++){
      if(ax[i]<=k) count[ax[i]]++;
      else samples++;
    }
    for(int c:count) std::cout<<c<<" ";
    std::cout<<"misses="<<samples<<"\n";
   }
   int nmiss=0;
   std::cout<<"  first 16 holes above 21:";
   for(int i=11;i<cap&&nmiss<16;i++)if(ax[i]>k){std::cout<<" "<<(2*i+1);nmiss++;}
   std::cout<<"\n";
 }
 void run(){
  V x; x.s[0][0][0]=1;
  for(b=0;b<(1u<<n);b++) if(__builtin_popcount(b)<=k) push(0,0,x,x);
  std::cout<<"Tapes per orientation="<<words<<", max H by arc "<<maxab<<" "<<maxba<<"\n";
  summary(ab,"A to B");summary(ba,"B to A");
  combined();summary(ab,"union");
  for(auto [h,p]:hits){
    if(literal(n,p[0],p[1],p[2])!=h) std::abort();
    std::cout<<"literal witness "<<h<<" masks A,B "<<p[0]<<" "<<p[1]<<" arc "<<p[2]<<"\n";
  }
 }
};
int main(int argc,char**argv){
 verify_small();
 if(argc==2) Census(atoi(argv[1]),3).run();
 else {Census(13,3).run();Census(14,3).run();}
}
