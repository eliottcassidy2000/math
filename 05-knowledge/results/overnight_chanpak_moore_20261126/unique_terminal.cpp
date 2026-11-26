// Exact census of opposite exits in the b=9 mixed-spine family.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using U=unsigned __int128;
int ix[16][16];
void prep(){for(int j=1;j<16;j++)for(int i=0;i<j;i++)ix[i][j]=ix[j][i]=j*(j-1)/2+i;}
array<int,16> adj(int n,U w){
 array<int,16>a{};
 for(int j=1;j<n;j++)for(int i=0;i<j;i++)
 if(w>>ix[i][j]&1)a[j]|=1<<i;else a[i]|=1<<j;
 return a;
}
struct Profile {uint64_t h,e,R;};
const int N=14;
uint64_t D[1<<N][N],C[1<<N][N],P[1<<N];
array<uint64_t,16> finish(int n,const array<int,16>& a){
 const int M=(1<<n)-1;
 memset(D,0,(M+1)*sizeof(D[0]));
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 for(int s=1;s<M;s++)for(int i=0;i<n;i++)if(D[s][i])
  for(int b=a[i]&~s;b;b&=b-1){int k=b&-b;D[s|k][__builtin_ctz(k)]+=D[s][i];}
 array<uint64_t,16> out{};copy(D[M],D[M]+n,out.begin());return out;
}
Profile profile(int n,U w){
 auto a=adj(n,w);
 const int M=(1<<n)-1,r=n-1,R=1<<r;
 memset(D,0,(M+1)*sizeof(D[0]));
 memset(C,0,(M+1)*sizeof(C[0]));
 for(int i=0;i<n;i++)D[1<<i][i]=1;
 C[R][r]=1;P[0]=1;
 for(int s=1;s<=M;s++){
  P[s]=0;
  for(int i=0;i<n;i++){
   auto x=D[s][i],y=C[s][i];P[s]+=x;
   if(x||y)for(int b=a[i]&~s;b;b&=b-1){
    const int k=b&-b,j=__builtin_ctz(k);D[s|k][j]+=x;C[s|k][j]+=y;
   }
  }
 }
 uint64_t total=0;
 for(int s=R;s<=M;s++)for(int i=0;i<n;i++)total+=C[s][i]*P[M^s];
 return {P[M],D[M][r],total};
}
void test(int n,const map<U,vector<int>>& bases){
 uint64_t roots=0,tested=0,forks=0;
 for(const auto& [w,good]:bases){
  auto a=adj(n-1,w);auto fin=finish(n-1,a);
  roots+=good.size();
  vector<Profile> old;old.reserve(good.size());
  for(auto g:good){
   const U orient=((1<<(n-1))-1)^g;
   old.push_back(profile(n,w|(U(orient)<<ix[0][n-1])));
   assert(old.back().e==1);
  }
  for(int j=1;j<n;j++)for(int i=0;i<j;i++){
   auto fin2=fin;
   if(j<n-1){a[i]^=1<<j;a[j]^=1<<i;fin2=finish(n-1,a);a[i]^=1<<j;a[j]^=1<<i;}
   for(int k=0;k<(int)good.size();k++){
    int g=good[k];if(j==n-1)g^=1<<i;
    uint64_t e=0;for(int x=g;x;x&=x-1)e+=fin2[__builtin_ctz(x)];
    if(e%9!=8)continue;
    tested++;
    const U orient=((1<<(n-1))-1)^good[k], x=w|(U(orient)<<ix[0][n-1]);
    auto y=profile(n,x^(U(1)<<ix[i][j]));
    assert(y.e==e);
    int64_t delta=9*(int64_t(y.h+y.R)-int64_t(old[k].h+old[k].R))-int64_t(y.e)+1;
    if(delta==2){
     if(forks++<100)cout<<"fork "<<n<<" "<<uint64_t(x>>64)<<":"<<uint64_t(x)<<" "<<i<<" "<<j<<" "
     <<old[k].h<<" "<<old[k].R<<" "<<y.h<<" "<<y.e<<" "<<y.R<<endl;
    }
   }
  }
 }
 cout<<"order "<<n<<" bases "<<bases.size()<<" eligible_roots "<<roots<<" profiled_edits "<<tested<<" forks "<<forks<<endl;
}

int main(int argc,char**argv){
 const int maxn=argc==1?11:atoi(argv[1]);
 assert(maxn>=2&&maxn<=N);prep();
 vector<vector<U>> strong(maxn+1);strong[1].push_back(0);
 uint64_t A=1,B=1;
 for(int n=2;n<=maxn;n++){
  map<U,vector<int>> bases;
  int r=n-1;
  for(int c=1;c<n;c++)for(U sub:strong[c]){
   const int q=r-c;
   U w=0;
   for(int i=0;i<c;i++)for(int j=i+1;j<c;j++)if(sub>>ix[i][j]&1)
    w|=U(1)<<ix[q+i][q+j];
   for(int m=0;m<1<<q;m++){
    int g=(1<<(r-1))|m;
    bases[w].push_back(g);
    if(q?!(m&1):c>1){
     int orient=((1<<r)-1)^g;
     strong[n].push_back(w|(U(orient)<<ix[0][r]));
    }
   }
  }
  if(n==2)A=0;else {A=A+B;B=B+A;}
  uint64_t roots=0;
  for(auto& [w,v]:bases)roots+=v.size();
  assert(strong[n].size()==A&&roots==B);
  test(n,bases);
 }
}
