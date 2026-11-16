#include <array>
#include <iostream>
#include <vector>
#include <cstdlib>
using namespace std;
using V=array<int,4>;using M=array<int,16>;
const int b[]={8,16,16,32};
const int shift[]={0,3,7,11};
int enc(V x){
 int v=0;for(int i=0;i<4;i++)v|=(x[i]&(b[i]-1))<<shift[i];return v;
}
V dec(int z){V a;for(int i=0;i<4;i++)a[i]=(z>>shift[i])&(b[i]-1);return a;}
int step(M A,int z){
 V x=dec(z),v{};
 for(int i=0;i<4;i++)for(int j=0;j<4;j++)v[i]+=A[4*i+j]*x[j];
 return enc(v);
}
int main(){
 int n;V x,f;cin>>n;for(auto&a:x)cin>>a;for(auto&a:f)cin>>a;
 vector<M> T(n);for(auto&M:T)for(auto&m:M){cin>>m;m&=31;}
 auto read=[&](int z){V v=dec(z);int t=0;for(int i=0;i<4;i++)t+=f[i]*v[i];return t&31;};
 vector<int> q{enc(x)};vector<bool> seen(65536);
 vector<unsigned> cov(65536);seen[q[0]]=true;int bad=0;
 for(size_t k=0;k<q.size();k++){
  auto a=q[k];for(auto M:T){auto z=step(M,a);auto t=read(z);if(!(t&1))return 1;
   cov[a]|=1u<<(t>>1);if(!seen[z]){seen[z]=true;q.push_back(z);}
  }bad+=(cov[a]!=65535);
 }
 if(q.size()!=32768||bad!=80)return 2;
 for(auto a:q){unsigned C=0;for(auto M:T){C|=cov[step(M,a)];if(C==65535)break;}if(C!=65535)return 3;}
 cout<<"32768 states; only 80 need a second return block; all pass the two-block test.\n";
}
