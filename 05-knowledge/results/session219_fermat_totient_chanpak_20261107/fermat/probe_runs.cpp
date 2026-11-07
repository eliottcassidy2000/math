#include <cstdio>
#include <array>
#include <cstdint>
using ll=long long;const ll P=641,Q=6700417,M=P*Q;
std::array<ll,65> powers;long long checked[5]={};long long hits[5]={};
struct S{ll D,E,U;};
ll md(ll x){x%=M;return x<0?x+M:x;}
ll mul(ll x,ll y){return (ll)((__int128)x*y%M);}
S run(S s,int op,int k){ll t=powers[k];if(op==0)return {s.D,0,mul(t,s.U)};
 return {md(s.D-s.E+s.U+mul(md(t-2),md(s.U-1))),1,md(mul(t,md(s.U-1))+1)};}
void visit(S s,int op,int depth){if(depth){checked[depth]++;ll h=md(s.D-s.E+s.U);if(h%Q==0&&h%P!=0){hits[depth]++;printf("escape depth%d hmod%lld\n",depth,h);}}
 if(depth==4)return;for(int k=1;k<=64;k++)visit(run(s,op,k),1-op,depth+1);}
int main(){powers[0]=1;for(int k=1;k<=64;k++)powers[k]=2*powers[k-1]%M;visit({1,0,2},0,0);visit({1,0,2},1,0);for(int d=1;d<=4;d++){printf("runs%d checked%lld escape%lld\n",d,checked[d],hits[d]);if(hits[d])return 2;long long expected=2;for(int k=0;k<d;k++)expected*=64;if(checked[d]!=expected)return 3;}return 0;}
