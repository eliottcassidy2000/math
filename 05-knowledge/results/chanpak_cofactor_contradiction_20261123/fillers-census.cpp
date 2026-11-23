#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
#include <set>
using namespace std;
bool prime(unsigned p){if(p<2)return false; for(unsigned i=2;i<=p/i;i++)if(p%i==0)return false;return true;}
template<int n> void run(){
 constexpr unsigned F=(1u<<n)-1;
 unsigned dp[F+1][n]={},h[F+1]={};
 h[0]=1;
 vector<pair<unsigned,unsigned>> E;
 for(unsigned i=0;i<n;i++)for(unsigned j=i+2;j<n;j++)E.emplace_back(i,j);
 map<unsigned,array<unsigned,n>> seeds;
 map<pair<unsigned,unsigned>,unsigned> profiles;
 unsigned maxp=0; array<unsigned,n> maxA{};
 for (unsigned mask=0;mask<1u<<E.size();mask++){
  array<unsigned,n>A={};
  for(unsigned i=0;i+1<n;i++)A[i]|=1u<<(i+1);
  for(unsigned i=0;i<E.size();i++){auto [u,v]=E[i];if(mask>>i&1)A[u]|=1u<<v;else A[v]|=1u<<u;}
  for(unsigned s=1;s<=F;s++){
    unsigned tot=0;
    for(unsigned u=s;u;u&=u-1){
     unsigned i=countr_zero(u),rest=s^(1u<<i);
     unsigned z=0;
     if(!rest)z=1;
     else for(unsigned v=A[i]&rest;v;v&=v-1)z+=dp[rest][countr_zero(v)];
     dp[s][i]=z;tot+=z;
    }
    h[s]=tot;
  }
  unsigned Q=0;
  for(unsigned s=0;s<=F;s++)Q+=h[s]*h[F^s];
  auto K=make_pair(h[F],Q);
  if(!profiles.emplace(K,mask).second)continue;
  unsigned p=40*Q+h[F];
  if(prime(p)){seeds.emplace(p,A);if(p>maxp)maxp=p,maxA=A;}
 }
 unsigned off=0;for(auto &[p,A]:seeds){unsigned v=p-1;
 bool on=(v%40==0)&&(v/40>=2)&&has_single_bit(v/40);if(!on)++off;}
 cout<<"order "<<n<<" profiles "<<profiles.size()<<" prime cofactors "<<seeds.size()
     <<" off dyadic "<<off<<" greatest "<<maxp<<" source count "<<7*maxp<<" masks";
 for(auto v:maxA)cout<<" "<<v;cout<<"\n";
 cout<<"prime list "<<n<<":";
 for(auto &[p,A]:seeds)cout<<" "<<p;
 cout<<"\n";
}
int main(int argc,char**argv){run<1>();run<2>();run<3>();run<4>();run<5>();run<6>();run<7>();if(argc>1)run<8>();}
