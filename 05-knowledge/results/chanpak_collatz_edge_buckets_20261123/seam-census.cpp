#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>
// Bit i in a row is the arc to vertex i. Fix a Hamilton path only in the original host.
using namespace std;
template<int N> struct Search {
 using A=array<int,N>; static constexpr int F=(1<<N)-1;
 array<array<unsigned short,N>,1<<N> d{};
 vector<pair<int,int>> E;
 vector<uint32_t> centers;
 unsigned strongest=0;
 Search(){for(int i=0;i<N;i++)for(int j=i+2;j<N;j++)E.emplace_back(i,j);}
 unsigned H(const A&a){
  for(int s=1;s<=F;s++){
   for(unsigned q=s;q;q&=q-1){
    int i=countr_zero(q), t=s^(1<<i);
    if(!t){d[s][i]=1;continue;}
    unsigned z=0;
    for(unsigned b=t&a[i];b;b&=b-1)z+=d[t][countr_zero(b)];
    d[s][i]=z;
   }
  }
  unsigned z=0;for(int i=0;i<N;i++)z+=d[F][i];
  return z;
 }
 A get(uint32_t m){
  A a{}; for(int i=0;i+1<N;i++)a[i]|=1<<(i+1);
  for(unsigned b=0;b<E.size();b++){
   auto [i,j]=E[b];if(m>>b&1)a[i]|=1<<j;else a[j]|=1<<i;
  }
  return a;
 }
 void run(uint32_t lo,uint32_t hi){
  unsigned num=0,n61=0,n65=0,c61=0,c63=0,c65=0;
  for(uint32_t m=lo;m<hi;m++){
   auto a=get(m);
   unsigned h=H(a);strongest=max(h,strongest);
   if(h==61)++c61;
   if(h==65)++c65;
   if(h!=63&&h!=61)continue;
   if(h==63)++c63;
   vector<pair<int,int>> L,R;
   for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){
    a[i]^=1<<j;a[j]^=1<<i;
    unsigned p=H(a);
    if(p==(h==63?61:63))L.emplace_back(i,j);
    if(p==65)R.emplace_back(i,j);
    a[i]^=1<<j;a[j]^=1<<i;
   }
   if(h==63){n61+=!L.empty();n65+=!R.empty();}
   if(h==61){
    if(!R.empty()){
      cout<<"DIRECT61 "<<N<<" "<<m;
      for(int x:a)cout<<" "<<x;
      for(auto [i,j]:R)cout<<" "<<i<<" "<<j;
      cout<<endl;
    }
    continue;
   }
   if(L.empty())continue;
   ++num;
   if(num<=200){
    cout<<"HIT order "<<N<<" mask "<<m<<" rows";
    for(int z:a)cout<<" "<<z;
    for(auto [i,j]:L)cout<<" L "<<i<<" "<<j;
    for(auto [i,j]:R)cout<<" R "<<i<<" "<<j;
    cout<<endl;
   }
  }
  cout<<"order "<<N<<" lo "<<lo<<" hi "<<hi<<" max "<<strongest<<" centers 61/63/65 "<<c61<<"/"<<c63<<"/"<<c65<<" centers adjacent to 61 "<<num<<" pair 61/65 "<<n61<<"/"<<n65<<endl;
 }
};
int main(int argc,char**argv){
 int n=atoi(argv[1]);
 uint32_t l=argc>2?strtoul(argv[2],nullptr,0):0;
 uint32_t r=argc>3?strtoul(argv[3],nullptr,0):(1u<<((n-1)*(n-2)/2));
 if(n==6)Search<6>().run(l,r);
 if(n==7)Search<7>().run(l,r);
 if(n==8)Search<8>().run(l,r);
}
