#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <tuple>
#include <vector>
#include <set>
using R=std::array<int,8>;
struct D {
 int a[256][8]={},h[256]={1};
 D(const R&r,int n){
  for(int i=0;i<n;i++)a[1<<i][i]=1;
  for(int s=1;s<1<<n;s++)for(int i=0;i<n;i++)if(int c=a[s][i]){
   h[s]+=c;
   for(int v=r[i]&~s;v;v&=v-1){int b=v&-v;a[s|b][__builtin_ctz(b)]+=c;}
  }
 }
};
int main(int argc,char**argv){
 if(argc!=2)return 1;
 std::ifstream f(argv[1]);std::string s;std::getline(f,s);
 std::map<std::tuple<int,int,int>,std::pair<R,int>> bank;
 int n=8,N=255;
 while(std::getline(f,s)){
  std::istringstream t(s);unsigned long long b;char c;int h;
  if(!(t>>b>>c>>h)||c!=',')return 2;
  R r{},rev{};
  for(int i=0;i<n;i++){r[i]=(b>>(8*i))&255;rev[i]=255^r[i]^(1<<i);}
  D d(r,n),e(rev,n);
  if(d.h[N]!=h)return 3;
  for(int i=0;i<n;i++){
   int q[2]={},g=d.h[N^(1<<i)];
   for(int s=1;s<=N;s++)if(s>>i&1)q[__builtin_popcount(s)&1]+=e.a[s][i]*d.h[N^s];
   if(q[0] && q[1]-q[0]==g)bank.emplace(std::tuple{h,g,q[0]},std::pair{r,i});
  }
 }
 int fval[31001]={1};
 for(int v=1;v<=31000;v++){int u=v/2;fval[v]=fval[u]+(v%2?0:u+u%2);}
 std::set<int> target={126,378,10008,30012,30090};
 int pairs=0;
 for(auto const& [kg,a]:bank){
  auto [h,g,r]=kg;auto b=bank.find({h+2,g,r});if(b==bank.end())continue;
  ++pairs;
  for(int z=1;h+1+2*r*z<=*target.rbegin();z++){
   int c=h+1+2*r*z+g*(fval[z]-1);
   if(!target.count(c))continue;
   std::cout<<c<<" h,g,r,z "<<h<<" "<<g<<" "<<r<<" "<<z;
   for(auto [a,j]:{a,b->second}){
    std::cout<<" root "<<j<<" rows ";for(auto v:a)std::cout<<v<<",";
   }
   std::cout<<"\n";
  }
 }
 std::cout<<bank.size()<<" stable centre keys; "<<pairs<<" matched pairs\n";
}
