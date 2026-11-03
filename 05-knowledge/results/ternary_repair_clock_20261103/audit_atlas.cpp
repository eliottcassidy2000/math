// Independent permutation / ordered-cut expansion over every source bit cube.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using I=long long;using V=std::array<I,3>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
struct Table{std::vector<I> h,e,c;};
Table build(int old,int tail){
 int n=old+tail,m=(old-1)*(old-2)/2,N=1<<m;int bit[10][10];
 for(auto& row:bit)for(int& x:row)x=-1;int k=0;
 for(int u=0;u<old;++u)for(int v=u+2;v<old;++v)bit[u][v]=bit[v][u]=k++;
 Table out{std::vector<I>(N),std::vector<I>(N),std::vector<I>(N)};
 std::vector<int> p(n);std::iota(p.begin(),p.end(),0);
 auto constraints=[&](int cut,unsigned&ones,unsigned&zero){
  ones=zero=0;
  for(int j=0;j+1<n;++j){if(j==cut)continue;int u=p[j],v=p[j+1],b=bit[u][v];
   if(b>=0){if(u<v)zero|=1u<<b;else ones|=1u<<b;}
   else if(!((std::abs(u-v)==1)?u<v:u>v))return false;
  }
  return !(ones&zero);
 };
 auto add=[&](unsigned ones,unsigned zero,std::vector<I>&a,std::vector<I>*b){
  unsigned free=(N-1)^(ones|zero),s=free;
  for(;;){++a[ones|s];if(b)++(*b)[ones|s];if(!s)break;s=(s-1)&free;}
 };
 do{unsigned a,b;if(constraints(-1,a,b))add(a,b,out.h,p.back()==0?&out.e:nullptr);}
 while(std::next_permutation(p.begin(),p.end()));
 std::iota(p.begin(),p.end(),0);
 do{for(int cut=0;cut<n;++cut){unsigned a,b;
   // Distinguished first component starts at0; cut=n-1 means empty second.
   if(constraints(cut,a,b))add(a,b,out.c,nullptr);
  }}while(std::next_permutation(p.begin()+1,p.end()));
 for(int c=0;c<N;++c)need(out.h[c]%2==1&&2*out.e[c]<out.h[c],"parity and nonsink endpoint bound");
 return out;
}
V at(const Table&t,int code){need((t.c[code]-t.e[code])%2==0,"parity");
 return {t.h[code],(t.c[code]-t.e[code])/2,(t.c[code]+t.e[code])/2};}
bool strong_cut(int n,unsigned code){
 for(int cut=1;cut<n;++cut){unsigned mask=0;int bit=0;
  for(int u=0;u<n;++u)for(int v=u+1;v<n;++v)if(v!=u+1){if(u<cut&&v>=cut)mask|=1u<<bit;++bit;}
  if((code&mask)==mask)return false;
 }
 return true;
}
int main(int argc,char**argv){
 need(argc==2||argc==3,"atlas7 tsv and optional eight-fringe tsv");
 std::array<Table,3> t={build(7,0),build(7,1),build(7,2)};
 std::ifstream in(argv[1]);int code,count=0;std::set<int> codes;
 while(in>>code){std::array<I,9>P;for(I&v:P)in>>v;need(codes.insert(code).second,"unique chart");
  V a=at(t[0],code),b=at(t[1],code),c=at(t[2],code);
  for(int j=0;j<3;++j){need(P[j]==a[j],"initial");need(P[3+j]==c[j]-b[j],"row1");need(P[6+j]==b[j]-a[j],"row2");}++count;
 }
 for(int c=0;c<(1<<15);++c)need(bool(codes.count(c))==strong_cut(7,c),"complete strong inventory");
 std::cout<<"independent7 matrices "<<count<<'\n';
 if(argc==3){std::ifstream target(argv[2]);V v;std::set<V>wanted;
  while(target>>v[0]>>v[1]>>v[2])wanted.insert(v);
  Table eight=build(8,0);int total=0,hits=0;
  for(int c=0;c<(1<<21);++c)if(strong_cut(8,c)){++total;if(wanted.count(at(eight,c)))++hits;}
  need(hits==0,"independent8 zero fringe hits");
  std::cout<<"independent8 strong_charts "<<total<<" fringe_hits "<<hits<<'\n';
 }
}
