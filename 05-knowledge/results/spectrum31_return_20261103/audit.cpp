// Independent literal Hamilton counts and descendant-support certificate.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>
using I=long long;using V=std::array<I,24>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
I hamilton(const std::vector<unsigned>&T){
 int n=T.size();need(n<=16,"literal order bound");std::vector<I>d((1u<<n)*n);
 for(int j=0;j<n;++j)d[(1u<<j)*n+j]=1;
 for(unsigned m=1;m<(1u<<n);++m)for(int j=0;j<n;++j)if(I z=d[m*n+j]){
  unsigned choices=T[j]&~m;while(choices){unsigned bit=choices&-choices;choices-=bit;
   d[(m|bit)*n+__builtin_ctz(bit)]+=z;}}
 return std::accumulate(d.end()-n,d.end(),I(0));
}
int main(int argc,char**argv){
 need(argc==3,"matrix and state paths");std::ifstream mf(argv[1]),sf(argv[2]);
 std::array<std::vector<std::tuple<int,int,I>>,6>M;V root,obs;
 for(auto&x:root)mf>>x;for(auto&x:obs)mf>>x;
 for(auto&A:M){int count;mf>>count;while(count--){int i,j;I x;mf>>i>>j>>x;A.emplace_back(i,j,x);}}
 need(bool(mf),"matrix input");
 std::vector<V> states;std::map<V,int> index;std::vector<I>H;
 std::vector<std::vector<unsigned>>graphs;std::vector<std::array<int,2>>ports;
 int parent,letter;V v;int maximum=0;
 while(sf>>parent>>letter){for(auto&x:v)sf>>x;need(bool(sf),"complete state record");
  int id=states.size();need(index.emplace(v,id).second,"unique state");
  std::vector<unsigned>T;std::array<int,2>p;
  if(id==0){need(v==root&&parent==-1&&letter==-1,"root");T={2,4,1};p={0,1};}
  else{
   need(parent>=0&&parent<id&&letter>=0&&letter<6,"parent address");
   V w{};for(auto[i,j,c]:M[letter])w[j]+=states[parent][i]*c;need(w==v,"parent transfer");
   T=graphs[parent];p=ports[parent];int n=T.size(),r=p[letter/3];
   T[r]|=1u<<n;T.push_back(((1u<<n)-1)^(1u<<r));if(letter%3<2)p[letter%3]=n;
  }
  for(int a=0;a<int(T.size());++a){need(!(T[a]>>a&1),"no loops");for(int b=0;b<a;++b)need(((T[a]>>b)^(T[b]>>a))&1,"one orientation");}
  I h=0;for(int j=0;j<24;++j)h+=v[j]*obs[j];
  need(h>=3&&h<=295&&h%2&&h!=7&&h!=21,"value universe");
  need(hamilton(T)==h,"independent literal Hamilton count");
  maximum=std::max(maximum,int(T.size()));states.push_back(v);H.push_back(h);graphs.push_back(T);ports.push_back(p);
 }
 need(!states.empty(),"nonempty state certificate");
 std::vector<std::vector<int>>edges(states.size());long long eligible=0;
 for(int id=0;id<int(states.size());++id)for(auto&A:M){V w{};for(auto[i,j,c]:A)w[j]+=states[id][i]*c;
  I h=0;for(int j=0;j<24;++j)h+=w[j]*obs[j];need(h>H[id],"strict actual growth");
  if(h<=295){auto it=index.find(w);need(it!=index.end(),"complete transition closure");edges[id].push_back(it->second);++eligible;}}
 std::vector<int>order(states.size());std::iota(order.begin(),order.end(),0);
 std::sort(order.begin(),order.end(),[&](int a,int b){return H[a]>H[b];});
 std::cout<<"literal sources "<<states.size()<<" maximum order "<<maximum<<" closed transitions "<<eligible<<'\n';
 for(int h:{31,95,97}){std::vector<unsigned>tag(states.size());std::map<unsigned,int>hist;
  for(int i:order){for(int j=0;j<3;++j)if(H[i]==3*h+2*j)tag[i]|=1u<<j;for(int k:edges[i])tag[i]|=tag[k];}
  for(int i=0;i<int(states.size());++i)if(H[i]==h)++hist[tag[i]];
  std::cout<<h;for(auto[t,n]:hist)std::cout<<' '<<t<<':'<<n;std::cout<<'\n';
 }
}
