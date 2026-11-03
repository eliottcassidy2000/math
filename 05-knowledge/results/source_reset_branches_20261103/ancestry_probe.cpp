// Exact bounded adaptive ancestry census; no port-observer quotient.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>
using Graph=std::vector<unsigned>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
long long H(const Graph&T){
 int n=T.size();need(n<=20,"order resource guard reached: census incomplete");
 std::vector<long long>d((1u<<n)*n);
 for(int i=0;i<n;++i)d[(1u<<i)*n+i]=1;
 for(unsigned s=1;s<(1u<<n);++s)for(int i=0;i<n;++i)if(auto x=d[s*n+i]){
  unsigned avail=T[i]&~s;while(avail){unsigned b=avail&-avail;avail-=b;d[(s|b)*n+__builtin_ctz(b)]+=x;}}
 return std::accumulate(d.end()-n,d.end(),0LL);
}
struct Node{Graph t;std::vector<int>roots,children;long long h;};
int main(int argc,char**argv){
 need(argc==3,"bound and output");int bound=std::stoi(argv[1]);std::ofstream out(argv[2]);
 std::vector<Node>nodes={{{2,4,1},{},{},3}};long long tried=0;
 for(int i=0;i<(int)nodes.size();++i){
  Graph t=nodes[i].t;auto roots=nodes[i].roots;int n=t.size();
  for(int r=0;r<(i==0?1:n);++r){Graph u=t;u[r]|=1u<<n;u.push_back(((1u<<n)-1)^(1u<<r));
   auto h=H(u);++tried;need(h>nodes[i].h,"strict count growth");if(h>bound)continue;
   auto w=roots;w.push_back(r);nodes[i].children.push_back(nodes.size());nodes.push_back({u,w,{},h});}
 }
 int maxn=0;std::map<int,int>orders;
 for(int i=0;i<(int)nodes.size();++i){auto&s=nodes[i];maxn=std::max(maxn,(int)s.t.size());++orders[s.t.size()];
  out<<s.h<<' '<<s.t.size();for(auto x:s.t)out<<' '<<x;out<<' '<<s.roots.size();for(auto x:s.roots)out<<' '<<x;out<<' '<<s.children.size();for(auto x:s.children)out<<' '<<x;out<<'\n';}
 std::cout<<"classes "<<nodes.size()<<" probes "<<tried<<" max_order "<<maxn<<" orders";for(auto[n,c]:orders)std::cout<<' '<<n<<':'<<c;std::cout<<'\n';
 for(int h:{31,95,97}){std::vector<int>mask(nodes.size());std::map<int,int>hist;
  for(int i=(int)nodes.size()-1;i>=0;--i){for(int j=0;j<3;++j)if(nodes[i].h==3*h+2*j)mask[i]|=1<<j;for(int c:nodes[i].children)mask[i]|=mask[c];if(nodes[i].h==h)++hist[mask[i]];}
  std::cout<<h;for(auto[m,c]:hist)std::cout<<' '<<m<<':'<<c;std::cout<<'\n';}
}
