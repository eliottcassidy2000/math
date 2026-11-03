// Exact labelled one-ear ternary branching, with the full ordered cut kernel.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <vector>
using I=long long;using G=std::vector<unsigned>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
std::vector<I> ends(const G&T){int n=T.size();need(n<=20,"order guard: incomplete");std::vector<I>d((1u<<n)*n);
 for(int i=0;i<n;++i)d[(1u<<i)*n+i]=1;
 for(unsigned s=1;s<(1u<<n);++s)for(int a=0;a<n;++a)if(I z=d[s*n+a]){unsigned r=T[a]&~s;while(r){unsigned b=r&-r;r-=b;d[(s|b)*n+__builtin_ctz(b)]+=z;}}
 return d;}
G converse(const G&T){G U=T;unsigned full=(1u<<T.size())-1;for(int i=0;i<(int)T.size();++i)U[i]=full^T[i]^(1u<<i);return U;}
I H(const G&T){auto d=ends(T);return std::accumulate(d.end()-T.size(),d.end(),I(0));}
G ear(G T,unsigned in){unsigned n=T.size(),full=(1u<<n)-1;for(int i=0;i<(int)n;++i)if(in>>i&1)T[i]|=1u<<n;T.push_back(full^in);return T;}
std::vector<I> allcuts(const G&T){
 int n=T.size();need(n<=18,"cut count integer guard");unsigned full=(1u<<n)-1;auto e=ends(T),s=ends(converse(T));std::vector<I>w(n*n);
 for(unsigned mask=1;mask<full;++mask)for(int a=0;a<n;++a)if(e[mask*n+a])for(int b=0;b<n;++b)w[a*n+b]+=e[mask*n+a]*s[(full^mask)*n+b];
 I h=std::accumulate(e.end()-n,e.end(),I(0));std::vector<I>values(1u<<n,h),linear(n),pair(n*n);
 for(int a=0;a<n;++a){linear[a]=e[full*n+a]-s[full*n+a];for(int b=0;b<n;++b){linear[a]+=w[a*n+b];pair[a*n+b]=w[a*n+b]+w[b*n+a];}}
 for(unsigned mask=1;mask<=full;++mask){unsigned bit=mask&-mask,old=mask^bit;int a=__builtin_ctz(bit);I v=values[old]+linear[a];for(int b=0;b<n;++b)if(old>>b&1)v-=pair[a*n+b];values[mask]=v;}
 need(values[0]==h&&values[full]==h,"padding cuts");return values;
}
struct Node{G t;I h;int parent,cut,depth;std::array<std::vector<int>,2>child;unsigned survive=1;};
#ifndef EAR_LIBRARY
int main(int argc,char**argv){need(argc>=4,"depth out rows");int depth=std::stoi(argv[1]);std::ofstream out(argv[2]);G source;for(int i=3;i<argc;++i)source.push_back(std::stoul(argv[i]));
 std::vector<Node>nodes={{source,H(source),-1,-1,0,{}}};I checked=0;
 for(int i=0;i<(int)nodes.size();++i){if(nodes[i].depth==depth)continue;G t=nodes[i].t;I h=nodes[i].h;int dep=nodes[i].depth;auto counts=allcuts(t);checked+=counts.size();
  for(unsigned cut=0;cut<counts.size();++cut)for(int j=0;j<2;++j)if(counts[cut]==3*h+2+2*j){G u=ear(t,cut);need(H(u)==counts[cut],"independent child count");nodes[i].child[j].push_back(nodes.size());nodes.push_back({u,counts[cut],i,(int)cut,dep+1,{}});}
 }
 std::map<std::pair<int,unsigned>,int>hist;
 for(int i=nodes.size()-1;i>=0;--i){unsigned left=0,right=0;for(int c:nodes[i].child[0])left|=nodes[c].survive;for(int c:nodes[i].child[1])right|=nodes[c].survive;nodes[i].survive=1|((left&right)<<1);++hist[{nodes[i].depth,nodes[i].survive}];}
 std::cout<<"nodes "<<nodes.size()<<" cut evaluations "<<checked<<" root survival "<<nodes[0].survive<<'\n';for(auto[k,c]:hist)std::cout<<k.first<<' '<<k.second<<' '<<c<<'\n';
 for(int i=0;i<(int)nodes.size();++i){auto&t=nodes[i];out<<i<<' '<<t.parent<<' '<<t.cut<<' '<<t.depth<<' '<<t.h<<' '<<t.survive<<' '<<t.t.size();for(auto r:t.t)out<<' '<<r;out<<'\n';}
}
#endif
