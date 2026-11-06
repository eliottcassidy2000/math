// Independent append-at-end path DP and physical construction controls.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using G=std::vector<unsigned>;using I=std::int64_t;
// h,s,e,L,R,d,C,b,T
using P=std::array<I,9>;
int profiles=0;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<=16,"order");for(int i=0;i<n;++i){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"range/loop");for(int j=0;j<i;++j)need((a[i]>>j&1)+(a[j]>>i&1)==1,"orientation");}}
P profile(const G&a,int root,bool full=true){valid(a);++profiles;int n=a.size(),N=1<<n,F=N-1;std::vector<I>all(std::size_t(N)*n),marked(std::size_t(N)*n),h(N),s(N),e(N);h[0]=1;for(int i=0;i<n;++i)all[std::size_t(1<<i)*n+i]=1;marked[std::size_t(1<<root)*n+root]=1;
 for(int mask=1;mask<N;++mask)for(unsigned ends=mask;ends;ends&=ends-1){int last=__builtin_ctz(ends);std::size_t at=std::size_t(mask)*n+last;I v=all[at],w=marked[at];h[mask]+=v;s[mask]+=w;if(last==root)e[mask]=v;if(!v&&!w)continue;for(unsigned next=a[last]&~mask;next;next&=next-1){int u=__builtin_ctz(next);std::size_t to=std::size_t(mask|(1<<u))*n+u;all[to]+=v;marked[to]+=w;}}
 I L=0,R=0,C=0,b=0,T=0;for(int mask=0;mask<N;++mask){L+=e[mask]*h[F^mask];R+=s[mask]*h[F^mask];C+=h[mask]*h[F^mask];}int body=F^(1<<root);for(int mask=body;;mask=(mask-1)&body){b+=h[mask]*h[body^mask];if(!mask)break;}
 if(full){std::vector<I>bins(N);for(int mask=body;;mask=(mask-1)&body){for(int cut=mask;;cut=(cut-1)&mask){bins[mask]+=h[cut]*h[mask^cut];if(!cut)break;}T+=s[F^mask]*bins[mask];if(!mask)break;}}
 return {h[F],s[F],e[F],L,R,h[body],C,b,T};}
G move(G a,int root,char op){int n=a.size();a.resize(n+1);for(int i=0;i<n;++i){if(op!='O'&&i==root)a[n]|=1u<<i;else a[i]|=1u<<n;}valid(a);return a;}
G child(const G&a,int root){G c={2,8,3,5};c.resize(a.size()+3);std::vector<int>map(a.size());map[root]=0;int next=4;for(int i=0;i<(int)a.size();++i)if(i!=root)map[i]=next++;for(int i=0;i<(int)a.size();++i)for(int j=0;j<(int)a.size();++j)if(a[i]>>j&1)c[map[i]]|=1u<<map[j];for(int i=1;i<4;++i)for(int j=0;j<(int)a.size();++j)if(j!=root)c[i]|=1u<<map[j];valid(c);return c;}
P transfer(P p,char op,int order){auto[h,s,e,L,R,d,C,b,T]=p;if(op=='O')return{h,s,0,L,2*R,d,2*C,2*b,3*T};if(op=='Z')return{h+R-e,order==1?0:s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b};return{h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C};}
std::array<I,7> seven(P p){return{p[0],p[1],p[2],p[3],p[4],p[5],p[7]};}
I base(P p,int j=0){return p[0]+p[3]-p[1]+2*j;}
std::array<I,4> row(P p,int j){return{base(p,j),p[2],p[4]+p[7],p[5]};}
bool strong(const G&a){int n=a.size();for(int root=0;root<n;++root){unsigned seen=1u<<root,old;do{old=seen;for(int i=0;i<n;++i)if(seen>>i&1)seen|=a[i];}while(old!=seen);if(seen!=(1u<<n)-1)return false;}return true;}
I literal_f(const G&a,int root){int n=a.size();std::vector<int>p;for(int i=0;i<n;++i)if(i!=root)p.push_back(i);I f=0;do{std::vector<int>w{root};w.insert(w.end(),p.begin(),p.end());for(int cut=2;cut<n;++cut){bool ok=bool(a[root]>>w.back()&1);for(int j=0;j<n-1;++j)if(j!=cut-1&&!(a[w[j]]>>w[j+1]&1))ok=false;f+=ok;}}while(std::next_permutation(p.begin(),p.end()));return f;}
G graph(int n,unsigned code){G a(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(code>>bit++&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
