// Independent bank closure using literal path enumeration and exact isomorphism search.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>
using G=std::vector<unsigned>;using I=std::uint64_t;using Pair=std::pair<int,int>;
I path_visits=0,iso_calls=0,iso_nodes=0;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<=15,"safe order");for(int i=0;i<n;++i){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"row");for(int j=0;j<i;++j)need((a[i]>>j&1)+(a[j]>>i&1)==1,"orientation");}}
bool strong(const G&a){int n=a.size();for(int r=0;r<n;++r){unsigned seen=1u<<r,old;do{old=seen;for(int i=0;i<n;++i)if(seen>>i&1)seen|=a[i];}while(old!=seen);if(seen!=(1u<<n)-1)return false;}return true;}
I hp(const G&a,I cap=1000000){I result=0;unsigned full=(1u<<a.size())-1;auto dfs=[&](auto&&self,int last,unsigned used)->void{++path_visits;if(result>cap)return;if(used==full){++result;return;}unsigned next=a[last]&~used;while(next&&result<=cap){unsigned bit=next&-next;next-=bit;self(self,__builtin_ctz(bit),used|bit);}};for(int r=0;r<(int)a.size()&&result<=cap;++r)dfs(dfs,r,1u<<r);return result;}
using Features=std::vector<std::vector<int>>;
Features features(const G&a){int n=a.size();std::vector<int>d(n);for(int i=0;i<n;++i)d[i]=__builtin_popcount(a[i]);Features f(n,std::vector<int>(n+1));for(int i=0;i<n;++i){f[i][0]=d[i];for(int j=0;j<n;++j)if(a[i]>>j&1)++f[i][1+d[j]];}return f;}
std::vector<int> signature(const G&a,int H){auto f=features(a);std::sort(f.begin(),f.end());std::vector<int>key{int(a.size()),H};for(auto&v:f)key.insert(key.end(),v.begin(),v.end());return key;}
bool isomorphic(const G&a,const G&b,bool fixed_root=false){++iso_calls;if(a.size()!=b.size())return false;int n=a.size();auto fa=features(a),fb=features(b);std::vector<int>map(n,-1);unsigned used=0;auto search=[&](auto&&self,int depth)->bool{++iso_nodes;if(depth==n)return true;int v=-1;std::vector<int>best;for(int i=0;i<n;++i)if(map[i]<0){std::vector<int>candidates;for(int j=0;j<n;++j)if(!(used>>j&1)&&fa[i]==fb[j]){bool ok=true;for(int k=0;k<n;++k)if(map[k]>=0&&bool(a[i]>>k&1)!=bool(b[j]>>map[k]&1)){ok=false;break;}if(ok)candidates.push_back(j);}if(candidates.empty())return false;if(v<0||candidates.size()<best.size()){v=i;best=std::move(candidates);}}for(int j:best){map[v]=j;used|=1u<<j;if(self(self,depth+1))return true;used^=1u<<j;map[v]=-1;}return false;};if(fixed_root){map[0]=0;used=1;return search(search,1);}return search(search,0);}
struct Entry{G a;int H;I code;};
struct Profile{I h,s,L;};
Profile profile(const G&a,int root){valid(a);int n=a.size(),N=1<<n,F=N-1;std::vector<I>h(N),e(N),s(N);h[0]=1;auto dfs=[&](auto&&self,int first,int last,unsigned used)->void{++path_visits;++h[used];if(last==root)++e[used];if(first==root)++s[used];unsigned next=a[last]&~used;while(next){unsigned bit=next&-next;next-=bit;self(self,first,__builtin_ctz(bit),used|bit);}};for(int r=0;r<n;++r)dfs(dfs,r,r,1u<<r);I L=0;for(int mask=1;mask<N;++mask)L+=e[mask]*h[F^mask];return{h[F],s[F],L};}
G deletion(const G&a,int skip,int root){std::vector<int>p{root};for(int i=0;i<(int)a.size();++i)if(i!=root&&i!=skip)p.push_back(i);G b(p.size());for(int i=0;i<(int)p.size();++i)for(int j=0;j<(int)p.size();++j)if(a[p[i]]>>p[j]&1)b[i]|=1u<<j;return b;}
G readout(const G&a){G out{2,8,3,5};out.resize(a.size()+3);for(int i=0;i<(int)a.size();++i)for(int j=0;j<(int)a.size();++j)if(a[i]>>j&1)out[i?i+3:0]|=1u<<(j?j+3:0);for(int c=1;c<4;++c)for(int j=1;j<(int)a.size();++j)out[c]|=1u<<(j+3);valid(out);return out;}
std::set<Pair> load_pairs(std::string file){std::ifstream in(file);need(bool(in),"pair file missing");std::set<Pair>r;std::string line;while(std::getline(in,line)){std::istringstream s(line);int B,start;s>>B>>start;need(bool(s),"pair parse");need(r.insert({B,start}).second,"duplicate profile key");}return r;}
