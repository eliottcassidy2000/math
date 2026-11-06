// Independent literal directed-prefix counting and physical root amalgamation.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using G=std::vector<unsigned>; using I=std::uint64_t;
I visits=0;int profile_calls=0,full_calls=0;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void valid(const G&a){int n=a.size();need(n<20,"order");for(int i=0;i<n;++i){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"row");for(int j=0;j<i;++j)need(((a[i]>>j&1)+(a[j]>>i&1))==1,"orientation");}}
G induced(const G&a,const std::vector<int>&p){G b(p.size());for(int i=0;i<(int)p.size();++i)for(int j=0;j<(int)p.size();++j)if(a[p[i]]>>p[j]&1)b[i]|=1u<<j;return b;}
struct Data{I h,s,e,L,R,d,beta;};
Data profile(const G&a,int r){valid(a);++profile_calls;int n=a.size(),full=(1<<n)-1;std::vector<I>h(1<<n),s(1<<n),e(1<<n);h[0]=1;
 auto dfs=[&](auto&&self,int first,int last,unsigned mask)->void{++visits;++h[mask];if(first==r)++s[mask];if(last==r)++e[mask];unsigned todo=a[last]&~mask;while(todo){int v=__builtin_ctz(todo);todo&=todo-1;self(self,first,v,mask|(1u<<v));}};
 for(int first=0;first<n;++first)dfs(dfs,first,first,1u<<first);
 I L=0,R=0,beta=0;for(int mask=0;mask<=full;++mask){L+=e[mask]*h[full^mask];R+=s[mask]*h[full^mask];}
 int body=full^(1<<r);for(int mask=body;;mask=(mask-1)&body){beta+=h[mask]*h[body^mask];if(!mask)break;}
 return{h[full],s[full],e[full],L,R,h[body],beta};}
I hp(const G&a){valid(a);++full_calls;if(a.empty())return 1;I total=0;unsigned full=(1u<<a.size())-1;auto dfs=[&](auto&&self,int v,unsigned mask)->void{++visits;if(mask==full){++total;return;}unsigned next=a[v]&~mask;while(next){int w=__builtin_ctz(next);next&=next-1;self(self,w,mask|(1u<<w));}};for(int v=0;v<(int)a.size();++v)dfs(dfs,v,1u<<v);return total;}
// Explicit identification: old labels first, then right labels except its root.
G amalgam(const G&a,int r,const G&b,int t){int n=a.size();G c(n+b.size()-1);std::vector<int>map(b.size());map[t]=r;int next=n;for(int j=0;j<(int)b.size();++j)if(j!=t)map[j]=next++;
 for(int i=0;i<n;++i)c[i]=a[i];
 for(int i=0;i<(int)b.size();++i)for(int j=0;j<(int)b.size();++j)if(b[i]>>j&1)c[map[i]]|=1u<<map[j];
 for(int i=0;i<n;++i)if(i!=r)for(int j=0;j<(int)b.size();++j)if(j!=t)c[i]|=1u<<map[j];
 valid(c);return c;}
G pad(G a,int r,bool z){int n=a.size();a.resize(n+1);for(int i=0;i<n;++i){if(z&&i==r)a[n]|=1u<<i;else a[i]|=1u<<n;}valid(a);return a;}
G edge(const G&a,int r){return amalgam(G{2,0},0,a,r);}
G child(const G&a,int r){return amalgam(G{2,8,3,5},0,a,r);}
bool strong(const G&a){for(int i=0;i<(int)a.size();++i){unsigned reach=1u<<i,old;do{old=reach;for(int j=0;j<(int)a.size();++j)if(reach>>j&1)reach|=a[j];}while(old!=reach);if(reach!=(1u<<a.size())-1)return false;}return true;}
using Row=std::array<I,4>;
Row row(const Data&q){return{q.h+q.L-q.s,q.e,q.R+q.beta,q.d};}
Row step(Row q,bool z){return z?Row{q[0]+q[2]-q[1],q[3],2*q[2]-q[3],q[3]}:Row{q[0],0,2*q[2],q[3]};}
void check_edge(const G&a,int r,const Data&q){G b=edge(a,r);Data p=profile(b,0);need(Row{p.h,p.e,p.R,p.d}==row(q),"literal edge retained profile");std::vector<int>old{0};for(int j=2;j<(int)b.size();++j)old.push_back(j);G back=induced(b,old);Data recovered=profile(back,0);need(recovered.s==q.s&&recovered.h==q.h,"marked deletion sidecar");}
void check_output(const G&a,int r,const Data&q){G c=child(a,r);need(hp(c)==3*(q.h+q.L-q.s)+2*q.s,"actual four-core readout");}

G inflate(const G&a,int v,const G&x){std::vector<std::vector<int>>slots(a.size());int n=0;for(int i=0;i<(int)a.size();++i)for(int j=0;j<(i==v?(int)x.size():1);++j)slots[i].push_back(n++);G out(n);for(int i=0;i<(int)a.size();++i)for(int j=0;j<(int)a.size();++j)if(a[i]>>j&1)for(int s:slots[i])for(int t:slots[j])out[s]|=1u<<t;for(int i=0;i<(int)x.size();++i)for(int j=0;j<(int)x.size();++j)if(x[i]>>j&1)out[slots[v][i]]|=1u<<slots[v][j];valid(out);return out;}
std::array<I,5> covers(const G&a){std::array<I,5>q{};int n=a.size();if(!n){q[0]=1;return q;}std::vector<int>p(n);std::iota(p.begin(),p.end(),0);do{for(unsigned cuts=0;cuts<(1u<<(n-1));++cuts){int bins=1+__builtin_popcount(cuts);if(bins>4)continue;bool ok=true;for(int i=0;i<n-1;++i)if(!(cuts>>i&1)&&!(a[p[i]]>>p[i+1]&1)){ok=false;break;}if(ok)++q[bins];}}while(std::next_permutation(p.begin(),p.end()));return q;}
I dot(std::array<I,5>p,std::array<I,5>c){I s=0;for(int i=0;i<5;++i)s+=p[i]*c[i];return s;}
G tournament(int n,unsigned code){G a(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(code>>bit++&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
int main(){G seeds[2]={{102,0,10,83,135,30,182,47},{142,0,34,22,71,155,47,94}};int replace[2]={5,3};std::vector<G>inputs;for(int n=0;n<=4;++n)for(unsigned code=0;code<(1u<<(n*(n-1)/2));++code)inputs.push_back(tournament(n,code));need(inputs.size()==76,"complete labeled module universe");for(unsigned code:{0u,1023u,345u,765u,500u,817u,101u,236u})inputs.push_back(tournament(5,code));inputs.push_back(tournament(6,32767));inputs.push_back(G{58,60,57,16,32,8});inputs.push_back(tournament(6,20321));inputs.push_back(tournament(6,11755));int graphs=0,padded=0,compiled=0,nonempty_strict=0;std::vector<I>boundary;
 for(const G&x:inputs){auto p=covers(x);int n=x.size();for(int j=1;j<=4;++j){I fact=1;for(int k=2;k<=j;++k)fact*=k;need(p[j]%fact==0,"labeled bin factorial divisibility");}I B=dot(p,{93,253,204,48,0}),V=dot(p,{300,1022,1178,520,72}),d=dot(p,{11,29,21,4,0}),s=dot(p,{5,16,15,4,0}),h=dot(p,{37,111,95,24,0}),L=dot(p,{61,158,124,28,0}),R=dot(p,{100,382,486,236,36}),beta=dot(p,{200,640,692,284,36});Data values[2];for(int i=0;i<2;++i){G a=inflate(seeds[i],replace[i],x);Data q=profile(a,0);I delta=i*2*(p[1]+p[2]);need(q.h==h+delta&&q.s==s+delta&&q.e==0&&q.L==L&&q.R==R&&q.d==d&&q.beta==beta,"all seven literal profile coordinates");need(row(q)==Row{B,0,V,d},"common padding row");need(a[1]==0&&!strong(a),"retained global sink");values[i]=q;++graphs;if(n<=3){check_output(a,0,q);++compiled;}
 if(n<=2){for(int len=0;len<=2;++len)for(int code=0;code<(1<<len);++code){G b=a;Row expected=row(q);for(int k=0;k<len;++k){bool z=code>>k&1;b=pad(b,0,z);expected=step(expected,z);}Data real=profile(b,0);need(row(real)==expected&&real.s==q.s,"actual common padding");++padded;}}
 }need(values[1].h-values[0].h==2*(p[1]+p[2])&&values[1].s-values[0].s==2*(p[1]+p[2]),"endpoint and Hamilton difference");if(n){need(p[1]>0,"nonempty gap positive");++nonempty_strict;}else boundary={values[0].h,values[0].s,values[0].e,values[0].L,values[0].R,values[0].d,values[0].beta};
 }
 std::cout<<"{\"status\":\"PASS\",\"complete_labeled_modules_through4\":76,\"extra_order5_modules\":8,\"extra_order6_modules\":4,\"literal_context_graphs\":"<<graphs<<",\"literal_padding_profiles\":"<<padded<<",\"literal_four_core_children\":"<<compiled<<",\"nonempty_positive_gap_modules\":"<<nonempty_strict<<",\"profile_calls\":"<<profile_calls<<",\"directed_prefix_visits\":"<<visits<<",\"empty_profile_h_s_e_L_R_d_beta\":[";for(int i=0;i<(int)boundary.size();++i)std::cout<<(i?",":"")<<boundary[i];std::cout<<"],\"coverage_claim\":false}\n";
}
