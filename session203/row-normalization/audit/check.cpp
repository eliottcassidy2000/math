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
int main(){G seeds[2]={{10,20,33,38,13,19},{10,12,17,52,35,7}};Row common{111,9,161,13};int words=0,right_graphs=0,right_roots=0,literal_pairs=0,small_padding=0;std::vector<I>seed_details;
 for(int i=0;i<2;++i){Data q=profile(seeds[i],0);need(strong(seeds[i]),"strong source");need(row(q)==common&&q.s==I(6+i),"initial common row");need(q.h==I(45-2*i)&&q.L==I(72+3*i),"initial h/s/L");check_edge(seeds[i],0,q);seed_details.insert(seed_details.end(),{q.h,q.s,q.e,q.L,q.R,q.d,q.beta});}
 for(int len=0;len<=4;++len)for(int code=0;code<(1<<len);++code){Row states[2];for(int i=0;i<2;++i){G a=seeds[i];Row expected=common;for(int j=0;j<len;++j){a=pad(a,0,code>>j&1);expected=step(expected,code>>j&1);}Data q=profile(a,0);need(row(q)==expected&&q.s==I(6+i),"all padding states");check_output(a,0,q);states[i]=row(q);
 if(code==(1<<len)-1){need(strong(a)&&strong(child(a,0)),"Z-only strong graph");if(len)need(q.h+q.L-q.s==148*(I(1)<<len)-33,"Z ray including first boundary");}
 if(len<=2)check_edge(a,0,q);
 }need(states[0]==states[1],"common word row");++words;}
 for(int n=1;n<=4;++n)for(int code=0;code<(1<<(n*(n-1)/2));++code){G u(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(code>>bit++&1)u[i]|=1u<<j;else u[j]|=1u<<i;}++right_graphs;
 for(int root=0;root<n;++root){Data w=profile(u,root);++right_roots;
 for(int i=0;i<2;++i){G a=amalgam(seeds[i],0,u,root);Data q=profile(a,0);std::int64_t B=13*w.h+148*w.L-4*w.s-46*w.d;need(B>0,"positive right base");Row expected{I(B),n==1?9:13*w.e,148*w.beta+13*w.R,13*w.d};
 need(row(q)==expected&&q.s==I(6+i)*w.d,"arbitrary right input and endpoint scale");check_output(a,0,q);++literal_pairs;
 G left=edge(a,0),right=amalgam(edge(seeds[i],0),0,u,root);need(left==right,"literal edge/amalgamation associativity");
 }
 if(n>=2){for(bool z:{false,true}){G a=pad(u,root,z);Data q=profile(a,root);need(row(q)==step(row(w),z)&&q.s==w.s,"generic O/Z state law");++small_padding;}}
 }}
 // The excluded singleton Z case is an explicit hostile control.
 Data one=profile(G{0},0),oz=profile(pad(G{0},0,true),0);need(one.s==1&&oz.s==0,"singleton Z exception");need(row(oz)!=step(row(one),true),"singleton row update is excluded");
 // Frozen first missing-base witnesses, built from their original seeds.
 G a={2,68,9,19,39,15,61};a=pad(a,6,false);a=pad(a,6,true);Data qa=profile(a,6);need(qa.h==485&&qa.s==16&&qa.L==520&&row(qa)[0]==989,"base999 residue1");need(hp(child(a,6))==2999,"target2999");
 G b={2,20,9,3,13};for(bool z:{false,false,true,false,true})b=pad(b,3,z);Data qb=profile(b,3);need(qb.h==341&&qb.s==2&&qb.L==660&&row(qb)[0]==999,"base999 residue2");need(hp(child(b,3))==3001,"target3001");
 std::cout<<"{\"status\":\"PASS\",\"padding_words\":"<<words<<",\"right_graphs\":"<<right_graphs<<",\"right_root_choices\":"<<right_roots<<",\"right_literal_pairs\":"<<literal_pairs<<",\"generic_small_padding_cases\":"<<small_padding<<",\"profile_calls\":"<<profile_calls<<",\"literal_full_count_calls\":"<<full_calls<<",\"directed_prefix_visits\":"<<visits<<",\"seed_profiles_h_s_e_L_R_d_beta\":[";for(int i=0;i<(int)seed_details.size();++i)std::cout<<(i?",":"")<<seed_details[i];std::cout<<"],\"base999_targets\":[2999,3001],\"unconditional_coverage_claim\":false}\n";
}
