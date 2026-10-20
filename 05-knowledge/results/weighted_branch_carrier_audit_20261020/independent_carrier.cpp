// Independent exhaustive switching-class audit. No repository mathematical imports.
// A state is the flip vector on 1 <= i < j < n after switching vertex 0 to dominate.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Perm = vector<int>;
struct Action { Perm p; vector<uint32_t> masks; uint32_t c=0; };
struct Rat {long long a=0,b=1; Rat(long long x=0,long long y=1):a(x),b(y){long long g=gcd(a,b);a/=g;b/=g;} Rat operator+(Rat z)const{return Rat(a*z.b+z.a*b,b*z.b);} bool operator==(Rat z)const{return a==z.a&&b==z.b;} string str()const{return to_string(a)+"/"+to_string(b);} };
Perm mul(const Perm&p,const Perm&q){Perm r(p.size());for(int i=0;i<(int)p.size();++i)r[i]=p[q[i]];return r;}
Perm inv(const Perm&p){Perm r(p.size());for(int i=0;i<(int)p.size();++i)r[p[i]]=i;return r;}
int order(const Perm&p){int o=1;vector<bool>seen(p.size());for(int i=0;i<(int)p.size();++i)if(!seen[i]){int j=i,l=0;do{seen[j]=true;j=p[j];++l;}while(j!=i);o=lcm(o,l);}return o;}
uint32_t act(uint32_t s,const Action&a){uint32_t r=a.c;for(int k=0;k<(int)a.masks.size();++k)r^=(uint32_t)(__builtin_parity(s&a.masks[k]))<<k;return r;}
void check(bool b,string why){if(!b)throw runtime_error(why);}
string plist(Perm p){string s="[";for(int i=0;i<(int)p.size();++i){if(i)s+=",";s+=to_string(p[i]);}return s+"]";}
int main(){cout<<"{\n\"method\":\"Exhaustive normalized affine switching states, all S_n relabellings, exact rational centralizer weights, centralizer-constrained decorated flips\",\n\"cases\":[\n";bool firstN=true;
for(int n: {2,4,6,8}){
vector<pair<int,int>>edges;int idx[8][8];for(auto& row:idx)fill(begin(row),end(row),-1);for(int i=1;i<n;++i)for(int j=i+1;j<n;++j){idx[i][j]=idx[j][i]=edges.size();edges.push_back({i,j});}int d=edges.size();
Perm p(n);iota(p.begin(),p.end(),0);vector<Action> acts;map<Perm,int> permidx;
do{Action a;a.p=p;Perm q=inv(p);int c=q[0];for(auto [i,j]:edges){int u=q[i],v=q[j];uint32_t m=0;for(auto [h,k]:vector<pair<int,int>>{{u,v},{c,u},{c,v}})if(idx[h][k]>=0)m^=1u<<idx[h][k];a.masks.push_back(m);int constant=(u>v)^(c>u)^(c>v);if(constant)a.c|=1u<<(a.masks.size()-1);}permidx[p]=acts.size();acts.push_back(a);}while(next_permutation(p.begin(),p.end()));
vector<int>orbit(1u<<d,-1);vector<uint32_t>reps;vector<vector<int>>groups;
for(uint32_t s=0;s<(1u<<d);++s)if(orbit[s]<0){int num=reps.size();reps.push_back(s);vector<int>G;for(int i=0;i<(int)acts.size();++i){uint32_t z=act(s,acts[i]);check(orbit[z]<0||orbit[z]==num,"overlapping orbits");orbit[z]=num;if(z==s)G.push_back(i);}groups.push_back(G);}
Rat total;map<string,int>weights;int evens=0,twocount=0,invcount=0,fixed=0,coarsefixed=0;map<string,int>strata_counts;map<string,Rat>strata_mass;
if(!firstN){cout<<",\n";}firstN=false;cout<<"{\"n\":"<<n<<",\"states\":"<<(1u<<d)<<",\"permutations\":"<<acts.size()<<",\"classes\":"<<reps.size()<<",\"rows\":[\n";bool firstrow=true;
for(int oi=0;oi<(int)reps.size();++oi){uint32_t s=reps[oi];auto G=groups[oi];if(G.size()%2)continue;++evens;int odd=0;map<int,int>orders;for(int h:G){int o=order(acts[h].p);++orders[o];odd+=(o%2);}
Rat contribution(G.size()-odd,G.size());total=total+contribution;set<int>todo;for(int h:G){int o=order(acts[h].p);if(o>1&&(o&(o-1))==0)todo.insert(h);}Rat sum;vector<int>xx,central_sizes,central_odd;vector<string>ww;vector<int>flags,coarseflags;vector<Perm>witnesses,flipmarks; vector<uint32_t>flipstates;
while(!todo.empty()){int x=*todo.begin();xx.push_back(x);auto px=acts[x].p;set<int>conj;for(int h:G){Perm z=mul(mul(acts[h].p,px),inv(acts[h].p));conj.insert(permidx.at(z));}for(int y:conj)todo.erase(y);int cz=0,co=0;for(int h:G)if(mul(acts[h].p,px)==mul(px,acts[h].p)){++cz;co+=(order(acts[h].p)%2);}Rat w(co,cz);sum=sum+w;weights[w.str()]++;ww.push_back(w.str());central_sizes.push_back(cz);central_odd.push_back(co);++twocount;
if(order(px)==2){++invcount; // Construct matching flip in full edge coordinates, then normalize root.
uint32_t t=s;auto matching=[&](int i,int j){return px[i]==j;};for(int k=0;k<d;++k){auto[i,j]=edges[k];if(matching(i,j)^matching(0,i)^matching(0,j))t^=1u<<k;}
check(act(t,acts[x])==t,"flip lost marked involution");bool cf=orbit[t]==oi;bool df=false;Perm wit;for(auto&a:acts)if(mul(a.p,px)==mul(px,a.p)&&act(s,a)==t){df=true;wit=a.p;break;}check(!df||cf,"decorated fixed not coarse fixed");flags.push_back(df);coarseflags.push_back(cf);witnesses.push_back(wit);
int target=orbit[t], transporter=-1;for(int h=0;h<(int)acts.size();++h)if(act(t,acts[h])==reps[target]){transporter=h;break;}check(transporter>=0,"no flip transporter");Perm pg=mul(mul(acts[transporter].p,px),inv(acts[transporter].p));int mark=permidx.at(pg);for(int h:groups[target])mark=min(mark,permidx.at(mul(mul(acts[h].p,pg),inv(acts[h].p))));flipstates.push_back(reps[target]);flipmarks.push_back(acts[mark].p);fixed+=df;coarsefixed+=cf;check([&](){uint32_t tt=t;for(int k=0;k<d;++k){auto[i,j]=edges[k];if(matching(i,j)^matching(0,i)^matching(0,j))tt^=1u<<k;}return tt==s;}(),"flip not involution");}}
check(sum==contribution,"centralizer identity failed");int syl=1;while(G.size()%(2*syl)==0)syl*=2;string ty;if(syl==2)ty="C2";else if(syl==4)ty=orders.count(4)?"C4":"V4";else if(syl==8){ty=orders.count(8)?"C8":"D8";if(ty=="D8"){bool found=false;for(int z:G)if(order(acts[z].p)==2){vector<int>CG;for(int h:G)if(mul(acts[h].p,acts[z].p)==mul(acts[z].p,acts[h].p))CG.push_back(h);if(CG.size()!=8)continue;int i2=0,i4=0;for(int h:CG){i2+=order(acts[h].p)==2;i4+=order(acts[h].p)==4;}if(i2==5&&i4==2)found=true;}check(found,"D8 Sylow not certified");}}else ty="other";strata_counts[ty]++;strata_mass[ty]=strata_mass[ty]+contribution;
if(!firstrow){cout<<",\n";}firstrow=false;cout<<"{\"state\":"<<s<<",\"aut_order\":"<<G.size()<<",\"sylow_order\":"<<syl<<",\"sylow_type\":\""<<ty<<"\",\"element_orders\":{";bool a=true;for(auto[k,v]:orders){if(!a)cout<<",";a=false;cout<<"\""<<k<<"\":"<<v;}cout<<"},\"contribution\":\""<<contribution.str()<<"\",\"two_element_classes\":[";
for(int k=0;k<(int)xx.size();++k){if(k)cout<<",";cout<<"{\"element\":"<<plist(acts[xx[k]].p)<<",\"order\":"<<order(acts[xx[k]].p)<<",\"centralizer_size\":"<<central_sizes[k]<<",\"centralizer_odd_elements\":"<<central_odd[k]<<",\"weight\":\""<<ww[k]<<"\"}";}cout<<"],\"decorated_fixed\":[";for(int k=0;k<(int)flags.size();++k){if(k)cout<<",";cout<<flags[k];}cout<<"],\"coarse_fixed\":[";for(int k=0;k<(int)coarseflags.size();++k){if(k)cout<<",";cout<<coarseflags[k];}cout<<"],\"decorated_fixed_witnesses\":[";for(int k=0;k<(int)witnesses.size();++k){if(k)cout<<",";cout<<plist(witnesses[k]);}cout<<"],\"flip_targets\":[";for(int k=0;k<(int)flipstates.size();++k){if(k)cout<<",";cout<<"{\"state\":"<<flipstates[k]<<",\"involution\":"<<plist(flipmarks[k])<<"}";}cout<<"]}";
}
cout<<"\n],\"even_aut_classes\":"<<evens<<",\"nontrivial_two_element_classes\":"<<twocount<<",\"involution_classes\":"<<invcount<<",\"decorated_fixed\":"<<fixed<<",\"coarse_fixed\":"<<coarsefixed<<",\"weight_counts\":{";bool a=true;for(auto[k,v]:weights){if(!a)cout<<",";a=false;cout<<"\""<<k<<"\":"<<v;}cout<<"},\"N_lev\":\""<<total.str()<<"\",\"strata\":{";a=true;for(auto[k,v]:strata_counts){if(!a)cout<<",";a=false;cout<<"\""<<k<<"\":{\"classes\":"<<v<<",\"mass\":\""<<strata_mass[k].str()<<"\"}";}cout<<"}}";
}
cout<<"\n]}\n";}
