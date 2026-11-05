#include <algorithm>
#include <array>
#include <deque>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using ll=long long;
using State=array<ll,3>;
void need(bool ok,const string& msg){if(!ok)throw runtime_error(msg);}
struct Data{int n,full;vector<int>out;vector<ll>h;vector<array<ll,12>>e,s;};
Data compute(const vector<int>&out){
 int n=out.size(),N=1<<n;Data d{n,N-1,out,vector<ll>(N),vector<array<ll,12>>(N),vector<array<ll,12>>(N)};d.h[0]=1;
 for(int m=1;m<N;++m)for(int v=0;v<n;++v)if(m>>v&1){int z=m^(1<<v);if(!z)d.e[m][v]=d.s[m][v]=1;else for(int u=0;u<n;++u)if(z>>u&1){if(out[u]>>v&1)d.e[m][v]+=d.e[z][u];if(out[v]>>u&1)d.s[m][v]+=d.s[z][u];}d.h[m]+=d.e[m][v];}
 return d;
}
int reachable(const vector<int>&out,int m,int v){int x=1<<v,old=0;while(x!=old){old=x;for(int u=0;u<(int)out.size();++u)if(x>>u&1)x|=out[u]&m;}return x;}
bool strong(const vector<int>&out,int m){if(!m)return false;for(int v=0;v<(int)out.size();++v)if(m>>v&1&&reachable(out,m,v)!=m)return false;return true;}
ll cover(const Data&d,int m,int v){ll z=0;for(int a=m;a;a=(a-1)&m)z+=d.e[a][v]*d.h[m^a];return z;}
State state(const Data&d,int m,int v){return {d.h[m],d.s[m][v],cover(d,m,v)};}
vector<int> bits(int a){vector<int>b;for(;a;a>>=1)b.push_back(a&1);reverse(b.begin(),b.end());return b;}
struct Word{ll g,A,C,q,eps;};
// Direct before/after-apex path partition, independent of the producer recurrence.
Word word(int a){if(!a)return {1,0,0,1,1};auto b=bits(a);int m=b.size();ll g=b.front()+(1-b.back()),C=0;
 for(int i=0;i<m;++i){if(!b[i])C+=1LL<<i;if(i+1<m&&!b[i]&&b[i+1])++g;for(int j=0;j<i;++j)if(b[j]&&!b[i])g+=1LL<<(i-j-1);}
 return {g,a,C,1LL<<m,1-b.back()};}
State transfer(State st,int a){auto[h,s,L]=st;auto w=word(a);return {w.g*h-w.eps*s+(w.A+1)*L,a?h:h-s,(1+w.C)*h+w.q*L};}
State alphabet(State st,int a){if(!a)return transfer(st,a);auto[h,s,L]=st;ll U=h,X=0,Y=0,Qh=h,QL=L,V=0,Ph=0,Ps=0;
 for(int b:bits(a)){ll nU=U+(1-b)*(X+Ph),nX=2*X+b*h,nY=2*Y+b*L,nV=V+(1-b)*Qh;U=nU;X=nX;Y=nY;V=nV;Qh*=2;QL*=2;Ph=b*h;Ps=b*s;}
 return {U-s+Ps+Y+L,h,h+V+QL};}
vector<int> construct(const vector<int>&C,int u,int a){int k=C.size();auto b=bits(a);int m=b.size(),v=k+m;vector<int>o=C;o.resize(v+1);
 for(int i=0;i<m;++i){for(int j=0;j<k;++j)o[k+i]|=1<<j;for(int j=i+1;j<m;++j)o[k+i]|=1<<(k+j);if(b[i])o[v]|=1<<(k+i);else o[k+i]|=1<<v;}
 for(int i=0;i<k;++i)if(i==u)o[i]|=1<<v;else o[v]|=1<<i;
 return o;}
void literal_check(const Data&d){int n=d.n;vector<int>p(n);iota(p.begin(),p.end(),0);ll h=0;array<ll,12>e{},s{},L{};
 do{bool hp=true;for(int i=1;i<n;++i)if(!(d.out[p[i-1]]>>p[i]&1))hp=false;if(hp){++h;++e[p.back()];++s[p.front()];}
 for(int k=1;k<=n;++k){bool ok=true;for(int i=1;i<n;++i)if(i!=k&&!(d.out[p[i-1]]>>p[i]&1))ok=false;if(ok)++L[p[k-1]];}
 }while(next_permutation(p.begin(),p.end()));need(h==d.h.back(),"literal H");for(int v=0;v<n;++v)need(e[v]==d.e.back()[v]&&s[v]==d.s.back()[v]&&L[v]==cover(d,d.full,v),"literal marked coordinates");}
#include <fstream>
#include <unordered_set>
using Profile=array<ll,8>;
using Signature=array<ll,7>;
Profile profile(const Data&d,int v){int m=d.full,del=m^(1<<v);ll L=cover(d,m,v),R=0,B=0;for(int a=m;a;a=(a-1)&m)R+=d.s[a][v]*d.h[m^a];for(int a=del;;a=(a-1)&del){B+=d.h[a]*d.h[del^a];if(!a)break;}return {d.h[m],d.s[m][v],L,d.e[m][v],R,B,d.n==1,d.h[del]};}
vector<Profile> literal_profiles(const vector<int>&o){int n=o.size();vector<int>p(n);iota(p.begin(),p.end(),0);ll h=0;array<ll,12>e{},s{},L{},R{};auto path=[&](const vector<int>&p,int a,int b){for(int i=a+1;i<b;++i)if(!(o[p[i-1]]>>p[i]&1))return false;return true;};
 do{if(path(p,0,n)){++h;++e[p.back()];++s[p.front()];}for(int k=0;k<=n;++k)if(path(p,0,k)&&path(p,k,n)){if(k)++L[p[k-1]];if(k<n)++R[p[k]];}}while(next_permutation(p.begin(),p.end()));
 vector<Profile>ps;for(int v=0;v<n;++v){vector<int>q;for(int u=0;u<n;++u)if(u!=v)q.push_back(u);ll B=0,d=0;do{d+=path(q,0,n-1);for(int k=0;k<n;++k)B+=path(q,0,k)&&path(q,k,n-1);}while(next_permutation(q.begin(),q.end()));ps.push_back({h,s[v],L[v],e[v],R[v],B,n==1,d});}return ps;
}
Signature sig(Profile p){Signature z;copy_n(p.begin(),7,z.begin());return z;}
State triple(Profile p){return {p[0],p[3],p[4]};}
State move_root(State x,Signature y){auto[h,e,R]=x;auto[g,s,L,E,Q,B,delta]=y;return {g*h-s*e+L*R,E*h-delta*e,Q*h+B*R};}
struct Mark{vector<int>out;int root;Profile p;};
vector<int> seam(const vector<int>&A,int a,const vector<int>&G,int b){int n=A.size(),m=G.size();vector<int>o;for(int row:A)o.push_back(row|(((1<<m)-1)<<n));for(int row:G)o.push_back(row<<n);o[a]^=1<<(n+b);o[n+b]|=1<<a;return o;}
struct Hash{size_t operator()(const State&x)const{size_t z=x[0];z=(z*1000003)^x[1];return(z*1000033)^x[2];}};
struct Census{ll states,edges;set<ll>values;};
Census census(int bound,const vector<Signature>&bank,State seed={1,1,1}){unordered_set<State,Hash>seen;seen.reserve(bound*60);seen.insert(seed);vector<State>todo{seed};set<ll>values{seed[0]};ll edges=0;for(size_t i=0;i<todo.size();++i){State x=todo[i];for(auto y:bank){auto z=move_root(x,y);need(z[0]>=max(x[0],x[2])&&z[2]>x[2],"state monotonicity");if(z[0]>bound)continue;++edges;values.insert(z[0]);if(z[2]<=bound&&seen.insert(z).second)todo.push_back(z);}}return{(ll)seen.size(),edges,values};}
int main(){try{
 map<Signature,Mark>inventory;vector<Mark>small;ll graphs=0,roots=0;set<State>equal_seeds;
 for(int n=1;n<=5;++n){int E=n*(n-1)/2;for(unsigned code=0;code<(1u<<E);++code){vector<int>o(n);int t=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++t)if(code>>t&1)o[i]|=1<<j;else o[j]|=1<<i;auto d=compute(o);auto brute=literal_profiles(o);++graphs;for(int v=0;v<n;++v){auto p=profile(d,v);need(p==brute[v],"literal full profile inventory");++roots;auto y=sig(p);if(inventory.count(y))need(inventory.at(y).p[7]==p[7],"finite bank deletion signature");else inventory.emplace(y,Mark{o,v,p});if(n<=3)small.push_back({o,v,p});auto st=triple(p);if(st[0]>=5&&st[2]==st[0]&&st[0]-st[1]>=2)equal_seeds.insert(st);}}}
 need(graphs==1099&&roots==5405&&inventory.size()==65,"complete inventory counts");ofstream out("inventory.tsv");vector<Signature>bank;for(auto&[y,m]:inventory){bank.push_back(y);for(ll k:y)out<<k<<'\t';out<<m.p[7]<<'\t'<<m.root<<'\t'<<m.out.size();for(int row:m.out)out<<'\t'<<row;out<<'\n';}out.close();
 ll pairs=0,literal_pairs=0,old_checks=0;for(auto&A:small)for(auto&G:small){auto o=seam(A.out,A.root,G.out,G.root);auto d=compute(o);auto p=profile(d,A.out.size()+G.root);need(triple(p)==move_root(triple(A.p),sig(G.p)),"universal join transfer");need(p[7]==A.p[0]*G.p[7],"new root deletion product");auto brute=literal_profiles(o);need(p==brute[A.out.size()+G.root],"literal physical output profile");++literal_pairs;++pairs;}
 vector<Mark>controls=small;for(vector<int>o: {vector<int>{2,12,9,17,7},vector<int>{24,9,3,20,6}}){auto d=compute(o);for(int v=0;v<(int)o.size();++v)controls.push_back({o,v,profile(d,v)});}
 for(auto&A:controls)for(auto&[y,G]:inventory){auto o=seam(A.out,A.root,G.out,G.root);auto d=compute(o);auto p=profile(d,A.out.size()+G.root);need(triple(p)==move_root(triple(A.p),y),"full bank graph transfer");need(p[7]==A.p[0]*G.p[7],"full bank deletion product");++pairs;}
 for(auto&A:controls){auto[h,e,R]=triple(A.p);ll del=A.p[7];auto o=seam(A.out,A.root,{0},0);auto p=profile(compute(o),A.root);need(triple(p)==State{h+R-e,del,2*R-del}&&p[7]==del,"stationary root K");o=A.out;for(int&row:o)row|=1<<A.out.size();o.push_back(0);p=profile(compute(o),A.root);need(triple(p)==State{h,0,2*R}&&p[7]==del,"stationary sink D");old_checks+=2;}
 cout<<"graphs="<<graphs<<" rooted_gadgets="<<roots<<" matrices="<<inventory.size()<<" physical_join_checks="<<pairs<<" literal_join_checks="<<literal_pairs<<" stationary_checks="<<old_checks<<"\n";
 array<array<ll,3>,3>M={array<ll,3>{1,1,3},array<ll,3>{1,3,5},array<ll,3>{3,5,9}};ll det=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);need(det==-4,"minimal continuation determinant");
 auto hostile=profile(compute({28,17,26,2,8}),3);need(triple(hostile)==State{9,3,9}&&hostile[7]==3,"gcd hostile");for(auto y:bank)for(ll k:move_root({9,3,9},y))need(k%3==0,"content preservation");
 auto H13=profile(compute({26,16,3,6,12}),2);need(triple(H13)==State{13,2,18},"H13 root witness");auto primitive=seam({26,16,3,6,12},2,{0,1},0);auto dp=compute(primitive);auto pr=profile(dp,5);need(strong(primitive,dp.full)&&triple(pr)==State{49,13,49},"primitive equal response hostile");equal_seeds.insert({49,13,49});
 for(auto st:equal_seeds){auto c=census(3*st[0]+4,bank,st);need(c.values==set<ll>{st[0],2*st[0]-st[1],3*st[0]},"equal-response bounded support");}cout<<"minimal_determinant="<<det<<" equal_response_seed_checks="<<equal_seeds.size()<<"\n";
 for(int X:{1,2,3,5,13,1000,20000}){auto c=census(X,bank);cout<<"X="<<X<<" states="<<c.states<<" values="<<c.values.size()<<" retained_edges="<<c.edges<<" missing=";for(int h=1;h<=X;h+=2)if(!c.values.count(h))cout<<h<<",";cout<<"\n";if(X==1000)need(c.states==12711&&c.values.size()==493,"1000 grammar census");if(X==20000)need(c.states==1085407&&c.values.size()==9993,"20000 grammar census");}
 struct Added{int target;vector<int>G;int b;vector<pair<vector<int>,int>>prefix;};
 vector<Added> added={
 {67,{52,17,2,7,12,30},3,{}},
 {163,{44,1,2,6,15,30},3,{{{0,1},0}}},
 {173,{56,1,3,6,14,30},1,{{{0},0}}},
 {247,{48,1,3,7,14,30},2,{{{0,1},0}}},
 {511,{42,16,3,6,13,30},3,{{{0,1,3},2}}}};
 for(auto&add:added){auto d=compute(add.G);auto p=profile(d,add.b);need(p==literal_profiles(add.G)[add.b],"literal order-six added profile");inventory.emplace(sig(p),Mark{add.G,add.b,p});vector<int>o={0};int r=0;State st={1,1,1};auto steps=add.prefix;steps.push_back({add.G,add.b});for(auto&[G,b]:steps){int nr=o.size()+b;auto gp=profile(compute(G),b);o=seam(o,r,G,b);r=nr;st=move_root(st,sig(gp));}auto dp=compute(o);auto final=profile(dp,r);need(triple(final)==st&&final[0]==add.target,"physical extra-gadget witness");cout<<"extra_target="<<add.target<<" n="<<o.size()<<" root="<<r<<" state="<<st[0]<<","<<st[1]<<","<<st[2]<<" strong="<<strong(o,dp.full)<<"\n";}
 need(inventory.size()==70,"extended inventory distinctness");vector<Signature>bank70;for(auto&[y,G]:inventory)bank70.push_back(y);auto ext=census(1000,bank70);need(ext.states==12863&&ext.values.size()==498,"extended1000 census");for(int h=1;h<=1000;h+=2)need(ext.values.count(h)==(h!=7&&h!=21),"extended exact misses");cout<<"extended_X=1000 matrices=70 states="<<ext.states<<" values="<<ext.values.size()<<" missing=7,21\n";
 for(auto st:equal_seeds){auto c=census(3*st[0]+4,bank70,st);need(c.values==set<ll>{st[0],2*st[0]-st[1],3*st[0]},"extended equal-response bounded support");}
 auto single_sig=sig(profile(compute({0}),0));auto terminal=census(5,{single_sig});need(terminal.states==3&&terminal.values==set<ll>{1,3,5},"R clipping terminal boundary");cout<<"singleton_bank_terminal_control: X=5 retained_states=3 values=1,3,5; terminal_state=5,2,6\n";
 cout<<"PASS\n";
}catch(exception&e){cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
