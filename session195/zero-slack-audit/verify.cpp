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
using Four=array<ll,4>;
Four four(Profile p){return {p[0],p[3],p[4],p[7]};}
ll slack(Four z){return z[2]+z[1]-z[0]-z[3];}
struct FourHash{size_t operator()(Four x)const{size_t z=0;for(ll v:x)z=z*1000033^v;return z;}};
bool path(const vector<int>&o,const vector<int>&p){for(int i=1;i<(int)p.size();++i)if(!(o[p[i-1]]>>p[i]&1))return false;return true;}
ll injection_check(const Data&d,int x,int y){vector<int>w;for(int i=0;i<d.n;++i)if(i!=x&&i!=y)w.push_back(i);set<pair<vector<int>,vector<int>>>images;ll words=0;
 do{if(!path(d.out,w))continue;++words;int cut=w.size(),a=y,b=x;for(int i=0;i<(int)w.size();++i)if((d.out[x]>>w[i]&1)||(d.out[y]>>w[i]&1)){cut=i;a=(d.out[x]>>w[i]&1)?x:y;b=(a==x)?y:x;break;}vector<int>P(w.begin(),w.begin()+cut);P.push_back(b);vector<int>Q{a};Q.insert(Q.end(),w.begin()+cut,w.end());need(path(d.out,P)&&path(d.out,Q),"endpoint extension paths");vector<int>back=P;back.pop_back();back.insert(back.end(),Q.begin()+1,Q.end());need(back==w,"endpoint extension recovery");need(images.insert({P,Q}).second,"endpoint extension injection");}while(next_permutation(w.begin(),w.end()));need(words==d.h[d.full^(1<<x)^(1<<y)],"literal residual Hamilton count");return words;
}
Four update(Four x,int op,const vector<Profile>&bank){auto[h,e,R,d]=x;if(op==0)return{h+R-e,d,2*R-d,d};if(op==1)return{h,0,2*R,d};auto p=bank.at(op-2);auto z=move_root({h,e,R},sig(p));return {z[0],z[1],z[2],p[7]*h};}
void envelope_trace(const vector<int>&word,Four final,const vector<Profile>&bank,ll C0,int D){vector<Four>st{{1,1,1,1}};for(int op:word)st.push_back(update(st.back(),op,bank));need(st.back()==final&&slack(final)==0&&final[3]<=D,"envelope trace target");int last=-1;for(int i=0;i<(int)word.size();++i)if(word[i]>=2)last=i;int start=last+1;Four z=st[start];
 if(last>=0){Four input=st[last];need(input[0]<=D&&z[3]==bank[word[last]-2][7]*input[0],"last switch bounded input");int last_non_D=-1;for(int i=0;i<last;++i)if(word[i]!=1)last_non_D=i;Four base=st[last_non_D+1];int k=last-(last_non_D+1);need(base[0]==input[0]&&base[2]<=C0*input[0]&&input[2]==base[2]*(1LL<<k),"bounded prefix dyadic coordinate");if(k)need(input[1]==0,"D tail erases endpoint");}
 for(int i=start;i<(int)st.size();++i)need(slack(st[i])==0,"zero suffix cannot erase positive slack");
 int sinks=0,where=-1;for(int i=start;i<(int)word.size();++i){need(word[i]<2,"last root switch classification");if(word[i]==1){++sinks;where=i;auto q=st[i];need(q[0]==q[1]&&q[0]==q[2]&&q[0]==q[3],"D requires all equal");}}need(sinks<=1,"only one zero-slack suffix D");
 if(sinks){int j=word.size()-where-1;ll d=final[3];need(final[0]==(j?d+d*(1LL<<j):d),"D suffix geometric value");}
 else {int j=word.size()-start;need(final[0]==(j?z[3]+(1LL<<j)*(z[0]-z[1]):z[0]),"K suffix geometric formula");}
}
int main(){try{
 ll graphs=0,rooted=0,zeroroots=0,unitzero=0,boundaries=0,injections=0,injectwords=0;set<ll>unitvalues;
 for(int n=1;n<=6;++n){int E=n*(n-1)/2;for(unsigned code=0;code<(1u<<E);++code){vector<int>o(n);int t=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++t)if(code>>t&1)o[i]|=1<<j;else o[j]|=1<<i;auto d=compute(o);++graphs;vector<Profile>literal;if(n<=5)literal=literal_profiles(o);
 for(int r=0;r<n;++r){auto p=profile(d,r);if(n<=5)need(p==literal[r],"literal rooted profile");Four z=four(p);ll f=slack(z),sum=0;int S=d.full^(1<<r),X=o[r];for(int a=S;a;a=(a-1)&S){int b=S^a;if(!b)continue;ll ends=0,starts=0;for(int v=0;v<n;++v)if(X>>v&1){ends+=d.e[a][v];starts+=d.s[b][v];}sum+=ends*starts;}need(f==sum,"exact ordered boundary slack sum");need((f==0)==(__builtin_popcount((unsigned)X)<=1),"zero slack iff root outdegree");++rooted;zeroroots+=(f==0);
 if(z[3]==1&&f==0){++unitzero;unitvalues.insert(z[0]);ll k=z[0]-1;need(z[0]==1||(k>=2&&(k&(k-1))==0),"unit deletion zero slack spectrum");}}
 for(int x=0;x<n;++x)for(int y=x+1;y<n;++y){ll B=0;for(int a=1;a<d.full;++a){int b=d.full^a;B+=d.e[a][x]*d.s[b][y]+d.e[a][y]*d.s[b][x];}need(B>=d.h[d.full^(1<<x)^(1<<y)],"two-endpoint positive bound");++boundaries;if(n<=5){injectwords+=injection_check(d,x,y);++injections;}}
 }}
 cout<<"all_labelled_orders=1..6 graphs="<<graphs<<" rooted_slack_checks="<<rooted<<" zero_slack_roots="<<zeroroots<<" unit_deletion_zero_roots="<<unitzero<<" positive_boundary_checks="<<boundaries<<" literal_injections="<<injections<<" injection_words="<<injectwords<<"\nunit_deletion_values=";for(ll h:unitvalues)cout<<h<<",";cout<<"\n";
 vector<int>host={10,12,9,16,7};auto h15=compute(host);need(strong(host,h15.full)&&h15.h.back()==15,"H15 positive witness");for(int r=0;r<5;++r){auto z=four(profile(h15,r));need(!(z[3]==1&&slack(z)==0),"H15 normalization hostile");}
 // Every intrinsic unit-deletion value has the claimed explicit strong fan.
 for(int m=1;m<=8;++m){int n=m+2,r=n-1;vector<int>o(n);for(int i=0;i<n-1;++i)for(int j=i+1;j<n-1;++j)o[i]|=1<<j;for(int i=0;i<n-1;++i)if(i==0)o[r]|=1<<i;else o[i]|=1<<r;auto d=compute(o);need(strong(o,d.full)&&four(profile(d,r))==Four{(1LL<<m)+1,1,(1LL<<m)+1,1},"strong dyadic fan");}
 vector<Mark>G={{ {0},0,{}},{{0,1},0,{}},{{0,1},1,{}},{{4,1,2},0,{}},{{8,1,3,6},1,{}},{{0,1,3},2,{}},{{52,17,2,7,12,30},3,{}}};vector<Profile>bank;ll C0=2;for(auto&g:G){g.p=profile(compute(g.out),g.root);need(g.p==literal_profiles(g.out)[g.root],"finite envelope bank literal profile");bank.push_back(g.p);C0=max(C0,g.p[4]+g.p[5]);}
 int X=127;vector<pair<Four,vector<int>>>queue{{{1,1,1,1},{}}};unordered_set<Four,FourHash>seen{{1,1,1,1}};ll envelopechecks=0;array<set<ll>,3>values;for(size_t i=0;i<queue.size();++i){auto state=queue[i].first;for(int op=0;op<(int)bank.size()+2;++op){Four z=update(state,op,bank);need(slack(z)>=0,"grammar f nonnegative");if(op!=1)need(z[0]>=max(state[0],state[2])&&z[2]<=C0*z[0],"non-D bounded ratio");if(z[0]>X)continue;auto w=queue[i].second;w.push_back(op);for(int a=0;a<3;++a){int D=array<int,3>{1,3,5}[a];if(z[3]<=D&&slack(z)==0){envelope_trace(w,z,bank,C0,D);values[a].insert(z[0]);++envelopechecks;}}if(z[2]<=X&&seen.insert(z).second)queue.push_back({z,w});}}
 cout<<"finite_bank_size="<<bank.size()<<" C0="<<C0<<" envelope_X="<<X<<" retained_states="<<seen.size()<<" certified_zero_slack_traces="<<envelopechecks<<"\n";for(int a=0;a<3;++a){cout<<"bounded_d="<<array<int,3>{1,3,5}[a]<<" zero_slack_values=";for(ll h:values[a])cout<<h<<",";cout<<"\n";}
 // Literal graph checks of the multiplicative macro, including strongness.
 ll products=0;for(int len=0;len<=3;++len){int count=1;for(int i=0;i<len;++i)count*=3;for(int code=0;code<count;++code){int v=code,cost=0;vector<int>ks;for(int i=0;i<len;++i){int k=array<int,3>{1,2,4}[v%3];v/=3;ks.push_back(k);cost+=k+1;}if(cost>10)continue;vector<int>o={0};int r=0;Four expected={1,1,1,1};for(int k:ks){ll oldh=expected[0];for(int j=0;j<k;++j){int n=o.size();for(int&row:o)row|=1<<n;o.push_back(0);}int oldn=o.size();o=seam(o,r,{0},0);r=oldn;expected={(1+(1LL<<k))*oldh,oldh,(1+(1LL<<k))*oldh,oldh};auto d=compute(o);need(strong(o,d.full)&&four(profile(d,r))==expected,"physical multiplicative returning macro");}++products;}}
 set<ll>distinct;for(int a=0;a<=4;++a)for(int b=0;b<=4;++b)for(int c=0;c<=4;++c){ll h=1;for(int i=0;i<a;++i)h*=3;for(int i=0;i<b;++i)h*=5;for(int i=0;i<c;++i)h*=17;need(distinct.insert(h).second,"prime exponent uniqueness");}cout<<"physical_product_words="<<products<<" distinct_exponent_products="<<distinct.size()<<"\nPASS\n";
}catch(exception&e){cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
