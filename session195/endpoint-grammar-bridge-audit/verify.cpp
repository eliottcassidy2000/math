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
#include <unordered_set>
struct Rooted{vector<int>out;int root;};
Rooted step(Rooted a,char c){int n=a.out.size();a.out.push_back(0);for(int i=0;i<n;++i)if(c!='D'&&i==a.root)a.out[n]|=1<<i;else a.out[i]|=1<<n;if(c=='A')a.root=n;return a;}
Rooted run(string w,Rooted a={{0},0}){for(char c:w)a=step(a,c);return a;}
vector<int>opposite(const vector<int>&o){int full=(1<<o.size())-1;vector<int>r;for(int i=0;i<(int)o.size();++i)r.push_back(full^o[i]^(1<<i));return r;}
vector<int>induced(const vector<int>&o,const vector<int>&names){vector<int>q(names.size());for(int i=0;i<(int)names.size();++i)for(int j=0;j<(int)names.size();++j)if(o[names[i]]>>names[j]&1)q[i]|=1<<j;return q;}
vector<int>chain_order(const vector<int>&o,int m){vector<int>v;for(int i=0;i<(int)o.size();++i)if(m>>i&1)v.push_back(i);sort(v.begin(),v.end(),[&](int a,int b){return __builtin_popcount((unsigned)(o[a]&m))>__builtin_popcount((unsigned)(o[b]&m));});for(int i=0;i<(int)v.size();++i)for(int j=i+1;j<(int)v.size();++j)need(o[v[i]]>>v[j]&1,"neutral vertices form chain");return v;}
int sink_component(const vector<int>&o,int mask){int C=mask;for(int u=0;u<(int)o.size();++u)if(mask>>u&1)C&=reachable(o,mask,u);return C;}
struct Compiled{string word;vector<int>names;};
Compiled compile_endpoint(const Data&d,int mask,int root){need(d.e[mask][root]==1,"compiler input unit endpoint");if((mask&(mask-1))==0)return {"",{root}};int C=sink_component(d.out,mask);if(C!=mask){need(C>>root&1,"root lies in sink SCC");int P=mask^C;need(d.h[P]==1,"weak neutral prefix");auto out=compile_endpoint(d,C,root);auto chain=chain_order(d.out,P);reverse(chain.begin(),chain.end());for(int z:chain){out.word+='D';out.names.push_back(z);}return out;}
 int rem=mask^(1<<root);C=sink_component(d.out,rem);int Z=rem^C,u=-1;for(int x=0;x<d.n;++x)if((C>>x&1)&&(d.out[x]>>root&1)){need(u<0,"strong unique predecessor");u=x;}need(u>=0&&d.h[Z]==1,"strong inverse");auto out=compile_endpoint(d,C,u);out.word+='A';out.names.push_back(root);auto chain=chain_order(d.out,Z);reverse(chain.begin(),chain.end());for(int z:chain){out.word+=(d.out[root]>>z&1)?'K':'D';out.names.push_back(z);}return out;}
vector<int>endpoint_macro(const vector<int>&o,int u,const vector<int>&b){int n=o.size(),m=b.size(),v=n+m;vector<int>q=o;q.resize(v+1);for(int i=0;i<n;++i)if(i==u)q[i]|=1<<v;else q[v]|=1<<i;for(int i=0;i<m;++i){q[n+i]|=(1<<n)-1;for(int j=i+1;j<m;++j)q[n+i]|=1<<(n+j);if(b[i])q[v]|=1<<(n+i);else q[n+i]|=1<<v;}return q;}
using Four=array<ll,4>;
struct FourHash{size_t operator()(Four x)const{size_t z=0;for(ll v:x)z=z*1000033^v;return z;}};
int main(){try{
 ll wordchecks=0,prefixchecks=0,literals=0;for(int len=0;len<=8;++len){int count=1;for(int i=0;i<len;++i)count*=3;for(int code=0;code<count;++code){int x=code;string w;for(int i=0;i<len;++i){w+="AKD"[x%3];x/=3;}auto a=run(w);int t=0;while(t<len&&w[t]=='K')++t;auto b=run(w.substr(t));vector<int>names{0};for(int i=t+1;i<=len;++i)names.push_back(i);need(induced(a.out,names)==b.out&&names[b.root]==a.root,"source prefix graph erasure");for(int i=1;i<=t;++i){for(int v:names)need(a.out[i]>>v&1,"permanent source prefix");for(int j=i+1;j<=t;++j)need(a.out[i]>>j&1,"source prefix ordering");}auto da=compute(a.out),db=compute(b.out);need(da.h.back()==db.h.back()&&db.s.back()[b.root]==1,"source stripping count and unit start");need(da.s.back()[a.root]==(t?0:1),"initial K exception only");if(len<=5){literal_check(da);++literals;}++wordchecks;prefixchecks+=(t>0);}}
 cout<<"all_words_through_length=8 word_checks="<<wordchecks<<" initial_K_words="<<prefixchecks<<" literal_word_checks="<<literals<<"\n";
 ll compiled=0,strong_compiled=0,weak_compiled=0,macrochecks=0,invariance=0;for(int n=1;n<=6;++n){int E=n*(n-1)/2;for(unsigned code=0;code<(1u<<E);++code){vector<int>o(n);int t=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++t)if(code>>t&1)o[i]|=1<<j;else o[j]|=1<<i;auto d=compute(o);for(int r=0;r<n;++r){if(d.e.back()[r]==1){auto c=compile_endpoint(d,d.full,r);auto got=run(c.word);need(got.out==opposite(induced(o,c.names)),"recursive compiler exact graph");need(c.names[got.root]==r,"recursive compiler exact root");need(c.word.empty()||c.word[0]!='K',"compiler uses no initial K");++compiled;if(strong(o,d.full))++strong_compiled;else ++weak_compiled;}
 if(n<=4){for(char op:string("AKD")){auto q=step({o,r},op);auto dq=compute(q.out);ll want=(n==1&&op=='K')?0:d.s.back()[r];need(dq.s.back()[q.root]==want,"arbitrary rooted start invariance");++invariance;}
 for(int m=0;m<=5;++m)for(int a=0;a<(1<<m);++a){vector<int>b(m);for(int i=0;i<m;++i)b[i]=(a>>(m-1-i))&1;auto target=endpoint_macro(o,r,b);string w="A";for(int i=m-1;i>=0;--i)w+=b[i]?'K':'D';auto got=run(w,{opposite(o),r});vector<int>names;for(int i=0;i<n;++i)names.push_back(i);for(int i=0;i<m;++i)names.push_back(n+m-i);names.push_back(n);need(induced(got.out,names)==opposite(target)&&names[n+m]==got.root,"arbitrary macro exact labelled compiler");++macrochecks;}}
 }}}
 cout<<"all_labelled_inputs_through_order=6 recursive_endpoint_compilations="<<compiled<<" strong="<<strong_compiled<<" weak="<<weak_compiled<<" arbitrary_rooted_macro_checks="<<macrochecks<<" rooted_start_invariance_checks="<<invariance<<"\n";
 array<string,4>inputs={"","A","AA","AAD"},contexts={"","A","AA","KA"};array<array<ll,4>,4>M{};for(int i=0;i<4;++i)for(int j=0;j<4;++j)M[i][j]=compute(run(inputs[j]+contexts[i]).out).h.back();ll det=0;array<int,4>p={0,1,2,3};do{int inversions=0;ll term=1;for(int i=0;i<4;++i){term*=M[i][p[i]];for(int j=i+1;j<4;++j)inversions+=(p[i]>p[j]);}det+=(inversions%2?-term:term);}while(next_permutation(p.begin(),p.end()));need(det==-16,"primitive four-coordinate rank minor");cout<<"physical_primitive_continuation_determinant="<<det<<"\n";
 int X=1000;Four seed={1,1,1,1};vector<Four>queue{seed};unordered_set<Four,FourHash>seen{seed};set<ll>values{1};for(size_t i=0;i<queue.size();++i){auto[h,e,R,d]=queue[i];for(Four z: {Four{h+R-e,h-e,h+R,h},Four{h+R-e,d,2*R-d,d},Four{h,0,2*R,d}}){need(0<=z[1]&&z[1]<=z[3]&&z[3]<=z[0]&&z[0]<=z[2],"physical four-state inequalities");if(z[0]>X)continue;values.insert(z[0]);if(z[2]<=X&&seen.insert(z).second)queue.push_back(z);}}
 vector<ll>missing;for(int h=1;h<=X;h+=2)if(!values.count(h))missing.push_back(h);need(missing==vector<ll>{7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821},"independent scalar support matches endpoint census");cout<<"elementary_X=1000 retained_four_states="<<seen.size()<<" scalar_values="<<values.size()<<" missing=";for(ll h:missing)cout<<h<<",";cout<<"\nPASS\n";
}catch(exception&e){cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
