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
struct Marked{vector<int>out;int root;State st;};
int main(){try{
 vector<Marked> small;set<State> observed;ll totalroot=0;
 for(int n=1;n<=6;++n){ll rooted=0,strong_t=0;int E=n*(n-1)/2;for(unsigned code=0;code<(1u<<E);++code){vector<int>o(n);int t=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++t)if(code>>t&1)o[i]|=1<<j;else o[j]|=1<<i;auto d=compute(o);if(n<=5)literal_check(d);if(!strong(o,d.full))continue;++strong_t;
 for(int v=0;v<n;++v)if(d.e.back()[v]==1){++rooted;auto st=state(d,d.full,v);observed.insert(st);if(n<=4)small.push_back({o,v,st});if(n==1){need(st==State{1,1,1},"singleton");continue;}
 int rem=d.full^(1<<v),C=rem;for(int x=0;x<n;++x)if(rem>>x&1)C&=reachable(o,rem,x);need(strong(o,C),"independent sink SCC");int Z=rem^C;vector<int>chain,core;for(int x=0;x<n;++x)if(Z>>x&1)chain.push_back(x);else if(C>>x&1)core.push_back(x);
 sort(chain.begin(),chain.end(),[&](int x,int y){return __builtin_popcount((unsigned)(o[x]&Z))>__builtin_popcount((unsigned)(o[y]&Z));});
 for(int i=0;i<(int)chain.size();++i){for(int j=i+1;j<(int)chain.size();++j)need(o[chain[i]]>>chain[j]&1,"prefix transitive");need((o[chain[i]]&C)==C,"prefix dominates core");}
 int u=-1;for(int x:core)if(o[x]>>v&1){need(u<0,"unique predecessor");u=x;}need(u>=0&&d.e[C][u]==1&&d.h[Z]==1,"inverse marked factorization");
 int a=0;for(int x:chain)a=2*a+((o[v]>>x)&1);need(chain.empty()?core.size()>1:(o[v]>>chain[0]&1),"strongness boundary");need((int)bits(a).size()==(int)chain.size(),"leading bit and length");need(transfer(state(d,C,u),a)==st,"inverse state transfer");
 vector<int>compressed(core.size());int uid=-1;for(int i=0;i<(int)core.size();++i){if(core[i]==u)uid=i;for(int j=0;j<(int)core.size();++j)if(o[core[i]]>>core[j]&1)compressed[i]|=1<<j;}
 auto rebuilt=construct(compressed,uid,a);vector<int>names=core;names.insert(names.end(),chain.begin(),chain.end());names.push_back(v);for(int i=0;i<n;++i)for(int j=0;j<n;++j)need(bool(rebuilt[i]>>j&1)==bool(o[names[i]]>>names[j]&1),"literal graph reconstruction");
 }}totalroot+=rooted;cout<<"n="<<n<<" labelled="<<(1u<<E)<<" strong="<<strong_t<<" endpoint_one_roots="<<rooted<<"\n";}
 ll forward=0,literal_forward=0;for(auto&x:small)for(int a=0;a<32;++a){if(!a&&x.st[0]==1)continue;auto out=construct(x.out,x.root,a);auto d=compute(out);int v=out.size()-1;need(strong(out,d.full)&&d.e.back()[v]==1,"forward strong endpoint-one");need(state(d,d.full,v)==transfer(x.st,a),"physical forward transfer");need(state(d,d.full,v)==alphabet(x.st,a),"physical alphabet transfer");if(d.n<=6){literal_check(d);++literal_forward;}++forward;}
 ll words=0;for(int a=1;a<4096;++a)for(State st:{State{1,0,0},State{0,1,0},State{0,0,1}}){need(transfer(st,a)==alphabet(st,a),"all linear basis word transfers");++words;}
 cout<<"total_endpoint_one_roots="<<totalroot<<" forward_graphs="<<forward<<" literal_forward_graphs="<<literal_forward<<" basis_word_checks="<<words<<"\n";
 for(int X:{13,1000,10000}){set<State>found{{1,1,1}};deque<State>queue{{1,1,1}};ll edges=0;while(!queue.empty()){auto st=queue.front();queue.pop_front();vector<int>params;if(st[0]>1)params.push_back(0);for(int a=1;a<=X/st[2]-1;++a)params.push_back(a);for(int a:params){auto ch=transfer(st,a);if(ch[0]>X)continue;need(ch[0]>st[0],"census strict growth");need(ch[1]>=1&&ch[1]<=ch[0]&&ch[2]>=ch[0],"census state invariants");++edges;if(found.insert(ch).second)queue.push_back(ch);}}
 set<ll>values;for(auto st:found)values.insert(st[0]);cout<<"X="<<X<<" states="<<found.size()<<" values="<<values.size()<<" retained_edges="<<edges<<" missing=";for(int h=1;h<=X;h+=2)if(!values.count(h))cout<<h<<",";cout<<"\n";
 if(X==13){for(auto st:found)cout<<"state13="<<st[0]<<","<<st[1]<<","<<st[2]<<"\n";need(found.size()==13&&!values.count(13),"H13 obstruction");}
 if(X==1000)need(found.size()==9293&&values.size()==475,"1000 census");
 if(X==10000)need(found.size()==284844&&values.size()==4971,"10000 census");
 }
 auto h13=compute({2,12,9,17,7});need(h13.h.back()==13,"H13 witness");for(int i=0;i<5;++i)need(h13.e.back()[i]==array<ll,5>{4,3,2,2,2}[i],"H13 endpoints");auto witness=construct({0},0,0b111010111);auto d511=compute(witness);need(strong(witness,d511.full)&&d511.e.back()[10]==1&&state(d511,d511.full,10)==State{511,1,553},"H511 witness");cout<<"H511_state=511,1,553; rows=";for(int row:witness)cout<<row<<",";cout<<"\n";
 State seq={3,1,3};vector<ll>vals;for(int j=0;j<8;++j){vals.push_back(seq[0]);seq=transfer(seq,0);}need(vals==vector<ll>{3,5,9,17,31,57,105,193},"empty-word sequence");
 array<array<ll,3>,3>M={array<ll,3>{1,-1,1},array<ll,3>{1,-1,0},array<ll,3>{1,0,1}},M2{},M3{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)M2[i][j]+=M[i][k]*M[k][j];for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)M3[i][j]+=M2[i][k]*M[k][j];for(int i=0;i<3;++i)for(int j=0;j<3;++j)need(M3[i][j]==M2[i][j]+M[i][j]+(i==j),"all-order matrix recurrence");cout<<"tribonacci_matrix_identity=PASS\n";
 cout<<"PASS\n";
 }catch(exception&e){cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
