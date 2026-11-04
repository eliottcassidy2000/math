#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;using ll=long long;
void need(bool b){if(!b)throw runtime_error("check failed");}
struct Profile{ll H,P2;vector<ll>s,e,F,tot,end;};
Profile profile(const vector<int>&r){int n=r.size(),N=1<<n;Profile p;p.end.assign(N*n,0);vector<ll>start(N*n);p.tot.assign(N,0);p.tot[0]=1;for(int v=0;v<n;v++)p.end[(1<<v)*n+v]=start[(1<<v)*n+v]=1;for(int M=1;M<N;M++){for(int v=0;v<n;v++)if(M>>v&1){int old=M^(1<<v);if(old){for(int u=0;u<n;u++)if(old>>u&1){if(r[u]>>v&1)p.end[M*n+v]+=p.end[old*n+u];if(r[v]>>u&1)start[M*n+v]+=start[old*n+u];}}p.tot[M]+=p.end[M*n+v];}}p.H=p.tot.back();p.s.resize(n);p.e.resize(n);p.F.assign(n,0);for(int v=0;v<n;v++){p.s[v]=start[(N-1)*n+v];p.e[v]=p.end[(N-1)*n+v];for(int M=1;M<N;M++)p.F[v]+=p.end[M*n+v]*p.tot[(N-1)^M];}p.P2=0;for(int M=0;M<N;M++)p.P2+=p.tot[M]*p.tot[(N-1)^M];return p;}
bool strong(const vector<int>&r){int n=r.size();for(int v=0;v<n;v++){int seen=1<<v;while(true){int next=seen;for(int u=0;u<n;u++)if(seen>>u&1)next|=r[u];if(next==seen)break;seen=next;}if(seen!=(1<<n)-1)return false;}return true;}
vector<ll> joint(const vector<int>&r,const Profile&p,int root){int n=r.size(),N=1<<n;vector<ll>J(n);for(int v=0;v<n;v++){vector<ll>d(N*n);d[(1<<v)*n+v]=1;for(int M=1;M<N;M++){for(int u=0;u<n;u++)if(M>>u&1){int old=M^(1<<u);if(old)for(int z=0;z<n;z++)if((old>>z&1)&&(r[z]>>u&1))d[M*n+u]+=d[old*n+z];}J[v]+=d[M*n+root]*p.tot[(N-1)^M];}}return J;}
template<class C>void array_out(const C&c){cout<<"[";bool first=true;for(auto x:c){if(!first)cout<<",";first=false;cout<<x;}cout<<"]";}

int main(){int n;map<int,set<ll>>E;map<int,ll>hits9;long long parents=0,extensions=0,inputs=0;auto check=[&](const vector<int>&r,const Profile&p){if(p.H!=89&&p.H!=91&&p.H!=95&&p.H!=97)return;if(!strong(r))return;if(r.size()==9)hits9[p.H]++;for(int root=0;root<(int)r.size();root++)if(p.F[root]==99)E[p.H].insert(p.e[root]);};while(cin>>n){vector<int>r(n);for(auto&v:r)cin>>v;inputs++;auto p=profile(r);check(r,p);if(n!=8||p.H>97)continue;parents++;for(int mask=0;mask<256;mask++){auto child=r;for(int v=0;v<8;v++)if(mask>>v&1)child[v]|=256;child.push_back(255^mask);extensions++;auto q=profile(child);check(child,q);}}cout<<"{\"input_rows\":"<<inputs<<",\"order8_parents\":"<<parents<<",\"order9_extensions\":"<<extensions<<",\"cycle_sets\":{";bool first=true;for(auto[h,e]:E){if(!first)cout<<",";first=false;cout<<"\""<<h<<"\":";array_out(e);}cout<<"},\"target_H_order9_records\":{";first=true;for(auto[h,c]:hits9){if(!first)cout<<",";first=false;cout<<"\""<<h<<"\":"<<c;}cout<<"}}\n";}
