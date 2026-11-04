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
int main(){int n;map<int,ll>counts,strongcounts;set<ll>F61,P61;vector<vector<int>>rows61;bool first=true;cout<<"{\"target_profiles\":[";ll total=0;while(cin>>n){need(n>=3&&n<=8);vector<int>r(n);for(auto&x:r)cin>>x;for(int i=0;i<n;i++){need(!(r[i]>>i&1));need((r[i]>>n)==0);for(int j=0;j<i;j++)need(((r[i]>>j&1)+(r[j]>>i&1))==1);}counts[n]++;total++;if(!strong(r))continue;strongcounts[n]++;auto p=profile(r);if(p.H==61){for(auto f:p.F)F61.insert(f);P61.insert(p.P2);rows61.push_back(r);}if(p.H!=53&&p.H!=55&&p.H!=59&&p.H!=61)continue;for(int root=0;root<n;root++)if(p.F[root]==63){auto J=joint(r,p,root);need(*min_element(J.begin(),J.end())>0);if(!first)cout<<",";first=false;cout<<"{\"H\":"<<p.H<<",\"root\":"<<root<<",\"rows\":";array_out(r);cout<<",\"J\":";array_out(J);cout<<",\"s\":";array_out(p.s);cout<<"}";}}
cout<<"],\"F61\":";array_out(F61);cout<<",\"P2_61\":";array_out(P61);cout<<",\"rows61\":[";for(size_t i=0;i<rows61.size();i++){if(i)cout<<",";array_out(rows61[i]);}cout<<"],\"total\":"<<total<<",\"counts\":{";first=true;for(auto [n,c]:counts){if(!first)cout<<",";first=false;cout<<"\""<<n<<"\":"<<c;}cout<<"},\"strong_counts\":{";first=true;for(auto [n,c]:strongcounts){if(!first)cout<<",";first=false;cout<<"\""<<n<<"\":"<<c;}cout<<"}}\n";
}
