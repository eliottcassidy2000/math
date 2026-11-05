#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std; using W=array<int,3>; using G=vector<unsigned>; using ll=long long;
void need(bool b,const char*s){if(!b)throw runtime_error(s);} bool arc(const G&a,int i,int j){return a[i]>>j&1;}
W key(W w){return min(w,min(W{w[1],w[2],w[0]},W{w[2],w[0],w[1]}));}
G graph(W w){vector<int>o;for(int i=0;i<3;i++)for(int j=0;j<w[i];j++)o.push_back(i);G a(o.size());for(int i=0;i<(int)a.size();i++)for(int j=i+1;j<(int)a.size();j++){bool f=o[i]==o[j]||(o[j]-o[i]+3)%3==1;a[f?i:j]|=1u<<(f?j:i);}return a;}
bool anti(const G&a){int n=a.size();vector<int>p(n,-1),used(n);function<bool(int)>go=[&](int d){if(d==n)return true;int u=-1;vector<int>best;for(int i=0;i<n;i++)if(p[i]<0){vector<int>opts;for(int j=0;j<n;j++)if(!used[j]&&__builtin_popcount(a[i])+__builtin_popcount(a[j])==n-1){bool ok=true;for(int k=0;k<n;k++)if(p[k]>=0&&arc(a,i,k)==arc(a,j,p[k])){ok=false;break;}if(ok)opts.push_back(j);}if(opts.empty())return false;if(u<0||opts.size()<best.size()){u=i;best=opts;}}for(int j:best){p[u]=j;used[j]=1;if(go(d+1))return true;used[j]=0;p[u]=-1;}return false;};return go(0);}
int main(){const int N=150;vector<ll>part(N+1),mom(N+1);part[0]=1;for(int s:{1,2,3})for(int n=s;n<=N;n++)part[n]+=part[n-s];mom[0]=1;for(int s:{1,2,3,3})for(int n=s;n<=N;n++)mom[n]+=mom[n-s];ll classes=0,physical=0,actualtests=0;for(int n=3;n<=N;n++){set<W>bank;for(int a=1;a<n;a++)for(int b=1;a+b<n;b++)bank.insert(key(W{a,b,n-a-b}));map<int,ll>inventory;ll sum=0;for(W w:bank){int d=min({abs(w[0]-w[1]),abs(w[1]-w[2]),abs(w[2]-w[0])});inventory[d]++;sum+=d;classes++;if(n<=12){int found=-1;for(int add=0;found<0;add++)for(int a=0;a<=add&&found<0;a++)for(int b=0;b<=add-a&&found<0;b++){W v{w[0]+a,w[1]+b,w[2]+add-a-b};actualtests++;if(anti(graph(v)))found=add;}need(found==d,"literal optimal completion");physical++;}}need((ll)bank.size()==((ll)(n-1)*(n-2)/2+2*(n%3==0))/3,"ordinary Burnside");need(inventory[0]==(n-1)/2,"exact zero cost");need(inventory.rbegin()->first==(n-3)/3,"maximum cost");for(int d=1;d<=n/3+1;d++){ll exact=n>=3+3*d?2*((n-1-3*d)/2):0;auto it=inventory.find(d);need((it==inventory.end()?0:it->second)==exact,"exact positive cost coefficient");ll tail=0;for(auto [j,c]:inventory)if(j>=d)tail+=c;need(tail==(n>=3+3*d?2*part[n-3-3*d]:0),"tail coefficient");}need(sum==(n>=6?2*mom[n-6]:0),"moment coefficient");if(n%25==0)cout<<"ORDER "<<n<<" classes "<<bank.size()<<" cost_sum "<<sum<<" zero_cost "<<inventory[0]<<'\n';}cout<<"PASS order_limit "<<N<<" ordinary_classes "<<classes<<" physical_optima "<<physical<<" literal_anti_tests "<<actualtests<<'\n';}
