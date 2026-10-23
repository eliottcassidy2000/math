#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
#include <string>
using namespace std;
vector<int>fib={0,1,1};
vector<int> z(int n){vector<int>x;for(int i=fib.size()-1;i>=2;i--)if(fib[i]<=n){n-=fib[i];x.push_back(i);}return x;}
vector<int>vf(int n){auto zz=z(n);int s=0;for(int x:zz)s+=x;return {n,__builtin_popcount((unsigned)n),(int)sqrt(n),(int)zz.size(),zz.back(),s};}
struct Carrier{int n;vector<vector<int>>p,e;};
int main(){while(fib.back()<10000)fib.push_back(fib.back()+fib[fib.size()-2]);vector<pair<int,int>>E;int ei[81][81]={};for(int b=1;b<=80;b++)for(int a=1;a<b;a++){int k=sqrt(a+b);if(k*k==a+b){ei[a][b]=ei[b][a]=E.size();E.push_back({a,b});}}
vector<Carrier>C;for(int n,N;cin>>n>>N;){Carrier c;c.n=n;while(N--){vector<int>p(n),e;for(int&i:p)cin>>i;if(p.front()>p.back())reverse(p.begin(),p.end());for(int i=1;i<n;i++)e.push_back(ei[p[i-1]][p[i]]);c.p.push_back(p);c.e.push_back(e);}C.push_back(c);}
bool first=true;cout<<"[\n";
for(int kind=0;kind<2;kind++){vector<vector<int>>fs;if(!kind)for(int n=1;n<=80;n++)fs.push_back(vf(n));else for(auto[a,b]:E){auto A=vf(a),B=vf(b);fs.push_back({a,b,(int)sqrt(a+b),A[1]+B[1],A[2]+B[2],A[3]+B[3],A[4]+B[4]});}int dim=fs[0].size();for(int q:{2,3}){set<string>seen;int lim=1;for(int i=0;i<dim;i++)lim*=q;for(int code=1;code<lim;code++){vector<int>co(dim);int x=code;for(int i=dim-1;i>=0;i--){co[i]=x%q;x/=q;}string colors;for(auto f:fs){int v=0;for(int i=0;i<dim;i++)v+=co[i]*f[i];colors+=char(v%q);}string canon;int ren[3]={-1,-1,-1},next=0;for(unsigned char v:colors){if(ren[v]<0)ren[v]=next++;canon+=char(ren[v]);}if(!seen.insert(canon).second)continue;
 vector<int>counts;int failureN=0,fa=-1,fb=-1;bool training=true;
 for(auto &c:C){unordered_map<string,int>labels;for(int j=0;j<(int)c.p.size();j++){string sig;sig+=char(c.p[j].front());sig+=char(c.p[j].back());if(!kind)for(int v:c.p[j])sig+=colors[v-1];else for(int v:c.e[j])sig+=colors[v];auto [it,inserted]=labels.emplace(sig,j);if(!inserted&&!failureN){failureN=c.n;fa=it->second;fb=j;}}
 counts.push_back(labels.size());if(c.n==32&&labels.size()!=c.p.size()){training=false;break;}}
 if(!training)continue;
 if(!first)cout<<",\n";first=false;cout<<"{\"kind\":\""<<(kind?"edge":"vertex")<<"\",\"q\":"<<q<<",\"coefficients\":[";for(int i=0;i<dim;i++){if(i)cout<<",";cout<<co[i];}cout<<"],\"counts32to36\":[";for(int i=0;i<(int)counts.size();i++){if(i)cout<<",";cout<<counts[i];}cout<<"],\"first_failure\":"<<failureN<<",\"witness_indices\":["<<fa<<","<<fb<<"]}";
 }}}
cout<<"\n]\n";
}
