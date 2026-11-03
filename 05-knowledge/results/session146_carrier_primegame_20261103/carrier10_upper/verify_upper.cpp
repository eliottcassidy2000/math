#include <bits/stdc++.h>
using namespace std;
constexpr int n=10,P=101; int D, rr;
struct Forest {int sign; uint64_t mask; vector<pair<int,int>> e;};
using Key=pair<vector<int>,vector<array<int,3>>>;
vector<pair<int,int>> pairs;
map<Key,Forest> generate(){
 map<Key,Forest> reps; int ix[n][n]={};for(int i=0;i<(int)pairs.size();i++)ix[pairs[i].first][pairs[i].second]=i;
 vector<pair<int,int>> matching; vector<int> U;for(int i=2*rr;i<n;i++)U.push_back(i); vector<int> lens(rr);
 function<void(unsigned)> match=[&](unsigned left){if(left){int x=__builtin_ctz(left);left^=1u<<x;for(int y=x+1;y<2*rr;y++)if(left>>y&1){matching.push_back({x,y});match(left^(1u<<y));matching.pop_back();}return;}
  auto perm=U;do {function<void(int,int)> alloc=[&](int i,int rem){if(i==rr-1){lens[i]=rem;vector<int>S;vector<array<int,3>>T;Forest f{1,0,{}};int pos=0;for(int k=0;k<rr;k++){auto [a,b]=matching[k];int last=a;for(int z=0;z<lens[k];z++){int v=perm[pos++],lo=min(last,v),hi=max(last,v);f.mask|=1ull<<ix[lo][hi];f.sign*=last<v?1:-1;last=v;}int lo=min(last,b),hi=max(last,b);f.mask|=1ull<<ix[lo][hi];f.sign*=last<b?1:-1;T.push_back({a,b,lens[k]%2});}sort(T.begin(),T.end());auto key=make_pair(S,T);auto it=reps.find(key);if(it==reps.end()||f.mask<it->second.mask){for(int j=0;j<(int)pairs.size();j++)if(f.mask>>j&1)f.e.push_back(pairs[j]);reps[key]=f;}return;}for(int x=0;x<=rem;x++){lens[i]=x;alloc(i+1,rem-x);}};alloc(0,U.size());}while(next_permutation(perm.begin(),perm.end()));
 };match((1u<<(2*rr))-1);return reps;
}

int pw(int a,int b){int r=1;for(;b;b>>=1,a=a*a%P)if(b&1)r=r*a%P;return r;}
int coeff(int c,int r){int f=1;for(int i=1;i<=c;i++)f=f*i%P;return f*pw(2,r)%P*pw(pw(2,c-1),P-2)%P;}
int vals[100][100];
int coef(const Forest &f,const vector<pair<int,int>>&g,int N){
 int adj[100][2],deg[100]={},seen[100]={};
 auto add=[&](int a,int b){if(deg[a]>=2||deg[b]>=2)return false;adj[a][deg[a]++]=b;adj[b][deg[b]++]=a;return true;};
 for(auto[a,b]:f.e)if(!add(a,b))return 0;for(auto[a,b]:g)if(!add(a,b))return 0;
 int c=0,r=0,s=1;
 for(int i=0;i<N;i++)if(!seen[i]&&deg[i]<=1){c++;seen[i]=1;if(!deg[i])continue;r++;int at=i,prev=-1,len=0;while(true){int nx=-1;for(int k=0;k<deg[at];k++)if(adj[at][k]!=prev)nx=adj[at][k];if(nx<0)break;if(seen[nx])return 0;seen[nx]=1;if(at>nx)s=-s;len++;prev=at;at=nx;}if(len%2)return 0;}
 for(int i=0;i<N;i++)if(!seen[i])return 0;
 int v=vals[c][r];return s*f.sign<0?(P-v)%P:v;
}
vector<vector<int>> components(const Forest&f){vector<vector<int>> adj(n);for(auto[a,b]:f.e){adj[a].push_back(b);adj[b].push_back(a);}vector<vector<int>>cs;vector<int>seen(n);for(int i=0;i<n;i++)if(adj[i].size()<=1&&!seen[i]){vector<int>v;int at=i,prev=-1;while(at>=0){v.push_back(at);seen[at]=1;int nx=-1;for(int b:adj[at])if(b!=prev)nx=b;prev=at;at=nx;}cs.push_back(v);}return cs;}
int main(){
 rr=4;D=735;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)pairs.emplace_back(i,j);auto reps=generate();vector<Forest>fs;for(auto &z:reps)fs.push_back(z.second);for(int c=1;c<100;c++)for(int r=0;r<=c;r++)vals[c][r]=coeff(c,r);
 int U[735][56];ifstream un("kernel.txt");for(auto &row:U)for(int &x:row)un>>x;if(!un)return 4;
 int V[735][56];for(int i=0;i<735;i++)for(int j=0;j<56;j++)V[i][j]=(U[i][j]%P+P)%P;int rk=0;for(int j=0;j<56;j++){int pivot=rk;while(pivot<735&&!V[pivot][j])pivot++;if(pivot==735)continue;for(int k=0;k<56;k++)swap(V[rk][k],V[pivot][k]);int inv=pw(V[rk][j],P-2);for(int k=j;k<56;k++)V[rk][k]=V[rk][k]*inv%P;for(int i=rk+1;i<735;i++)if(V[i][j]){int t=V[i][j];for(int k=j;k<56;k++)V[i][k]=(V[i][k]-t*V[rk][k]%P+P)%P;}rk++;}cout<<"kernel candidate rank "<<rk<<endl;if(rk!=56)return 5;
 long count=0, zero=0;
 function<bool(unsigned,int,vector<pair<int,int>>&)> rec=[&](unsigned left,int N,vector<pair<int,int>>&g){if(!left){count++;int sums[56]={},first=0;for(int j=0;j<D;j++){int z=coef(fs[j],g,N);if(!z)continue;if(!first)first=z;int t=z*pw(first,P-2)%P;if(t!=1&&t!=100){cerr<<"scale failure "<<t<<endl;return false;}int sg=t==1?1:-1;for(int q=0;q<56;q++)sums[q]+=sg*U[j][q];}if(!first)zero++;for(int q=0;q<56;q++)if(sums[q]){cerr<<"violation diagram="<<count<<" col="<<q<<" sum="<<sums[q]<<" N="<<N<<" k="<<g.size();for(auto [a,b]:g)cerr<<' '<<a<<' '<<b;cerr<<endl;return false;}if(count%5000==0)cerr<<"checked "<<count<<endl;return true;}
 int x=__builtin_ctz(left);left^=1u<<x;
 for(int len:{1,2}){int at=x;int start=N;auto old=g.size();for(int z=0;z<len;z++){g.push_back({at,start});at=start++;}if(!rec(left,start,g))return false;g.resize(old);}
 for(int y=x+1;y<8;y++)if(left>>y&1)for(int len:{2,3}){int at=x;int start=N;auto old=g.size();for(int z=1;z<len;z++){g.push_back({at,start});at=start++;}g.push_back({at,y});if(!rec(left^(1u<<y),start,g))return false;g.resize(old);}
 return true;}; vector<pair<int,int>>g;bool ok=rec(255,n,g);ok=ok&&count==32400&&zero==17042;cout<<"checked "<<count<<" zero "<<zero<<" result "<<(ok?"PASS":"FAIL")<<endl;return !ok;
}
