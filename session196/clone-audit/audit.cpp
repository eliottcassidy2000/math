#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>
using namespace std;
using ll=long long; using G=vector<unsigned>;
void ck(bool b,const char*s){if(!b)throw runtime_error(s);}
bool arc(const G&g,int u,int v){return (g[u]>>v)&1;}
G tour(int n,unsigned bits){G g(n);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(bits>>k&1)g[i]|=1u<<j;else g[j]|=1u<<i;return g;}
void valid(const G&g){for(int i=0;i<(int)g.size();i++){ck(!(g[i]>>i&1),"loop");for(int j=0;j<i;j++)ck(arc(g,i,j)!=arc(g,j,i),"not tournament");}}
// Delete the first vertex: F[S,v] counts words on S that begin at v.
ll hp(const G&g){int n=g.size();if(!n)return 1;vector<ll>d((1<<n)*n);for(int S=1;S<(1<<n);S++)for(int v=0;v<n;v++)if(S>>v&1){int R=S^(1<<v);ll z=0;if(!R)z=1;else for(int w=0;w<n;w++)if((R>>w&1)&&arc(g,v,w))z+=d[R*n+w];d[S*n+v]=z;}ll z=0;for(int v=0;v<n;v++)z+=d[((1<<n)-1)*n+v];return z;}
G induced(const G&g,const vector<int>&vs){G h(vs.size());for(int i=0;i<(int)vs.size();i++)for(int j=0;j<(int)vs.size();j++)if(arc(g,vs[i],vs[j]))h[i]|=1u<<j;return h;}
// Construct in a fresh ordered vertex list, not by iterative cloning.
G blow(const G&g,int r,int m){vector<int>owner;for(int i=0;i<(int)g.size();i++)for(int k=0;k<(i==r?m:1);k++)owner.push_back(i);G h(owner.size());for(int i=0;i<(int)owner.size();i++)for(int j=0;j<(int)owner.size();j++)if(owner[i]==owner[j]?i<j:arc(g,owner[i],owner[j]))h[i]|=1u<<j;valid(h);return h;}
int nu_rec(const G&g,int r,const vector<int>&O,int pos,unsigned used){if(pos==(int)O.size())return 0;int best=nu_rec(g,r,O,pos+1,used);for(int y=0;y<(int)g.size();y++)if(y!=r&&arc(g,y,r)&&!(used>>y&1)&&arc(g,O[pos],y))best=max(best,1+nu_rec(g,r,O,pos+1,used|(1<<y)));return best;}
int nu(const G&g,int r){vector<int>O;for(int x=0;x<(int)g.size();x++)if(arc(g,r,x))O.push_back(x);return nu_rec(g,r,O,0,0);}
// Marked-word state DP: root is reusable but no consecutive root symbols.
vector<ll> coeff(const G&g,int r){int n=g.size();vector<int>V;for(int i=0;i<n;i++)if(i!=r)V.push_back(i);int full=(1<<(n-1))-1;vector<vector<vector<ll>>>d(1<<(n-1),vector<vector<ll>>(n,vector<ll>(n+1)));
 for(int j=0;j<n-1;j++){d[1<<j][j][0]=1;}
 d[0][n-1][1]=1;
 for(int S=0;S<=full;S++)for(int k=0;k<=n;k++)for(int u=0;u<n;u++){ll z=d[S][u][k];if(!z)continue;int old=u==n-1?r:V[u];for(int j=0;j<n-1;j++)if(!(S>>j&1)&&arc(g,old,V[j]))d[S|(1<<j)][j][k]+=z;if(u!=n-1&&k<n&&arc(g,old,r))d[S][n-1][k+1]+=z;}
 vector<ll>a(n+1);for(int u=0;u<n;u++)for(int k=0;k<=n;k++)a[k]+=d[full][u][k];if(n==1)a[0]=1;return a;}
ll evaluate(const vector<ll>&a,int m){ll S[16][16]={};S[0][0]=1;for(int i=1;i<=m;i++)for(int j=1;j<=i;j++)S[i][j]=j*S[i-1][j]+S[i-1][j-1];ll f=1,z=0;for(int k=0;k<(int)a.size();k++){if(k)f*=k;z+=a[k]*f*S[m][k];}return z;}
struct P{G g;int r;ll h,e,R,L,d;};
P prof(const G&g,int r){int n=g.size();vector<int>order(n);iota(order.begin(),order.end(),0);ll h=0,e=0,R=0,L=0;auto legal=[&](int lo,int hi){for(int j=lo;j+1<hi;j++)if(!arc(g,order[j],order[j+1]))return false;return true;};do{int p=find(order.begin(),order.end(),r)-order.begin();if(legal(0,n)){h++;e+=p==n-1;}R+=legal(0,p)&&legal(p,n);L+=legal(0,p+1)&&legal(p+1,n);}while(next_permutation(order.begin(),order.end()));vector<int>V;for(int i=0;i<n;i++)if(i!=r)V.push_back(i);return {g,r,h,e,R,L,hp(induced(g,V))};}
G synth(const P&x,const P&y,int m){int a=x.g.size()-1,b=y.g.size()-1,n=a+b+1+m;vector<int>A,B;for(int i=0;i<(int)x.g.size();i++)if(i!=x.r)A.push_back(i);for(int j=0;j<(int)y.g.size();j++)if(j!=y.r)B.push_back(j);G g(n);auto add=[&](int i,int j,bool dir){if(dir)g[i]|=1u<<j;else g[j]|=1u<<i;};for(int i=0;i<a;i++){add(0,1+i,arc(x.g,x.r,A[i]));for(int j=0;j<i;j++)add(1+j,1+i,arc(x.g,A[j],A[i]));}for(int i=0;i<b;i++){add(0,1+a+i,arc(y.g,y.r,B[i]));for(int j=0;j<i;j++)add(1+a+j,1+a+i,arc(y.g,B[j],B[i]));}for(int i=0;i<a;i++)for(int j=0;j<b;j++)add(1+i,1+a+j,true);for(int z=1+a+b;z<n;z++){add(z,0,true);for(int i=0;i<a;i++)add(1+i,z,true);for(int j=0;j<b;j++)add(z,1+a+j,true);for(int w=1+a+b;w<z;w++)add(w,z,true);}valid(g);return g;}
int main(){long long gs=0,rs=0,cs=0,ps=0,zeros=0;vector<P>bank;for(int n=1;n<=6;n++){for(unsigned bits=0;bits<(1u<<(n*(n-1)/2));bits++){G g=tour(n,bits);gs++;for(int r=0;r<n;r++){rs++;auto a=coeff(g,r);int top=0;for(int k=0;k<(int)a.size();k++)if(a[k])top=k;ck(top==1+nu(g,r),"matching degree");ck(a[1]==hp(g),"base m1");ck(a[0]==hp(blow(g,r,0)),"base m0");for(int m=1;m<=4;m++){ck(hp(blow(g,r,m))==evaluate(a,m),"physical clone law");cs++;}if(n<=4)bank.push_back(prof(g,r));}}cerr<<"cloning n="<<n<<" done\n";}
for(const P&x:bank){for(const P&y:bank){ll U=x.R-x.d,sig=x.R-x.h-x.d+x.e,alpha=U*y.L,beta=x.d*y.h-sig*y.d;for(int m=1;m<=3;m++){ck(hp(synth(x,y,m))==alpha*(1ll<<m)+beta,"physical synthesis law");ps++;}if(hp(synth(x,y,0))!=alpha+beta)zeros++;}}
int higher=0;for(int rank: {3,4}){int n=2*rank+1;G g(n);auto add=[&](int u,int v,bool dir){if(dir)g[u]|=1u<<v;else g[v]|=1u<<u;};for(int i=0;i<rank;i++){add(0,1+i,true);add(0,1+rank+i,false);for(int j=0;j<i;j++){add(1+j,1+i,true);add(1+rank+j,1+rank+i,true);}for(int j=0;j<rank;j++)add(1+i,1+rank+j,i==j);}valid(g);auto a=coeff(g,0);ck(nu(g,0)==rank,"higher rank matching");ck(a[rank+1]>0,"higher rank coefficient");for(int m=1;m<=5;m++){ck(hp(blow(g,0,m))==evaluate(a,m),"higher rank physical cloning");higher++;}}
cout<<"PASS graphs="<<gs<<" rooted="<<rs<<" actual_clone_checks="<<cs<<" synthesis_pairs="<<bank.size()*bank.size()<<" synthesis_checks="<<ps<<" higher_rank_physical_checks="<<higher<<" synthesis_m0_exponential_mismatches="<<zeros<<"\n";
}
