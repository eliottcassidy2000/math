#include <array>
#include <vector>
#include <iostream>
#include <cstdint>
#include <chrono>
using namespace std;
constexpr int N=12,F=(1<<N)-1;
using A=array<int,N>;
uint16_t dp[1<<N], sp[1<<N][N];
int ends(const A&a){
 A in{};for(int v=0;v<N;v++)in[v]=F^a[v]^(1<<v);
 dp[0]=0;for(int m=1;m<=F;m++){int z=0;for(int x=m;x;x&=x-1){int v=__builtin_ctz((unsigned)x),r=m^(1<<v);if(!r||__builtin_parity((unsigned)(dp[r]&in[v])))z|=1<<v;}dp[m]=z;}return dp[F];
}
A rev(const A&a){A b;for(int i=0;i<N;i++)b[i]=F^a[i]^(1<<i);return b;}
bool prime(const A&a){for(int s=1;s<F;s++){if(!(s&(s-1)))continue;bool mod=true;for(int x=F^s;x;x&=x-1){int u=__builtin_ctz((unsigned)x),q=a[u]&s;if(q&&q!=s){mod=false;break;}}if(mod)return false;}return true;}
A oddskew(const A&a){A in{},g{};for(int v=0;v<N;v++)in[v]=F^a[v]^(1<<v);for(int m=1;m<=F;m++){for(int x=m;x;x&=x-1){int v=__builtin_ctz((unsigned)x),r=m^(1<<v),z=0;if(!r)z=1<<v;else for(int y=r&in[v];y;y&=y-1)z^=sp[r][__builtin_ctz((unsigned)y)];sp[m][v]=z;if(__builtin_parity((unsigned)m)&&m!=(1<<v))g[v]^=z;}}for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){int b=((g[u]>>v)^(g[v]>>u))&1;g[u]=(g[u]&~(1<<v))|(b<<v);g[v]=(g[v]&~(1<<u))|(b<<u);}return g;}
int main(int argc,char**argv){int depth=argc>1?atoi(argv[1]):3;vector<pair<int,int>> ed;for(int u=0;u<N;u++)for(int v=u+1;v<N;v++)ed.emplace_back(u,v);
 vector<A>bases={A{1284,1297,1322,1411,1325,1355,1311,1399,1536,1279,2048,1023},A{3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002}};
 for(int bi=0;bi<2;bi++){A a=bases[bi];uint64_t tested=0,balanced=0,p=0,skewzero=0;auto visit=[&](){tested++;if(ends(a)!=ends(rev(a)))return;balanced++;bool pr=prime(a);p+=pr;A g=oddskew(a);bool zero=true;for(int x:g)if(x)zero=false;if(zero){skewzero++;cout<<"SURVIVOR ";for(int x:a)cout<<x<<",";cout<<" prime "<<pr<<endl;}};
 auto dfs=[&](auto&&self,int start,int left)->void{if(!left){visit();return;}for(int k=start;k<=(int)ed.size()-left;k++){auto[u,v]=ed[k];a[u]^=1<<v;a[v]^=1<<u;self(self,k+1,left-1);a[u]^=1<<v;a[v]^=1<<u;}};
 for(int d=0;d<=depth;d++){auto t=chrono::steady_clock::now();dfs(dfs,0,d);cout<<"base "<<bi<<" through_depth "<<d<<" tested "<<tested<<" endpoint_balanced "<<balanced<<" prime_balanced "<<p<<" joint_survivors "<<skewzero<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;}}
}
