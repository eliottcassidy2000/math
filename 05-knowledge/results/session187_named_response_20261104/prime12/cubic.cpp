#include <vector>
#include <array>
#include <iostream>
#include <map>
#include <cstdint>
using namespace std;
using I=long long;
I hp(const vector<int>&a){int n=a.size(),f=(1<<n)-1;vector<I>d((1<<n)*n);vector<int>in(n);for(int v=0;v<n;v++)in[v]=f^a[v]^(1<<v);for(int m=1;m<=f;m++)for(int z=m;z;z&=z-1){int v=__builtin_ctz((unsigned)z),r=m^(1<<v);I t=0;if(!r)t=1;else for(int y=r&in[v];y;y&=y-1)t+=d[r*n+__builtin_ctz((unsigned)y)];d[m*n+v]=t;}I z=0;for(int v=0;v<n;v++)z+=d[f*n+v];return z;}
vector<int>op(vector<int>a){int f=(1<<a.size())-1;for(int i=0;i<(int)a.size();i++)a[i]^=f^(1<<i);return a;}
vector<int>ex(vector<int>a,int x,int c){int n=a.size();a.resize(n+2);a[n]=1<<(n+1);for(int v=0;v<n;v++){if(x>>v&1)a[v]|=1<<n;else a[n]|=1<<v;if(c>>v&1)a[v]|=1<<(n+1);else a[n+1]|=1<<v;}return a;}
int main(int argc,char**argv){vector<int>a={1284,273,1322,1411,1325,1355,1311,1399,1536,1279,2050,1023};if(argc>2){int order=stoi(argv[1]);unsigned long long code=stoull(argv[2]);if(order<2||order>6||order%2)return 2;a.assign(order,0);int k=0;for(int u=0;u<order;u++)for(int v=u+1;v<order;v++,k++)if(code>>k&1)a[u]|=1<<v;else a[v]|=1<<u;}int n=a.size(),f=(1<<n)-1;auto b=op(a);map<pair<int,int>,I>vals;auto val=[&](int x,int y){auto key=make_pair(x,y);if(!vals.count(key))vals[key]=hp(ex(a,x,y))-hp(ex(b,x,y));return vals[key];};vector<array<I,3>> coeff;vector<int>sm{0};for(int i=0;i<n;i++)sm.push_back(1<<i);for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)sm.push_back((1<<i)|(1<<j));
 for(int x:sm)for(int y:sm){if(__builtin_popcount((unsigned)x)+__builtin_popcount((unsigned)y)>3)continue;I c=0;for(int xx=x;;xx=(xx-1)&x){for(int yy=y;;yy=(yy-1)&y){int s=__builtin_popcount((unsigned)(x^xx))+__builtin_popcount((unsigned)(y^yy));c+=(s&1?-1:1)*val(xx,yy);if(!yy)break;}if(!xx)break;}coeff.push_back({x,y,c});}
 // Independently compute all AAC cubic parities by truncated disjoint path supports.
 vector<vector<int>> G(1<<n,vector<int>(n));vector<int>S(1<<n);for(int m=1;m<=f;m++){for(int v=0;v<n;v++)if(m>>v&1){int r=m^(1<<v);if(!r)G[m][v]=1<<v;else for(int z=r&a[v];z;z&=z-1)G[m][v]^=G[r][__builtin_ctz((unsigned)z)];if(__builtin_parity((unsigned)G[m][v]))S[m]|=1<<v;}}
 vector<vector<int>> K(1<<n,vector<int>(n));long long partitions=0;for(int A=1;A<=f;A++){int rem=f^A;for(int B=rem;B;B=(B-1)&rem){int sz=__builtin_popcount((unsigned)(A|B));if(sz%2==0||sz>n-3)continue;partitions++;for(int us=S[A];us;us&=us-1){int u=__builtin_ctz((unsigned)us);for(int vs=B;vs;vs&=vs-1){int v=__builtin_ctz((unsigned)vs);for(int ws=G[B][v];ws;ws&=ws-1)K[(1<<u)|(1<<v)][__builtin_ctz((unsigned)ws)]^=1;}}}}
 int nonzero=0,cn=0;for(auto[x,y,c]:coeff){int sx=__builtin_popcount((unsigned)x),sy=__builtin_popcount((unsigned)y);if(c%2){cerr<<"odd difference coefficient";return 1;}if(sx+sy==3){int k=sx==2?K[x][__builtin_ctz((unsigned)y)]:K[y][__builtin_ctz((unsigned)x)];if(((c/2)%2+2)%2!=k){cerr<<"cubic mismatch "<<x<<","<<y;return 1;}cn+=k;}nonzero+=c!=0;}
 uint32_t r=193;for(int k=0;k<24;k++){r=1664525*r+1013904223;int x=r&f;r=1664525*r+1013904223;int y=r&f;I z=0;for(auto[i,j,c]:coeff)if((i&x)==i&&(j&y)==j)z+=c;if(z!=val(x,y)){cerr<<"evaluation mismatch";return 1;}}
 cout<<"{\"rows\":[";for(int i=0;i<n;i++)cout<<(i?",":"")<<a[i];cout<<"],\"H\":"<<hp(a)<<",\"literal_query_assignments\":"<<vals.size()<<",\"independent_random_assignments\":24,\"nonzero_exact_coefficients\":"<<nonzero<<",\"nonzero_modular_cubic_coefficients\":"<<cn<<",\"ordered_disjoint_support_pairs\":"<<partitions<<",\"coefficients\":[";bool first=true;for(auto[x,y,c]:coeff){cout<<(first?"":",")<<"["<<x<<","<<y<<","<<c<<"]";first=false;}cout<<"]}\n";
}
