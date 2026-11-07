// Exact all-vertex count via a transitive-twin block generating function.
#include <cstdio>
#include <vector>
#include <cstdint>
int main(){int n=16,closing=15;std::vector<unsigned> row={2,0};int root=0;const char* w="OOZJOZJJOJJOZJ";
for(int k=0;w[k];k++){int v=row.size();row.push_back(0);for(int u=0;u<v;u++){if(w[k]!='O'&&u==root)row[v]|=1u<<u;else row[u]|=1u<<v;}if(w[k]=='J')root=v;}
if(row.size()!=16)return 1;unsigned N=1u<<n;
auto ix=[n](unsigned mask,int v,int k){return ((size_t)mask*n+v)*3+k;};std::vector<unsigned long long> dp((size_t)N*n*3);
for(int v=0;v<n;v++)dp[ix(1u<<v,v,0)]=1;dp[ix(1u<<closing,closing,1)]=1;
for(unsigned mask=1;mask<N;mask++)for(int v=0;v<n;v++)for(int k=0;k<=2;k++){auto h=dp[ix(mask,v,k)];if(!h)continue;
for(int u=0;u<n;u++)if(!(mask>>u&1)){unsigned next=mask|(1u<<u);if(row[v]>>u&1)dp[ix(next,u,k)]+=h;if(k<2&&v!=closing&&u==closing)dp[ix(next,u,k+1)]+=h;}}
unsigned long long B[3]={};for(int v=0;v<n;v++)for(int k=0;k<=2;k++){B[k]+=dp[ix(N-1,v,k)];if(k<2&&v!=closing)B[k+1]+=dp[ix(N-1,v,k)];}
auto full=B[1]+((1ull<<20)-2)*B[2];printf("B0=%llu B1=%llu B2=%llu full36=%llu\n",B[0],B[1],B[2],full);
if(B[0]!=10063||B[1]!=16935||B[2]!=6013||full!=6305092397ull)return 2;return 0;}
