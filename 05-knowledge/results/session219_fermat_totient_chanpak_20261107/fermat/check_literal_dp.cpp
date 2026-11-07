#include <cstdio>
#include <vector>
#include <cstdint>
int main(){int cases; if(scanf("%d",&cases)!=1)return 1;int maximum=0;
 for(int c=0;c<cases;c++){int n;unsigned long long expected;scanf("%d%llu",&n,&expected);maximum=maximum<n?n:maximum;
 std::vector<unsigned int> row(n);for(int i=0;i<n;i++)scanf("%u",&row[i]);
 unsigned N=1u<<n;std::vector<unsigned long long> dp((size_t)N*n);for(int v=0;v<n;v++)dp[((size_t)1<<v)*n+v]=1;
 for(unsigned mask=1;mask<N;mask++)for(int v=0;v<n;v++){auto h=dp[(size_t)mask*n+v];if(!h)continue;
 unsigned next=row[v]&~mask;while(next){int w=__builtin_ctz(next);next&=next-1;dp[(size_t)(mask|(1u<<w))*n+w]+=h;}}
 unsigned long long H=0;for(int v=0;v<n;v++)H+=dp[(size_t)(N-1)*n+v];if(H!=expected){printf("FAIL n%d expected%llu got%llu\n",n,expected,H);return 2;}}
 printf("PASS %d literal Hamilton-path subset-DP checks, maximum order%d; no source-profile recurrence used by DP.\n",cases,maximum);}
