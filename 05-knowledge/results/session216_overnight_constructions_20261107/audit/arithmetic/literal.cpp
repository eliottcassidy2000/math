#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
using namespace std;
int main(){uint64_t want;int n,tests=0,maxn=0;while(cin>>want>>n){if(n>13||n<1)return 1;vector<unsigned>a(n),in(n);for(auto&x:a)cin>>x;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[j]>>i&1)in[i]|=1u<<j;vector<uint64_t>dp((1<<n)*n);for(unsigned S=1;S<(1u<<n);S++)for(int j=0;j<n;j++)if(S>>j&1){unsigned R=S^(1u<<j);if(!R){dp[S*n+j]=1;continue;}unsigned pred=R&in[j];while(pred){int i=__builtin_ctz(pred);pred&=pred-1;dp[S*n+j]+=dp[R*n+i];}}uint64_t got=0;for(int j=0;j<n;j++)got+=dp[((1<<n)-1)*n+j];if(got!=want){cerr<<tests<<" "<<want<<" != "<<got<<endl;return 2;}
// Check actual strongness independently using bitset transitive closure.
vector<unsigned>reach=a;for(int i=0;i<n;i++)reach[i]|=1u<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];for(auto r:reach)if(r!=(1u<<n)-1)return 3;tests++;maxn=max(n,maxn);}
if(tests!=891)return 4;cout<<"{\"literal_checks\":"<<tests<<",\"strongness_checks\":"<<tests<<",\"max_order\":"<<maxn<<",\"integer_bound\":\"13! < 2^64\"}"<<endl;}
