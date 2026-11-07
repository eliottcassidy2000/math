#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
int main(){uint64_t want;int n,strong,tests=0,sc=0;while(cin>>want>>strong>>n){if(n>16||n<1)return 1;vector<unsigned>a(n),in(n);for(auto&x:a)cin>>x;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[j]>>i&1)in[i]|=1u<<j;vector<uint64_t>dp((1<<n)*n);for(unsigned S=1;S<(1u<<n);S++)for(int j=0;j<n;j++)if(S>>j&1){unsigned R=S^(1u<<j);if(!R){dp[S*n+j]=1;continue;}unsigned pred=R&in[j];while(pred){int i=__builtin_ctz(pred);pred&=pred-1;dp[S*n+j]+=dp[R*n+i];}}uint64_t got=0;for(int j=0;j<n;j++)got+=dp[((1<<n)-1)*n+j];if(got!=want){cerr<<tests<<" "<<want<<" != "<<got<<endl;return 2;}if(strong){vector<unsigned>reach=a;for(int i=0;i<n;i++)reach[i]|=1u<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];for(auto r:reach)if(r!=(1u<<n)-1)return 3;sc++;}tests++;}cout<<"{\"literal_checks\":"<<tests<<",\"strongness_checks\":"<<sc<<"}"<<endl;}
