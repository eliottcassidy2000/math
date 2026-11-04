#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
int main(){int Q;cin>>Q;for(int z=0;z<Q;z++){int n;cin>>n;vector<unsigned>r(n+2);for(int i=0;i<n;i++)cin>>r[i];unsigned a,c;cin>>a>>c;for(int i=0;i<n;i++){if(a>>i&1)r[i]|=1u<<n;else r[n]|=1u<<i;if(c>>i&1)r[i]|=1u<<(n+1);else r[n+1]|=1u<<i;}r[n]|=1u<<(n+1);int N=n+2,M=1<<N;vector<uint64_t>dp((size_t)M*N);for(int i=0;i<N;i++)dp[((size_t)1<<i)*N+i]=1;for(unsigned s=1;s<(unsigned)M;s++)for(int j=0;j<N;j++)if(s>>j&1){unsigned prev=s^(1u<<j);if(!prev)continue;uint64_t v=0;for(unsigned t=prev;t;t&=t-1){int i=__builtin_ctz(t);if(r[i]>>j&1)v+=dp[(size_t)prev*N+i];}dp[(size_t)s*N+j]=v;}uint64_t H=0;for(int i=0;i<N;i++)H+=dp[(size_t)(M-1)*N+i];cout<<H<<'\n';}}
