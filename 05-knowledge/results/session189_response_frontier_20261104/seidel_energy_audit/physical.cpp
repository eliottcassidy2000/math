#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(){int n,H;while(cin>>n>>H){if(n<1||n>16||H<1||H>5)throw runtime_error("bound");vector<int>out(n);for(int&i:out)cin>>i;size_t N=1<<n;vector<unsigned long long>dp(N*n*H);auto ix=[&](int mask,int last,int h){return ((size_t)mask*n+last)*H+h;};for(int j=0;j<n;j++)dp[ix(1<<j,j,0)]=1;for(int mask=1;mask<(int)N;mask++)for(int j=0;j<n;j++)if(mask>>j&1)for(int h=0;h<H;h++){auto v=dp[ix(mask,j,h)];if(!v)continue;int rem=((int)N-1)^mask;while(rem){int k=__builtin_ctz(rem);rem&=rem-1;int m=mask|(1<<k);if(out[j]>>k&1)dp[ix(m,k,h)]+=v;if(h+1<H)dp[ix(m,k,h+1)]+=v;}}for(int h=0;h<H;h++){unsigned long long ans=0;for(int j=0;j<n;j++)ans+=dp[ix(N-1,j,h)];if(h)cout<<' ';cout<<ans;}cout<<'\n';}}
