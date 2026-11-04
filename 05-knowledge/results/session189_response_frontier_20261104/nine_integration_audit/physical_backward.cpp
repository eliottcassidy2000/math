#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;if(n>15)throw std::runtime_error("order bound");std::vector<unsigned>r(n),incoming(n);for(auto&x:r)std::cin>>x;for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(r[i]>>j&1)incoming[j]|=1u<<i;unsigned full=(1u<<n)-1;std::vector<uint64_t>dp((1u<<n)*n);for(unsigned mask=1;mask<=full;++mask){unsigned ends=mask;while(ends){int j=__builtin_ctz(ends);ends&=ends-1;unsigned prior=mask^(1u<<j);if(!prior){dp[mask*n+j]=1;continue;}unsigned prev=prior&incoming[j];uint64_t v=0;while(prev){int i=__builtin_ctz(prev);prev&=prev-1;v+=dp[prior*n+i];}dp[mask*n+j]=v;}}uint64_t answer=0;for(int j=0;j<n;++j)answer+=dp[full*n+j];std::cout<<answer<<'\n';}}
