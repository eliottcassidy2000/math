#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using V=std::array<uint64_t,5>;
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;if(n>13)throw std::runtime_error("bound13");std::vector<unsigned>out(n);for(auto&r:out)std::cin>>r;unsigned full=(1u<<n)-1;std::vector<V>dp((1u<<n)*n);for(int i=0;i<n;++i)dp[(1u<<i)*n+i][0]=1;for(unsigned m=1;m<=full;++m)for(int i=0;i<n;++i)if(m>>i&1){const auto &v=dp[m*n+i];for(int j=0;j<n;++j)if(!(m>>j&1)){auto&w=dp[(m|(1u<<j))*n+j];bool e=out[i]>>j&1;for(int h=0;h<5;++h)w[h]+=(e?v[h]:0)+(h?v[h-1]:0);}}V ans{};for(int i=0;i<n;++i)for(int h=0;h<5;++h)ans[h]+=dp[full*n+i][h];for(auto v:ans)std::cout<<v<<' ';std::cout<<'\n';}}
