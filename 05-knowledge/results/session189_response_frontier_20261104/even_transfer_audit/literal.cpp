#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using V=std::array<uint64_t,5>;
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;if(n<1||n>18)throw std::runtime_error("vertex bound");std::vector<unsigned> rows(n);for(auto&r:rows)std::cin>>r;unsigned end=1u<<n;std::vector<V> at(end*n),sum(end);for(unsigned S=1;S<end;++S){unsigned left=S;while(left){unsigned b=left&-left;left-=b;int v=__builtin_ctz(b);unsigned M=S^b;auto &z=at[S*n+v];if(!M){z[0]=1;}else{for(int h=1;h<5;++h)z[h]=sum[M][h-1];unsigned rem=M;while(rem){unsigned a=rem&-rem;rem-=a;int u=__builtin_ctz(a);if(rows[u]>>v&1)for(int h=0;h<5;++h)z[h]+=at[M*n+u][h];}}for(int h=0;h<5;++h)sum[S][h]+=z[h];}}for(auto x:sum[end-1])std::cout<<x<<' ';std::cout<<'\n';}}
