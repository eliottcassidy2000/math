#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;if(n>20)throw std::runtime_error("order bound20");std::vector<unsigned>r(n);for(auto&x:r)std::cin>>x;unsigned full=(1u<<n)-1;std::vector<uint64_t>d((1u<<n)*n);for(int i=0;i<n;++i)d[(1u<<i)*n+i]=1;for(unsigned m=1;m<=full;++m)for(int i=0;i<n;++i)if(m>>i&1){uint64_t v=d[m*n+i];unsigned todo=r[i]&~m;while(todo){int j=__builtin_ctz(todo);todo&=todo-1;d[(m|(1u<<j))*n+j]+=v;}}uint64_t h=0;for(int i=0;i<n;++i)h+=d[full*n+i];std::cout<<h<<'\n';}}
