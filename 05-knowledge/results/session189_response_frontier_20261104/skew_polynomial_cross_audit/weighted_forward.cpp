#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using W=__int128_t;
std::string dec(W x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;std::string s;while(x){s+=char('0'+x%10);x/=10;}if(neg)s+='-';std::reverse(s.begin(),s.end());return s;}
std::vector<uint64_t> forward(const std::vector<unsigned>&R){int n=R.size();unsigned full=(1u<<n)-1;std::vector<uint64_t>d(size_t(full+1)*n*n);for(int v=0;v<n;++v)d[(size_t(1u<<v)*n+v)*n]=1;for(unsigned S=1;S<=full;++S){int degree=__builtin_popcount(S)-1;unsigned ends=S;while(ends){unsigned b=ends&-ends;ends-=b;int u=__builtin_ctz(b);const uint64_t*from=&d[(size_t(S)*n+u)*n];unsigned unused=full^S;while(unused){unsigned c=unused&-unused;unused-=c;int v=__builtin_ctz(c);uint64_t*to=&d[(size_t(S|c)*n+v)*n];bool edge=R[u]>>v&1;for(int k=0;k<=degree;++k){if(edge)to[k]+=from[k];to[k+1]+=from[k];}}}}return d;}
int main(){int count;std::cin>>count;while(count--){int n,z,ng;std::cin>>n>>z>>ng;if(n>17)throw std::runtime_error("order bound");std::vector<unsigned>R(n),T(n);std::vector<int>groups(n);for(auto&r:R)std::cin>>r;for(auto&g:groups)std::cin>>g;unsigned full=(1u<<n)-1;for(int i=0;i<n;++i)T[i]=full^(1u<<i)^R[i];auto E=forward(R);auto A=forward(T);std::vector<W>out(ng*(n-1));for(unsigned S=1;S<full;++S)if(S>>z&1){unsigned U=full^S;int a=__builtin_popcount(S)-1,b=__builtin_popcount(U)-1;const auto*ez=&E[(size_t(S)*n+z)*n];const auto*az=&A[(size_t(S)*n+z)*n];unsigned left=U;while(left){unsigned bit=left&-left;left-=bit;int i=__builtin_ctz(bit);const auto*ei=&E[(size_t(U)*n+i)*n];const auto*ai=&A[(size_t(U)*n+i)*n];for(int p=0;p<=a;++p)for(int q=0;q<=b;++q)out[groups[i]*(n-1)+p+q]+=W(ez[p])*ai[q]-W(az[p])*ei[q];}}for(auto x:out)std::cout<<dec(x)<<' ';std::cout<<'\n';}}
