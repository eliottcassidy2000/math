#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
U choose(int n,int r){if(r<0||r>n)return 0;U z=1;for(int i=1;i<=r;++i)z=z*U(n-i+1)/U(i);return z;}
int main(){int n,cases=0;U total=0;while(std::cin>>n){if(n<1||n>9)throw std::runtime_error("documented range");std::vector<unsigned>a(n);std::vector<U>want(n),hist(n),f(n);for(auto&v:a)std::cin>>v;for(auto&v:want)std::cin>>v;std::vector<int>p(n);std::iota(p.begin(),p.end(),0);U words=0;do{int d=0;for(int i=0;i+1<n;++i)d+=!((a[p[i]]>>p[i+1])&1u);++hist[d];++words;}while(std::next_permutation(p.begin(),p.end()));for(int r=1;r<=n;++r)for(int d=0;d<r;++d)f[r-1]+=hist[d]*choose(n-1-d,r-1-d);if(f!=want)throw std::runtime_error("literal profile mismatch");std::cout<<"ROW "<<cases++<<" n "<<n<<" words "<<words<<" profile";for(auto v:f)std::cout<<' '<<v;std::cout<<'\n';total+=words;}std::cout<<"PASS cases "<<cases<<" permutations "<<total<<'\n';}
