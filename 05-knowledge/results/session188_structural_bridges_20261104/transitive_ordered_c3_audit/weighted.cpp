#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using U=unsigned __int128;
std::string decimal(U z){std::string s;do{s.push_back(char('0'+z%10));z/=10;}while(z);std::reverse(s.begin(),s.end());return s;}
U binom(unsigned n,unsigned k){U z=1;for(unsigned i=1;i<=k;++i)z=z*(n+1-i)/i;return z;}
U surj(unsigned n,unsigned k){__int128 z=0;for(unsigned j=0;j<=k;++j){U p=1;for(unsigned t=0;t<n;++t)p*=j;__int128 v=__int128(binom(k,j)*p);z+=((k-j)&1)?-v:v;}return U(z);}
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n,a,b;std::cin>>n>>a>>b;if(n>11||a>8||b>8||!a||!b)return 2;unsigned m=n+2,N=1u<<n;std::vector<unsigned> rows(m),pred(m);for(auto&x:rows)std::cin>>x;for(unsigned u=0;u<m;++u)for(unsigned v=0;v<m;++v)if(rows[u]>>v&1)pred[v]|=1u<<u;std::vector<U>D(size_t(N)*(a+1)*(b+1)*m);auto at=[&](unsigned M,unsigned i,unsigned j,unsigned v)->U&{return D[(((size_t(M)*(a+1)+i)*(b+1)+j)*m)+v];};for(unsigned M=0;M<N;++M)for(unsigned i=0;i<=a;++i)for(unsigned j=0;j<=b;++j){unsigned length=__builtin_popcount(M)+i+j;if(!length)continue;for(unsigned v=0;v<m;++v){unsigned R=M,k=i,l=j;if(v<n){if(!(M>>v&1))continue;R^=1u<<v;}else if(v==n){if(!i)continue;--k;}else{if(!j)continue;--l;}if(length==1){at(M,i,j,v)=1;continue;}unsigned options=pred[v];U z=0;while(options){unsigned u=__builtin_ctz(options);options&=options-1;if((u<n&&(R>>u&1))||(u==n&&k)||(u==n+1&&l))z+=at(R,k,l,u);}at(M,i,j,v)=z;}}U answer=0;for(unsigned i=1;i<=a;++i)for(unsigned j=1;j<=b;++j){U endpoint=0;for(unsigned v=0;v<m;++v)endpoint+=at(N-1,i,j,v);answer+=endpoint*surj(a,i)*surj(b,j);}std::cout<<decimal(answer)<<'\n';}}
