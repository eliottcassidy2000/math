#include <array>
#include <vector>
#include <iostream>
#include <cstdint>
using V=std::array<uint64_t,9>;
V count(const std::vector<unsigned>&Q,int distinguished,int copies){int n=Q.size(),states=1;std::vector<int>stride(n),cap(n,1);cap[distinguished]=copies;for(int i=0;i<n;++i){stride[i]=states;states*=cap[i]+1;}std::vector<V>dp(states*n);for(int i=0;i<n;++i)dp[stride[i]*n+i][0]=1;for(int st=1;st<states;++st)for(int i=0;i<n;++i){const auto&v=dp[st*n+i];for(int j=0;j<n;++j)if(st/stride[j]%(cap[j]+1)<cap[j]){auto&w=dp[(st+stride[j])*n+j];for(int h=0;h<9;++h)w[h]+=((Q[i]>>j&1)?v[h]:0)+(h?v[h-1]:0);}}V out{};for(int i=0;i<n;++i)for(int h=0;h<9;++h)out[h]+=dp[(states-1)*n+i][h];return out;}
int main(){for(int n=2;n<=5;++n)for(unsigned code=0;code<(1u<<(n*(n-1)/2));++code){std::vector<unsigned>Q(n);int bit=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)if(code>>bit&1)Q[i]|=1<<j;else Q[j]|=1<<i;for(int i=0;i<n;++i)for(int k:{2,3,5}){auto a=count(Q,i,k);std::cout<<n<<' '<<code<<' '<<i<<' '<<k;for(auto x:a)std::cout<<' '<<x;std::cout<<'\n';}}}
