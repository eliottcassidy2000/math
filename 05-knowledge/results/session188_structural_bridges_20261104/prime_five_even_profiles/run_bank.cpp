#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <functional>
#include <stdexcept>
const int N=5, B=6, STATES=7776; int pow6[5]={1,6,36,216,1296};
struct Counter{std::array<unsigned,5>Q;unsigned a,c;std::vector<int8_t> memo;Counter(std::array<unsigned,5>q,unsigned aa,unsigned cc):Q(q),a(aa),c(cc),memo(STATES*4*7,-1){}
 bool arc(int u,int v){if(u==v)return false;if(u<5&&v<5)return Q[u]>>v&1;if(u==5)return v==6||(v<5&&!(a>>v&1));if(u==6)return v<5&&!(c>>v&1);return (v==5?a:c)>>u&1;}
 int f(int code,int ports,int last){int key=(code*4+ports)*7+last;if(memo[key]>=0)return memo[key];int next=code,pm=ports;if(last<5){if(!(code/pow6[last]%6))return memo[key]=0;next-=pow6[last];}else{int bit=1<<(last-5);if(!(ports&bit))return memo[key]=0;pm^=bit;}if(next==0&&pm==0)return memo[key]=1;int z=0;for(int u=0;u<7;++u)if(arc(u,last))z+=f(next,pm,u);return memo[key]=z%16;}
 int word(const std::array<int,5>&h){int code=0;for(int i=0;i<5;++i)code+=pow6[i]*h[i];int z=0;for(int u=0;u<7;++u)z+=f(code,3,u);return z%16;}
 std::array<int,16> bank(){std::array<int,16>out{};std::array<int,5>h;h.fill(1);const int p[5]={0,1,14,4,8};std::function<void(int,int,int)>go=[&](int i,int cost,int weight){if(i==5){out[0]=(out[0]+weight*word(h))%16;return;}for(int k=1;k<=4&&cost+k-1<=3;++k){h[i]=k;go(i+1,cost+k-1,weight*p[k]%16);}h[i]=1;};go(0,0,1);int at=1;for(int i=0;i<5;++i)for(int k:{2,3,5}){h.fill(1);h[i]=k;out[at++]=word(h)%2;}return out;}}
;
int main(){int cases;std::cin>>cases;std::vector<unsigned>masks;for(unsigned x=0;x<32;++x)masks.push_back(x);while(cases--){std::array<unsigned,5>Q;for(auto&x:Q)std::cin>>x;for(unsigned a:masks)for(unsigned c:masks)for(int rev=0;rev<2;++rev){auto R=Q;if(rev)for(int i=0;i<5;++i)R[i]=31u^(1u<<i)^R[i];Counter z(R,a,c);auto v=z.bank();for(auto x:v)std::cout<<x<<' ';std::cout<<'\n';}}}
