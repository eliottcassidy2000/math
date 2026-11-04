#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
#include <functional>
const int STATES=7776;int pw[5]={1,6,36,216,1296};
struct Count{std::array<unsigned,5>Q;unsigned a,c;std::vector<int>memo;Count(std::array<unsigned,5>q,unsigned aa,unsigned cc):Q(q),a(aa),c(cc),memo(STATES*4*7,-1){}
 bool arc(int u,int v){if(u==v)return false;if(u<5&&v<5)return Q[u]>>v&1;if(u==5)return v==6||(v<5&&!(a>>v&1));if(u==6)return v<5&&!(c>>v&1);return (v==5?a:c)>>u&1;}
 int f(int code,int ports,int last){int key=(code*4+ports)*7+last;if(memo[key]>=0)return memo[key];int cc=code,mm=ports;if(last<5){if(!(code/pw[last]%6))return memo[key]=0;cc-=pw[last];}else{int b=1<<(last-5);if(!(ports&b))return memo[key]=0;mm^=b;}if(!cc&&!mm)return memo[key]=1;int z[3]={0,0,0};for(int u=0;u<7;++u){int p=f(cc,mm,u);bool edge=arc(u,last);for(int k=0;k<3;++k){int x=(p>>(4*k))&15;if(edge)z[k]+=x;if(k<2)z[k+1]+=x;}}return memo[key]=(z[0]&15)|((z[1]&15)<<4)|((z[2]&15)<<8);}
 std::array<int,3>word(const std::array<int,5>&h){int code=0;for(int i=0;i<5;++i)code+=pw[i]*h[i];std::array<int,3>z{};for(int u=0;u<7;++u){int p=f(code,3,u);for(int k=0;k<3;++k)z[k]=(z[k]+((p>>(4*k))&15))%16;}return z;}
 std::array<int,48>bank(){std::array<int,48>out{};std::array<int,5>h;h.fill(1);int q[5]={0,1,14,4,8};std::function<void(int,int,int)>go=[&](int i,int cost,int w){if(i==5){auto v=word(h);for(int k=0;k<3;++k)out[k]=(out[k]+w*v[k])%16;return;}for(int j=1;j<=4&&cost+j-1<=3;++j){h[i]=j;go(i+1,cost+j-1,w*q[j]%16);}h[i]=1;};go(0,0,1);int at=3;for(int i=0;i<5;++i)for(int j:{2,3,5}){h.fill(1);h[i]=j;auto v=word(h);for(int k=0;k<3;++k)out[at++]=v[k]%2;}return out;}}
;
int main(){int cases;std::cin>>cases;while(cases--){std::array<unsigned,5>Q;for(auto&x:Q)std::cin>>x;int queries;std::cin>>queries;while(queries--){unsigned a,c;std::cin>>a>>c;for(int rev=0;rev<2;++rev){auto R=Q;if(rev)for(int i=0;i<5;++i)R[i]=31^(1<<i)^R[i];Count z(R,a,c);auto v=z.bank();for(auto x:v)std::cout<<x<<' ';std::cout<<'\n';}}}}
