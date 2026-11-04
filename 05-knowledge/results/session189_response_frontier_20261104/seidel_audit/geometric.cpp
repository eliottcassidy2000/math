#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>
using V=std::array<int,5>;const int MOD=256,BASE=9;
struct Engine{int n,N;std::vector<int>pw,codes,id,weight;std::vector<std::vector<int>>prev;std::vector<unsigned>Q;unsigned a,c;std::vector<V>memo;std::vector<bool>done;std::vector<std::vector<bool>>arc;
 Engine(std::vector<unsigned>q,unsigned aa,unsigned cc):n(q.size()),N(n+2),Q(q),a(aa),c(cc){int cap=1;for(int i=0;i<n;++i){pw.push_back(cap);cap*=BASE;}id.assign(cap,-1);for(int code=0;code<cap;++code){int cost=0,x=code,w=1;for(int i=0;i<n;++i){int k=x%BASE;x/=BASE;if(k==0)w=0;else for(int j=1;j<k;++j)w=w*(-2)%MOD;cost+=k>0?k-1:0;}if(cost<8){id[code]=codes.size();codes.push_back(code);weight.push_back((w+MOD)%MOD);}}
 for(auto code:codes){std::vector<int>p(n);for(int i=0;i<n;++i)p[i]=(code/pw[i]%BASE)?id[code-pw[i]]:-1;prev.push_back(p);}memo.resize(codes.size()*4*N);done.resize(codes.size()*4*N);arc.assign(N,std::vector<bool>(N));for(int u=0;u<N;++u)for(int v=0;v<N;++v){bool b=false;if(u!=v){if(u<n&&v<n)b=Q[u]>>v&1;else if(u==n)b=v==n+1||(v<n&&!(a>>v&1));else if(u==n+1)b=v<n&&!(c>>v&1);else b=(v==n?a:c)>>u&1;}arc[u][v]=b;}}
 const V& f(int z,int ports,int last){int key=(z*4+ports)*N+last;if(done[key])return memo[key];done[key]=true;auto&out=memo[key];int pz=z,pp=ports;if(last<n){pz=prev[z][last];if(pz<0)return out;}else{int b=1<<(last-n);if(!(ports&b))return out;pp^=b;}if(codes[pz]==0&&pp==0){out[0]=1;return out;}for(int u=0;u<N;++u){const auto&v=f(pz,pp,u);for(int h=0;h<5;++h)out[h]=(out[h]+(arc[u][last]?v[h]:0)+(h?v[h-1]:0))%MOD;}return out;}
 V sum(int ports){V out{};for(size_t z=0;z<codes.size();++z)if(weight[z])for(int u=0;u<N;++u){const auto&v=f(z,ports,u);for(int h=0;h<5;++h)out[h]=(out[h]+weight[z]*v[h])%MOD;}return out;}
};
int main(){int cases;std::cin>>cases;while(cases--){int n;std::cin>>n;if(n>4)throw std::runtime_error("audit bound4");std::vector<unsigned>Q(n);for(auto&r:Q)std::cin>>r;for(int side=0;side<2;++side){auto R=Q;if(side)for(int i=0;i<n;++i)R[i]=((1u<<n)-1)^(1u<<i)^R[i];Engine e(R,0,0);auto v=e.sum(0);std::cout<<side<<" -1 -1";for(int z:v)std::cout<<' '<<z;std::cout<<'\n';for(unsigned a=0;a<(1u<<n);++a)for(unsigned c=0;c<(1u<<n);++c){Engine g(R,a,c);auto w=g.sum(3);std::cout<<side<<' '<<a<<' '<<c;for(int z:w)std::cout<<' '<<z;std::cout<<'\n';}}}}
