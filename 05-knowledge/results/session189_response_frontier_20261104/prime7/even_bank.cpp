#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdint>
using V=std::array<unsigned char,5>;const int N=7, L=9, CODES=279936;std::array<int,N> pw={1,6,36,216,1296,7776,46656};
std::vector<int> codes,id(CODES,-1);std::vector<std::array<int,N>> counts;std::vector<int> baseline;std::vector<std::array<int,N>> prevs;
void setup(){for(int code=0;code<CODES;++code){int c=code,cost=0;std::array<int,N>ct;for(int i=0;i<N;++i){ct[i]=c%6;c/=6;cost+=std::max(0,ct[i]-1);}
 if(cost<=3 || (cost==4 && std::count(ct.begin(),ct.end(),5)==1 && std::all_of(ct.begin(),ct.end(),[](int k){return k<=1||k==5;}))){id[code]=codes.size();codes.push_back(code);counts.push_back(ct);}}
 for(size_t z=0;z<codes.size();++z){std::array<int,N>pr;int weight=1;for(int i=0;i<N;++i){pr[i]=counts[z][i]?id[codes[z]-pw[i]]:-1;int r[]={0,1,14,4,8,0};weight=weight*r[counts[z][i]]%16;}prevs.push_back(pr);baseline.push_back(weight);}}
struct Counter{std::array<unsigned,N> Q;unsigned a,c;std::vector<V> memo;std::vector<unsigned char> done;bool arc[L][L];
 Counter(std::array<unsigned,N>q,unsigned aa,unsigned cc):Q(q),a(aa),c(cc),memo(codes.size()*4*L),done(codes.size()*4*L){for(int u=0;u<L;++u)for(int v=0;v<L;++v){bool z=false;if(u!=v){if(u<N&&v<N)z=Q[u]>>v&1;else if(u==N)z=v==N+1||(v<N&&!(a>>v&1));else if(u==N+1)z=v<N&&!(c>>v&1);else z=(v==N?a:c)>>u&1;}arc[u][v]=z;}}
 const V& f(int z,int ports,int last){int key=(z*4+ports)*L+last;if(done[key])return memo[key];done[key]=1;V&out=memo[key];int pz=z,pp=ports;if(last<N){pz=prevs[z][last];if(pz<0)return out;}else {int bit=1<<(last-N);if(!(ports&bit))return out;pp^=bit;}
 if(codes[pz]==0&&pp==0){out[0]=1;return out;}for(int u=0;u<L;++u){const V&v=f(pz,pp,u);for(int h=0;h<5;++h){int d=(arc[u][last]?v[h]:0)+(h?v[h-1]:0);out[h]=(out[h]+d)%16;}}return out;}
 V word(int code){int z=id[code];V out{};for(int u=0;u<L;++u){const V&v=f(z,3,u);for(int h=0;h<5;++h)out[h]=(out[h]+v[h])%16;}return out;}
 void output(){V b{};for(size_t z=0;z<codes.size();++z)if(baseline[z]){V v=word(codes[z]);for(int h=0;h<5;++h)b[h]=(b[h]+baseline[z]*v[h])%16;}for(int h=0;h<5;++h)std::cout<<int(b[h])<<' ';int one=0;for(auto p:pw)one+=p;for(int k:{2,3,5})for(int i=0;i<N;++i){V v=word(one+(k-1)*pw[i]);for(int h=0;h<5;++h)std::cout<<int(v[h]%2)<<' ';}std::cout<<'\n';}
};
int main(){setup();int cases,nq;std::cin>>cases>>nq;std::vector<std::array<unsigned,2>> queries(nq);for(auto&v:queries)std::cin>>v[0]>>v[1];while(cases--){std::array<unsigned,N>Q;for(auto&r:Q)std::cin>>r;for(auto v:queries)for(int side=0;side<2;++side){auto R=Q;if(side)for(int i=0;i<N;++i)R[i]=127u^(1u<<i)^R[i];Counter x(R,v[0],v[1]);x.output();}}}
