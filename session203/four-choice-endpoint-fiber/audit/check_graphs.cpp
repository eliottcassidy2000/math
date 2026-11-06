// Independent append-at-end subset DP, arbitrary-root profiles and physical colored body construction.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using G=std::vector<unsigned>;using I=std::uint64_t;using S=std::int64_t;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void valid(const G&a){int n=a.size();need(n>0&&n<=18,"order bound");for(int i=0;i<n;++i){need(a[i]<(1u<<n)&&!(a[i]>>i&1),"self or outside arc");for(int j=0;j<i;++j)need((a[i]>>j&1)+(a[j]>>i&1)==1,"tournament pair");}}
using Profile=std::array<I,9>; // h,s,e,L,R,d,beta,B,V
Profile profile(const G&a,int root){valid(a);int n=a.size(),N=1<<n,F=N-1;std::vector<I>all(std::size_t(N)*n),marked(std::size_t(N)*n),h(N),s(N),e(N);h[0]=1;
 for(int i=0;i<n;++i)all[std::size_t(1<<i)*n+i]=1;
 marked[std::size_t(1<<root)*n+root]=1;
 for(int mask=1;mask<N;++mask){for(unsigned ends=mask;ends;ends&=ends-1){int last=__builtin_ctz(ends);std::size_t index=std::size_t(mask)*n+last;I v=all[index],w=marked[index];h[mask]+=v;s[mask]+=w;if(last==root)e[mask]=v;if(!v&&!w)continue;for(unsigned next=a[last]&~mask;next;next&=next-1){int u=__builtin_ctz(next);std::size_t target=std::size_t(mask|(1<<u))*n+u;all[target]+=v;marked[target]+=w;}}}
 I L=0,R=0,beta=0;for(int mask=1;mask<N;++mask){L+=e[mask]*h[F^mask];R+=s[mask]*h[F^mask];}int body=F^(1<<root);for(int mask=body;;mask=(mask-1)&body){beta+=h[mask]*h[body^mask];if(!mask)break;}return{h[F],s[F],e[F],L,R,h[body],beta,h[F]+L-s[F],R+beta};}
G context(const G&body,unsigned colors){G a={102,0,10,83,135,30,182,47};std::vector<int>core;for(int i=0;i<8;++i)if(i!=5)core.push_back(i);G out(7+body.size());for(int i=0;i<7;++i)for(int j=0;j<7;++j)if(a[core[i]]>>core[j]&1)out[i]|=1u<<j;for(int j=0;j<(int)body.size();++j){for(int i=0;i<7;++i){bool win=core[i]==3?bool(colors>>j&1):bool(a[core[i]]>>5&1);if(win)out[i]|=1u<<(7+j);else out[7+j]|=1u<<i;}for(int k=0;k<(int)body.size();++k)if(body[j]>>k&1)out[7+j]|=1u<<(7+k);}valid(out);return out;}
I power(I x,int k){I r=1;while(k--)r*=x;return r;}
std::array<int,8> eigen={0,1,2,4,6,9,12,16};
std::array<std::array<S,8>,9>coeff={{
 {{1,9,-341,2480,-6600,11664,0,0}},{{0,8,-110,456,-1100,1944,0,0}},{{0,0,0,0,0,0,0,0}},
 {{0,16,-462,3352,-7700,13608,0,0}},{{8,-166,132,4832,-26400,71244,-91224,110592}},
 {{0,11,-132,648,-1100,1944,0,0}},{{8,-156,-528,10240,-39600,94572,-91224,110592}},
 {{1,17,-693,5376,-13200,23328,0,0}},{{16,-322,-396,15072,-66000,165816,-182448,221184}}
}};
int main(){int graphs=0;std::cout<<"{\"status\":\"PASS\",\"graphs\":[";for(int k=0;k<=3;++k){int m=5+2*k;G body(m);for(int i=0;i<m;++i)for(int j=i+1;j<m;++j)body[i]|=1u<<j;int seed[4]={10,14,17,21},offset[4]={0,8,14,22};for(int t=0;t<4;++t){unsigned word=(seed[t]<<k)|((1u<<k)-1);G a=context(body,word);Profile q=profile(a,0);for(int i=0;i<9;++i){S want=0;for(int j=0;j<8;++j)want+=coeff[i][j]*S(power(eigen[j],k));if(i<=1)want+=offset[t]*S(power(2,k));need(want>=0&&q[i]==I(want),"literal graph versus full coefficient table");}need(a[1]==0&&q[2]==0,"global sink retained");std::cout<<(graphs?",":"")<<"{\"k\":"<<k<<",\"seed\":"<<seed[t]<<",\"vertices\":"<<a.size()<<",\"profile_h_s_e_L_R_d_beta_B_V\":[";for(int i=0;i<9;++i)std::cout<<(i?",":"")<<q[i];std::cout<<"]}";++graphs;}}
 std::cout<<"],\"literal_graphs\":"<<graphs<<",\"method\":\"Append-at-end DP with independent construction, separate root-start DP, and cover sums; includes all nine coordinates.\",\"coverage_claim\":false}\n";
}
