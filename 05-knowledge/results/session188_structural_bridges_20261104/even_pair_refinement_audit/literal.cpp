#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int paths(const vector<unsigned>&G){unsigned n=G.size(),N=1u<<n;vector<uint8_t>d(size_t(N)*n);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;
for(unsigned s=1;s<N;s++){unsigned tails=s;while(tails){unsigned b=tails&-tails;tails-=b;unsigned v=__builtin_ctz(b),ways=d[size_t(s)*n+v];if(!ways)continue;unsigned next=G[v]&~s;while(next){unsigned t=next&-next;next-=t;auto &z=d[size_t(s|t)*n+__builtin_ctz(t)];z=(z+ways)&7;}}}
unsigned total=0;for(unsigned i=0;i<n;i++)total+=d[size_t(N-1)*n+i];return total&7;}
int main(){vector<unsigned>Q={598,508,344,865,40,69,656,317,113,438};vector<int>v;
for(int flipped: {0,683})for(int converse=0;converse<2;converse++){vector<unsigned>G(22);for(unsigned i=0;i<20;i++)for(unsigned j=i+1;j<20;j++){bool edge;if(i/2==j/2)edge=!(flipped>>(i/2)&1);else edge=bool(Q[i/2]>>(j/2)&1)^bool(converse);G[edge?i:j]|=1u<<(edge?j:i);}for(unsigned i=0;i<20;i++){if(i<2)G[i]|=1u<<20;else G[20]|=1u<<i;G[21]|=1u<<i;}G[20]|=1u<<21;v.push_back(paths(G));}
if(v!=vector<int>({3,7,3,7})){throw runtime_error("independent physical witness");}
cout<<"{\"augmented_order\":22,\"physical_a_mask\":3,\"physical_c_mask\":0,\"pair_flip_masks\":[0,683],\"values_mod8\":[3,7,3,7]}\n";}
