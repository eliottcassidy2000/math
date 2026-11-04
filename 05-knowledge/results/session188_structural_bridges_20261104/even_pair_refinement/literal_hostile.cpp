#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
int count_mod8(vector<unsigned>R){unsigned n=R.size(),N=1u<<n;vector<unsigned>incoming(n);for(unsigned i=0;i<n;i++)for(unsigned j=0;j<n;j++)if(R[i]>>j&1)incoming[j]|=1u<<i;vector<unsigned char>d((size_t)N*n);for(unsigned i=0;i<n;i++)d[(size_t)(1u<<i)*n+i]=1;
for(unsigned s=1;s<N;s++){if(!(s&(s-1)))continue;unsigned lastbits=s;while(lastbits){unsigned b=lastbits&-lastbits;lastbits-=b;unsigned j=__builtin_ctz(b),r=s^b,preds=incoming[j]&r;unsigned value=0;while(preds){unsigned t=preds&-preds;preds-=t;value+=d[(size_t)r*n+__builtin_ctz(t)];}d[(size_t)s*n+j]=value&7;}}
unsigned value=0;for(unsigned j=0;j<n;j++)value+=d[(size_t)(N-1)*n+j];return value&7;}
int main(){vector<unsigned>Q={598,508,344,865,40,69,656,317,113,438};vector<int>values;for(int theta=0;theta<2;theta++){vector<unsigned>R(22);for(unsigned i=0;i<10;i++){R[2*i]|=1u<<(2*i+1);unsigned row=theta?(1023^Q[i]^(1u<<i)):Q[i];for(unsigned j=0;j<10;j++)if(row>>j&1){R[2*i]|=3u<<(2*j);R[2*i+1]|=3u<<(2*j);}}R[20]=1u<<21;for(unsigned i=0;i<20;i++){if(i<2)R[i]|=1u<<20;else R[20]|=1u<<i;R[21]|=1u<<i;}values.push_back(count_mod8(R));}if(values!=vector<int>({3,7}))throw runtime_error("literal physical twenty-name separator");cout<<"{\"status\":\"PASS\",\"old_names\":20,\"augmented_names\":22,\"a_mask\":3,\"c_mask\":0,\"source_mod8\":"<<values[0]<<",\"target_mod8\":"<<values[1]<<",\"method\":\"Literal predecessor subset DP on all twenty-two physical names\"}\n";}
