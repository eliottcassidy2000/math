// Independent subset endpoint convolution, all arithmetic checked in chosen controls.
#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>
using namespace std;
int main(){int n;unsigned sel;while(cin>>n>>sel){if(n<1||n>20)throw runtime_error("literal control order");vector<unsigned> row(n),pred(n);for(auto &x:row)cin>>x;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(row[j]>>i&1)pred[i]|=1u<<j;unsigned N=1u<<n;vector<long long> e((size_t)N*n),s((size_t)N*n),h(N),ee(N),ss(N);for(int v=0;v<n;v++)e[(size_t)(1u<<v)*n+v]=s[(size_t)(1u<<v)*n+v]=1;for(unsigned mask=1;mask<N;mask++){for(unsigned bits=mask;bits;bits&=bits-1){int v=__builtin_ctz(bits);unsigned r=mask^(1u<<v);if(r){for(unsigned b=r&pred[v];b;b&=b-1)e[(size_t)mask*n+v]+=e[(size_t)r*n+__builtin_ctz(b)];for(unsigned b=r&row[v];b;b&=b-1)s[(size_t)mask*n+v]+=s[(size_t)r*n+__builtin_ctz(b)];}h[mask]+=e[(size_t)mask*n+v];if(sel>>v&1){ee[mask]+=e[(size_t)mask*n+v];ss[mask]+=s[(size_t)mask*n+v];}}}__int128 W=0,U=ss[N-1];for(unsigned mask=1;mask<N-1;mask++){W+=(__int128)ee[mask]*ss[(N-1)^mask];U+=(__int128)h[mask]*ss[(N-1)^mask];}if(W>numeric_limits<long long>::max()||U>numeric_limits<long long>::max())throw runtime_error("overflow");cout<<h[N-1]<<" "<<(long long)(h[N-1]-W)<<" "<<(long long)U<<"\n";}}
