#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(){unsigned n,a,c;if(!(cin>>n>>a>>c))return 1;if(n+2>17)throw runtime_error("control order exceeds exact declared range");vector<unsigned> R(n+2);for(unsigned i=0;i<n;i++)cin>>R[i];R[n]=1u<<(n+1);for(unsigned i=0;i<n;i++){if((a>>i)&1u)R[i]|=1u<<n;else R[n]|=1u<<i;if((c>>i)&1u)R[i]|=1u<<(n+1);else R[n+1]|=1u<<i;}unsigned m=n+2,N=1u<<m;vector<uint64_t>D(size_t(N)*m);for(unsigned i=0;i<m;i++)D[size_t(1u<<i)*m+i]=1;for(unsigned S=1;S<N;S++)for(unsigned u=0;u<m;u++){auto v=D[size_t(S)*m+u];if(!v)continue;unsigned E=R[u]&~S;while(E){unsigned w=__builtin_ctz(E);E&=E-1;D[size_t(S|(1u<<w))*m+w]+=v;}}uint64_t h=0;for(unsigned u=0;u<m;u++)h+=D[size_t(N-1)*m+u];cout<<h<<'\n';}
