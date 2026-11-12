#include <vector>
#include <cstdint>
extern "C" long long defect(const long long* r,int n,int O){int N=1<<n;std::vector<long long>s(N*n),e(N*n);for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){int t=m^(1<<v);if(!t)s[m*n+v]=e[m*n+v]=1;else for(int u=0;u<n;u++){if(r[v]>>u&1)s[m*n+v]+=s[t*n+u];if(r[u]>>v&1)e[m*n+v]+=e[t*n+u];}}long long h=0,w=0;for(int v=0;v<n;v++)h+=e[(N-1)*n+v];for(int m=1;m<N-1;m++){long long a=0,b=0;for(int v=0;v<n;v++)if(O>>v&1){a+=e[m*n+v];b+=s[((N-1)^m)*n+v];}w+=a*b;}return h-w;}
