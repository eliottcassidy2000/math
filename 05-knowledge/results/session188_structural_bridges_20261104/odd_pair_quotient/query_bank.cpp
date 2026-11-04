#include <iostream>
#include <vector>
using namespace std;
unsigned long long H(const vector<unsigned>&R){unsigned n=R.size(),N=1u<<n;vector<unsigned long long>d(N*n);for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;for(unsigned S=1;S<N;S++)for(unsigned u=0;u<n;u++){auto z=d[S*n+u];if(!z)continue;unsigned left=R[u]&~S;while(left){unsigned v=__builtin_ctz(left),b=1u<<v;left^=b;d[(S|b)*n+v]+=z;}}unsigned long long h=0;for(unsigned v=0;v<n;v++)h+=d[(N-1)*n+v];return h;}
int main(){unsigned n,m;if(!(cin>>n>>m))return 1;vector<unsigned>T(n);for(auto&r:T)cin>>r;cout<<"[";for(unsigned k=0;k<m;k++){unsigned a,c;cin>>a>>c;auto R=T;R.push_back(1u<<(n+1));R.push_back(0);for(unsigned u=0;u<n;u++){if(a>>u&1)R[u]|=1u<<n;else R[n]|=1u<<u;if(c>>u&1)R[u]|=1u<<(n+1);else R[n+1]|=1u<<u;}if(k)cout<<",";cout<<H(R);}cout<<"]\n";}
