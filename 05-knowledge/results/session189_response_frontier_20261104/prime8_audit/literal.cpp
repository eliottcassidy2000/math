#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
int main(){int count;cin>>count;while(count--){int n;cin>>n;if(n>20)return 2;vector<unsigned>r(n);for(auto&x:r)cin>>x;int all=(1<<n)-1;vector<uint64_t>f((1<<n)*n);for(int i=0;i<n;i++)f[(1<<i)*n+i]=1;for(int S=1;S<=all;S++){if(!(S&(S-1)))continue;unsigned vv=S;while(vv){int v=__builtin_ctz(vv);vv&=vv-1;int R=S^(1<<v);unsigned w=r[v]&R;uint64_t z=0;while(w){int j=__builtin_ctz(w);w&=w-1;z+=f[R*n+j];}f[S*n+v]=z;}}uint64_t H=0;for(int v=0;v<n;v++)H+=f[all*n+v];cout<<H<<"\n";}}
