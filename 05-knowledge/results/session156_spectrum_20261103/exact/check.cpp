#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
using namespace std;
uint64_t dp_count(int d){string w="10011"+string(d,'0'); int m=w.size(), n=m+1, lim=1<<n; vector<unsigned> out(n); for(int i=0;i<m;i++){for(int j=i+1;j<m;j++)out[i]|=1u<<j; if(w[i]=='1')out[m]|=1u<<i; else out[i]|=1u<<m;} vector<uint64_t> dp((size_t)lim*n,0);for(int i=0;i<n;i++)dp[((1<<i)*n)+i]=1;for(int s=1;s<lim;s++)for(int v=0;v<n;v++)if((s>>v)&1){auto val=dp[(size_t)s*n+v]; if(!val)continue; unsigned avail=out[v]&(~(unsigned)s)&(lim-1);while(avail){int u=__builtin_ctz(avail);avail&=avail-1;dp[((size_t)(s|1<<u))*n+u]+=val;}}uint64_t t=0;for(int v=0;v<n;v++)t+=dp[((size_t)lim-1)*n+v];return t;}
int main(){for(int d: {2,5,8,11})cout<<"d="<<d<<" H="<<dp_count(d)<<"\n";}
