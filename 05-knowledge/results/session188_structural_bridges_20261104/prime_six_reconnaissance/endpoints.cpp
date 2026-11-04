#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
vector<unsigned> ends(const vector<unsigned>&Q){unsigned n=Q.size(),N=1u<<n;vector<unsigned>d(size_t(N)*n);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;for(unsigned S=1;S<N;S++)for(unsigned i=0;i<n;i++){unsigned v=d[size_t(S)*n+i];if(!v)continue;unsigned js=Q[i]&~S;while(js){unsigned b=js&-js;js-=b;d[size_t(S|b)*n+__builtin_ctz(b)]+=v;}}vector<unsigned>out(n);for(unsigned i=0;i<n;i++)out[i]=d[size_t(N-1)*n+i];return out;}
int main(){for(unsigned code=0;code<32768;code++){vector<unsigned>Q(6);int b=0;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,b++)if(code>>b&1)Q[i]|=1u<<j;else Q[j]|=1u<<i;unsigned witness=0;for(unsigned M=1;M<63;M++){if(__builtin_popcount(M)<2)continue;bool ok=true;for(int i=0;i<6;i++)if(!(M>>i&1)){unsigned a=Q[i]&M;if(a&&a!=M){ok=false;break;}}if(ok){witness=M;break;}}auto O=Q;for(int i=0;i<6;i++)O[i]=63^(1<<i)^Q[i];auto S=ends(O),E=ends(Q);unsigned p=0;for(int i=0;i<6;i++)if((S[i]+E[i])%2)p|=1u<<i;cout<<code<<' '<<witness<<' '<<p;for(auto v:S)cout<<' '<<v;for(auto v:E)cout<<' '<<v;cout<<'\n';}}
