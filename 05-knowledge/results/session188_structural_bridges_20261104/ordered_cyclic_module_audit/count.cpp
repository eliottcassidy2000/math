#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using Count=uint64_t;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
Count count(vector<unsigned>old,unsigned a,unsigned c){
 unsigned n=old.size(),m=n+2;need(m<=17,"exact uint64 boundary");old.push_back(1u<<(n+1));old.push_back(0);
 for(unsigned i=0;i<n;i++){if(a>>i&1)old[i]|=1u<<n;else old[n]|=1u<<i;if(c>>i&1)old[i]|=1u<<(n+1);else old[n+1]|=1u<<i;}
 vector<unsigned>incoming(m);for(unsigned i=0;i<m;i++)for(unsigned j=0;j<m;j++)if(old[i]>>j&1)incoming[j]|=1u<<i;
 unsigned N=1u<<m;vector<Count>d(size_t(N)*m);for(unsigned i=0;i<m;i++)d[size_t(1u<<i)*m+i]=1;
 for(unsigned s=1;s<N;s++){if(!(s&(s-1)))continue;unsigned last=s;while(last){unsigned b=last&-last;last-=b;unsigned i=__builtin_ctz(b),rest=s^b,pred=incoming[i]&rest;Count z=0;while(pred){unsigned t=pred&-pred;pred-=t;z+=d[size_t(rest)*m+__builtin_ctz(t)];}d[size_t(s)*m+i]=z;}}
 Count answer=0;for(unsigned i=0;i<m;i++)answer+=d[size_t(N-1)*m+i];return answer;
}
int main(){unsigned cases;cin>>cases;while(cases--){unsigned n,k;cin>>n>>k;vector<unsigned>T(n),U(n);for(auto&x:T)cin>>x;for(auto&x:U)cin>>x;while(k--){unsigned a,c;cin>>a>>c;cout<<count(T,a,c)<<' '<<count(U,a,c)<<'\n';}}}
