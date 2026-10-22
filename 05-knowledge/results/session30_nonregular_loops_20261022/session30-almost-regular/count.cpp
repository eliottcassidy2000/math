#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std; using U=unsigned __int128;
string decimal(U x){string s;do{s+=char('0'+x%10);x/=10;}while(x);reverse(s.begin(),s.end());return s;}
int main(){int n;string bits;while(cin>>n>>bits){vector<unsigned> in(n,0);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){if(bits[k++]=='1')in[j]|=1u<<i;else in[i]|=1u<<j;}size_t sz=1u<<n;vector<U> dp(sz*n,0);for(int i=0;i<n;i++)dp[(1u<<i)*n+i]=1;for(unsigned m=1;m<sz;m++){if(!(m&(m-1)))continue;unsigned js=m;while(js){int j=__builtin_ctz(js);js&=js-1;unsigned p=m^(1u<<j),ks=p&in[j];U v=0;while(ks){int k=__builtin_ctz(ks);ks&=ks-1;v+=dp[(size_t)p*n+k];}dp[(size_t)m*n+j]=v;}}U total=0;for(int j=0;j<n;j++)total+=dp[(sz-1)*n+j];cout<<decimal(total)<<endl;}}
