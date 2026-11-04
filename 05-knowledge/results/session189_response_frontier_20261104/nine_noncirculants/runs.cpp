#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
uint64_t count(const vector<unsigned>&R,const vector<int>&k){int n=R.size(),states=1;vector<int>stride(n);for(int i=0;i<n;i++){stride[i]=states;states*=k[i]+1;}vector<uint64_t>dp(states*n);for(int i=0;i<n;i++)dp[stride[i]*n+i]=1;for(int st=1;st<states;st++){unsigned avail=0;for(int j=0;j<n;j++)if(st/stride[j]%(k[j]+1)<k[j])avail|=1<<j;for(int i=0;i<n;i++){auto v=dp[st*n+i];if(!v)continue;unsigned todo=R[i]&avail;while(todo){int j=__builtin_ctz(todo);todo&=todo-1;dp[(st+stride[j])*n+j]+=v;}}}uint64_t z=0;for(int i=0;i<n;i++)z+=dp[(states-1)*n+i];return z;}
uint64_t cycle(const vector<unsigned>&R){int n=R.size(),root=n-2;unsigned full=(1<<n)-1;vector<uint64_t>d((1<<n)*n);d[(1<<root)*n+root]=1;for(unsigned s=1;s<=full;s++)for(int i=0;i<n;i++)if(s>>i&1){auto v=d[s*n+i];if(!v)continue;unsigned t=R[i]&~s;while(t){int j=__builtin_ctz(t);t&=t-1;d[(s|(1<<j))*n+j]+=v;}}uint64_t z=0;for(int i=0;i<n;i++)if(R[i]>>root&1)z+=d[full*n+i];return z;}
int main(){int cases;cin>>cases;while(cases--){int m;cin>>m;vector<unsigned>R(m+2);for(auto&r:R)cin>>r;vector<int>k(m+2,1);cout<<count(R,k);for(int i=0;i<m;i++){k[i]=2;cout<<' '<<count(R,k);k[i]=1;}for(int i=0;i<m;i++){k[i]=3;cout<<' '<<count(R,k);k[i]=1;}for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){k[i]=k[j]=2;cout<<' '<<count(R,k);k[i]=k[j]=1;}cout<<' '<<cycle(R)<<'\n';}}
