#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
int main(){int cases;cin>>cases;for(int t=0;t<cases;t++){int n;cin>>n;vector<unsigned>r(n);for(auto &x:r)cin>>x;int N=1<<n;vector<uint64_t>dp((size_t)N*n);for(int mask=1;mask<N;mask++){unsigned todo=mask;while(todo){int j=__builtin_ctz(todo);todo&=todo-1;int prev=mask^(1<<j);uint64_t z=prev?0:1;unsigned p=prev;while(p){int i=__builtin_ctz(p);p&=p-1;if(r[i]>>j&1)z+=dp[(size_t)prev*n+i];}dp[(size_t)mask*n+j]=z;}}uint64_t ans=0;for(int j=0;j<n;j++)ans+=dp[(size_t)(N-1)*n+j];cout<<ans<<'\n';}}
