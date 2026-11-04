#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
using I=long long;
I H(const vector<int>& old,int a,int c){
 int n=old.size(),N=n+2,full=(1<<N)-1; vector<int> rows=old; rows.resize(N); rows[n]=1<<(n+1);
 for(int i=0;i<n;++i){if(a>>i&1)rows[i]|=1<<n;else rows[n]|=1<<i;if(c>>i&1)rows[i]|=1<<(n+1);else rows[n+1]|=1<<i;}
 vector<I>d((1<<N)*N,0);for(int i=0;i<N;++i)d[(1<<i)*N+i]=1;
 for(int s=1;s<=full;++s){int v=s;while(v){int i=__builtin_ctz(v);v&=v-1;I z=d[s*N+i];if(!z)continue;int w=rows[i]&~s;while(w){int j=__builtin_ctz(w);w&=w-1;d[(s|(1<<j))*N+j]+=z;}}}
 I ans=0;for(int i=0;i<N;++i)ans+=d[full*N+i];return ans;
}
template<class T>void array_out(ostream&o,const vector<T>&a){o<<"[";for(size_t i=0;i<a.size();++i){if(i)o<<",";o<<a[i];}o<<"]";}
int main(){int count;cin>>count;for(int k=0;k<count;++k){int n,a,c;cin>>n;vector<int>r(n);for(int&i:r)cin>>i;cin>>a>>c;cout<<H(r,a,c)<<"\n";}}
