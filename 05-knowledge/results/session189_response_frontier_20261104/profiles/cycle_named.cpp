#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
uint64_t C(const vector<unsigned>&r){int n=r.size(),x=n-2,full=(1<<n)-1;vector<uint64_t>d((1<<n)*n);d[(1<<x)*n+x]=1;for(int m=1;m<=full;m++){unsigned mm=m;while(mm){int i=__builtin_ctz(mm);mm&=mm-1;auto v=d[m*n+i];if(!v)continue;unsigned z=r[i]&~m;while(z){int j=__builtin_ctz(z);z&=z-1;d[(m|(1<<j))*n+j]+=v;}}}uint64_t o=0;for(int i=0;i<n;i++)if(r[i]>>x&1)o+=d[full*n+i];return o;}
int main(int argc,char**argv){(void)argv;bool all=argc>1;int cases;cin>>cases;for(int k=0;k<cases;k++){int n;cin>>n;vector<unsigned>t(n);for(auto &v:t)cin>>v;vector<int> ms;for(int a=0;a<(1<<n);a++)if(all||__builtin_popcount((unsigned)a)<=2)ms.push_back(a);cout<<n<<" "<<ms.size()<<"\n";for(auto a:ms)for(auto c:ms){auto r=t;r.push_back((1<<(n+1))|(((1<<n)-1)^a));r.push_back(((1<<n)-1)^c);for(int i=0;i<n;i++){if(a>>i&1)r[i]|=1<<n;if(c>>i&1)r[i]|=1<<(n+1);}cout<<C(r)<<"\n";}}}
