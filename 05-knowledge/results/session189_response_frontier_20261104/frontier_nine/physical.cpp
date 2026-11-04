#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
uint64_t H(const vector<unsigned>& r){int n=r.size();int full=(1<<n)-1;vector<uint64_t>d((1<<n)*n);for(int i=0;i<n;i++)d[(1<<i)*n+i]=1;for(int m=1;m<=full;m++){unsigned mm=m;while(mm){int i=__builtin_ctz(mm);mm&=mm-1;auto v=d[m*n+i];if(!v)continue;unsigned z=r[i]&~m;while(z){int j=__builtin_ctz(z);z&=z-1;d[(m|(1<<j))*n+j]+=v;}}}uint64_t out=0;for(int i=0;i<n;i++)out+=d[full*n+i];return out;}
int main(int argc,char**argv){(void)argv;bool all=argc>1;int cases;cin>>cases;for(int k=0;k<cases;k++){int n;cin>>n;vector<unsigned> t(n);for(auto &v:t)cin>>v;vector<int> masks;for(int a=0;a<(1<<n);a++)if(all||__builtin_popcount((unsigned)a)<=2)masks.push_back(a);cout<<n<<" "<<masks.size()<<"\n";for(auto a:masks)for(auto c:masks){vector<unsigned> r=t;r.push_back((1<<(n+1))|(((1<<n)-1)^a));r.push_back(((1<<n)-1)^c);for(int i=0;i<n;i++){if(a>>i&1)r[i]|=1<<n;if(c>>i&1)r[i]|=1<<(n+1);}cout<<H(r)<<"\n";}cerr<<k+1<<"/"<<cases<<"\n";}}
