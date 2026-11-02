#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
using I=std::int64_t;using U=std::uint64_t;
U choose(int n,int r){if(r<0||r>n)return 0;U z=1;for(int j=1;j<=r;++j)z=z*U(n-j+1)/U(j);return z;}
int main(){U totalT=0,totalWords=0,totalPairs=0,totalProfiles=0;for(int n=1;n<=6;++n){unsigned N=1u<<(n*(n-1)/2);std::set<std::vector<I>>profiles;for(unsigned code=0;code<N;++code){std::vector<unsigned>a(n);unsigned c=code;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(c&1)a[i]|=1u<<j;else a[j]|=1u<<i;c>>=1;}std::vector<int>p(n);std::iota(p.begin(),p.end(),0);std::vector<U>hist(n);do{int d=0;for(int j=0;j+1<n;++j)d+=!((a[p[j]]>>p[j+1])&1u);++hist[d];++totalWords;}while(std::next_permutation(p.begin(),p.end()));std::vector<I>f(n);for(int r=1;r<=n;++r)for(int d=0;d<r;++d)f[r-1]+=I(hist[d]*choose(n-1-d,r-1-d));profiles.insert(f);}U pairs=0;for(auto it=profiles.begin();it!=profiles.end();++it){const auto&f=*it;for(int i=1;i+1<n;++i)if(f[i]*f[i]<f[i-1]*f[i+1])throw std::runtime_error("non-log-concave actual small profile");auto jt=it;for(++jt;jt!=profiles.end();++jt){std::vector<I>d(n+3);for(int i=0;i<n;++i)d[i+1]=(*jt)[i]-f[i];bool pos=false,neg=false;for(int i=1;i<=n+1;++i){I z=d[i-1]+3*d[i]+d[i+1];pos|=z>0;neg|=z<0;if(i<=n){z=d[i]+d[i+1];pos|=z>0;neg|=z<0;}}if(pos&&neg)throw std::runtime_error("mixed-sign actual small kernel");++pairs;}}totalT+=N;totalPairs+=pairs;totalProfiles+=profiles.size();std::cout<<"ROW n "<<n<<" tournaments "<<N<<" profiles "<<profiles.size()<<" pairs "<<pairs<<" all_kernels_one_signed 1 all_profiles_log_concave 1\n";}std::cout<<"PASS tournaments "<<totalT<<" permutations "<<totalWords<<" profiles "<<totalProfiles<<" pairs "<<totalPairs<<'\n';}
