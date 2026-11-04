#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
using R=std::vector<unsigned>;void need(bool x,const char*s){if(!x)throw std::runtime_error(s);}
std::uint64_t H(const R&a){unsigned n=a.size(),N=1u<<n;std::vector<std::uint64_t>d(size_t(N)*n);for(unsigned u=0;u<n;u++)d[size_t(1u<<u)*n+u]=1;for(unsigned mask=1;mask<N;mask++)for(unsigned u=0;u<n;u++){auto z=d[size_t(mask)*n+u];if(!z)continue;unsigned next=a[u]&~mask;while(next){unsigned b=next&-next;next-=b;d[size_t(mask|b)*n+__builtin_ctz(b)]+=z;}}std::uint64_t out=0;for(unsigned u=0;u<n;u++)out+=d[size_t(N-1)*n+u];return out;}
R addports(R a,unsigned aa,unsigned cc){unsigned n=a.size(),x=n,y=n+1;a.push_back(1u<<y);a.push_back(0);for(unsigned u=0;u<n;u++){if(aa&(1u<<u))a[u]|=1u<<x;else a[x]|=1u<<u;if(cc&(1u<<u))a[u]|=1u<<y;else a[y]|=1u<<u;}return a;}
R blow(unsigned n,bool reversed){R a(3*n);for(unsigned i=0;i<n;i++){for(unsigned u=0;u<3;u++)a[3*i+u]|=1u<<(3*i+(u+1)%3);for(unsigned j=0;j<n;j++)if(i!=j&&((((j+n-i)%n)<=n/2)^reversed))for(unsigned u=0;u<3;u++)a[3*i+u]|=7u<<(3*j);}return a;}
R local(R a){for(unsigned u=0;u<3;u++)a[u]^=7u^(1u<<u);return a;}
std::array<std::uint64_t,16> query(const R&a,std::array<unsigned,4> marks){std::array<std::uint64_t,16> out{};for(unsigned mask=0;mask<16;mask++){unsigned aa=0,cc=0;for(unsigned i=0;i<2;i++)if(mask&(1u<<i))aa|=1u<<marks[i];for(unsigned i=2;i<4;i++)if(mask&(1u<<i))cc|=1u<<marks[i];out[mask]=H(addports(a,aa,cc));}return out;}
std::int64_t coefficient(const std::array<std::uint64_t,16>&a){std::int64_t c=0;for(unsigned i=0;i<16;i++)c+=(__builtin_popcount(i)&1)?-std::int64_t(a[i]):std::int64_t(a[i]);return c;}
void array(const std::array<std::uint64_t,16>&a){std::cout<<"[";for(unsigned i=0;i<16;i++){if(i)std::cout<<",";std::cout<<a[i];}std::cout<<"]";}
int main(){std::cout<<"{\"status\":\"PASS\",\"records\":[";bool first=true;unsigned count=0;for(unsigned n:{3u,5u}){R T=blow(n,false),W=local(T),U=blow(n,true);for(unsigned same=0;same<2;same++){std::array<unsigned,4> marks=same?std::array<unsigned,4>{0,1,3,6}:std::array<unsigned,4>{0,3,1,6};auto t=query(T,marks),w=query(W,marks),u=query(U,marks);auto dl=coefficient(t)-coefficient(w),dr=coefficient(t)-coefficient(u);need(dl==dr,"balanced root restoration");need(dl==(same?0:(n==3?60:-109584)),"literal expected quartic");if(!first)std::cout<<",";first=false;std::cout<<"{\"quotient_order\":"<<n<<",\"same_port_double_mark\":"<<(same?"true":"false")<<",\"local_contrast\":"<<dl<<",\"root_contrast\":"<<dr<<",\"T_values\":";array(t);std::cout<<",\"local_values\":";array(w);std::cout<<",\"root_values\":";array(u);std::cout<<"}";count+=48;}}
std::cout<<"],\"physical_query_count\":"<<count<<"}\n";}
