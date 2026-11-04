#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
unsigned hp(unsigned n,unsigned code){std::array<unsigned,6>a{};unsigned k=0;for(unsigned u=0;u<n;u++)for(unsigned v=u+1;v<n;v++,k++)if(code&(1u<<k))a[u]|=1u<<v;else a[v]|=1u<<u;std::array<std::array<unsigned,6>,64>d{};for(unsigned u=0;u<n;u++)d[1u<<u][u]=1;for(unsigned mask=1;mask<(1u<<n);mask++)for(unsigned u=0;u<n;u++)if(d[mask][u]){unsigned next=a[u]&~mask;while(next){unsigned b=next&-next;next-=b;d[mask|b][__builtin_ctz(b)]+=d[mask][u];}}unsigned h=0;for(unsigned u=0;u<n;u++)h+=d[(1u<<n)-1][u];return h;}
int main(){std::cout<<"[";bool first=true;for(unsigned n:{4u,6u}){std::array<unsigned,721>hist{};unsigned max=0,firstcode=0,total=1u<<(n*(n-1)/2);for(unsigned c=0;c<total;c++){unsigned h=hp(n,c);hist[h]++;if(h>max){max=h;firstcode=c;}}if(max!=(n==4?5u:45u))throw std::runtime_error("small maximum");if(!first)std::cout<<",";first=false;std::cout<<"{\"n\":"<<n<<",\"labelled_tournaments\":"<<total<<",\"maximum\":"<<max<<",\"first_maximizer_code\":"<<firstcode<<",\"histogram\":[";bool f=true;for(unsigned h=0;h<hist.size();h++)if(hist[h]){if(!f)std::cout<<",";f=false;std::cout<<"["<<h<<","<<hist[h]<<"]";}std::cout<<"]}";}std::cout<<"]\n";}
