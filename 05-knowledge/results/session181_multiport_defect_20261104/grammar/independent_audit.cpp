#include <array>
#include <vector>
#include <iostream>
using A=std::vector<unsigned>;
struct Fields { long long h=0; std::array<long long,8> f{}; };
void words(const A&a,unsigned mask,int last,int bad,Fields& z){
 if(mask==((1u<<a.size())-1)){if(bad<0){++z.h;for(unsigned i=0;i<a.size();++i)++z.f[i];}else ++z.f[bad];return;}
 for(unsigned v=0;v<a.size();++v)if(!(mask>>v&1)){
  bool reverse=last>=0 && !(a[last]>>v&1);
  if(reverse && bad>=0)continue;
  words(a,mask|(1u<<v),v,reverse?last:bad,z);
 }
}
Fields calc(const A&a){Fields z;words(a,0,-1,-1,z);return z;}
A extend(A a,unsigned c){unsigned n=a.size();for(unsigned i=0;i<n;++i)if(c>>i&1)a[i]|=1u<<n;a.push_back(((1u<<n)-1)^c);return a;}
bool strong(const A&a){bool r[8][8]{};for(unsigned i=0;i<a.size();++i)for(unsigned j=0;j<a.size();++j)r[i][j]=(i==j)||(a[i]>>j&1);for(unsigned k=0;k<a.size();++k)for(unsigned i=0;i<a.size();++i)for(unsigned j=0;j<a.size();++j)r[i][j]=r[i][j]||(r[i][k]&&r[k][j]);for(unsigned i=0;i<a.size();++i)for(unsigned j=0;j<a.size();++j)if(!r[i][j])return false;return true;}
int main(){for(A a: {A{56,33,51,54,2,16},A{56,1,19,22,2,30}}){auto old=calc(a);std::cout<<"SOURCE "<<old.h<<" "<<old.f[0]<<"\n";int counts[2][2]{};int false0=0;for(unsigned x=0;x<64;++x)for(unsigned y=0;y<128;++y){A b=extend(extend(a,x),y);auto z=calc(b);bool s=strong(b);for(int d=2;d<=4;d+=2)if(z.h==3*old.h+d)for(int w=0;w<8;++w)if(z.f[w]==3*old.f[0]){if(s){++counts[1][d/2-1];if(w==0)++counts[0][d/2-1];std::cout<<"HIT "<<d<<" "<<x<<" "<<y<<" MARK "<<w<<"\n";}else if(w==0){++false0;std::cout<<"NONSTRONG "<<d<<" "<<x<<" "<<y<<" SINKS";for(int v=0;v<8;++v)if(!b[v])std::cout<<" "<<v;std::cout<<"\n";}}}std::cout<<"FIXED "<<counts[0][0]<<" "<<counts[0][1]<<" ALL "<<counts[1][0]<<" "<<counts[1][1]<<" NONSTRONG_FIXED "<<false0<<"\n";}}
