#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;
vector<unsigned> graph(unsigned code){vector<unsigned>a(8);unsigned k=0;for(unsigned i=0;i<8;i++)for(unsigned j=i+1;j<8;j++,k++){if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
vector<unsigned> endbits(vector<unsigned>a){vector<unsigned>d(256);for(unsigned i=0;i<8;i++)d[1u<<i]=1u<<i;for(unsigned s=1;s<256;s++){unsigned bits=d[s];while(bits){unsigned b=bits&-bits;bits-=b;unsigned j=__builtin_ctz(b),next=a[j]&~s;while(next){unsigned t=next&-next;next-=t;d[s|t]^=t;}}}return d;}
int main(){ifstream file("blind8.bin",ios::binary);map<unsigned,unsigned>hist;unsigned count=0;char raw[4];while(file.read(raw,4)){unsigned code=0;for(unsigned i=0;i<4;i++)code|=unsigned((unsigned char)raw[i])<<(8*i);auto a=graph(code);auto e=endbits(a);for(unsigned i=0;i<8;i++)a[i]^=255^(1u<<i);auto s=endbits(a);unsigned h[8]={};for(unsigned mask=1;mask<256;mask++)if(!(__builtin_popcount(mask)&1)){for(unsigned j=0;j<8;j++)if(!(mask>>j&1))h[j]^=e[mask]^s[mask];}for(unsigned j=0;j<8;j++)if(h[j])throw runtime_error("false shared cubic blind");hist[e[255]^s[255]]++;count++;}cout<<"{\"verified_blind_quotients\":"<<count<<",\"endpoint_xor_histogram\":[";bool first=true;for(auto [b,n]:hist){if(!first)cout<<",";first=false;cout<<"["<<b<<","<<n<<"]";}cout<<"]}\n";}
