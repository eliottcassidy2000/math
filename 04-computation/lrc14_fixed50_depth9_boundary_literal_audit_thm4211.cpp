// Literal joint-wall boundary control for THM-4211's q=5681/5682 edge.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

using i64=std::int64_t; using i128=__int128_t;
constexpr std::array<int,30> P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr std::array<int,9> R={40,63,80,120,126,145,176,190,252};
void req(bool x,const char*m){if(!x)throw std::runtime_error(m);} 
i64 xlcm(i64 a,i64 b){i128 x=static_cast<i128>(a/std::gcd(a,b))*b;req(x<=std::numeric_limits<i64>::max(),"overflow");return static_cast<i64>(x);} 
bool safe(int s,i64 l,i64 r,i64 d){i128 z=static_cast<i128>(s)*(l+r)%(static_cast<i128>(2)*d);return static_cast<i128>(7)*z>=d&&static_cast<i128>(7)*z<=static_cast<i128>(13)*d;}
std::string dec(i128 x){if(!x)return"0";bool n=x<0;if(n)x=-x;std::string s;while(x){s.push_back(char('0'+x%10));x/=10;}if(n)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
void run(int q){std::vector<int>S;for(int p:P)if(std::find(R.begin(),R.end(),p)==R.end())S.push_back(p);S.push_back(50);S.push_back(q);i64 d=1;for(int s:S)d=xlcm(d,14LL*s);std::vector<i64>w={0,d};for(int s:S){i64 u=d/(14LL*s);for(int k=0;k<s;k++){w.push_back((14LL*k+1)*u);w.push_back((14LL*k+13)*u);}}std::sort(w.begin(),w.end());w.erase(std::unique(w.begin(),w.end()),w.end());i64 m=0;for(std::size_t j=0;j+1<w.size();j++){bool ok=true;for(int s:S)if(!safe(s,w[j],w[j+1],d)){ok=false;break;}if(ok)m+=w[j+1]-w[j];}i128 delta=static_cast<i128>(63)*m-static_cast<i128>(4)*d;std::cout<<"Q "<<q<<" DEN "<<d<<" CELLS "<<w.size()-1<<" MASS "<<m<<" DELTA "<<dec(delta)<<'\n';}
int main(){run(5681);run(5682);}
