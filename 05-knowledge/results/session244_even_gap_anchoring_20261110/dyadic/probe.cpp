#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using I=std::uint64_t;
int main(){
std::set<I> both,plus,minus;std::map<I,std::vector<int>> examples;
for(int n=2;n<=6;++n){int full=(1<<n)-1;std::vector<std::pair<int,int>>ed;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)ed.emplace_back(i,j);
for(I code=0;code<(I(1)<<ed.size());++code){std::array<int,6>out{};for(int k=0;k<int(ed.size());k++){auto[i,j]=ed[k];if(code>>k&1)out[i]|=1<<j;else out[j]|=1<<i;}
I end[64][6]{},start[64][6]{},h[64]{};h[0]=1;
for(int m=1;m<=full;m++)for(int v=0;v<n;v++)if(m>>v&1){int q=m^(1<<v);if(!q)end[m][v]=start[m][v]=1;else for(int w=0;w<n;w++)if(q>>w&1){if(out[w]>>v&1)end[m][v]+=end[q][w];if(out[v]>>w&1)start[m][v]+=start[q][w];}h[m]+=end[m][v];}
bool p=false,m=false;for(int r=0;r<n;r++){I u=0;for(int s=1;s<=full;s++)if(s>>r&1)u+=start[s][r]*h[full^s];I d=h[full^(1<<r)];if(end[full][r]!=d)continue;if(u==h[full]+d+1)p=true;if(u+1==h[full]+d)m=true;}
if(p){plus.insert(h[full]);}
if(m){minus.insert(h[full]);}
if(p&&m){both.insert(h[full]);if(!examples.count(h[full]))examples[h[full]]=std::vector<int>(out.begin(),out.begin()+n);}
}}
auto print=[](const char*name,const std::set<I>&s){std::cout<<name<<":";for(auto h:s)std::cout<<" "<<h;std::cout<<"\n";};print("plus",plus);print("minus",minus);print("both",both);for(const auto&[h,a]:examples){std::cout<<h<<":";for(int x:a)std::cout<<" "<<x;std::cout<<"\n";}
}
