#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=std::uint64_t;
void need(bool b){if(!b) throw std::runtime_error("certificate failure");}
int main(int argc,char**argv){
 int maxn=6;if(argc>1 && std::string(argv[1])=="--fail")need(false);
 I graphs=0,roots=0,pairs=0,saturated=0,equality=0;std::array<I,7> mins{};mins.fill(~I(0));
 for(int n=1;n<=maxn;++n){
  int full=(1<<n)-1;std::vector<std::pair<int,int>> ed;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)ed.push_back({i,j});
  for(I code=0;code<(I(1)<<ed.size());++code){
   std::array<int,6> out{};for(int t=0;t<int(ed.size());++t){auto [i,j]=ed[t];if(code>>t&1)out[i]|=1<<j;else out[j]|=1<<i;}
   I end[64][6]{},start[64][6]{},h[64]{};h[0]=1;
   for(int m=1;m<=full;++m){for(int v=0;v<n;++v)if(m>>v&1){int q=m^(1<<v);if(!q)end[m][v]=start[m][v]=1;else for(int w=0;w<n;++w)if(q>>w&1){if(out[w]>>v&1)end[m][v]+=end[q][w];if(out[v]>>w&1)start[m][v]+=start[q][w];}h[m]+=end[m][v];}need(h[m]%2==1);}
   ++graphs;
   for(int r=0;r<n;++r){
    I U=0,eO=0;int k=__builtin_popcount(unsigned(out[r]));int comp=full^(1<<r);
    for(int m=1;m<=full;++m)if(m>>r&1)U+=start[m][r]*h[full^m];
    for(int a=0;a<n;++a)if(out[r]>>a&1)eO+=end[comp][a];
    I cs=0;
    for(int a=0;a<n;++a)if(out[r]>>a&1)for(int b=a+1;b<n;++b)if(out[r]>>b&1){I K=0;for(int m=comp;;m=(m-1)&comp){K+=start[m][a]*end[comp^m][b]+start[m][b]*end[comp^m][a];if(!m)break;}need(K>=2 && K%2==0);cs+=K;++pairs;}
    need(U-h[full]==eO+cs);need(U-h[full]>=I(k*(k-1))+eO);++roots;
    if(U-h[full]==I(k*(k-1))+eO)++equality;
    bool strong=true;for(int s=0;s<n;++s){int seen=1<<s,old;do{old=seen;for(int x=0;x<n;++x)if(seen>>x&1)seen|=out[x];}while(old!=seen);if(seen!=full)strong=false;}
    if(strong && n>1){mins[k]=std::min(mins[k],U-h[full]);if(U==h[full]){need(k==1);need(end[full][r]==h[comp]);++saturated;}}
   }
  }
 }
 std::cout<<"{\n  \"max_order\": "<<maxn<<",\n  \"labelled_tournaments\": "<<graphs<<",\n  \"rooted_cases\": "<<roots<<",\n  \"outneighbor_pair_checks\": "<<pairs<<",\n  \"bound_equalities\": "<<equality<<",\n  \"nontrivial_strong_saturated_roots\": "<<saturated<<",\n  \"minimum_strong_defect_by_outdegree\": [";
 bool first=true;for(int k=1;k<6;++k)if(mins[k]!=~I(0)){if(!first)std::cout<<", ";first=false;std::cout<<"["<<k<<", "<<mins[k]<<"]";}std::cout<<"]\n}\n";
}
