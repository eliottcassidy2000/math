#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
#include <filesystem>
#include <stdexcept>
using namespace std;
using St=array<int,4>; // h,F,q,t
constexpr int LIMIT=1000000;
vector<bool> support(LIMIT+1); uint64_t visits=0;
// No global state deduplication: enumerate the word tree. Only an exact
// self-loop is deleted. Emit F even at q>LIMIT, then stop there.
void walk(St a){
 ++visits; auto [h,f,q,t]=a;
 if(f>LIMIT)return;
 support[f]=true;
 if(q>LIMIT)return;
 St zero={h,f+q+h-t,2*q,h},one={h,f,2*q+h,0};
 if(zero!=a)walk(zero);
 if(one!=a)walk(one);
}
int main(int argc,char**argv){
 if(argc!=2){cerr<<"Usage: independent OUTPUT_DIRECTORY\n";return 2;}
 const filesystem::path dst(argv[1]);filesystem::create_directories(dst);
 set<St> seeds; ofstream sf(dst/"seeds.tsv");
 if(!sf)throw runtime_error("Cannot open seed output");
 for(int n=1;n<=5;++n)for(int code=0;code<(1<<(n*(n-1)/2));++code){
  bool a[6][6]={};int k=0;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k){a[i][j]=(code>>k)&1;a[j][i]=!a[i][j];}
  for(int cut=0;cut<(1<<n);++cut){
   int h=0,f=0,q=0,t=0; vector<int> p(n);iota(p.begin(),p.end(),0);
   do{
    bool path=true;for(int i=1;i<n;++i)path&=a[p[i-1]][p[i]];
    if(path){++h;if(!(cut>>p.back()&1))++t;}
    // Each permutation and split is a unique ordered path cover.
    for(int split=0;split<=n;++split){
     bool valid=true;for(int i=1;i<n;++i)if(i!=split)valid&=a[p[i-1]][p[i]];
     if(!valid)continue;
     bool left=split==0||!(cut>>p[split-1]&1);
     bool right=split==n||(cut>>p[split]&1);
     if(split<n&&right)++q;
     if(left&&right)++f;
    }
   }while(next_permutation(p.begin(),p.end()));
   seeds.insert({h,f,q,t});
  }
 }
 for(auto s:seeds){for(size_t i=0;i<s.size();++i){if(i)sf<<' ';sf<<s[i];}sf<<'\n';walk(s);}
 ofstream out(dst/"support.bin",ios::binary);
 if(!out)throw runtime_error("Cannot open support output");
 int count=0;for(int i=0;i<=LIMIT;++i){out.put(support[i]?1:0);count+=support[i];if(support[i]&&!(i%2))throw runtime_error("even support");}
 cout<<"independent literal permutation seeds "<<seeds.size()<<"\nword-tree visits "<<visits<<"\nsupport "<<count<<"\nallowed misses";
 for(int i=1;i<=LIMIT;i+=2)if(i!=7&&i!=21&&!support[i])cout<<' '<<i;
 cout<<'\n';if(support[7]||support[21])throw runtime_error("forbidden found");
}
