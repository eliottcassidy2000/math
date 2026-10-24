// Independent audit: numeric Gosper enumeration, 14+14 direct charge tables,
// and deletion-list BFS. No producer grouping or union-find code is used.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
struct Half {array<int,16384>x{},y{},w{};};
int main(int argc,char**argv){
 if(argc!=5){cerr<<"input first last budget\n";return 1;}
 array<int,28>x,y,w;ifstream in(argv[1]);int s,xx;
 for(int i=0;i<28;i++)if(!(in>>s>>xx>>x[i]>>y[i]>>w[i]))return 2;
 Half a,b;
 for(int m=1;m<16384;m++){int t=__builtin_ctz((unsigned)m),r=m&(m-1);a.x[m]=a.x[r]+x[t];a.y[m]=a.y[r]^y[t];a.w[m]=a.w[r]+w[t];b.x[m]=b.x[r]+x[t+14];b.y[m]=b.y[r]^y[t+14];b.w[m]=b.w[r]+w[t+14];}
 int budget=atoi(argv[4]);
 for(int k=atoi(argv[2]);k<=atoi(argv[3]);k++){
  array<vector<uint32_t>,2048> g;uint64_t total=0,checks=0;int nf=0,bad=0,depth_max=0;map<int,int> hist;
  uint32_t m=(1u<<k)-1;
  do{int lo=m&16383,hi=m>>14;g[(((a.x[lo]+b.x[hi])&127)<<4)|(a.y[lo]^b.y[hi])].push_back(m);total++;if(!k)break;uint32_t u=m&-m,v=u+m;m=v+(((v^m)/u)>>2);}while(m<(1u<<28));
  for(auto&v:g){if(v.empty())continue;nf++;int components=0;
   while(!v.empty()){
    components++;vector<uint32_t> todo;vector<int> depth;todo.push_back(v.back());depth.push_back(0);v.pop_back();
    for(size_t at=0;at<todo.size();at++){
     uint32_t current=todo[at];size_t j=0;
     while(j<v.size()){
      checks++;uint32_t d=current^v[j];
      if(__builtin_popcount(d)<=8 && a.w[d&16383]+b.w[d>>14]<=budget){todo.push_back(v[j]);depth.push_back(depth[at]+1);depth_max=max(depth_max,depth[at]+1);v[j]=v.back();v.pop_back();}else j++;
     }
    }
   }
   hist[components]++;if(components>1)bad++;
  }
  cout<<"{\"margin\":"<<k<<",\"budget\":"<<budget<<",\"support\":8,\"states\":"<<total<<",\"fibres\":"<<nf<<",\"disconnected_fibres\":"<<bad<<",\"edge_tests\":"<<checks<<",\"maximum_root_depth\":"<<depth_max<<",\"component_counts\":{";bool first=true;for(auto [n,c]:hist){if(!first)cout<<",";first=false;cout<<"\""<<n<<"\":"<<c;}cout<<"}}"<<endl;
 }
}
