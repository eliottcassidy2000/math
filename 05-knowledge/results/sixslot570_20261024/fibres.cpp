// Exact complete binary-fibre census, no search cutoffs or random sampling.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>
using namespace std;
int main(int argc,char**argv){
 if(argc!=5){cerr<<"labels.txt first-k last-k budget\n";return 2;}
 ifstream in(argv[1]);array<int,28>x{},y{};for(int i=0;i<28;i++)if(!(in>>y[i]>>x[i]))return 3;
 int lo=stoi(argv[2]),hi=stoi(argv[3]),budget=stoi(argv[4]);
 if(lo<0||hi>14||lo>hi||budget<2||budget%2)return 4;
 for(int k=lo;k<=hi;k++){
  auto t=chrono::steady_clock::now();array<vector<uint32_t>,2048>bank;
  uint32_t m=(1u<<k)-1,limit=1u<<28;uint64_t total=0,comparisons=0;
  do{int a=0,b=0;for(uint32_t z=m;z;z&=z-1){int j=countr_zero(z);a^=y[j];b+=x[j];}
   bank[(a<<7)|(b&127)].push_back(m);total++;
   if(k==0)break;uint32_t c=m&-m,r=m+c;m=(((r^m)>>2)/c)|r;
  }while(m<limit);
  uint64_t fibres=0,bad=0,edges=0,checksum=0;int maxdepth=0;vector<array<uint64_t,4>>rows;
  for(int g=0;g<2048;g++){
   auto todo=std::move(bank[g]);if(todo.empty())continue;fibres++;uint64_t count=todo.size(),components=0;
   while(!todo.empty()){
    components++;vector<pair<uint32_t,int>>queue{{todo.back(),0}};todo.pop_back();
    for(size_t qi=0;qi<queue.size()&&!todo.empty();qi++){
     auto [a,depth]=queue[qi];
     size_t j=0;
     while(j<todo.size()){
      auto b=todo[j];comparisons++;
      if(popcount(a^b)<=budget){
       queue.emplace_back(b,depth+1);edges++;maxdepth=max(maxdepth,depth+1);
       checksum^=uint64_t(a)*0x9e3779b185ebca87ULL+uint64_t(b)*0xc2b2ae3d27d4eb4fULL;
       todo[j]=todo.back();todo.pop_back();
      }else j++;
     }
    }
   }
   if(components>1)bad++;rows.push_back({uint64_t(g),count,components,0});
  }
  double sec=chrono::duration<double>(chrono::steady_clock::now()-t).count();
  cout<<"{\"k\":"<<k<<",\"budget\":"<<budget<<",\"states\":"<<total<<",\"fibres\":"<<fibres<<",\"disconnected\":"<<bad<<",\"forest_edges\":"<<edges<<",\"max_root_depth\":"<<maxdepth<<",\"forest_checksum\":"<<checksum<<",\"comparisons\":"<<comparisons<<",\"seconds\":"<<sec<<",\"rows\":[";
  for(size_t i=0;i<rows.size();i++){if(i)cout<<',';auto r=rows[i];cout<<'['<<r[0]<<','<<r[1]<<','<<r[2]<<']';}
  cout<<"]}"<<endl;
 }
}
