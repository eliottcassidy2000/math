// Second same-agent implementation: lexicographic subsets, intersection test,
// reversed owner order, and least-mask roots. No import of producer code.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=3)return 2;ifstream in(argv[1]);array<int,28>x{},y{};
 for(int i=0;i<28;i++)if(!(in>>y[i]>>x[i]))return 3;
 int budget=stoi(argv[2]);if(budget!=6&&budget!=8)return 4;
 for(int k=0;k<=14;k++){
  if(budget==8&&(k<4||k>5))continue;
  array<vector<uint32_t>,2048>bank;vector<int>sel(k);for(int j=0;j<k;j++)sel[j]=j;
  uint64_t count=0;
  while(true){
   int p=0,q=0;uint32_t m=0;for(int j:sel){p^=y[j];q+=x[j];m|=1u<<j;}bank[128*p+(q%128)].push_back(m);count++;
   int j=k-1;while(j>=0&&sel[j]==28-k+j)j--;if(j<0)break;
   sel[j]++;for(int h=j+1;h<k;h++)sel[h]=sel[h-1]+1;
  }
  cout<<"{\"k\":"<<k<<",\"states\":"<<count<<",\"rows\":[";bool first=true;
  for(int g=0;g<2048;g++){
   vector<uint32_t>unseen=std::move(bank[g]);if(unseen.empty())continue;
   size_t size=unseen.size();int cc=0;
   while(!unseen.empty()){
    auto it=min_element(unseen.begin(),unseen.end());uint32_t root=*it;*it=unseen.back();unseen.pop_back();cc++;
    vector<uint32_t>todo{root};
    for(size_t head=0;head<todo.size()&&!unseen.empty();head++){
     uint32_t a=todo[head];size_t j=unseen.size();
     while(j){j--;uint32_t b=unseen[j];if(popcount(a&b)>=k-budget/2){todo.push_back(b);unseen[j]=unseen.back();unseen.pop_back();}}
    }
   }
   if(!first)cout<<',';first=false;cout<<'['<<g<<','<<size<<','<<cc<<']';
  }
  cout<<"]}"<<endl;
 }
}
