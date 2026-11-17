#include <iostream>
#include <array>
#include <algorithm>
#include <vector>
#include <utility>
#include <cstdlib>

int main(){
 constexpr int n=7, full=(1<<n)-1;
 std::vector<std::pair<int,int>> edges;
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)edges.push_back({i,j});
 int witnesses=0, transfers=0;
 for(int word=0;word<1<<edges.size();++word){
  int out[n]={},D[full+1][n]={};
  for(unsigned b=0;b<edges.size();++b){
   auto [i,j]=edges[b];
   if(word>>b&1)out[i]|=1<<j;
   else out[j]|=1<<i;
  }
  for(int i=0;i<n;++i)D[1<<i][i]=1;
  for(int S=1;S<=full;++S)for(int i=0;i<n;++i)if(D[S][i]){
   int todo=out[i]&~S;while(todo){
    int bit=todo&-todo;todo^=bit;
    D[S|bit][__builtin_ctz((unsigned)bit)]+=D[S][i];
   }
  }
  int H=0;for(int i=0;i<n;++i)H+=D[full][i];
  if(H!=35)continue;
  ++witnesses;
  for(int S=1;S<=full;++S){
   H=0;for(int i=0;i<n;++i)H+=D[S][i];
   if(H==13)++transfers;
  }
 }
 std::cout<<"order "<<n<<", labelled H35 witnesses "<<witnesses
          <<", their induced H13 occurrences "<<transfers<<"\n";
 if(witnesses!=10080||transfers)return EXIT_FAILURE;
}
