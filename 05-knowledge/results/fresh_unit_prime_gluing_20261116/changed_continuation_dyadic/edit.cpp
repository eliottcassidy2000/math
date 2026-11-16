// Minimize the old arc reversals when embedding the marked 493 host.
// Old root 7 is matched to the last vertex of each new rooted witness.
#include <array>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
using Rows=vector<int>;
int main(){
 Rows src{118,124,504,369,224,192,128,11,243};
 vector<pair<int,Rows>> bank{
   {1479,{510,1004,1016,496,482,448,896,256,512,185}},
   {1481,{2014,4092,4056,4048,4064,4045,1920,1792,3584,1024,2048,705}},
   {1483,{2042,4092,4089,4016,4064,1984,1928,1792,3584,1024,2048,737}}
 };
 constexpr array<int,9> order{7,2,3,0,1,8,4,5,6};
 for(auto [q,dst]:bank){
  const int n=dst.size(), F=(1<<n)-1;
  int best=37;
  array<int,9> map{},witness{};
  auto dfs=[&](auto&& self,int t,int used,int v)->void{
   if(v>=best) return;
   if(t==9){best=v;witness=map;return;}
   int old=order[t];
   int rest=t?(F^used):(1<<(n-1));
   for(int b=rest;b;b&=b-1){
    int id=__builtin_ctz(unsigned(b));
    int k=v;
    for(int i=0;i<t;i++){int pre=order[i];k+=((src[old]>>pre&1)!=(dst[id]>>map[pre]&1));}
    map[old]=id;self(self,t+1,used|(1<<id),k);
   }
  };
  dfs(dfs,0,0,0);
  if(best!=4)throw runtime_error("unexpected edit minimum");
  cout<<q<<" edits "<<best<<" map";
  for(int x:witness)cout<<" "<<x;
  cout<<" arcs";
  for(int i=0;i<9;i++)for(int j=i+1;j<9;j++)if((src[i]>>j&1)!=(dst[witness[i]]>>witness[j]&1))cout<<" "<<i<<"-"<<j;
  cout<<"\n";
 }
}
