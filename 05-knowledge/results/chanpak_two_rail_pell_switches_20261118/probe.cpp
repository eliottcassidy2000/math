#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
using namespace std;
int main(){
 for(int n=3;n<=7;n++){
  vector<pair<int,int>> es;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)es.push_back({i,j});
  const int lim=1<<es.size(), full=(1<<n)-1;
  vector<int> hs(lim);
  for(int x=0;x<lim;x++){
   int in[7]{};
   for(unsigned i=0;i<es.size();i++){
    auto [a,b]=es[i];if((x>>i)&1)swap(a,b);
    in[b]|=1<<a;
   }
   int p[128][7]{};
   for(int s=1;s<=full;s++) for(int v=s;v;v&=v-1){
    int j=__builtin_ctz(unsigned(v)), t=s^(1<<j);
    if(!t)p[s][j]=1;
    else for(int u=t&in[j];u;u&=u-1)p[s][j]+=p[t][__builtin_ctz(unsigned(u))];
   }
   for(int v=0;v<n;v++)hs[x]+=p[full][v];
  }
  for(int value:{35,49,63}){
   auto it=find(hs.begin(),hs.end(),value);
   if(it!=hs.end())cout<<"fixture "<<n<<" "<<value<<" "<<it-hs.begin()<<'\n';
  }
  map<int,pair<int,int>> first;
  for(int x=0;x<lim;x++) for(unsigned e=0;e<es.size();e++){
   int y=x^(1<<e);
   if(hs[y]==7*hs[x])first.try_emplace(hs[x],x,y);
  }
  if(n==7){
   bool seen=false;
   for(int x=0;x<lim && !seen;x++) if(hs[x]==5)
    for(unsigned i=0;i<es.size() && !seen;i++)
     for(unsigned j=i+1;j<es.size() && !seen;j++){
      int y=x^(1<<i)^(1<<j);
      if(hs[y]==35) {
       cout<<"two flips 5 -> 35: "<<x<<" "<<y<<" "<<i<<" "<<j<<'\n';
       seen=true;
      }
     }
  }
  cout<<n<<":";
  for(auto [h,xy]:first)cout<<" "<<h<<"("<<xy.first<<","<<xy.second<<")";
  cout<<endl;
 }
}
