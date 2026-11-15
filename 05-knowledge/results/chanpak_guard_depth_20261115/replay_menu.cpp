#include <array>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
using namespace std;
int main(){
 for(int n=3;n<=6;n++){
  set<pair<int,int>> seen;
  vector<pair<int,int>> edges;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.emplace_back(i,j);
  int full=(1<<n)-1;
  for(int k=0;k<1<<(int)edges.size();k++){
   array<int,6> row{};
   for(int j=0;j<(int)edges.size();j++){
    auto [a,b]=edges[j]; if((k>>j)&1)swap(a,b); row[a]|=1<<b;
   }
   // Exclude all nontrivial source cuts, exactly the strong tournaments.
   bool strong=true;
   for(int m=1;m<full && strong;m++){
    int out=0;
    for(int i=0;i<n;i++)if((m>>i)&1)out|=row[i];
    if((out&(~m))==0)strong=false;
   }
   if(!strong)continue;
   array<array<int,6>,64> D{};
   array<int,64> H{}; H[0]=1;
   for(int m=1;m<=full;m++){
    for(int i=0;i<n;i++)if((m>>i)&1){
     int u=m^(1<<i);
     if(!u)D[m][i]=1;
     else for(int j=0;j<n;j++)if(((u>>j)&1)&&((row[j]>>i)&1))D[m][i]+=D[u][j];
     H[m]+=D[m][i];
    }
   }
   int q=0;for(int m=0;m<=full;m++)q+=H[m]*H[full^m];
   seen.emplace(H[full],q);
  }
  int want[]={1,1,5,16}; if((int)seen.size()!=want[n-3])return EXIT_FAILURE;
  cout<<n<<" "<<seen.size()<<" profiles:";
  for(auto [h,q]:seen)cout<<" "<<h<<":"<<q;
  cout<<"\n";
 }
}
