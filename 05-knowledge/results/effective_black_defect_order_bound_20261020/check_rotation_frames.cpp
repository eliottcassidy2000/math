#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
using namespace std;
int main(){
 vector<pair<int,int>> cases={{1,3},{1,5},{1,7},{1,9},{1,11},{1,13},{1,15},{2,3},{2,5},{2,7},{3,3},{4,3}};
 uint64_t totalframes=0,totalchecks=0;
 for(auto [q,m]:cases){
  int n=q*m, id[15][15],sgn[15][15]={};
  for(auto &row:id)fill(begin(row),end(row),-1);
  int vars=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(id[i][j]<0){
   for(int k=0;k<m;k++){
    int u=(i+k*q)%n,v=(j+k*q)%n;
    assert(id[u][v]<0 || (id[u][v]==vars && sgn[u][v]==1));
    id[u][v]=id[v][u]=vars;sgn[u][v]=1;sgn[v][u]=-1;
   }vars++;
  }
  vector<int> val(vars,-1),freevars;
  for(int i=0;i<n;i++){
   int j=(i+1)%n,k=id[i][j],v=sgn[i][j]>0;
   assert(val[k]<0 || val[k]==v);val[k]=v;
  }
  for(int i=0;i<vars;i++)if(val[i]<0)freevars.push_back(i);
  int lower=q==1?(m+2)/3:(m+1)/2,full=(1<<n)-1;
  vector<int> keeps;
  for(int d=0;d<=full;d++)if(__builtin_popcount((unsigned)d)<lower)keeps.push_back(full^d);
  uint64_t frames=1ULL<<freevars.size(),checks=0;
  for(uint64_t z=0;z<frames;z++){
   for(int i=0;i<(int)freevars.size();i++)val[freevars[i]]=(z>>i)&1;
   int out[15]={};
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j && val[id[i][j]]==(sgn[i][j]>0))out[i]|=1<<j;
   for(int i=0;i<n;i++){
    assert(out[i]>>(i+1)%n&1);
    int internal=0;for(int j=i%q;j<n;j+=q)internal+=out[i]>>j&1;
    assert(internal==(m-1)/2);
   }
   for(int keep:keeps){
    checks++;
    for(int rev=0;rev<2;rev++){
     int seen=keep&-keep;
     while(true){
      int next=seen;
      for(int u=0;u<n;u++)if(seen>>u&1)next|=(rev?(full^(out[u]|(1<<u))):out[u])&keep;
      if(next==seen)break;seen=next;
     }
     if(seen!=keep){cerr<<"FAIL q="<<q<<" m="<<m<<" z="<<z<<" keep="<<keep<<"\n";return 1;}
    }
   }
  }
  totalframes+=frames;totalchecks+=checks;
  cout<<"q="<<q<<" m="<<m<<" n="<<n<<" frames="<<frames<<" required_kappa="<<lower<<" deletion_checks="<<checks<<" PASS\n";
 }
 cout<<"TOTAL frames="<<totalframes<<" deletion_checks="<<totalchecks<<" PASS\n";
}
