#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using I=long long;
using Profile=array<I,7>; // H,E,p,q,s,t,F
struct Witness {vector<int> rows;int mask;};
map<Profile,Witness> profiles;
bool strong(const vector<int>& rows){
 int n=rows.size(),full=(1<<n)-1;
 for(int v=0;v<n;v++){
  int seen=1<<v,old=0;
  while(old!=seen){old=seen;for(int u=0;u<n;u++)if(old>>u&1)seen|=rows[u];}
  if(seen!=full)return false;
 }
 return true;
}
int main(int argc,char**argv){
 int maxn=argc>1?stoi(argv[1]):7,cap=argc>2?stoi(argv[2]):31;
 bool fixed=argc<4 || string(argv[3])!="labelled";
 cout<<"{\"cap_H\":"<<cap<<",\"max_order\":"<<maxn<<",\"gauge\":\""<<(fixed?"fixed_path":"all_labelled")<<"\",\"levels\":[";
 for(int n=1;n<=maxn;n++){
  uint64_t count=0,small=0;I minimum=INT64_MAX;
  int edges=n*(n-1)/2-(fixed?n-1:0),full=(1<<n)-1;
  for(uint64_t bits=0;bits<(1ULL<<edges);bits++){
   vector<int> rows(n);int pos=0;
   for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
    if(fixed && v==u+1){rows[u]|=1<<v;continue;}
    if(bits>>pos&1)rows[u]|=1<<v;else rows[v]|=1<<u;pos++;
   }
   if(!strong(rows))continue;count++;
   vector<vector<I>> end(1<<n,vector<I>(n)),start=end;
   vector<I> H(1<<n);H[0]=1;
   for(int v=0;v<n;v++)end[1<<v][v]=start[1<<v][v]=1;
   for(int m=1;m<=full;m++)for(int v=0;v<n;v++)if(m>>v&1){
    int a=m^(1<<v);
    if(a)for(int u=0;u<n;u++)if(a>>u&1){
     if(rows[u]>>v&1)end[m][v]+=end[a][u];
     if(rows[v]>>u&1)start[m][v]+=start[a][u];
    }
    H[m]+=end[m][v];
   }
   minimum=min(minimum,H[full]);if(H[full]>cap)continue;small++;
   I E=0;vector<I> pv(n),qv(n);vector<vector<I>> W(n,vector<I>(n));
   for(int m=0;m<=full;m++){
    E+=H[m]*H[full^m];
    for(int v=0;v<n;v++)if(m>>v&1){pv[v]+=end[m][v]*H[full^m];qv[v]+=start[m][v]*H[full^m];}
    if(m && m!=full)for(int u=0;u<n;u++)if(m>>u&1)for(int v=0;v<n;v++)if(!(m>>v&1))W[u][v]+=end[m][u]*start[full^m][v];
   }
   for(int mask=0;mask<=full;mask++){
    I p=0,q=0,s=0,t=0,C=0;
    for(int v=0;v<n;v++){
     if(mask>>v&1){q+=qv[v];s+=start[full][v];}
     else {p+=pv[v];t+=end[full][v];}
     if(mask>>v&1)for(int u=0;u<n;u++)if(!(mask>>u&1) && (rows[v]>>u&1)){
      if((W[u][v]+W[v][u])%2)throw 1;
      C+=(W[u][v]+W[v][u])/2;
     }
    }
    Profile b={H[full],E,p,q,s,t,H[full]+2*C};
    if(!profiles.count(b))profiles[b]={rows,mask};
   }
  }
  if(n>1)cout<<",";
  cout<<"{\"n\":"<<n<<",\"strong_presentations\":"<<count<<",\"small_presentations\":"<<small<<",\"minimum_H\":"<<(count?minimum:0)<<"}";
  cerr<<"n="<<n<<" small="<<small<<" profiles="<<profiles.size()<<"\n";
 }
 cout<<"],\"profiles\":[";bool first=true;
 for(auto [b,w]:profiles){if(!first)cout<<",";first=false;
  cout<<"{\"profile\":[";for(int i=0;i<7;i++){if(i)cout<<",";cout<<b[i];}
  cout<<"],\"mask\":"<<w.mask<<",\"rows\":[";for(int i=0;i<w.rows.size();i++){if(i)cout<<",";cout<<w.rows[i];}cout<<"]}";
 }
 cout<<"]}\n";
 return 0;
}
