#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <queue>
#include <vector>
using namespace std;
int H(int mask){
 int r[6]={},bit=0;
 for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,bit++) if(mask>>bit&1)r[i]|=1<<j;else r[j]|=1<<i;
 int dp[64][6]={};for(int i=0;i<6;i++)dp[1<<i][i]=1;
 for(int s=1;s<64;s++)for(int i=0;i<6;i++)if(dp[s][i])for(int j=0;j<6;j++)if((r[i]>>j&1)&&!(s>>j&1)) dp[s|1<<j][j]+=dp[s][i];
 int x=0;for(int i=0;i<6;i++)x+=dp[63][i];return x;
}
int common(int x,int y){
 int r[6]={},bit=0;
 for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,bit++)if(((x^y)>>bit&1)==0){if(x>>bit&1)r[i]|=1<<j;else r[j]|=1<<i;}
 int dp[64][6]={};for(int i=0;i<6;i++)dp[1<<i][i]=1;
 for(int S=1;S<64;S++)for(int u=0;u<6;u++)if(dp[S][u])for(int v=0;v<6;v++)if((r[u]>>v&1)&&!(S>>v&1))dp[S|1<<v][v]+=dp[S][u];
 int z=0;for(int u=0;u<6;u++)z+=dp[63][u];return z;
}
int main(){
 constexpr int N=32768;vector<int> h(N);for(int i=0;i<N;i++)h[i]=H(i);
 vector<int> isolated,strict;int exact=0;map<int,int> hd;
 for(int x=0;x<N;x++){
  int e=abs(h[x]-17);if(!e){exact++;continue;}
  bool a=true,b=true;for(int i=0;i<15;i++){int f=abs(h[x^(1<<i)]-17);if(f<e)a=false;if(f<=e)b=false;}
  if(a)isolated.push_back(x);if(b)strict.push_back(x);
 }
 vector<int> dist(N,-1);queue<int>q;for(int x=0;x<N;x++)if(h[x]==17){dist[x]=0;q.push(x);}
 while(!q.empty()){int x=q.front();q.pop();for(int i=0;i<15;i++){int y=x^(1<<i);if(dist[y]<0){dist[y]=dist[x]+1;q.push(y);}}}
 map<int,int> mindist;for(int x:strict)mindist[dist[x]]++;
 int sublevel_isolated=0;map<int,int> barriers;
 // Minimax absolute-error level needed to connect source to any target,
 // via Dijkstra on all labelled sources and single-edge arcs.
 vector<int> barrier(N,1000);priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> pq;
 for(int x=0;x<N;x++)if(h[x]==17){barrier[x]=0;pq.push({0,x});}
 while(!pq.empty()){
  auto [c,x]=pq.top();pq.pop();if(c!=barrier[x])continue;
  for(int i=0;i<15;i++){int y=x^(1<<i),b=max(c,abs(h[y]-17));if(b<barrier[y]){barrier[y]=b;pq.push({b,y});}}
 }
 for(int x=0;x<N;x++)if(abs(h[x]-17)<=2 && h[x]!=17){bool iso=true;for(int i=0;i<15;i++)if(abs(h[x^(1<<i)]-17)<=2)iso=false;if(iso){sublevel_isolated++;barriers[barrier[x]]++;}}
 vector<int> lowdist(N,-1),next(N,-1);queue<int> lq;for(int x=0;x<N;x++)if(h[x]==17){lowdist[x]=0;lq.push(x);}
 while(!lq.empty()){int x=lq.front();lq.pop();for(int i=0;i<15;i++){int y=x^(1<<i);if(abs(h[y]-17)<=6&&lowdist[y]<0){lowdist[y]=lowdist[x]+1;next[y]=x;lq.push(y);}}}
 map<int,int> lowhist,guardhist,complementhist;vector<array<int,4>> reps;
 for(int x:strict){
  lowhist[lowdist[x]]++;int best=100,besty=-1,shared=0;
  for(int i=0;i<15;i++)for(int j=i+1;j<15;j++){
   int y=x^(1<<i)^(1<<j);if(h[y]==17&&common(x,y)>0){best=2;besty=y;shared=common(x,y);}
  }
  if(best>2)for(int i=0;i<15;i++)for(int j=i+1;j<15;j++)for(int k=j+1;k<15;k++){
   int y=x^(1<<i)^(1<<j)^(1<<k);if(h[y]==17&&common(x,y)>0){best=3;besty=y;shared=common(x,y);}
  }
  guardhist[best]++;
  if(reps.empty()|| (best==3&&reps.size()==1))reps.push_back({x,besty,shared,best});
  int e[6][6]={},b=0;for(int i=0;i<6;i++)for(int j=i+1;j<6;j++)e[i][j]=e[j][i]=b++;
  array<int,6> perm={0,1,2,3,4,5};int hits=0;
  do{bool ok=true;int used=0;for(int j=0;j<5;j++){int a=perm[j],c=perm[j+1],k=e[a][c];if(((x>>k&1)!=0)!=(a<c)){ok=false;break;}used|=1<<k;}
   if(ok&&h[x^(32767^used)]==17)hits++;
  }while(next_permutation(perm.begin(),perm.end()));complementhist[hits]++;
 }
 cout<<"{\n\"labelled_sources\":32768,\"target\":17,\"target_sources\":"<<exact<<",\"weak_local_minima\":"<<isolated.size()<<",\"strict_local_minima\":"<<strict.size()<<",\"sublevel_isolated\":"<<sublevel_isolated<<",\n\"strict_minimum_distance_histogram\":{";
 bool first=true;for(auto [d,c]:mindist){if(!first)cout<<',';first=false;cout<<'"'<<d<<"\":"<<c;}cout<<"},\n\"isolated_minimum_barrier_histogram\":{";
 first=true;for(auto[d,c]:barriers){if(!first)cout<<',';first=false;cout<<'"'<<d<<"\":"<<c;}cout<<"},\n\"mask72\":{\"H\":"<<h[72]<<",\"distance\":"<<dist[72]<<",\"barrier_absolute\":"<<barrier[72]<<",\"one_flip_counts\":[";
 for(int i=0;i<15;i++){if(i)cout<<',';cout<<h[72^(1<<i)];}cout<<"],\"two_flip_targets\":[";
 first=true;for(int i=0;i<15;i++)for(int j=i+1;j<15;j++)if(h[72^(1<<i)^(1<<j)]==17){if(!first)cout<<',';first=false;cout<<'['<<i<<','<<j<<']';}cout<<"]},\n\"local_minima_summary\":[";
 first=true;map<array<int,3>,int> summary;for(int x:strict)summary[{h[x],dist[x],barrier[x]}]++;
 for(auto [k,c]:summary){if(!first)cout<<',';first=false;cout<<"{\"H\":"<<k[0]<<",\"distance\":"<<k[1]<<",\"barrier\":"<<k[2]<<",\"sources\":"<<c<<'}';}cout<<"],\n\"low_barrier_distance_histogram\":{";first=true;for(auto[d,c]:lowhist){if(!first)cout<<',';first=false;cout<<'"'<<d<<"\":"<<c;}
 cout<<"},\n\"shared_path_packet_distance_histogram\":{";first=true;for(auto[d,c]:guardhist){if(!first)cout<<',';first=false;cout<<'"'<<d<<"\":"<<c;}
 cout<<"},\n\"legal_path_complement_target_hits_histogram\":{";first=true;for(auto[d,c]:complementhist){if(!first)cout<<',';first=false;cout<<'"'<<d<<"\":"<<c;}
 cout<<"},\n\"examples\":[";first=true;for(auto a:reps){if(!first)cout<<',';first=false;cout<<"{\"source\":"<<a[0]<<",\"target\":"<<a[1]<<",\"shared_paths\":"<<a[2]<<",\"packet_size\":"<<a[3]<<",\"minimum_barrier_walk\":[";int x=a[0];bool f=true;while(x>=0){if(!f)cout<<',';f=false;cout<<'['<<x<<','<<h[x]<<']';x=next[x];}cout<<"]}";}
 cout<<"]\n}\n";
}
