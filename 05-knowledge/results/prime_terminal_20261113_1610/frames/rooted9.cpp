// Complete one-vertex extension of the supplied ordinary order-eight bank.
// Test the first three exponent charts for the same two target rays.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>
using namespace std;
using Rows=vector<unsigned>;
long long hp(const Rows &a){
  int n=a.size(); if(n>11)throw runtime_error("DP order limit");
  long long dp[2048][11]={};
  for(int v=0;v<n;v++)dp[1<<v][v]=1;
  for(unsigned s=1;s<(1u<<n);s++)for(int v=0;v<n;v++)if(dp[s][v]){
    unsigned z=a[v]&~s;
    while(z){int w=__builtin_ctz(z);z&=z-1;dp[s|(1u<<w)][w]+=dp[s][v];}
  }
  long long h=0;for(int v=0;v<n;v++)h+=dp[(1<<n)-1][v];return h;
}
bool strong(const Rows&a){
  int n=a.size();for(int v=0;v<n;v++){
    unsigned s=1u<<v,t=0;
    while(t!=s){t=s;for(int u=0;u<n;u++)if(s>>u&1)s|=a[u];}
    if(s!=(1u<<n)-1)return false;
  }return true;
}
bool two_run(const Rows&a,int r){
  vector<pair<int,int>>e;
  for(int o=0;o<int(a.size());o++)if(a[r]>>o&1)
    for(int i=0;i<int(a.size());i++)if((a[i]>>r&1)&&(a[o]>>i&1))e.emplace_back(o,i);
  for(auto x:e)for(auto y:e)if(x.first!=y.first&&x.second!=y.second)return false;
  return true;
}
Rows clone(const Rows&a,int r){
  int n=a.size();Rows b=a;
  for(int u=0;u<n;u++)if(u==r || (a[u]>>r&1))b[u]|=1u<<n;
  b.push_back(a[r]);return b;
}
int main(){
  int count;cin>>count;if(count!=6880)throw runtime_error("ordinary order-eight cardinality");
  long long retained=0,extensions=0,targets=0,roots=0,hits=0;
  map<long long,long long> wanted={{95,11},{113,81},{115,81},{275,162},{277,162},{599,324},{601,324}};
  map<long long,long long> counts;
  vector<tuple<Rows,int,long long>>found;
  vector<Rows> audit;
  for(int k=0;k<count;k++){
    Rows a(8);for(auto &x:a)cin>>x;
    if(hp(a)>601)continue;retained++;
    for(unsigned mask=0;mask<256;mask++){
      Rows b=a;b.push_back(0);
      for(int v=0;v<8;v++)if(mask>>v&1)b[v]|=1u<<8;else b[8]|=1u<<v;
      extensions++;
      if(!strong(b))continue;
      auto h=hp(b);if(!wanted.count(h))continue;targets++;counts[h]++;
      if(h==95||h==113||h==115)audit.push_back(b);
      for(int r=0;r<9;r++)if(two_run(b,r)){
        roots++;auto h2=hp(clone(b,r));
        if(h2-h==2*wanted[h]){hits++;found.emplace_back(b,r,h);}
      }
    }
  }
  cout<<"{\"ordinary_parents\":"<<count<<",\"parents_H_le601\":"<<retained
      <<",\"extensions\":"<<extensions<<",\"strong_target_presentations\":"<<targets
      <<",\"two_run_root_presentations\":"<<roots<<",\"hit_presentations\":"<<hits<<",\"hits\":[";
  bool first=true;
  for(auto &[b,r,h]:found){if(!first)cout<<",";first=false;
    cout<<"{\"H\":"<<h<<",\"root\":"<<r<<",\"rows\":[";
    for(int i=0;i<9;i++){if(i)cout<<",";cout<<b[i];}cout<<"]}";
  }
  cout<<"],\"target_counts\":{";first=true;
  for(auto &[h,c]:counts){if(!first)cout<<",";first=false;cout<<"\""<<h<<"\":"<<c;}
  cout<<"},\"audit_targets\":[";first=true;
  for(auto &b:audit){if(!first)cout<<",";first=false;cout<<"[";
    for(int i=0;i<9;i++){if(i)cout<<",";cout<<b[i];}cout<<"]";}
  cout<<"]}\n";
  return 0;
}
