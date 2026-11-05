// Independent census: degree classes only, no color refinement; all remaining
// permutations are explicit. Counts use predecessor-removal dynamic programs.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
using Rows=vector<unsigned>;

uint64_t canonical(const Rows&a){
 int n=a.size();map<int,vector<int>>groups;
 for(int i=0;i<n;i++)groups[__builtin_popcount(a[i])].push_back(i);
 vector<vector<int>>cells;for(auto [d,c]:groups)cells.push_back(c);
 vector<int>p(n);uint64_t best=UINT64_MAX;
 function<void(int,int)>rec=[&](int c,int k){
  if(c==(int)cells.size()){
   uint64_t key=0;
   for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[p[i]]&(1U<<p[j]))key|=1ULL<<(8*i+j);
   best=min(best,key);return;
  }
  auto q=cells[c];do{copy(q.begin(),q.end(),p.begin()+k);rec(c+1,k+q.size());}while(next_permutation(q.begin(),q.end()));
 };
 rec(0,0);return best;
}

Rows unpack(uint64_t key,int n){Rows a(n);for(int i=0;i<n;i++)a[i]=(key>>(8*i))&255;return a;}

pair<long long,long long>counts(const Rows&a){
 int n=a.size(),N=1<<n;array<array<long long,8>,256>dp{},cycle{};
 vector<unsigned>pred(n);
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[j]&(1<<i))pred[i]|=1<<j;
 for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++)if(mask&(1<<v)){
  int rest=mask^(1<<v);
  if(!rest){dp[mask][v]=1;if(v==0)cycle[mask][v]=1;continue;}
  unsigned prev=rest&pred[v];while(prev){int u=__builtin_ctz(prev);prev&=prev-1;dp[mask][v]+=dp[rest][u];cycle[mask][v]+=cycle[rest][u];}
 }
 long long H=0,C=0;for(int v=0;v<n;v++){H+=dp[N-1][v];if(n>=3&&(a[v]&1))C+=cycle[N-1][v];}
 return {H,C};
}

bool strong(const Rows&a){
 int n=a.size();for(int s=0;s<n;s++){unsigned seen=1<<s,old=0;while(old!=seen){old=seen;for(int v=0;v<n;v++)if(seen&(1<<v))seen|=a[v];}if(seen!=(1U<<n)-1)return false;}return true;
}

void jsonrows(const Rows&a){cout<<"[";for(size_t i=0;i<a.size();i++){if(i)cout<<",";cout<<a[i];}cout<<"]";}

int main(){
 set<uint64_t>classes={0};vector<int>expected={1,1,2,4,12,56,456,6880};
 cout<<"{\"inventory\":[";
 for(int n=1;n<=8;n++){
  if(n>1){set<uint64_t>next;for(auto key:classes){Rows a=unpack(key,n-1);for(int mask=0;mask<(1<<(n-1));mask++){Rows b=a;b.push_back(0);for(int v=0;v<n-1;v++){if(mask&(1<<v))b[v]|=1<<(n-1);else b[n-1]|=1<<v;}next.insert(canonical(b));}}classes.swap(next);}
  if((int)classes.size()!=expected[n-1])return 2;
  ofstream all("classes"+to_string(n)+".txt");
  if(n>1)cout<<",";cout<<"{\"order\":"<<n<<",\"classes\":"<<classes.size()<<",\"H67\":[";
  bool first=true;int stronghits=0;map<int,int>lowcycles;
  for(auto key:classes){Rows a=unpack(key,n);auto [H,C]=counts(a);bool st=strong(a);if((C>0)!=(st&&n>=3))return 3;
   all<<key<<" "<<H<<" "<<C<<" "<<st<<"\n";
   if(C==1||C==2)lowcycles[C]++;
   if(H==67){if(!first)cout<<",";first=false;cout<<"{\"rows\":";jsonrows(a);cout<<",\"H\":"<<H<<",\"C\":"<<C<<",\"strong\":"<<(st?"true":"false")<<"}";stronghits+=st;}
  }
  cout<<"],\"strong_H67\":"<<stronghits<<",\"C1_classes\":"<<lowcycles[1]<<",\"C2_classes\":"<<lowcycles[2]<<"}"<<flush;
  cerr<<"order "<<n<<": "<<classes.size()<<" classes, "<<stronghits<<" strong H67\n";
 }
 cout<<"]}\n";
}
