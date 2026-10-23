// Enumerate cut sizes in increasing order, always containing the queried edge.
// A tight cut has c(G-S)=|S|+1, hence |S|<=floor((n-1)/2).
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>
using namespace std;
int main(int argc,char**argv){
 if(argc!=2)return 2;ofstream out(argv[1]);out<<"{\"queries\":[";bool first=true;
 vector<tuple<int,int,int,int>> queries={{15,1,3,7},{23,1,15,11},{23,3,6,11},{23,3,13,11},{23,4,5,11},{25,1,15,12},{26,1,3,12},{28,7,9,13},{29,2,7,14},{29,5,20,14},{29,7,9,14},{30,2,7,14},{30,5,20,14},{30,7,9,14},{31,7,9,15}};
 for(auto [n,a,b,limit]:queries){
  vector<uint64_t>adj(n);for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)for(int s=2;s*s<=2*n;s++)if(u+v+2==s*s){adj[u]|=1ULL<<v;adj[v]|=1ULL<<u;}
  uint64_t full=(1ULL<<n)-1,base=(1ULL<<(a-1))|(1ULL<<(b-1)),witness=0,checked=0;vector<int>available;for(int v=0;v<n;v++)if(!(base>>v&1))available.push_back(v);
  auto count=[&](uint64_t cut){uint64_t rem=full^cut;int components=0;while(rem){components++;uint64_t front=rem&-rem;rem^=front;while(front){int v=__builtin_ctzll(front);front&=front-1;uint64_t next=adj[v]&rem;rem^=next;front|=next;}}return components;};
  int tested=1;for(int k=2;k<=min(limit,(n-1)/2)&&!witness;k++){
   function<void(int,int,uint64_t)>choose=[&](int start,int left,uint64_t cut){if(witness)return;if(!left){checked++;if(count(cut)==k+1)witness=cut;return;}for(int i=start;i<=(int)available.size()-left&&!witness;i++)choose(i+1,left-1,cut|(1ULL<<available[i]));};
   choose(0,k-2,base);tested=k;
  }
  if(!first)out<<',';first=false;out<<"{\"n\":"<<n<<",\"edge\":["<<a<<','<<b<<"],\"searched_through_size\":"<<tested<<",\"all_possible_sizes\":"<<(tested==(n-1)/2?"true":"false")<<",\"cuts_tested\":"<<checked<<",\"minimal_tight_cut\":[";bool sep=false;for(int v=0;v<n;v++)if(witness>>v&1){if(sep)out<<',';sep=true;out<<v+1;}out<<"]}";
  cerr<<n<<' '<<a<<'-'<<b<<" size "<<__builtin_popcountll(witness)<<" tested "<<checked<<'\n';
 }out<<"]}\n";
}
