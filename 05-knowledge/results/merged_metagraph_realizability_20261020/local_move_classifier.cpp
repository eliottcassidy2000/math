// Standalone classifier certified by the independent Burnside class counts.
// No canonical codes, stored atlas, or enumeration of vertex permutations.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
struct DSU {
 vector<uint32_t> p, sz; uint32_t count;
 DSU(uint32_t n):p(n),sz(n,1),count(n){iota(p.begin(),p.end(),0);}
 uint32_t root(uint32_t a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;}
 bool join(uint32_t a,uint32_t b){a=root(a);b=root(b);if(a==b)return false;if(sz[a]<sz[b])swap(a,b);p[b]=a;sz[a]+=sz[b];--count;return true;}
};
struct Move {int from,to;uint32_t ones=0,zeros=0,constant=0;vector<pair<int,int>> bits;};
int main(int argc,char**argv){
 int top=argc>1?stoi(argv[1]):8;
 string output_dir=argc>2?argv[2]:".";
 if(top<3||top>8){cerr<<"audit range is 3..8"<<endl;return 2;}
 for(int n=3;n<=top;++n){
  int index[12][12];for(auto &row:index)for(int &entry:row)entry=-1;vector<pair<int,int>>free;
  for(int i=0;i<n;++i)for(int j=i+2;j<n;++j){index[i][j]=free.size();free.push_back({i,j});}
  int m=free.size();uint32_t N=1u<<m;vector<Move>moves;
  for(int a=0;a<n;++a)for(int b=a+2;b<n;++b){
   vector<int> p(n);iota(p.begin(),p.end(),0);int v=p[a];p.erase(p.begin()+a);p.insert(p.begin()+b,v);
   Move move;move.from=a;move.to=b;bool valid=true;
   for(int i=0;i<n-1;++i){int s=p[i],t=p[i+1];if(s+1==t)continue;if(t+1==s){valid=false;break;}if(s<t)move.ones|=1u<<index[s][t];else move.zeros|=1u<<index[t][s];}
   if(!valid)continue;
   for(int e=0;e<m;++e){auto[i,j]=free[e];int s=p[i],t=p[j];if(s+1==t)move.constant|=1u<<e;else if(t+1==s){}else{if(s>t){swap(s,t);move.constant|=1u<<e;}move.bits.push_back({index[s][t],e});}}
   moves.push_back(move);
  }
  DSU d(N);uint64_t valid_moves=0,unions=0;
  for(uint32_t x=0;x<N;++x)for(const auto&v:moves){
   if((x&v.ones)!=v.ones||(x&v.zeros))continue;
   ++valid_moves;uint32_t y=v.constant;for(auto[a,b]:v.bits)y^=((x>>a)&1u)<<b;
   if(d.join(x,y))++unions;
  }
  uint32_t ordinary=d.count;map<uint32_t,uint32_t> histord;
  for(uint32_t x=0;x<N;++x)if(d.root(x)==x)++histord[d.sz[x]];
  vector<int>ref;
  for(auto[i,j]:free)ref.push_back(index[n-1-j][n-1-i]);
  uint64_t runions=0;uint32_t sc=0,ordinary_even=0;
  for(uint32_t x=0;x<N;++x)if(d.root(x)==x){uint32_t y=0;for(int a=0;a<m;++a)y|=((x>>a)&1u)<<ref[a];if(d.root(y)==x)++sc;if(d.sz[x]%2==0)++ordinary_even;}
  for(uint32_t x=0;x<N;++x){uint32_t y=0;for(int a=0;a<m;++a)y|=((x>>a)&1u)<<ref[a];if(d.join(x,y))++runions;}
  uint32_t merged=d.count;map<uint32_t,uint32_t>hist;
  for(uint32_t x=0;x<N;++x)if(d.root(x)==x)++hist[d.sz[x]];
  cout<<"{\"n\":"<<n<<",\"tilings\":"<<N<<",\"move_templates\":"<<moves.size()<<",\"valid_directed_moves\":"<<valid_moves<<",\"path_move_unions\":"<<unions<<",\"ordinary_components\":"<<ordinary<<",\"ordinary_even_size_components\":"<<ordinary_even<<",\"self_converse_components\":"<<sc<<",\"reflection_unions\":"<<runions<<",\"merged_components\":"<<merged<<",\"merged_size_histogram\":{";
  bool first=true;for(auto[size,count]:hist){if(!first)cout<<",";first=false;cout<<"\""<<size<<"\":"<<count;}cout<<"}}"<<endl;
  if(n<=7){ofstream f(output_dir+"/local_assignment_n"+to_string(n)+".txt");if(!f){cerr<<"cannot write assignment file"<<endl;return 3;}map<uint32_t,uint32_t> ranks;for(uint32_t x=0;x<N;++x){auto r=d.root(x);if(!ranks.count(r))ranks[r]=ranks.size();f<<ranks[r]<<'\n';}}
 }
}
