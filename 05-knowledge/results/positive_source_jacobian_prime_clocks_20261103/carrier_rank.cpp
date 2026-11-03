// Exact signed block contraction; modular row bases give rational lower bounds.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using I=std::int64_t;
struct Block{int a,b,sg;};
std::vector<std::vector<Block>> types;
I coefficient(const std::vector<Block>& old,const std::vector<int>& word){
 std::vector<unsigned> adj(6+word.size(),0);int n=6;
 for(int op:word){for(int i=0;i<n;i++)if(op<6&&op==i)adj[n]|=1u<<i;else adj[i]|=1u<<n;n++;}
 std::vector<Block>b=old;for(int i=6;i<n;i++)b.push_back({i,i,1});
 int k=int(b.size()),states=2*k;
 std::vector<int> starts(states),ends(states),sg(states),valid(states);
 for(int i=0;i<k;i++){
  starts[2*i]=b[i].a;ends[2*i]=b[i].b;sg[2*i]=1;valid[2*i]=1;
  starts[2*i+1]=b[i].b;ends[2*i+1]=b[i].a;sg[2*i+1]=b[i].sg;valid[2*i+1]=b[i].a!=b[i].b;
 }
 std::vector<int> trans(states*states);
 for(int a=0;a<states;a++)for(int c=0;c<states;c++){
  int u=ends[a],v=starts[c];
  trans[a*states+c]=(u<6&&v<6)?1:2*int((adj[u]>>v)&1);
 }
 std::vector<I>dp(std::size_t(1u<<k)*states,0);
 for(int a=0;a<states;a++)if(valid[a])dp[(1u<<(a/2))*states+a]=sg[a];
 for(unsigned mask=1;mask<(1u<<k);mask++)for(int a=0;a<states;a++){
  I v=dp[mask*states+a];if(!v)continue;
  unsigned free=((1u<<k)-1)^mask;
  while(free){int j=__builtin_ctz(free);free&=free-1;
   for(int o=0;o<2;o++){int c=2*j+o;if(valid[c])
    dp[(mask|(1u<<j))*states+c]+=v*trans[a*states+c]*sg[c];}
  }
 }
 I total=0;for(int a=0;a<states;a++)total+=dp[((1u<<k)-1)*states+a];return total;
}
struct Basis{
 int p,rank=0;std::vector<std::vector<int>> row;
 explicit Basis(int q,int size):p(q),row(size){}
 bool insert(const std::vector<I>&v){
  std::vector<int>a(v.size());for(std::size_t i=0;i<v.size();i++)a[i]=int((v[i]%p+p)%p);
  for(std::size_t j=0;j<a.size();j++)if(a[j]){
   if(row[j].empty()){
    int inv=1;while(inv*a[j]%p!=1)inv++;
    for(std::size_t z=j;z<a.size();z++)a[z]=a[z]*inv%p;
    row[j]=std::move(a);rank++;return true;
   }
   int c=a[j];for(std::size_t z=j;z<a.size();z++)a[z]=(a[z]-c*row[j][z]%p+p)%p;
  }return false;
 }
};
int main(){
 int d;if(!(std::cin>>d)||d!=571)return 2;types.resize(d);
 for(auto&t:types){int k;std::cin>>k;t.resize(k);for(auto &b:t)std::cin>>b.a>>b.b>>b.sg;}
 Basis a(101,d),b(103,d);std::vector<std::vector<int>>selected;std::vector<int>layers;
 int examined=0;bool full=false;
 for(int depth=0;depth<=5&&!full;depth++){
  int choices=1;for(int i=0;i<depth;i++)choices*=7;
  int cap=depth<5?choices:4000;
  for(int ix=0;ix<cap;ix++){
   int code=(ix*137)%choices;std::vector<int>w(depth);
   for(int &z:w){z=code%7;code/=7;}
   std::vector<I>v;for(const auto&t:types)v.push_back(coefficient(t,w));examined++;
   bool fresh=a.insert(v);b.insert(v);if(fresh)selected.push_back(w);
   if(a.rank==d&&b.rank==d){full=true;break;}
  }
  layers.push_back(a.rank);std::cerr<<"depth "<<depth<<" rank "<<a.rank<<"/"<<b.rank<<" queries "<<examined<<"\n";
 }
 std::cout<<"{\"rank_mod101\":"<<a.rank<<",\"rank_mod103\":"<<b.rank<<",\"queries_examined\":"<<examined<<",\"ranks_by_depth\":[";
 for(std::size_t i=0;i<layers.size();i++){if(i)std::cout<<",";std::cout<<layers[i];}
 std::cout<<"],\"independent_words_mod101\":[";
 for(std::size_t i=0;i<selected.size();i++){
  if(i)std::cout<<",";std::cout<<"[";for(std::size_t j=0;j<selected[i].size();j++){if(j)std::cout<<",";std::cout<<selected[i][j];}std::cout<<"]";
 }std::cout<<"],\"operations\":\"0..5=B_old_port; 6=sink\"}\n";
}
