// Independent word-tree traversal: literal block profiles, no state deduplication.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
using I=int64_t;
struct State{I F,q,g,h;bool operator==(const State&b)const{return F==b.F&&q==b.q&&g==b.g&&h==b.h;}};
struct Block{I s=0,t=0,p=0,q=0,F=0;};
array<Block,16> blocks;
vector<unsigned char> support;
I ceiling;uint64_t nodes=0,loops=0;
bool path(const vector<int>&a,const vector<int>&r,int lo,int hi){for(int i=lo+1;i<hi;++i)if(!(r[a[i-1]]>>a[i]&1))return false;return true;}
void profiles(){
 vector<int> r={8,1,3,6};
 for(int mask=0;mask<16;++mask){
  auto&b=blocks[mask];int h=0,E=0;vector<int>a={0,1,2,3};
  do{
   if(path(a,r,0,4)){++h;b.s+=mask>>a[0]&1;b.t+=!(mask>>a[3]&1);}
   for(int k=0;k<=4;++k)if(path(a,r,0,k)&&path(a,r,k,4)){
    ++E;if(k)b.p+=!(mask>>a[k-1]&1);if(k<4)b.q+=mask>>a[k]&1;
   }
  }while(next_permutation(a.begin(),a.end()));
  if(h!=5||E!=32)throw runtime_error("block base counts");
  vector<int> augmented=r;augmented.push_back(mask);
  for(int i=0;i<4;++i)if(!(mask>>i&1))augmented[i]|=16;
  a={0,1,2,3,4};do{b.F+=path(a,augmented,0,5);}while(next_permutation(a.begin(),a.end()));
 }
}
void visit(State a){
 ++nodes;if(a.F>ceiling)return;
 for(I h=a.F;h<=ceiling;h*=5)support[h]=1;
 if(a.q>ceiling)return;
 vector<State> children={{a.F+a.q+a.g,2*a.q,0,a.h},{a.F,2*a.q+a.h,a.h,a.h}};
 for(auto b:blocks)children.push_back({5*a.F+b.p*a.q+(5-b.s)*a.g+(b.F-5)*a.h,
                                      32*a.q+b.q*a.h,(5-b.t)*a.h,5*a.h});
 for(auto b:children){
  if(b==a){++loops;continue;}
  if(b.F<a.F||b.q<a.q||b.h<1||b.F<b.h||b.g<0||b.g>b.h)throw runtime_error("state invariant");
  visit(b);
 }
}
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("usage: audit ceiling bitmap-output");
 ceiling=stoll(argv[1]);if(ceiling<1||ceiling>1000000)throw runtime_error("ceiling outside audited integer range");
 profiles();support.resize(ceiling+1);visit({1,0,0,1});visit({1,1,1,1});
 I count=0;for(I h=0;h<=ceiling;++h)if(support[h]){if(!(h&1)||h==7||h==21)throw runtime_error("odd spectrum hostile");++count;}
 ofstream out(argv[2],ios::binary);out.write(reinterpret_cast<const char*>(support.data()),support.size());
 if(!out)throw runtime_error("bitmap output failed");
 cout<<"{\"status\":\"passed\",\"limit\":"<<ceiling<<",\"word_tree_nodes\":"<<nodes
     <<",\"exact_self_loops_skipped\":"<<loops<<",\"distinct_counts\":"<<count<<"}\n";
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
