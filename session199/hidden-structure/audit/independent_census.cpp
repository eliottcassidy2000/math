#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using U=uint32_t; using K=uint64_t;
struct V { U a,f,e,p; };
constexpr K MASK=(1ULL<<20)-1;
K encode(V v) {
 U tag,x,y,z;
 if(v.e==0){tag=0;x=v.a;y=v.f;z=v.p;}
 else if(v.p==0){tag=1;x=v.a;y=v.f;z=v.e;}
 else {if(v.f)throw std::runtime_error("nonphysical face");tag=2;x=v.a;y=v.e;z=v.p;}
 if(x>MASK||y>MASK||z>MASK)throw std::runtime_error("encoding overflow");
 return (K(tag)<<60)|(K(x)<<40)|(K(y)<<20)|z;
}
V decode(K k){U t=k>>60,x=(k>>40)&MASK,y=(k>>20)&MASK,z=k&MASK;
 if(t==0)return{x,y,0,z};if(t==1)return{x,y,z,0};if(t==2)return{x,0,y,z};throw std::runtime_error("tag");}
V O(V v){return{v.a+v.e,v.a+2*v.f+v.e+2*v.p,0,v.e+v.p};}
V Z(V v){return{2*v.a+v.f+v.e+v.p,v.f,v.e+v.p,0};}
V J(V v){return{v.f+v.p,0,v.a+v.e+v.p,v.e};}
U G(V v){return 3*v.a+2*v.f+6*v.e+5*v.p;}
int main(int argc,char**argv){
 U cap=argc>1?std::stoul(argv[1]):1000001;
 if(cap>1000001||cap<5)throw std::runtime_error("declared cap");
 std::unordered_set<K> seen;seen.max_load_factor(0.8);seen.reserve(cap==1000001?21000000:600000);
 std::vector<K> previous,next,stack;
 std::vector<uint8_t> scalar(cap/2+1,255);scalar[1]=0;
 std::array<size_t,4> counts{};
 for(U cost=0;cost<=3;++cost){
  auto close=[&](V v){
   if(G(v)>cap)return;
   stack.push_back(encode(v));
   while(!stack.empty()){
    K key=stack.back();stack.pop_back();V s=decode(key);U g=G(s);
    if(g>cap)continue;
    if(!seen.insert(key).second)continue;
    if(encode(s)!=key)throw std::runtime_error("encoding roundtrip");
    next.push_back(key);++counts[cost];
    if(!(g&1))throw std::runtime_error("nonodd target");
    if(scalar[g/2]>cost)scalar[g/2]=cost;
    V a=O(s),b=Z(s);
    if(G(a)<=g||G(b)<=g)throw std::runtime_error("nonincreasing transfer");
    if(G(a)<=cap)stack.push_back(encode(a));
    if(G(b)<=cap)stack.push_back(encode(b));
   }
  };
  if(cost==0)close({0,0,0,1});
  else for(K key:previous){V s=decode(key),v=J(s);if(G(v)<=G(s))throw std::runtime_error("J growth");close(v);}
  std::cerr<<"cost "<<cost<<": "<<counts[cost]<<" new, "<<seen.size()<<" total\n";
  previous.swap(next);next.clear();
 }
 std::array<size_t,4> hist{};std::vector<U> gaps,two;
 for(U n=3;n<=cap;n+=2){auto c=scalar[n/2];if(c==255)gaps.push_back(n);else ++hist[c];if(c>2)two.push_back(n);}
 auto arr=[](const auto&v){std::cout<<'[';bool comma=false;for(auto x:v){if(comma)std::cout<<',';comma=true;std::cout<<x;}std::cout<<']';};
 std::cout<<"{\"cap\":"<<cap<<",\"source_count\":"<<seen.size()<<",\"source_histogram\":";arr(counts);
 std::cout<<",\"scalar_histogram\":";arr(hist);std::cout<<",\"gaps\":";arr(gaps);std::cout<<",\"two_J_gaps\":";arr(two);std::cout<<"}\n";
}
