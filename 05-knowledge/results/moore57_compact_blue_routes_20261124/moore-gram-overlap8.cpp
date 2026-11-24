#include <array>
#include <cstdint>
#include <iostream>
int main(){
  // Every tournament has a Hamilton path. Choose its order as labels,
  // leaving 28-7=21 unfixed arcs. No isomorphism pruning is needed.
  constexpr unsigned n=8;
  unsigned low=0, tri_free=0, square_free=0, both_free=0;
  for(uint32_t m=0;m<(1u<<21);m++){
    std::array<unsigned,n> r{}, f{};
    uint32_t c=m;
    for(unsigned i=0;i<n;i++) for(unsigned j=i+1;j<n;j++){
      if(j==i+1) r[i]|=1u<<j;
      else {
        if(c&1)r[j]|=1u<<i;else r[i]|=1u<<j;
        c>>=1;
      }
    }
    bool valid=true;
    for(unsigned i=0;i<n;i++) for(unsigned j=i+1;j<n;j++){
      int x=__builtin_popcount(r[i]&r[j]);
      if(x>2)valid=false;
      if(x==2){f[i]|=1u<<j;f[j]|=1u<<i;}
    }
    if(!valid)continue;
    low++;
    bool triangles=true,squares=true;
    for(unsigned i=0;i<n;i++) for(unsigned j=i+1;j<n;j++) {
      int x=__builtin_popcount(f[i]&f[j]);
      if(x>1)squares=false;
      if(x&&((f[i]>>j)&1))triangles=false;
    }
    tri_free+=triangles;
    square_free+=squares;
    both_free+=triangles&&squares;
  }
  std::cout<<"bounded codegrees "<<low<<", forced adjacency triangle free "
           <<tri_free<<", square free "<<square_free<<", both free "<<both_free<<"\n";
  return low>0 && both_free==0 && tri_free==0 && square_free==0 ? 0 : 1;
}
