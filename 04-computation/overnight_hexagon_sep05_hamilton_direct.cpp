// Independent K8 Hamilton audit by literal cycle parity, no transform/imports.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

void check(bool ok,const char*msg){if(!ok)throw std::runtime_error(msg);}

int main(){
  const int n=8;
  int index[n][n]={},dim=0;
  for(int u=1;u<n;u++)for(int v=u+1;v<n;v++)index[u][v]=index[v][u]=dim++;
  std::vector<int> word={0,1,2,3,4,5,6,7};
  std::set<uint32_t> unique;
  do{
    uint32_t mask=0;
    for(int i=0;i<n;i++){
      int u=word[i],v=word[(i+1)%n];
      if(u&&v)mask^=uint32_t(1)<<index[u][v];
    }
    unique.insert(mask);
  }while(std::next_permutation(word.begin(),word.end()));
  check(unique.size()==2520,"full-permutation cycle universe");
  std::vector<uint32_t> cycles(unique.begin(),unique.end());
  const uint32_t end=uint32_t(1)<<dim;
  int minimum=2520;
  std::set<uint32_t> ties;
  uint64_t semantic=1469598103934665603ULL;
  for(uint32_t h=0;h<end;h++){
    int negative=0;
    for(auto cycle:cycles)negative+=__builtin_parity(h&cycle);
    semantic^=uint16_t(2520-2*negative);semantic*=1099511628211ULL;
    if(h==0||h==end-1){check(negative==0,"two zero classes");continue;}
    if(negative<minimum){minimum=negative;ties.clear();}
    if(negative==minimum)ties.insert(h);
  }
  std::set<uint32_t> expected;
  for(int j=0;j<dim;j++)expected.insert(uint32_t(1)<<j);
  for(int v=1;v<n;v++){
    uint32_t mask=0;
    for(int u=1;u<n;u++)if(u!=v)mask|=uint32_t(1)<<index[u][v];
    expected.insert(mask);
  }
  auto positive=expected;for(auto h:positive)expected.insert((end-1)^h);
  check(minimum==720,"K8 minimum");check(ties==expected,"K8 equality classes");
  std::cout<<"n=8 classes="<<end<<" cycles="<<cycles.size()<<" minimum="<<minimum
           <<" ties="<<ties.size()<<" equality_exact=1\n";
  std::cout<<"semantic_fnv64="<<std::hex<<semantic<<std::dec<<"\nRESULT=PASS\n";
}
