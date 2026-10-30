// Exhaustive normalized B2 DAG states: complements are absorbed in gate tables.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>
using U=uint16_t;
U normalize(U a){return (a&1)?U(a^65535):a;}
int main(){
 std::array<U,4> inputs={0xaaaa,0xcccc,0xf0f0,0xff00};
 U target=normalize(U(0xaaaa^0xcccc^0xf0f0^0xff00^0x8000));
 std::unordered_set<uint64_t> states={0};
 for(int d=0;d<=4;++d){
  std::unordered_set<uint64_t> next;
  bool found=false;
  for(uint64_t key:states){
   std::array<U,8> sig{};std::copy(inputs.begin(),inputs.end(),sig.begin());
   for(int j=0;j<d;++j)sig[4+j]=U(key>>(16*j));
   for(int i=0;i<4+d;++i)for(int j=i+1;j<4+d;++j){
    U a=sig[i],b=sig[j];
    std::array<U,5> opts={U(a&b),U(a&U(~b)),U(U(~a)&b),U(a|b),U(a^b)};
    for(U z:opts){
     z=normalize(z);
     if(z==target)found=true;
     if(d==4||z==0||std::find(sig.begin(),sig.begin()+4+d,z)!=sig.begin()+4+d)continue;
     std::array<U,4> added{};
     for(int h=0;h<d;++h)added[h]=sig[4+h];
     added[d]=z;std::sort(added.begin(),added.begin()+d+1);
     uint64_t nk=0;for(int h=0;h<=d;++h)nk|=uint64_t(added[h])<<(16*h);
     next.insert(nk);
    }
   }
  }
  std::cout<<"available_gates="<<d<<" states="<<states.size()<<" target_reachable_with_next_gate="<<found<<"\n";
  if(found){std::cerr<<"Target reached; claimed six-gate lower bound fails\n";return 1;}
  states.swap(next);
 }
 // Six-gate upper bound is the independently checked formula in Python.
 std::cout<<"FINITE-EXACT: four-input parity plus point defect requires six B2 gates.\n";
}
