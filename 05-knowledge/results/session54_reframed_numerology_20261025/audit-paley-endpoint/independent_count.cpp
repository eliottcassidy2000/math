#include <iostream>
#include <vector>
#include <cstdint>
// Independently written full vertex-mask predecessor recurrence. Start fixed at 0.
int main(){for(int q: {3,7,11,19}){
 unsigned lim=1u<<q; std::vector<uint64_t> d(uint64_t(lim)*q,0);d[q]=1;
 std::vector<int> residue(q,0);for(int k=1;k<q;k++)residue[k*k%q]=1;
 for(unsigned mask=3;mask<lim;mask+=2){
  for(int end=1;end<q;end++)if(mask&(1u<<end)){
   unsigned before=mask^(1u<<end);uint64_t value=0;
   for(int pred=0;pred<q;pred++)if((before&(1u<<pred)) && residue[(end-pred+q)%q])value+=d[uint64_t(before)*q+pred];
   d[uint64_t(mask)*q+end]=value;
  }
 }
 std::cout<<q;for(int end=0;end<q;end++)std::cout<<" "<<d[uint64_t(lim-1)*q+end];std::cout<<"\n";
}}
