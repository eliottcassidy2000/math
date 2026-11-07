#include <cstdint>
#include <iostream>
#include <array>
#include <stdexcept>
using u64=uint64_t;
constexpr u64 p=641,q=6700417,M=p*q;
struct S{u64 D,E,U,d;};
S step(S s,bool z){
 if(!z) return {s.D,0,(2*s.U)%M,s.d};
 return {(s.D+M-s.E+s.U)%M,s.d,(2*s.U+M-s.d)%M,s.d};
}
u64 F(S s){return (s.D+M-s.E+s.U)%M;}
std::array<u64,5> tested{}, qzero{}, pzero{}, bad{};
void visit(S initial,bool z,int depth){
 S s=initial;
 for(int k=1;k<=64;++k){
  s=step(s,z); u64 f=F(s);
  ++tested[depth]; if(f%p==0)++pzero[depth];
  if(f%q==0){++qzero[depth];if(f%p!=0)++bad[depth];}
  if(depth<4)visit(s,!z,depth+1);
 }
}
int main(){
 visit({1,0,2,1},false,1);visit({1,0,2,1},true,1);
 u64 expect=2;
 std::cout<<"{\"status\":\"PASS\",\"method\":\"repeated single-letter recurrence modulo 641*6700417; lengths 1..64\",\"rows\":[";
 for(int r=1;r<=4;++r){expect*=64;if(tested[r]!=expect||bad[r])throw std::runtime_error("certificate failed");
 if(r>1)std::cout<<",";
 std::cout<<"{\"runs\":"<<r<<",\"tested\":"<<tested[r]<<",\"q_zero\":"<<qzero[r]<<",\"p_zero\":"<<pzero[r]<<",\"q_only\":"<<bad[r]<<"}";}
 std::cout<<"]}\n";
}
