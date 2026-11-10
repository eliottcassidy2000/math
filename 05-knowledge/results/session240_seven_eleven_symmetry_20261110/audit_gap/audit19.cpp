// Independent direct Hamilton-path DP; no substitution/path-cover formula.
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(int argc,char**){
 const unsigned n=19, size=1u<<n; std::vector<unsigned> out(n);
 for(unsigned i=0;i<n;i++) for(unsigned j=0;j<n;j++) if(i!=j){
  bool edge=false;
  if(i<7&&j<7) { unsigned d=(j+7-i)%7; edge=d>=1&&d<=3; }
  else if(i>=7&&i<18&&j>=7&&j<18){unsigned d=(j+11-i)%11; edge=d==2||d==3||d==4||d==5||d==10;}
  else edge=(i<7&&j>=7&&j<18)||(i>=7&&i<18&&j==18)||(i==18&&j<7);
  if(edge) out[i]|=1u<<j;
 }
 std::vector<uint64_t> dp(static_cast<size_t>(size)*n);
 for(unsigned s=1;s<size;s++) for(unsigned vs=s;vs;vs&=vs-1){
  unsigned v=__builtin_ctz(vs), rest=s^(1u<<v); uint64_t val=0;
  if(!rest) val=1;
  else for(unsigned ws=rest&out[v];ws;ws&=ws-1) val+=dp[static_cast<size_t>(rest)*n+__builtin_ctz(ws)];
  dp[static_cast<size_t>(s)*n+v]=val;
 }
 uint64_t h=0; for(unsigned v=0;v<n;v++)h+=dp[static_cast<size_t>(size-1)*n+v];
 if(h!=(argc>1?3049061632ULL:3049061631ULL)) throw std::runtime_error("wrong expected H");
 std::cout<<"PASS: direct 19-vertex H = "<<h<<"\n";
}
