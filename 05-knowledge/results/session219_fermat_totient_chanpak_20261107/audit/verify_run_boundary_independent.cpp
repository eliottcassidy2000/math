#include <array>
#include <cstdint>
#include <iostream>
using I=int64_t;
constexpr I p=641,q=6700417;
struct State {I D,E,U;};
I mod(I x,I m){x%=m;return x<0?x+m:x;}
// Formula derived by summing U_i-1 after first Z; d remains1.
State run(State s,bool z,I t,I m){
 if(!z)return {s.D,0,s.U*t%m};
 I u=mod(s.U-1,m);
 return {mod(s.D-s.E+1+(t-1)*u,m),1,(t*u+1)%m};
}
I read(State s,I m){return mod(s.D-s.E+s.U,m);}
std::array<I,65> pp,qq;
std::array<int,4> lengths;
I checked[5]={},qhits[5]={},escapes[5]={};
void search(State state,bool z,int depth,bool initial){
 if(depth){
  ++checked[depth];
  if(read(state,q)==0){
   ++qhits[depth];State small{1,0,2};
   for(int j=0;j<depth;j++)small=run(small,initial^(j%2),pp[lengths[j]],p);
   if(read(small,p)!=0)++escapes[depth];
  }
 }
 if(depth==4)return;
 for(int k=1;k<=64;k++){lengths[depth]=k;search(run(state,z,qq[k],q),!z,depth+1,initial);}
}
int main(){
 pp[0]=qq[0]=1;
 for(int k=1;k<=64;k++){pp[k]=pp[k-1]*2%p;qq[k]=qq[k-1]*2%q;}
 for(bool z:{false,true})search({1,0,2},z,0,z);
 for(int d=1;d<=4;d++){
  std::cout<<d<<" "<<checked[d]<<" "<<qhits[d]<<" "<<escapes[d]<<"\n";
  if(escapes[d])return 1;
 }
 State sp{1,0,2},sq=sp;
 int witness[]={31,15,9,5,9};
 for(int i=0;i<5;i++){sp=run(sp,i%2,pp[witness[i]],p);sq=run(sq,i%2,qq[witness[i]],q);}
 if(read(sp,p)!=275 || read(sq,q)!=0)return 2;
 std::cout<<"PASS five-run witness: mod641="<<read(sp,p)<<", mod6700417="<<read(sq,q)<<"\n";
}
