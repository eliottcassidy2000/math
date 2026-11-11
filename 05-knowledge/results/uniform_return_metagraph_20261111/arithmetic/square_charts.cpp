// Exact first-three-core inverse, recording minimum short flank for r^2+1.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>
using I=std::int64_t;
int bits(I x){int k=0;while(x){x>>=1;k++;}return k;}
struct Hit{I A=-1,B=-1;int ci=-1,cut=-1,order=0,mask=0;};
int main(int argc,char**argv){
 if(argc!=4)throw std::runtime_error("usage: raw-signatures bound output.tsv");
 I N=std::stoll(argv[2]);
 if(N<5||N>16777215)throw std::runtime_error("declared bound exceeded");
 std::vector<int> g(2*N+1,1);
 for(I x=2;x<=2*N;x++)g[x]=g[x/2]+((x&1)?0:x/2+((x/2)&1));
 std::vector<std::int16_t> target(N+1,-1);
 std::vector<Hit> hits;
 for(I r=2;r*r+1<=N;r+=2){target[r*r+1]=int(hits.size());hits.emplace_back();}
 I checked=0;
 auto mark=[&](I v,int ci,int cut,I A,I B){
  checked++;
  if(v>N)return;
  if(v<1||!(v&1))throw std::runtime_error("invalid count");
  int t=target[v];if(t<0)return;
  auto &h=hits[t];h.mask|=1<<ci;
  int core_order=ci==0?1:ci==1?3:4;
  int order=core_order+1+bits(A)+bits(B);
  if(h.A<0||std::make_tuple(std::min(A,B),order,ci,cut,A,B)<std::make_tuple(std::min(h.A,h.B),h.order,h.ci,h.cut,h.A,h.B)){
   h.A=A;h.B=B;h.ci=ci;h.cut=cut;h.order=order;
  }
 };
 std::ifstream f(argv[1]);int count;f>>count;
 for(int row=0;row<count;row++){
  int ci,cut;I h,d,L,i,o,s,e;f>>ci>>cut>>h>>d>>L>>i>>o>>s>>e;
  if(!f)throw std::runtime_error("bad observer input");
  if(ci>2||h>N)continue;
  if(i==0||o==0){
   if(h!=d||(i==0&&s!=0)||(o==0&&e!=0))throw std::runtime_error("axis guard");
   for(I q=0;q<N/d;q++){
    if(i==0)mark(d*g[2*q],ci,cut,2*q,0);
    if(o==0)mark(d*g[2*q],ci,cut,0,2*q);
   }
  }
  I cap=N-h,maxB=o?cap/o:cap/L;
  for(I B=0;B<=maxB;B++){
   I den=L*B+i;if(!den)continue;
   I base=h+d*(g[B]-1)+o*B+e*(B&1);if(base>N)continue;
   I maxA=(N-base)/den,minA=o?0:1;
   for(I A=minA;A<=maxA;A++)mark(base+d*(g[A]-1)+den*A+s*(A&1),ci,cut,A,B);
  }
 }
 std::ofstream out(argv[3]);out<<"root\tvalue\tcore\tcut\tA\tB\tmin_flank\torder\tcore_mask\n";
 for(unsigned t=0;t<hits.size();t++){
  I r=2*(t+1);auto h=hits[t];
  out<<r<<'\t'<<r*r+1<<'\t'<<h.ci<<'\t'<<h.cut<<'\t'<<h.A<<'\t'<<h.B<<'\t'<<std::min(h.A,h.B)<<'\t'<<h.order<<'\t'<<h.mask<<'\n';
 }
 if(!out)throw std::runtime_error("output failure");
 std::cout<<"pairs="<<checked<<" targets="<<hits.size()<<'\n';
}
