#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using I=std::int64_t;
int main(int argc,char**argv){
 if(argc!=4)throw std::runtime_error("signatures targets output");
 std::ifstream ts(argv[2]); std::set<I> targets;I t;
 while(ts>>t)targets.insert(t);
 I N=*targets.rbegin();std::vector<int> g(2*N+1,1);
 for(I v=2;v<=2*N;v++){I q=v/2;g[v]=g[q]+((v&1)?0:q+(q&1));}
 std::ifstream in(argv[1]);int ns;in>>ns;
 std::ofstream out(argv[3]);out<<"[\n";bool first=true;I tested=0,found=0;
 for(int z=0;z<ns;z++){
  int ci;I d,kp,km,ip,im,op,om,jp,jm;
  in>>ci>>d>>kp>>km>>ip>>im>>op>>om>>jp>>jm;
  if(ci<3)continue;
  I h=d+jp+jm,cap=N-h;
  auto mark=[&](I v,I a,I b){tested++;if(targets.count(v)){
   if(!first)out<<",\n";first=false;found++;
   out<<"["<<v<<","<<z<<","<<ci<<","<<a<<","<<b<<"]";
  }};
  if(!op||!ip)for(I q=0;q<N/d;q++)mark(d*g[2*q],0,2*q);
  I maxa=ip?cap/ip:cap/(d+kp);
  for(I a=0;a<=maxa;a++){
   I den=(d+kp)*a+op;if(!den)continue;
   I maxb=(cap-ip*a)/den;
   for(I b=0;b<=maxb;b++){
    if((!ip&&!b)||(!op&&!a))continue;
    I v=d*(I(g[a])+g[b]-1)+(d+kp)*a*b+ip*a-im*(a&1)+op*b-om*(b&1)+jp+jm;
    if(v<=N)mark(v,a,b);
   }
  }
 }
 out<<"\n]\n";std::cout<<"pairs="<<tested<<" matches="<<found<<"\n";
}
