#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using I = std::int64_t;

int main(int argc, char** argv) {
  if(argc != 4 && argc != 5) throw std::runtime_error("usage: support signatures.txt bound output.json [mask.bin]");
  const I N=std::stoll(argv[2]);
  if(N<1 || N>1073741823) throw std::runtime_error("bound outside exact int range");
  std::ifstream in(argv[1]); int ns; in>>ns;
  if(!in || ns<1) throw std::runtime_error("missing signatures");
  std::vector<int> g(2*N+1,1);
  for(I v=2;v<=2*N;v++) {
    I q=v/2; g[v]=g[q]+((v&1)?0:q+(q&1));
  }
  std::vector<unsigned char> cover(N/2+1,0);
  I evaluations=0;
  auto mark=[&](I value,int ci) {
    if(value<=N) {
      if(value<1 || !(value&1)) throw std::runtime_error("nonpositive/even count");
      cover[value/2]|=(1<<ci);
    }
    evaluations++;
  };
  for(int z=0;z<ns;z++) {
    int ci; I d,kp,km,ip,im,op,om,jp,jm;
    in>>ci>>d>>kp>>km>>ip>>im>>op>>om>>jp>>jm;
    if(!in || ci<0 || ci>4) throw std::runtime_error("bad signature");
    I m=std::numeric_limits<I>::max();
    for(I x=0;x<2;x++) for(I y=0;y<2;y++)
      m=std::min(m,(d+km)*x*y-im*x-om*y+jp+jm);
    I cap=N-d-m;
    if(cap<0) continue;
    if(!op || !ip) for(I q=0;q<N/d;q++) mark(d*g[2*q],ci);
    I maxa=ip?cap/ip:cap/(d+kp);
    for(I a=0;a<=maxa;a++) {
      I denom=(d+kp)*a+op;
      if(!denom) continue;
      I maxb=(cap-ip*a)/denom;
      for(I b=0;b<=maxb;b++) {
        if((!ip && !b)||(!op && !a)) continue;
        I x=a&1,y=b&1;
        I value=d*(I(g[a])+g[b]-1)+(d+kp)*a*b+(d+km)*x*y
          +ip*a-im*x+op*b-om*y+jp+jm;
        mark(value,ci);
      }
    }
  }
  std::ofstream out(argv[3]);
  out<<"{\n\"bound\":"<<N<<",\n\"signature_count\":"<<ns
     <<",\n\"evaluated_parameter_pairs\":"<<evaluations<<",\n";
  auto array=[&](const char* label, auto predicate) {
    out<<'"'<<label<<"\":["; bool first=true;
    for(I h=1;h<=N;h+=2) if(predicate(cover[h/2],h)) {
      if(!first) out<<','; first=false; out<<h;
    } out<<"]";
  };
  array("five_core_missing",[](int c,I){return !c;}); out<<",\n";
  array("small_parent_extra_holes",[](int c,I h){return !(c&7)&&h!=7&&h!=21;}); out<<",\n";
  array("without_S5a_extra_holes",[](int c,I h){return !(c&23)&&h!=7&&h!=21;}); out<<",\n";
  array("without_S5b_extra_holes",[](int c,I h){return !(c&15)&&h!=7&&h!=21;}); out<<",\n";
  out<<"\"support_mask_counts\":{"; std::array<I,32> counts{};
  for(I h=1;h<=N;h+=2) counts[cover[h/2]]++;
  for(int i=0;i<32;i++){if(i)out<<',';out<<'"'<<i<<"\":"<<counts[i];}
  out<<"}\n}\n";
  std::cout<<"bound="<<N<<" signatures="<<ns<<" pairs="<<evaluations
           <<" missing="<<counts[0]<<"\n";
  if(argc==5) {
    std::ofstream raw(argv[4],std::ios::binary);
    raw.write(reinterpret_cast<const char*>(cover.data()),static_cast<std::streamsize>(cover.size()));
    if(!raw) throw std::runtime_error("binary support write failed");
  }
}
