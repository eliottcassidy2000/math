// Independent raw-cut positive compiler and B-first target traversal.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=std::int64_t;

int direct_g(unsigned n) {
  // Direct weighted binary inversion sum: each zero contributes its
  // preceding binary prefix plus the immediately preceding bit.
  if(!n) return 1;
  unsigned top=1;
  while(top<=n/2) top*=2;
  I ans=1,prefix=0;int previous=0;
  for(unsigned bit=top;bit;bit>>=1) {
    int current=(n&bit)?1:0;
    if(!current) ans+=prefix+previous;
    prefix=2*prefix+current;previous=current;
  }
  return static_cast<int>(ans);
}

int main(int argc,char**argv) {
  if(argc!=5) throw std::runtime_error("usage: support_positive raw.txt bound output.json mask.bin");
  I N=std::stoll(argv[2]);
  if(N<1 || N>1073741823) throw std::runtime_error("bound outside exact integer range");
  std::ifstream input(argv[1]);int rows;input>>rows;
  if(!input || rows!=90) throw std::runtime_error("all90 raw cuts required");
  std::vector<int> g(2*N+1);
  for(I n=0;n<=2*N;n++) g[n]=direct_g(static_cast<unsigned>(n));
  std::vector<unsigned char> support(N/2+1,0);
  I examined=0;
  auto mark=[&](I value,int core) {
    examined++;
    if(value>N) return;
    if(value<1 || !(value&1)) throw std::runtime_error("invalid Hamilton count");
    support[value/2]|=1<<core;
  };
  for(int row=0;row<rows;row++) {
    int core,cut;I h,d,L,i,o,s,e;
    input>>core>>cut>>h>>d>>L>>i>>o>>s>>e;
    if(!input || core<0 || core>4 || d<=0 || L<=0 || i<0 || o<0 || s<0 || e<0)
      throw std::runtime_error("invalid positive response");
    if(h>N) continue;
    // A zero slope is a uniform core cut. Reduce the sole free-axis word
    // by removing trailing ones; its even representative is2q with q<N/d.
    if(i==0 || o==0) {
      if(h!=d || (i==0 && s!=0) || (o==0 && e!=0))
        throw std::runtime_error("uniform-axis identity failed");
      for(I q=0;q<N/d;q++) mark(d*g[2*q],core);
    }
    // FixB first. Except on the normalized axis, oB or LAB boundsB.
    I cap=N-h,maxB=o?cap/o:cap/L;
    for(I B=0;B<=maxB;B++) {
      I denominator=L*B+i;
      if(!denominator) continue; // the A-only uniform axis handled above
      I base=h+d*(g[B]-1)+o*B+e*(B&1);
      if(base>N) continue;
      I maxA=(N-base)/denominator;
      I minA=o?0:1; // the B-only uniform axis handled above
      for(I A=minA;A<=maxA;A++) {
        I value=base+d*(g[A]-1)+denominator*A+s*(A&1);
        mark(value,core);
      }
    }
  }
  std::array<I,32> counts{};
  for(unsigned char mask:support) counts[mask]++;
  std::ofstream output(argv[3]);
  output<<"{\n\"bound\":"<<N<<",\n\"raw_cut_count\":"<<rows
        <<",\n\"evaluated_parameter_pairs\":"<<examined<<",\n\"support_mask_counts\":{";
  for(int k=0;k<32;k++){if(k)output<<',';output<<'\"'<<k<<"\":"<<counts[k];}
  output<<"}\n}\n";
  std::ofstream binary(argv[4],std::ios::binary);
  binary.write(reinterpret_cast<const char*>(support.data()),static_cast<std::streamsize>(support.size()));
  if(!output || !binary) throw std::runtime_error("output write failed");
  std::cout<<"bound="<<N<<" rawcuts="<<rows<<" pairs="<<examined<<" missing="<<counts[0]<<"\n";
}
