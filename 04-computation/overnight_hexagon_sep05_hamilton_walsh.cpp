// Exact Hamilton negative-cycle gap on every labelled switching class.
// Memory: one int16_t array; K9 uses exactly 512 MiB, no array copies.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>

void check(bool ok, const char* msg) { if (!ok) throw std::runtime_error(msg); }
uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}

int main(int argc, char** argv) {
  check(argc<=3,"arguments n and optional spectrum file");
  int n = argc >= 2 ? std::stoi(argv[1]) : 8;
  check(n==8 || n==9,"universe n=8 or9");
  int index[9][9]={},dim=0;
  for(int u=1;u<n;u++)for(int v=u+1;v<n;v++)index[u][v]=index[v][u]=dim++;
  uint32_t end=uint32_t(1)<<dim;
  std::vector<int16_t> spectrum(end,0);
  std::vector<int> path;
  for(int v=1;v<n;v++)path.push_back(v);
  int cycles=0;
  do {
    if(path.front()>path.back())continue;
    uint32_t mask=0;
    for(int i=0;i<n-2;i++)mask|=uint32_t(1)<<index[path[i]][path[i+1]];
    check(spectrum[mask]==0,"Hamilton-cycle mask unique");
    spectrum[mask]=1;
    cycles++;
  }while(std::next_permutation(path.begin(),path.end()));
  int factorial=1;
  for(int j=1;j<=n-1;j++)factorial*=j;
  check(cycles==factorial/2 && cycles<=std::numeric_limits<int16_t>::max(),"cycle count and int16 bound");
  // Each intermediate coefficient is a signed sum of a subset of cycles.
  // Hence its absolute value is <=cycles<2^15, including during this transform.
  for(uint32_t width=1;width<end;width*=2)
    for(uint32_t offset=0;offset<end;offset+=2*width)
      for(uint32_t j=0;j<width;j++) {
        int a=spectrum[offset+j],b=spectrum[offset+j+width];
        spectrum[offset+j]=int16_t(a+b);
        spectrum[offset+j+width]=int16_t(a-b);
      }
  check(spectrum[0]==cycles,"balanced spectrum");
  check(spectrum[end-1]==(n%2 ? -cycles : cycles),"antibalanced spectrum");
  int maximum=-cycles;
  std::set<uint32_t> ties;
  uint64_t semantic=1469598103934665603ULL;
  uint64_t keyed_sum=0;
  for(uint32_t h=0;h<end;h++) {
    check((cycles-int(spectrum[h]))%2==0,"integral parity count");
    semantic^=uint16_t(spectrum[h]); semantic*=1099511628211ULL;
    keyed_sum+=mix(h)*(uint64_t(int(spectrum[h])+cycles+1));
    if(h==0 || (n%2==0 && h==end-1))continue;
    if(spectrum[h]>maximum){maximum=spectrum[h];ties.clear();}
    if(spectrum[h]==maximum)ties.insert(h);
  }
  std::set<uint32_t> expected;
  for(int j=0;j<dim;j++)expected.insert(uint32_t(1)<<j);
  for(int v=1;v<n;v++) {
    uint32_t mask=0;
    for(int u=1;u<n;u++)if(u!=v)mask|=uint32_t(1)<<index[u][v];
    expected.insert(mask);
  }
  if(n%2==0){auto positive=expected;for(auto h:positive)expected.insert((end-1)^h);}
  int single=factorial/(n-1);
  std::cout<<"n="<<n<<" classes="<<end<<" cycles="<<cycles
           <<" minimum="<<(cycles-maximum)/2<<" expected="<<single
           <<" ties="<<ties.size()<<" equality_exact="<<(ties==expected)
           <<" memory_bytes="<<end*sizeof(int16_t)<<"\n";
  std::cout<<"semantic_fnv64="<<std::hex<<semantic<<std::dec<<"\n";
  std::cout<<"keyed_spectrum_sum64="<<std::hex<<keyed_sum<<std::dec<<"\n";
  check((cycles-maximum)/2==single,"Hamilton single-edge minimum");
  check(ties==expected,"Hamilton exact equality classes");
  if(argc==3){
    std::ofstream output(argv[2],std::ios::binary);
    output.write(reinterpret_cast<const char*>(spectrum.data()),std::streamsize(spectrum.size()*sizeof(int16_t)));
    check(bool(output),"spectrum artifact write");
  }
  std::cout<<"RESULT=PASS\n";
}
