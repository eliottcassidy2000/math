// Independent Hamilton-cycle audit by deletion to Hamilton paths and attachment
// quadratic forms. Does not build the full graph's cycle frequency/spectrum.
// At n=9 stores 28 arrays of 2^21 int16_t =112 MiB, plus a 128-entry buffer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

void check(bool ok,const char*msg){if(!ok)throw std::runtime_error(msg);}
uint64_t mix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);}

int main(int argc,char**argv){
  check(argc<=3,"arguments n and optional reference spectrum");
  int n=argc>=2?std::stoi(argv[1]):9;
  check(n==8||n==9,"universe n8 or9");
  const int m=n-1;
  int parent_index[9][9]={},dim=0;
  for(int u=1;u<m;u++)for(int v=u+1;v<m;v++)parent_index[u][v]=parent_index[v][u]=dim++;
  int full_index[9][9]={},full_dim=0;
  for(int u=1;u<n;u++)for(int v=u+1;v<n;v++)full_index[u][v]=full_index[v][u]=full_dim++;
  const uint32_t parents=uint32_t(1)<<dim, attachments=uint32_t(1)<<(m-1),end=uint32_t(1)<<full_dim;
  std::vector<int16_t> reference;
  if(argc==3){
    reference.resize(end);std::ifstream input(argv[2],std::ios::binary);
    input.read(reinterpret_cast<char*>(reference.data()),std::streamsize(reference.size()*sizeof(int16_t)));
    check(bool(input)&&input.peek()==EOF,"reference spectrum exact byte length");
  }
  std::vector<std::pair<int,int>> pairs;
  for(int u=0;u<m;u++)for(int v=u+1;v<m;v++)pairs.push_back({u,v});
  std::vector<int16_t> coefficients(uint64_t(parents)*pairs.size(),0);
  int path_count=1;for(int j=1;j<=m-2;j++)path_count*=j;
  for(size_t edge=0;edge<pairs.size();edge++){
    auto[u,v]=pairs[edge];
    auto* data=&coefficients[uint64_t(edge)*parents];
    std::vector<int> interior;
    for(int w=0;w<m;w++)if(w!=u&&w!=v)interior.push_back(w);
    int generated=0;
    do{
      int previous=u;uint32_t mask=0;
      for(int w:interior){if(previous&&w)mask^=uint32_t(1)<<parent_index[previous][w];previous=w;}
      if(previous&&v)mask^=uint32_t(1)<<parent_index[previous][v];
      data[mask]++;generated++;
    }while(std::next_permutation(interior.begin(),interior.end()));
    check(generated==path_count,"ordered-endpoint Hamilton path count");
    for(uint32_t width=1;width<parents;width*=2)
      for(uint32_t start=0;start<parents;start+=2*width)
        for(uint32_t j=0;j<width;j++){
          int a=data[start+j],b=data[start+j+width];
          data[start+j]=int16_t(a+b);data[start+j+width]=int16_t(a-b);
        }
    check(data[0]==path_count,"positive parent path sum");
  }
  const int total_cycles=path_count*int(pairs.size());
  std::vector<uint32_t> parent_full_bit(dim),attachment_full(attachments,0);
  for(int u=1;u<m;u++)for(int v=u+1;v<m;v++)parent_full_bit[parent_index[u][v]]=uint32_t(1)<<full_index[u][v];
  for(uint32_t t=1;t<attachments;t++){
    uint32_t bit=t&-t;int v=__builtin_ctz(bit)+1;
    attachment_full[t]=attachment_full[t^bit]|(uint32_t(1)<<full_index[v][m]);
  }
  int maximum=-total_cycles;std::set<uint32_t> ties;
  uint64_t keyed_sum=0,checked=0;
  uint32_t parent_full=0,previous_h=0;
  for(uint32_t h=0;h<parents;h++){
    uint32_t changes=h^previous_h;
    while(changes){uint32_t bit=changes&-changes;parent_full^=parent_full_bit[__builtin_ctz(bit)];changes^=bit;}
    previous_h=h;
    std::array<int32_t,128> q={};
    for(size_t edge=0;edge<pairs.size();edge++){
      auto[u,v]=pairs[edge];
      uint32_t monomial=(u?uint32_t(1)<<(u-1):0)|(uint32_t(1)<<(v-1));
      q[monomial]=coefficients[uint64_t(edge)*parents+h];
    }
    // Each output evaluates the exact quadratic attachment polynomial.
    for(uint32_t width=1;width<attachments;width*=2)
      for(uint32_t start=0;start<attachments;start+=2*width)
        for(uint32_t j=0;j<width;j++){
          int32_t a=q[start+j],b=q[start+j+width];
          q[start+j]=a+b;q[start+j+width]=a-b;
        }
    for(uint32_t t=0;t<attachments;t++){
      uint32_t full=parent_full|attachment_full[t];
      int value=q[t];checked++;
      if(!reference.empty())check(value==reference[full],"independent per-character spectrum equality");
      check(value>=-total_cycles&&value<=total_cycles&&((total_cycles-value)%2==0),"exact attachment parity range");
      keyed_sum+=mix(full)*uint64_t(value+total_cycles+1);
      if(full==0||(n%2==0&&full==end-1)){check(value==total_cycles,"zero class");continue;}
      if(value>maximum){maximum=value;ties.clear();}
      if(value==maximum)ties.insert(full);
    }
  }
  std::set<uint32_t> expected;
  for(int j=0;j<full_dim;j++)expected.insert(uint32_t(1)<<j);
  for(int v=1;v<n;v++){
    uint32_t mask=0;
    for(int u=1;u<n;u++)if(u!=v)mask|=uint32_t(1)<<full_index[u][v];
    expected.insert(mask);
  }
  if(n%2==0){auto edges=expected;for(auto h:edges)expected.insert((end-1)^h);}
  const int single=2*total_cycles/(n-1);
  check(checked==end,"complete parent times attachment universe");
  std::cout<<"n="<<n<<" classes="<<checked<<" path_channels="<<pairs.size()
           <<" paths_per_channel="<<path_count<<" minimum="<<(total_cycles-maximum)/2
           <<" expected="<<single<<" ties="<<ties.size()<<" equality_exact="<<(ties==expected)
           <<" memory_bytes="<<(coefficients.size()+reference.size())*sizeof(int16_t)
           <<" full_reference_compared="<<(!reference.empty())<<"\n";
  std::cout<<"keyed_spectrum_sum64="<<std::hex<<keyed_sum<<std::dec<<"\n";
  check((total_cycles-maximum)/2==single,"extension minimum");check(ties==expected,"extension equality");
  std::cout<<"RESULT=PASS\n";
}
