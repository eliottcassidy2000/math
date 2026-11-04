#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
int main(int argc,char**argv){
 if(argc!=2)throw std::runtime_error("owner path");
 std::ifstream in(argv[1],std::ios::binary);std::vector<unsigned char>bytes((std::istreambuf_iterator<char>(in)),{});
 if(bytes.size()!=2u*(1u<<21))throw std::runtime_error("owner byte count");
 auto owner=[&](unsigned c){return unsigned(bytes[2*c])|(unsigned(bytes[2*c+1])<<8);};
 int nr;std::cin>>nr;unsigned orbitchecks=0;
 for(int ri=1;ri<=nr;++ri){std::array<unsigned,7>r;for(auto&v:r)std::cin>>v;std::array<int,7>p={0,1,2,3,4,5,6};do{
  unsigned code=0,bit=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j,++bit)if(r[p[i]]>>p[j]&1)code|=1u<<bit;
  if(owner(code)!=unsigned(ri))throw std::runtime_error("orbit ownership");
  ++orbitchecks;
 }while(std::next_permutation(p.begin(),p.end()));}
 std::array<unsigned,7>rows{};std::array<int,21>left{},right{};int k=0;
 for(int i=0;i<7;++i)for(int j=i+1;j<7;++j){left[k]=i;right[k]=j;++k;rows[j]|=1u<<i;}
 unsigned primes=0;
 for(unsigned code=0;code<(1u<<21);++code){
  if(code){unsigned changed=code^(code-1);while(changed){unsigned bit=__builtin_ctz(changed);changed&=changed-1;rows[left[bit]]^=1u<<right[bit];rows[right[bit]]^=1u<<left[bit];}}
  bool pr=true;
  for(unsigned mask=1;mask<127&&pr;++mask){if(__builtin_popcount(mask)<2)continue;bool module=true;for(int v=0;v<7&&module;++v)if(!(mask>>v&1)){unsigned z=rows[v]&mask;if(z&&z!=mask)module=false;}if(module)pr=false;}
  if(pr!=(owner(code)!=0))throw std::runtime_error("prime coverage");
  if(owner(code)>unsigned(nr))throw std::runtime_error("owner range");
  primes+=pr;
 }
 std::cout<<"{\"labeled_tournaments\":2097152,\"prime_labels\":"<<primes<<",\"representatives\":"<<nr<<",\"permutation_ownership_checks\":"<<orbitchecks<<"}\n";
}
