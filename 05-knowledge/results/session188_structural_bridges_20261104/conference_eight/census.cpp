#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using namespace std;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
vector<unsigned> graph(unsigned n,unsigned code){vector<unsigned>a(n);unsigned k=0;for(unsigned i=0;i<n;i++)for(unsigned j=i+1;j<n;j++,k++){if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;}return a;}
bool normalized(const vector<unsigned>&a){
 for(unsigned i=0;i<7;i++)if(__builtin_popcount(a[i])!=3)return false;
 for(unsigned i=0;i<7;i++)for(unsigned j=i+1;j<7;j++){
  int dot=0;for(unsigned k=0;k<7;k++)if(k!=i&&k!=j)dot+=((a[i]>>k&1)==(a[j]>>k&1))?1:-1;
  if(dot!=-1)return false;
 }
 return true;
}
bool conference(const vector<unsigned>&a){for(unsigned i=0;i<8;i++)for(unsigned j=i+1;j<8;j++){int dot=0;for(unsigned k=0;k<8;k++)if(k!=i&&k!=j)dot+=((a[i]>>k&1)==(a[j]>>k&1))?1:-1;if(dot)return false;}return true;}
int main(){
 vector<unsigned>bases;unsigned regular=0;
 for(unsigned code=0;code<(1u<<21);code++){
  auto a=graph(7,code);bool reg=true;for(auto r:a)if(__builtin_popcount(r)!=3){reg=false;break;}
  if(!reg){continue;}regular++;if(normalized(a))bases.push_back(code);
 }
 set<unsigned>codes;
 for(unsigned base:bases)for(unsigned mask=0;mask<128;mask++){
  auto a=graph(8,(base<<7)|127);unsigned code=0,k=0;
  for(unsigned i=0;i<8;i++)for(unsigned j=i+1;j<8;j++,k++){
   unsigned di=i?((mask>>(i-1))&1):0,dj=(mask>>(j-1))&1;
   code|=(((a[i]>>j)&1)^di^dj)<<k;
  }
  need(conference(graph(8,code)),"generated conference identity");need(codes.insert(code).second,"nonunique normalized switching");
 }
 need(bases.size()==240&&codes.size()==30720,"conference counts");
 ifstream input("../eight_pair_quotient/blind8.bin",ios::binary);need(bool(input),"missing h-zero bank");vector<unsigned>old;char raw[4];while(input.read(raw,4)){unsigned c=0;for(unsigned k=0;k<4;k++)c|=unsigned((unsigned char)raw[k])<<(8*k);old.push_back(c);}need(input.eof()&&input.gcount()==0,"incomplete h-zero bank word");need(vector<unsigned>(codes.begin(),codes.end())==old,"full exact conference and h-zero set equality");
 cout<<"{\"status\":\"PASS\",\"all_seven_tournaments\":2097152,\"regular_seven_tournaments\":"<<regular<<",\"normalized_blocks\":"<<bases.size()<<",\"switchings_per_block\":128,\"complete_eight_conference_count\":"<<codes.size()<<",\"all_generated_matrices_checked\":"<<codes.size()<<",\"full_h_zero_set_equality\":true,\"normalized_seven_codes\":[";
 for(unsigned i=0;i<bases.size();i++){if(i)cout<<",";cout<<bases[i];}cout<<"]}\n";
}
