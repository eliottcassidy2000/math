// Exhaust all four-letter two-apex presentations through a stated chain length.
// Independent state encoding from verify.py: counts indexed by (order,slot statuses).
// Stores one compact witness word for every scalar H. No isomorphism filtering.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <limits>
using State=std::array<uint64_t,36>;
std::vector<int> transitions[4][36];std::vector<int> final_states;
struct Witness {uint32_t code;uint8_t length;};
std::map<uint64_t,Witness> witnesses;
std::vector<uint64_t> counts;
uint64_t visited=0;
uint64_t value(const State&s){uint64_t h=0;for(int i:final_states)h+=s[i];return h;}
void search(const State&s,int depth,int maximum,uint32_t word){
 ++visited;uint64_t h=value(s);auto it=witnesses.find(h);
 if(it==witnesses.end()||it->second.length>depth)witnesses[h]={word,(uint8_t)depth};
 if(depth==maximum)return;
 for(int letter=0;letter<4;letter++){
  State next{};
  for(int i=0;i<36;i++)if(s[i])for(int j:transitions[letter][i])next[j]+=s[i];
  search(next,depth+1,maximum,(word<<2)|letter);
 }
}
std::string word(Witness w){std::string r(w.length,'0');for(int i=w.length-1;i>=0;i--){r[i]='0'+(w.code&3);w.code>>=2;}return r;}
int main(int argc,char**argv){
 int maximum=argc>1?std::stoi(argv[1]):12;
 if(maximum<0||maximum>15){std::cerr<<"Supported chain lengths are 0..15.\n";return 2;}
 for(int order=0;order<2;order++)for(int before=0;before<3;before++)for(int between=0;between<3;between++)for(int after=0;after<2;after++){
  int i=order*18+before*6+between*2+after;
  if(before<2&&(between==1||(between==0&&order==0)))final_states.push_back(i);
  for(int letter=0;letter<4;letter++){
   int first=(letter>>order)&1,last=(letter>>(1-order))&1;
   transitions[letter][i].push_back(order*18+(first?2:1)*6+between*2+after);
   if(between!=0||first)transitions[letter][i].push_back(order*18+before*6+(last?2:1)*2+after);
   if(after!=0||last)transitions[letter][i].push_back(order*18+before*6+between*2+1);
  }
 }
 State initial{};initial[0]=initial[18]=1;search(initial,0,maximum,0);
 uint64_t expected=0,power=1;for(int i=0;i<=maximum;i++){expected+=power;power*=4;}
 if(visited!=expected){std::cerr<<"Universe mismatch\n";return 3;}
 uint64_t gap=1;while(gap==7||gap==21||witnesses.count(gap))gap+=2;
 std::cout<<"{\"chain_length_max\":"<<maximum<<",\"presentations\":"<<visited
 <<",\"distinct_values\":"<<witnesses.size()<<",\"maximum_H\":"<<witnesses.rbegin()->first
 <<",\"odd_prefix_except_7_21\":"<<gap-2<<",\"first_additional_finite_gap\":"<<gap<<"}\n";
 // Every attained value gets a witness. The release Python driver verifies
 // a chosen bounded prefix independently by full vertex-subset dynamic programming.
 std::ofstream out("scan_witnesses.tsv");
 for(const auto& [h,w]:witnesses)out<<h<<"\t"<<word(w)<<"\n";
 return 0;
}
