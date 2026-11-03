#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <cstdint>
#include <set>
#include <string>
#include <vector>
int main(int argc,char**argv){
 if(argc<2||argc>3)return 2;
 const int born=argc==3?std::stoi(argv[2]):2;
 std::ifstream in(argv[1]);int n,id;bool first=true;std::cout<<"[";
 while(in>>id>>n){
  std::vector<int> old(n);for(int &r:old)in>>r;
  int N=n+born,edge[16][16],d=0;
  for(int i=0;i<N;++i)for(int j=i+1;j<N;++j)if(j>=n)edge[i][j]=d++;
  if(d>24||N>15)return 3;
  std::vector<long long> coeff(1<<d,0),walsh(1<<d,0);std::vector<int> p(N);std::iota(p.begin(),p.end(),0);
  long long permutations=0,compatible=0;
  do{
   ++permutations;int on=0,off=0;bool legal=true;
   for(int k=0;k<N-1;++k){
    int u=p[k],v=p[k+1];
    if(u<n&&v<n){if(!(old[u]>>v&1)){legal=false;break;}}
    else {int bit=1<<edge[std::min(u,v)][std::max(u,v)];if(u<v)on|=bit;else off|=bit;}
   }
   if(!legal)continue;++compatible;
   for(int sub=off;;sub=(sub-1)&off){coeff[on|sub]+=(__builtin_popcount((unsigned)sub)&1)?-1:1;if(!sub)break;}
   int used=on|off,scale=1<<(2*born-__builtin_popcount((unsigned)used));
   for(int sub=used;;sub=(sub-1)&used){walsh[sub]+=((__builtin_popcount((unsigned)(sub&on))&1)?-scale:scale);if(!sub)break;}
  }while(std::next_permutation(p.begin(),p.end()));
  auto values=coeff;
  for(int bit=0;bit<d;++bit)for(int mask=0;mask<(1<<d);++mask)if(mask>>bit&1)values[mask]+=values[mask^(1<<bit)];
  for(int stride=1;stride<(1<<d);stride*=2)for(int start=0;start<(1<<d);start+=2*stride)
   for(int k=0;k<stride;++k){auto a=walsh[start+k],b=walsh[start+k+stride];walsh[start+k]=a+b;walsh[start+k+stride]=a-b;}
  std::map<long long,long long> hist;std::vector<int> hits;
  for(int mask=0;mask<(1<<d);++mask){auto h=values[mask];if(walsh[mask]!=(1LL<<(2*born))*h)return 6;
   if(h<=0||h%2==0)return 4;++hist[h];if(h==111)hits.push_back(mask);}
  if(!first)std::cout<<",";first=false;
  std::cout<<"{\"id\":"<<id<<",\"old_vertices\":"<<n<<",\"variables\":"<<d
           <<",\"born\":"<<born<<",\"transforms_agree_every_mask\":true,\"permutations\":"<<permutations<<",\"compatible_permutations\":"<<compatible
           <<",\"H111_count\":"<<hits.size()<<",\"H111_first_masks\":[";
  for(size_t i=0;i<std::min<size_t>(hits.size(),8);++i){if(i)std::cout<<",";std::cout<<hits[i];}
  std::cout<<"],\"histogram\":{";bool hf=true;
  for(auto [h,c]:hist){if(!hf)std::cout<<",";hf=false;std::cout<<"\""<<h<<"\":"<<c;}
  std::cout<<"},\"probe_H\":{";
  std::set<int> probes{0,(1<<d)-1};
  for(int bit=0;bit<d;++bit){probes.insert(1<<bit);probes.insert(((1<<d)-1)^(1<<bit));}
  for(size_t i=0;i<std::min<size_t>(hits.size(),8);++i)probes.insert(hits[i]);
  uint32_t random_state=55+100*born+id;
  for(int j=0;j<20;++j){random_state^=random_state<<13;random_state^=random_state>>17;random_state^=random_state<<5;
   probes.insert(random_state&((1<<d)-1));}
  bool pf=true;for(int mask:probes){if(!pf)std::cout<<",";pf=false;std::cout<<"\""<<mask<<"\":"<<values[mask];}
  std::cout<<"},\"coefficients\":[";bool cf=true;
  for(int mask=0;mask<(1<<d);++mask)if(coeff[mask]){
   if(__builtin_popcount((unsigned)mask)>2*born)return 5;
   if(born>2)continue;
   if(!cf)std::cout<<",";cf=false;std::cout<<"["<<mask<<","<<coeff[mask]<<"]";
  }
  std::cout<<"]}";
 }
 std::cout<<"]\n";
}
