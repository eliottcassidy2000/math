// Independent finite break census. Native 128-bit descent and packed-bit control.
// Usage: binary parity26.bits output.json
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=uint64_t; using I=int64_t; using Wide=__uint128_t;
void need(bool ok,const char* s){if(!ok)throw std::runtime_error(s);}
int W(U n){return n%3==0?0:n%3==1?-1:1;}
int main(int argc,char**argv){try{
 need(argc==3,"bits output");const U N=1ULL<<27,lo=1ULL<<22;
 std::vector<unsigned char> p(N+1,0);
 for(U n=2;n<=N;++n){Wide x=n;U steps=0;
  do {need(x<(Wide(1)<<120),"orbit overflow");x=(x&1)?(3*x+1)/2:x/2;need(++steps<1000000,"orbit cap");}while(x>=n);
  p[n]=(steps+p[U(x)])&1;
 }
 const U saved=1ULL<<26;std::vector<U> bits((saved+64)/64);
 std::ifstream in(argv[1],std::ios::binary);in.read(reinterpret_cast<char*>(bits.data()),bits.size()*8);need(bool(in),"read bits");
 for(U n=1;n<=saved;++n)need(p[n]==((bits[n/64]>>(n%64))&1),"full packed-source cross-check");
 U B=0,gaps=0,large=0,last=0,maxgap=0,gap_length=0,long_length=0,clock_square=0,long_square=0;I J=0;
 std::array<U,64> shortgap{};std::array<U,3> residues{},classes{};
 std::array<std::array<U,3>,3> transitions{};std::array<I,13> direct{},twist{};std::array<U,13> count{};
 std::array<int,13> histW{},histV{};int cursor=0,nhist=0;
 std::array<I,27> incoming{},current_sum{};
 std::array<U,27> current_square{},source_energy{};
 I H=0;int octave=0;
 for(U n=1;n<N;++n){
  if(n==(1ULL<<(octave+1)))++octave;
  if(n==(1ULL<<octave))incoming[octave]=J;
  int label=1-2*int(p[n]);H+=(n%3==2?2:-1)*label;
  need(H==W(n)*label-2*J,"full source jump identity");
  current_sum[octave]+=J;
  Wide j2=Wide(J<0?-J:J)*Wide(J<0?-J:J),h2=Wide(H<0?-H:H)*Wide(H<0?-H:H);
  need(current_square[octave]+j2<(Wide(1)<<62) && source_energy[octave]+h2<(Wide(1)<<62),"octave moment guard");
  current_square[octave]+=U(j2);source_energy[octave]+=U(h2);
  if(p[n]==p[n+1])continue;
  if(n<lo){J+=W(n)*(1-2*int(p[n+1]));continue;}
  int c=n%3,w=W(n),v=w*(1-2*int(p[n+1]));++B;++classes[c];
  if(last){U gap=n-last;++gaps;++residues[gap%3];++transitions[last%3][c];maxgap=std::max(maxgap,gap);
   Wide energy=Wide(gap)*Wide(J<0?-J:J)*Wide(J<0?-J:J);need(energy+clock_square<(Wide(1)<<62),"clock moment guard");
   gap_length+=gap;clock_square+=U(energy);
   if(gap<64)++shortgap[gap];else{++large;long_length+=gap;long_square+=U(energy);}}
  for(int l=1;l<=12;++l)if(nhist>=l){int j=(cursor-l+13)%13;direct[l]+=v*histV[j];twist[l]+=(l%2?-1:1)*w*histW[j];++count[l];}
  histW[cursor]=w;histV[cursor]=v;cursor=(cursor+1)%13;nhist=std::min(13,nhist+1);last=n;J+=v;
 }
 need(gaps+1==B,"gap count");for(int l=1;l<=12;++l)need(direct[l]==twist[l],"twisted residue identity");
 std::ofstream out(argv[2]);out<<"{\"status\":\"FINITE-EXACT, no stationary or infinite-memory premise\",\"source_limit\":"<<N<<",\"window_start\":"<<lo<<",\"window_end_exclusive\":"<<N<<",\"packed_sources_checked\":"<<saved<<",\"breaks\":"<<B<<",\"gaps\":"<<gaps<<",\"gaps_at_least_64\":"<<large<<",\"max_gap\":"<<maxgap<<",\"gap_residue_counts\":[";
 for(int j=0;j<3;++j)out<<(j?",":"")<<residues[j];out<<"],\"break_class_counts\":[";
 for(int j=0;j<3;++j)out<<(j?",":"")<<classes[j];out<<"],\"transition_counts\":[";
 for(int a=0;a<3;++a){out<<(a?",":"")<<'[';for(int b=0;b<3;++b)out<<(b?",":"")<<transitions[a][b];out<<']';}out<<"],\"lag_moments\":[";
 for(int l=1;l<=12;++l)out<<(l>1?",":"")<<"{\"lag\":"<<l<<",\"pairs\":"<<count[l]<<",\"sum\":"<<direct[l]<<'}';
 out<<"],\"short_gap_counts\":[";for(int j=0;j<64;++j)out<<(j?",":"")<<shortgap[j];
 out<<"],\"interior_gap_lengths\":"<<gap_length<<",\"long_gap_lengths\":"<<long_length<<",\"interior_gap_J_square\":"<<clock_square<<",\"long_gap_J_square\":"<<long_square<<",\"octaves\":[";
 for(int k=0;k<27;++k)out<<(k?",":"")<<"{\"k\":"<<k<<",\"X\":"<<(1ULL<<k)<<",\"incoming_current\":"<<incoming[k]<<",\"current_sum\":"<<current_sum[k]<<",\"current_square_sum\":"<<current_square[k]<<",\"H_energy\":"<<source_energy[k]<<'}';
 out<<"]}\n";need(bool(out),"output");
 std::cerr<<"breaks "<<B<<", gaps >=64 "<<large<<", max gap "<<maxgap<<"; all packed sources agree\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
