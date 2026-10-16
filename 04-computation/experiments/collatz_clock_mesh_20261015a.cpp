// Exact integer-clock and Haar census. No probabilistic or convergence premise.
// c++ -O3 -std=c++17 this.cpp -o /tmp/clock-mesh
// /tmp/clock-mesh 29 output.json
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using U=uint64_t; using I=int64_t; using Wide=__uint128_t;
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
std::string dec(Wide x){if(!x)return "0";std::string s;while(x){s.insert(s.begin(),char('0'+x%10));x/=10;}return s;}
Wide sq(I x){Wide y=x<0?Wide(-x):Wide(x);return y*y;}
int W(U n){return n%3==0?0:n%3==1?-1:1;}
int main(int argc,char**argv){try{
 need(argc==3,"exponent output.json");int K=std::stoi(argv[1]);need(1<=K&&K<=30,"exponent outside 1..30");U N=U(1)<<K;
 std::vector<U> bits((N+64)/64,0);
 auto parity=[&](U n){return int((bits[n/64]>>(n%64))&1);};
 U max_descent=0;
 for(U n=2;n<=N;++n){Wide x=n;U steps=0;
  do{need(x<(Wide(1)<<120),"native orbit overflow guard");x=(x&1)?(3*x+1)/2:x/2;need(++steps<1000000,"native descent cap");}while(x>=n);
  need(x>=1,"zero descent");if((steps+parity(U(x)))&1)bits[n/64]|=U(1)<<(n%64);
  if(steps>max_descent)max_descent=steps;
  if(!(n&(n-1)))std::cerr<<"grounded sources through "<<n<<'\n';
 }
 std::ofstream out(argv[2]);need(bool(out),"open output");
 out<<"{\"status\":\"FINITE-EXACT; all source labels independently grounded by descent\",\"source_limit\":"<<N<<",\"max_descent_steps\":"<<max_descent<<",\"octaves\":[";
 I J=0,H=0;
 for(int k=0;k<K;++k){U X=U(1)<<k;I incoming=J,sum=0;Wide S=0,E=0;
  std::array<I,31> stack{};std::array<Wide,31> innovations{};
  std::array<I,5> block_sum{};std::array<Wide,5> coarse_num{};std::array<int,5> mesh_log{};
  for(int q=0;q<5;++q)mesh_log[q]=std::max(0,k/2-q);
  for(U n=X;n<2*X;++n){int t=1-2*parity(n);H+=(n%3==2?2:-1)*t;
   need(H==W(n)*t-2*J,"every-source jump identity");sum+=J;S+=sq(J);E+=sq(H);
   I carry=J;U occupied=n-X;int level=0;
   while(occupied&1){innovations[level+1]+=sq(stack[level]-carry);carry+=stack[level];occupied>>=1;++level;}
   stack[level]=carry;
   for(int q=0;q<5;++q){block_sum[q]+=J;U b=U(1)<<mesh_log[q];if((n-X+1)%b==0){coarse_num[q]+=sq(block_sum[q]);block_sum[q]=0;}}
   if(parity(n)!=parity(n+1))J+=W(n)*(1-2*parity(n+1));
  }
  need(stack[k]==sum,"root sum");Wide rhs=sq(sum);for(int i=1;i<=k;++i)rhs+=innovations[i]*(X>>i);
  need(rhs==S*X,"exact Haar decomposition");
  out<<(k?",":"")<<"{\"k\":"<<k<<",\"X\":"<<X<<",\"incoming_current\":"<<incoming<<",\"current_sum\":"<<sum<<",\"current_square_sum\":\""<<dec(S)<<"\",\"H_energy\":\""<<dec(E)<<"\",\"H_at_last\":"<<H<<",\"haar_numerators\":[";
  for(int i=1;i<=k;++i)out<<(i>1?",":"")<<'"'<<dec(innovations[i])<<'"';out<<"],\"meshes\":[";
  for(int q=0;q<5;++q){U b=U(1)<<mesh_log[q];Wide coarse=sq(sum);for(int i=mesh_log[q]+1;i<=k;++i)coarse+=innovations[i]*(X>>i);
   need(coarse==coarse_num[q]*(X/b),"independent mesh sum versus Haar");
   out<<(q?",":"")<<"{\"j\":"<<q<<",\"b\":"<<b<<",\"coarse_numerator\":\""<<dec(coarse_num[q])<<"\",\"coarse_denominator\":"<<b<<'}';}
  out<<"]}";std::cerr<<"clock octave "<<k<<" checked\n";
 }
 out<<"]}\n";need(bool(out),"write output");
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
