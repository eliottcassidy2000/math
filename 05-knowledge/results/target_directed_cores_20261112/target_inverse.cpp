// Exact prescribed-target inversion of inherited binary-shell observers.
// Compile: c++ -O3 -std=c++17 target_inverse.cpp -o target_inverse
// Input: n rank h e s2 p2 e2 s3, one complete core per line.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using I = int64_t;
using W = __int128_t;
struct Core { I n,rank,h,e,s,p,f,t; };
I shell(I n) {
  I ans=1;
  while(n>1) { I m=n/2; if(!(n&1)) ans+=m+(m&1); n=m; }
  return ans;
}
uint64_t visits;
bool inverse(I N,I a,I b,I c,I d,I place,I prefix,I &hit) {
  ++visits;
  I numer=N-d-a-c;
  I lo=numer<=0?0:(numer+a+b-1)/(a+b);
  I hi=(N-d-a)/b;
  if(N<d+a || hi<lo) return false;
  if(hi-lo<=3) {
    for(I n=lo;n<=hi;++n)
      if(W(a)*shell(n)+W(b)*n+c*(n&1)+d==N) {
        hit=prefix+place*n; return true;
      }
    return false;
  }
  return inverse(N,a,2*b+a,a,d,2*place,prefix,hit)
      || inverse(N,a,2*b,0,d+b+c,2*place,prefix+place,hit);
}
int main(int argc,char**argv) {
  if(argc!=4) { std::cerr<<"usage: inverse core_bank.txt max_k output.json\n";return 2; }
  std::ifstream in(argv[1]); std::vector<Core> bank; Core c;
  while(in>>c.n>>c.rank>>c.h>>c.e>>c.s>>c.p>>c.f>>c.t) bank.push_back(c);
  if(bank.empty()) return 3;
  int maxk=std::stoi(argv[2]); if(maxk<1||maxk>48)return 4;
  std::ofstream out(argv[3]);out<<"{\n  \"scope\": \"First ready source found in the supplied complete finite core bank, with arbitrary outer word\",\n  \"core_count\": "<<bank.size()<<",\n  \"targets\": [\n";
  for(int k=1;k<=maxk;++k) {
    I N=405*(I(1)<<k)+71,hit=-1; visits=0; size_t tests=0;
    auto start=std::chrono::steady_clock::now(); const Core* chosen=nullptr;
    for(const auto &v:bank) {
      ++tests;I G=v.s+v.p,d=G-v.e-v.f-1;
      if(inverse(N,v.h,G,v.e,d,1,0,hit) && (hit>0||v.n>1)) {chosen=&v;break;}
    }
    double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    out<<"    {\"k\": "<<k<<", \"target\": "<<N<<", \"cores_tested\": "<<tests<<", \"visits\": "<<visits;
    if(chosen)out<<", \"core_order\": "<<chosen->n<<", \"core_rank\": "<<chosen->rank<<", \"word_integer\": "<<hit;
    else out<<", \"found\": false";
    out<<"}"<<(k<maxk?",":"")<<"\n";out.flush();
    std::cout<<"k="<<k<<" target="<<N<<" source="<<(chosen?chosen->n:-1)<<":"<<(chosen?chosen->rank:-1)<<" word="<<hit<<" tests="<<tests<<" visits="<<visits<<" seconds="<<secs<<std::endl;
  }
  out<<"  ]\n}\n";
  return 0;
}
