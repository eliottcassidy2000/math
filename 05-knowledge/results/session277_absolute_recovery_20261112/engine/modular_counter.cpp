// Adapted from Session270 counter.cpp: raw response modulo 2^b. Compile: g++ -O3 -std=c++17 modular_counter.cpp -lgmpxx -lgmp
// Algorithm and arithmetic invariants are proved in PROOF.md. No approximate arithmetic.
#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
struct Key {int d; mpz_class r; bool operator==(Key const& x)const{return d==x.d && r==x.r;}};
struct Hash {size_t operator()(Key const& k)const {
 size_t h=k.d+0x9e3779b97f4a7c15ULL;
 for(size_t i=0;i<mpz_size(k.r.get_mpz_t());++i)
   h^=mpz_getlimbn(k.r.get_mpz_t(),i)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);
 return h;
}};
using States=std::unordered_map<Key,mpz_class,Hash>;
struct Transfer {int k,nd,base;mpz_class factor,mask,offset;};
static int charge(mpz_class const& r,int shift,unsigned long digit){
 mpz_class t;r.get_mpz_t();mpz_tdiv_q_2exp(t.get_mpz_t(),r.get_mpz_t(),shift);return t.get_ui()&digit;
}
int main(int argc,char**argv){
 if(argc<3){std::cerr<<"Usage: modular_counter n bits\n";return 1;}
 int n=std::stoi(argv[1]), bits=std::stoi(argv[2]), sw=-1;
 if(n<0 || n>10000)throw std::invalid_argument("n must be in 0..10000; large values are not a runtime promise");
 if(bits<0)throw std::invalid_argument("bits must be nonnegative");
 if(bits==0){std::cout<<0<<"\n";return 0;}
 mpz_class modulus=mpz_class(1)<<bits;
 auto begin=std::chrono::steady_clock::now();
 int width=1;while((1UL<<width)<=unsigned(n))++width;unsigned long digit=(1UL<<width)-1;
 std::vector<mpz_class> odds(n+1),masks(n/2+2);std::vector<int> twos(n+1);
 odds[0]=1;for(int i=1;i<=n;++i){int x=i,v=0;while(x%2==0){x/=2;++v;}odds[i]=odds[i-1]*x;twos[i]=twos[i-1]+v;}
 for(int i=0;i<int(masks.size());++i){mpz_setbit(masks[i].get_mpz_t(),width*i);--masks[i];}
 States states;states.emplace(Key{0,0},1);bool gauge=false;int start=1;long long transitions=0;size_t peak=1;
 for(int ell=n%2?n:n-1;ell>=3;ell-=2){
  mpz_class row=0;for(int j=start;j<ell;j+=2){mpz_class t=std::gcd(ell,j);t<<=width*((j-start)/2);row+=t;}
  std::vector<std::vector<Transfer>> tables(n+1);std::vector<bool> done(n+1,false);
  for(auto const&kv:states){int d=kv.first.d;if(done[d])continue;done[d]=true;mpz_class p=1;
   for(int k=0;k<=(n-d)/ell;++k){int nd=d+ell*k,R=n-nd;mpz_class denom=odds[d]*p*odds[k],factor;
    if(!mpz_divisible_p(odds[nd].get_mpz_t(),denom.get_mpz_t()))throw std::runtime_error("nonintegral transfer");
    mpz_divexact(factor.get_mpz_t(),odds[nd].get_mpz_t(),denom.get_mpz_t());
    int base=k*(ell-1)/2+ell*k*(k-1)/2-twos[k];
    int size=std::min(ell-2,R)<start?0:(std::min(ell-2,R)-start)/2+1;
    tables[d].push_back(Transfer{k,nd,base,std::move(factor),masks[size],row*k});p*=ell;
   }
  }
  States next;next.reserve(states.size()*2);
  for(auto const&kv:states){int q=charge(kv.first.r,width*((ell-start)/2),digit);
   for(auto const&t:tables[kv.first.d]){
    int shift=t.base+t.k*q;if(shift>=bits)continue;
    mpz_class val=kv.second*t.factor;
    if(shift<0)throw std::runtime_error("raw negative shift");val<<=shift;
    mpz_fdiv_r_2exp(val.get_mpz_t(),val.get_mpz_t(),bits);if(val==0)continue;
    mpz_class r=(kv.first.r+t.offset)&t.mask;auto &dest=next[Key{t.nd,std::move(r)}];dest+=val;mpz_fdiv_r_2exp(dest.get_mpz_t(),dest.get_mpz_t(),bits);++transitions;
   }
  }
  states=std::move(next);peak=std::max(peak,states.size());
 }
 mpz_class total=0;
 for(auto const&kv:states){int d=kv.first.d,R=n-d,q=charge(kv.first.r,0,digit);mpz_class den=odds[d]*odds[R],fac;
  if(!mpz_divisible_p(odds[n].get_mpz_t(),den.get_mpz_t()))throw std::runtime_error("unit tail nonintegral");
  mpz_divexact(fac.get_mpz_t(),odds[n].get_mpz_t(),den.get_mpz_t());mpz_class val=kv.second*fac;
  int shift=R*(R-1)/2+R*q-twos[R];if(shift>=bits)continue;
  val<<=shift;total+=val;mpz_fdiv_r_2exp(total.get_mpz_t(),total.get_mpz_t(),bits);
 }
 mpz_class inverse,answer;
 if(!mpz_invert(inverse.get_mpz_t(),odds[n].get_mpz_t(),modulus.get_mpz_t()))throw std::runtime_error("odd denominator not a unit");
 answer=total*inverse;mpz_fdiv_r_2exp(answer.get_mpz_t(),answer.get_mpz_t(),bits);
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
 std::cout<<answer<<"\n";std::cerr<<"{\"seconds\":"<<seconds<<",\"transitions\":"<<transitions<<",\"peak_states\":"<<peak<<"}\n";
}
