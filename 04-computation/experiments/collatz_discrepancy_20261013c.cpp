// Exact all-prefix mod-3 discrepancy census; Terras map, ROOT parity.
// Usage: prog E K output.json prefix.bits; even jump K (4..16).
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using U=uint64_t;using W=__uint128_t;using I=int64_t;
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
std::string dec(W x){std::string s;do{s+='0'+x%10;x/=10;}while(x);std::reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){try{
 check(argc==5||argc==6,"E K output bits [block moments]");int E=std::stoi(argv[1]),K=std::stoi(argv[2]);
 check(E>=10&&E<=34&&K>=4&&K<=16&&K%2==0,"scope");U N=1ULL<<E,mask=(1ULL<<K)-1;
 struct J{U p,b;};std::vector<J>jump(1ULL<<K);
 for(U r=0;r<=mask;r++){U z=r,p=1;for(int k=0;k<K;k++){if(z&1){z=(3*z+1)/2;p*=3;}else z/=2;}jump[r]={p,z};}
 std::vector<U>bits((N+64)/64);auto bit=[&](U n){return(bits[n/64]>>(n%64))&1;};
 U maxloops=0;W maxstate=0;auto st=std::chrono::steady_clock::now();
 for(U n=2;n<=N;n++){U m=n,flip=0,l=0;
  if(!(n&1)){m=n/2;flip=1;}else do{J j=jump[m&mask];W z=W(m>>K)*j.p+j.b;maxstate=std::max(maxstate,z);++l;
   while(z>std::numeric_limits<U>::max()){check(z<(W(1)<<120),"overflow");for(int k=0;k<K;k++){check(z<(W(1)<<126),"step overflow");z=(z&1)?(3*z+1)/2:z/2;}++l;check(l<1000000,"descent guard");}
   m=U(z);check(l<1000000,"descent guard");
  }while(m>=n);check(m&&m<n,"earlier source");if(bit(m)^flip)bits[n/64]|=1ULL<<(n%64);maxloops=std::max(maxloops,l);
  if(n%(1ULL<<30)==0)std::cerr<<"generated "<<n<<" seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count()<<'\n';
 }
 U controls=0;auto control=[&](U n){W z=n;U t=0;while(z!=1){check(z<(W(1)<<120),"control overflow");z=(z&1)?(3*z+1)/2:z/2;check(++t<1000000,"control cap");}check(bit(n)==t%2,"strict control");++controls;};
 for(U n=1;n<=std::min(N,1ULL<<18);n++)control(n);U rng=20261013;for(int i=0;i<16384;i++){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;control(1+rng%N);}for(int e=0;e<=E;e++)control(1ULL<<e);
 U saved=std::min(N,1ULL<<26),words=(saved+64)/64;std::ofstream bf(std::string(argv[4])+".bits",std::ios::binary);bf.write(reinterpret_cast<char*>(bits.data()),words*8);check(bool(bf),"bits output");bf.close();
 std::ofstream out(argv[3]);out<<"{\"N\":"<<N<<",\"K\":"<<K<<",\"strict_controls\":"<<controls<<",\"saved_prefix_N\":"<<saved<<",\"max_jump_count\":"<<maxloops<<",\"max_jump_endpoint\":\""<<dec(maxstate)<<"\",\"rows\":[\n";
 if(argc==6){
  // For every [4096*j,4096*(j+1)), save sum t and sum (n-4096*j)t.
  // n=0 is forced to zero. The final n=N is not in these half-open blocks.
  const U B=4096;const U masks[6]={0xaaaaaaaaaaaaaaaaULL,0xccccccccccccccccULL,0xf0f0f0f0f0f0f0f0ULL,0xff00ff00ff00ff00ULL,0xffff0000ffff0000ULL,0xffffffff00000000ULL};
  std::ofstream mo(argv[5],std::ios::binary);
  for(U lo=0;lo<N;lo+=B){I sm=B,wm=B*(B-1)/2;
   for(U j=0;j<B/64;j++){U v=bits[lo/64+j],cnt=__builtin_popcountll(v),w=0;for(int b=0;b<6;b++)w+=(1ULL<<b)*__builtin_popcountll(v&masks[b]);sm-=2*cnt;wm-=2*(64*j*cnt+w);}
   if(lo==0)--sm;I record[2]={sm,wm};mo.write(reinterpret_cast<char*>(record),sizeof(record));
  }check(bool(mo),"moments output");
 }
 I H[3]={0,0,0},s=0;U mx[3]={0,0,0},arg[3]={0,0,0},ratio_h[3]={0,0,0},ratio_n[3]={1,1,1};W square[3]={0,0,0};int octave=0;bool comma=false;
 // Rows cover (2^(k-1),2^k], except row 0={1}. All prefixes are inspected.
 for(U n=1;n<=N;n++){I t=bit(n)?-1:1;s+=t;for(int r=0;r<3;r++){H[r]+=((n%3==U(r))?2:-1)*t;U h=H[r]<0?-H[r]:H[r];square[r]+=W(h)*h;if(h>mx[r]){mx[r]=h;arg[r]=n;}if(W(h)*h*ratio_n[r]>W(ratio_h[r])*ratio_h[r]*n){ratio_h[r]=h;ratio_n[r]=n;}}
  if((n&(n-1))==0){check(H[0]+H[1]+H[2]==0,"zero sum residues");if(comma)out<<",\n";comma=true;out<<"{\"k\":"<<octave++<<",\"N\":"<<n<<",\"S\":"<<s<<",\"residues\":[";
   for(int r=0;r<3;r++){if(r)out<<',';out<<"["<<H[r]<<','<<mx[r]<<','<<arg[r]<<",\""<<dec(square[r])<<"\","<<ratio_h[r]<<','<<ratio_n[r]<<"]";mx[r]=arg[r]=0;square[r]=0;ratio_h[r]=0;ratio_n[r]=1;}out<<"]}";
  }
 }out<<"]}\n";check(bool(out),"json output");std::cerr<<"complete seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count()<<'\n';
 }catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}}
