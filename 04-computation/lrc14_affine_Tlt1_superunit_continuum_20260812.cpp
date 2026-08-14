// Exact continuum census for the T<1 affine-resonance bank.
// Reads lines L j e f on stdin.  All arithmetic governing comparisons is
// integral; __int128 is used for cross-products.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
using i64=long long; using i128=__int128_t;
struct Ctx{i64 L,j,e,f;};

static i64 pos2(i64 x){return x>0?x*x:0;}
static i64 cdf2(i64 s,i64 T,i64 U){
  return pos2(s)-pos2(s-T)-pos2(s-U)+pos2(s-T-U);
}
static i64 mod(i64 x,i64 H){x%=H;return x<0?x+H:x;}
static i64 area_num(i64 y,i64 T,i64 U,i64 H){
  y=mod(y,H); i64 ans=0, h=H/14;
  for(i64 k=-3;k<=4;k++)
    ans+=cdf2(k*H+h-y,T,U)-cdf2(k*H-h-y,T,U);
  return ans;
}
static i64 overlap_num(i64 y,i64 U,i64 H){
  y=mod(y,H); i64 ans=0,h=H/14;
  for(i64 k=-3;k<=4;k++){
    i64 lo=std::max<i64>(0,k*H-h-y), hi=std::min(U,k*H+h-y);
    if(hi>lo) ans+=hi-lo;
  }
  return ans;
}
int main(){
  std::vector<Ctx> cs;Ctx x;while(std::cin>>x.L>>x.j>>x.e>>x.f)cs.push_back(x);
  if(cs.size()!=2530){std::cerr<<"bad context count "<<cs.size()<<"\n";return 2;}
  i64 checks=0,zero=0,equalities=0,boundary=0,boundary_bad=0,superunit=0,triples=0;
  i64 bn=1,bd=1; bool have=false; std::tuple<int,int,int,Ctx,i64> bw{};
  for(int d=1;d<=8;d++) for(int a=0;a<=7*d;a++){
    int D=d+a,g=std::gcd(a,d);
    for(int c=-13;c<=13;c++){
      if(c==0 || c%g || (a==0&&c<0) || (a==7*d&&c>0))continue;
      triples++;
      for(auto C:cs){
        i64 A=C.L*c+C.e*D-d*C.f;
        // Every finite T<1 point with p>=679 has
        // T_infinity=|A|/(Ld)<p/(p-r)<=679/672.  This deliberately
        // overcounts triples which have no physical residue lane.
        if((i128)std::llabs(A)*672>(i128)C.L*d*679)continue;
        i64 R=C.e*C.j%C.L,S=C.f*C.j%C.L,H=14*C.L*d;
        i64 X=14*(R*D-S*d),U=2*C.L*D,num,den;
        if(A==0){
          num=0;for(int s=0;s<d;s++)num+=overlap_num(X+14*C.L*a*s-U/2,U,H);
          den=(i64)D*H;zero++;
        }else{
          i64 T=14*std::llabs(A),start=X-(A<0?T:0);num=0;
          for(int s=0;s<d;s++)num+=area_num(start+14*C.L*a*s-U/2,T,U,H);
          den=392*C.L*d*(i64)D*std::llabs(A);
        }
        checks++;
        if(std::llabs(A)>C.L*d)superunit++;
        if(std::llabs(A)==C.L*d){
          boundary++;
          if((i128)num*49!=(i128)den)boundary_bad++;
        }
        if((i128)num*48048 < (i128)709*den){
          std::cerr<<"FAIL "<<d<<' '<<a<<' '<<c<<' '<<C.L<<' '<<C.j<<' '<<C.e<<' '<<C.f<<' '<<A<<' '<<num<<'/'<<den<<"\n";return 1;
        }
        if((i128)num*48048==(i128)709*den)equalities++;
        if(!have || (i128)num*bd<(i128)bn*den){have=true;bn=num;bd=den;bw={d,a,c,C,A};}
      }
    }
  }
  i64 gg=std::gcd(bn,bd);bn/=gg;bd/=gg;
  auto [d,a,c,C,A]=bw;
  std::cout<<"contexts="<<cs.size()<<";compatible_triples="<<triples<<";Tinf_le_679_over_672_checks="<<checks<<";Tinf_gt_1_checks="<<superunit<<";A0_checks="<<zero
           <<";Tinf_eq_1_checks="<<boundary<<";Tinf_eq_1_nonmean="<<boundary_bad<<"\n";
  std::cout<<"minimum="<<bn<<'/'<<bd<<";equalities="<<equalities<<";witness="
           <<d<<','<<a<<','<<c<<','<<C.L<<','<<C.j<<','<<C.e<<','<<C.f<<",A="<<A<<"\n";
  std::cout<<"status=FINITE-EXACT centered continuum floor 709/48048 on every compatible nonzero-c |c|<=13, T_infinity<=679/672 context; equality face exactly 1/49\n";
}
