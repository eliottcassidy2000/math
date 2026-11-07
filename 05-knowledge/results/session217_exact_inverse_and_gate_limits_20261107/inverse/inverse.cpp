// Exact complete three-core inverse. Domain: integer 1<=N<=10^15.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
using I=long long; using W=__int128_t;
I ff(I p){I out=0;while(p){I s=p/2;if(p&1)out+=s+(s&1);p=s;}return out;}
I rev(I p,int b){I z=0;for(int i=0;i<b;i++){z=2*z+(p&1);p/=2;}return z;}
I floorq(I x,I y){I v=x/y;return v-(x%y<0);}
I ceilq(I x,I y){return -floorq(-x,y);}
struct Node{I k,s,c,l;};
struct Solver{
 std::vector<std::array<I,7>> cores;I tests=0,branches=0; int ci,ansb; I ansp,ansr;
 bool test(int c,I P,I R,int b,I N){
  tests++;auto [h,q,A,G,T,S,E]=cores[c];
  if(P<0||R<0||R>=(1LL<<b)||(P==0&&A==0))return false;
  I a0=q*P+A;if((W)a0*(1LL<<b)+h>N)return false;
  W val=(W)a0*((2LL<<b)-rev(R,b))+h*(2*R+1-ff(R))-P*G-h*ff(P)-T-(P&1)*S-E*(b?R>>(b-1):0);
  if(val!=N)return false;ci=c;ansp=P;ansr=R;ansb=b;return true;
 }
 bool solve(I N){
 tests=branches=0;
 for(int b=0;(1LL<<b)<=N/5;b++){
  I M=1LL<<b;
  for(int c=0;c<(int)cores.size();c++){
   auto [h,q,A,G,T,S,E]=cores[c];
   if((W)M*M<N){ // short suffix: enumerate R, invert P
    int t=0;W den=(W)q*q*M*M;while(den<=2*(W)h*N){den*=4;t++;}
    for(I R=0;R<M;R++){
     I L=2*M-rev(R,b),K=q*L-G,C=A*L+h*(2*R+1-ff(R))-T-E*(b?R>>(b-1):0);
     std::function<bool(Node,int)> visit=[&](Node v,int i){
      if(i<t){
       if(visit({2*v.k,0,v.c,v.l},i+1))return true;
       return visit({2*v.k-h,h,v.c+v.k-v.s,v.l+(1LL<<i)},i+1);
      }
      branches++; if(test(c,v.l,R,b,N))return true;
      I x=N-v.c,lo=std::max<I>(1,ceilq(x,v.k)),hi=floorq(x-h+v.s,v.k-h);
      if(hi-lo>0){std::cerr<<"prefix uniqueness failure\n";exit(2);}
      return lo<=hi&&test(c,v.l+(lo<<t),R,b,N);
     };
     if(visit({K,S,C,0},0))return true;
    }
   }else{ // long suffix: enumerate P, invert reversed R
    I pmax=floorq((N-h)/M-A,q);
    for(I P=0;P<=pmax;P++){
     if(P==0&&A==0)continue;
     I a0=q*P+A,tau=P*G+h*ff(P)+T+(P&1)*S,D=2*M*a0+h-tau,X=D-N;
     int t=std::min(2,b);W factor=(W)a0*(1LL<<(2*t));while(t<b&&factor<24*(W)h*M){t++;factor*=4;}
     std::function<bool(Node,int)> visit=[&](Node v,int i){
      if(i<t){
       if(visit({2*v.k,0,v.c,v.l},i+1))return true;
       return visit({2*v.k+h,h,v.c+v.k+v.s-2*h*(1LL<<(b-i-1)),v.l+(1LL<<i)},i+1);
      }
      branches++;int m=b-t;I x=X-v.c,lo=std::max<I>(0,ceilq(x-h,v.k+h)),hi=std::min<I>((1LL<<m)-1,floorq(x+2*h*((1LL<<m)-1),v.k));
      if(hi-lo>0){std::cerr<<"suffix uniqueness failure\n";exit(2);}
      if(lo<=hi){I z=v.l+(lo<<t);return test(c,P,rev(z,b),b,N);}
      return false;
     };
     if(visit({a0,E,0,0},0))return true;
    }
   }
  }
 }
 return false;
 }
};
int main(int argc,char**argv){if(argc<3){std::cerr<<"usage: inverse CORE_STATES.txt N...\n";return 2;}Solver s;std::ifstream in(argv[1]);std::array<I,7>a;while(in>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5]>>a[6])s.cores.push_back(a);if(s.cores.size()!=80)return 2;
 for(int i=2;i<argc;i++){I N=std::stoll(argv[i]);if(N<1||N>1000000000000000LL)return 2;bool hit=s.solve(N);std::cout<<"{\"N\":"<<N<<",\"status\":\""<<(hit?"FOUND":"ABSENT_FROM_GRAMMAR")<<"\",\"branches\":"<<s.branches<<",\"tests\":"<<s.tests;if(hit)std::cout<<",\"core\":"<<s.ci<<",\"P\":"<<s.ansp<<",\"R\":"<<s.ansr<<",\"b\":"<<s.ansb;std::cout<<"}\n";}
}
