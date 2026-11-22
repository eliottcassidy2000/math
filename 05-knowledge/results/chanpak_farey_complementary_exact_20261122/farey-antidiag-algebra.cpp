#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdint>
#include <cassert>
using namespace std; using I=__int128; using U=uint64_t;
const I N=(I(1)<<121)-(I(1)<<60)-1;
struct Chart{long long h,q,c,t,e,r,k;};
I rootfloor(I R, I q, I B){
 if(R<0)return -1;
 long double b=B, Q=q;
 I x=(I)((sqrtl(b*b+4*Q*(long double)R)-b)/(2*Q));
 x=max(I(0),x);
 while(q*x*x+B*x>R)--x;
 while(q*(x+1)*(x+1)+B*(x+1)<=R)++x;
 return x;
}
I f(I x){I s=0;while(x){I y=x/2;if(!(x&1))s+=y+(y&1);x=y;}return s;}
int md(I a,int p){int x=a%p;return x<0?x+p:x;}
int main(int argc,char**argv){
 assert(argc>=4);
 int ell=stoi(argv[1]),rat=stoi(argv[2]),want=stoi(argv[3]);
 U B0=U(1)<<ell,rr=argc>=5?stoull(argv[4]):B0-1;
 int n;Chart z;
 while(cin>>n>>z.h>>z.q>>z.c>>z.t>>z.e>>z.r>>z.k){
  if(n!=want)continue;
  auto [h,q,c,t,e,r,k]=z;
  int H=h;
  vector<int> primes;
  for(int d=2;d<=H;++d)if(H%d==0){H/=d;primes.push_back(d);}
  struct Cell{int p;vector<array<bool,2>> entry;};
  vector<Cell> sieve;
  for(int p:primes){
   Cell T{p,vector<array<bool,2>>(2*p)};
   for(int s=0;s<2*p;++s)for(int a=0;a<2*p;++a){
    int b=(1&a)^(1&s);
    I v=I(k)+I(c)*s+I(q)*a*(s-a)+I(t-c)*a+I(e)*b+I(r)*(a&1);
    if(md(v-N,p)==0)T.entry[s][a&1]=true;
   }
   sieve.push_back(move(T));
  }
  auto upper=[&](I s){return k+I(q)*(s*s/4)+(max(c,t)+h)*s+r+e;};
  I L=0,R=I(1)<<62;
  // Only compare sums near the actual 121-bit root, limiting above formula to safe widths.
  R=2*rootfloor((N-k),q,0)+3;
  while(L<R){I m=(L+R)/2;if(upper(m)>=N)R=m;else L=m+1;}
  // x:y lies in [1:rat,rat:1], hence xy >= rat*S^2/(rat+1)^2.
  I high=rootfloor((N-k)*I(rat+1)*(rat+1),I(q)*rat,0);
  U start=U(L)+(rr+B0-U(L)%B0)%B0;
  U totals=0,masks=0,cands=0;
  for(U S=start;S<=high;S+=B0){
   ++totals;
   int mask=3;
   for(const auto&T:sieve){
    const auto &a=T.entry[S%(2*T.p)];
    mask&=int(a[0])|(int(a[1])<<1);
    if(!mask)break;
   }
   if(!mask)continue;
   ++masks;
   I center=(I(q)*S+t-c+q)/(2*q);
   I beta=I(q)*S+t-c-2*I(q)*center;
   I D=k+I(q)*center*(S-center)+I(t)*center+I(c)*(S-center)-N;
   I E=I(h)*S+r+e;
   I xmin=(I(S)+rat)/(rat+1),xmax=S-xmin;
   for(int sign:{-1,1}){
    I off= -sign*beta;
    I a=rootfloor(D-1,q,off)+1,b=rootfloor(D+E,q,off);
    // For the positive side with a negative linear term, the nonnegative roots
    // above a negative lower threshold are not necessarily an initial interval.
    // Offset zero is handled explicitly and the positive indices 1,2,... are monotone.
    a=max(a,I(sign==1));
    if(sign<0){a=max(a,center-xmax);b=min(b,center-xmin);}
    else{a=max(a,xmin-center);b=min(b,xmax-center);}
    for(I j=a;j<=b;++j){
     I x=center+sign*j,y=S-x;
     if(!(mask>>(x&1)&1))continue;
     I res=-D+I(q)*j*j+off*j-r*(x&1)-e*(y&1);
     if(res<0||res%(2*h))continue;
     ++cands;
     if(res==I(h)*(f(x)+f(y)))cout<<"HIT "<<(U)x<<" "<<(U)y<<"\n";
    }
   }
  }
  cerr<<"spine "<<n<<" ell "<<ell<<" aspect "<<rat<<" residue "<<rr<<" sums "<<totals<<" crt-live "<<masks<<" digital "<<cands<<"\n";
 }
}
