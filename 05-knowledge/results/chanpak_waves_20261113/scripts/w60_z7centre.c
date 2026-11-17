/* One-centre quotient Q = TT_m + u (u -> spine vertex i iff w_i = 1), centre u blown up into B in {I7, P7}.
   a_r(w) = walks R_0 u R_1 ... u R_r through Q (u visited r times, spine once): R_i increasing spine runs,
   R_1..R_{r-1} nonempty, last(R_i) -> u (w = 0) for i < r, u -> first(R_i) (w = 1) for i >= 1.
   Orbit count o(w) = sum_r a_r c_r(B)/7 with c_r/7 = (25|27, 234, 954, 2160, 2880, 2160, 720). Check a_1 = one-centre H. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static const u64 g[8]={0,25,234,954,2160,2880,2160,720};
static int pw3[10];
static u64 ar(int m,const int*w,int r){ /* slots 0..r; state digit per slot: 0 empty, 1 open-ok(end ok), 2 open-bad */
  int S=pw3[r+1]; static u64 dp[6561], nd[6561]; memset(dp,0,sizeof(u64)*S); dp[0]=1;
  for(int i=0;i<m;i++){ memset(nd,0,sizeof(u64)*S);
    for(int st=0;st<S;st++){ u64 c=dp[st]; if(!c) continue;
      for(int k=0;k<=r;k++){ int dig=(st/pw3[k])%3;
        if(dig==0 && k>=1 && w[i]!=1) continue;                 /* start condition u -> first */
        int nd_dig = (k<r) ? ((w[i]==0)?1:2) : 1;               /* end condition last -> u for k < r */
        nd[st+(nd_dig-dig)*pw3[k]]+=c; } }
    memcpy(dp,nd,sizeof(u64)*S); }
  u64 tot=0; for(int st=0;st<S;st++){ u64 c=dp[st]; if(!c) continue; int ok=1;
    for(int k=0;k<=r&&ok;k++){ int dig=(st/pw3[k])%3; if(k==0){ if(dig==2) ok=0; } else if(k<r){ if(dig!=1) ok=0; } }
    if(ok) tot+=c; }
  return tot; }
static int isprime(u64 q){ if(q<2) return 0; if(q%2==0) return q==2; for(u64 d=3;d*d<=q;d+=2) if(q%d==0) return 0; return 1; }
static u64 onecentre(int m,const int*w){ u64 acc=0,tot=0; int n10=0; for(int a=0;a<m;a++){ if(w[a]==0) tot+=acc; acc=2*acc+w[a]; if(a+1<m&&w[a]==1&&w[a+1]==0) n10++; } return 1+n10+tot; }
int main(int argc,char**argv){ int M=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10); pw3[0]=1; for(int i=1;i<10;i++) pw3[i]=pw3[i-1]*3;
  unsigned char*hitI=calloc(LIM+1,1),*hitP=calloc(LIM+1,1); long bad=0, words=0;
  for(int m=1;m<=M;m++) for(u64 W=0;W<(1ULL<<m);W++){ int w[32]; for(int i=0;i<m;i++) w[i]=W>>i&1; words++;
      u64 a[8]={0}; for(int r=1;r<=7&&r<=m+1;r++) a[r]=ar(m,w,r);
      if(a[1]!=onecentre(m,w)) bad++;
      u64 oI=0; for(int r=1;r<=7;r++) oI+=a[r]*g[r]; u64 oP=oI+2*a[1];
      if(oI<=LIM) hitI[oI]=1; if(oP<=LIM) hitP[oP]=1; }
  long pI=0,pP=0,pAny=0,primes=0; u64 firstmiss=0;
  for(u64 q=3;q<=LIM;q+=2){ if(!isprime(q)) continue; primes++; if(hitI[q]) pI++; if(hitP[q]) pP++; if(hitI[q]||hitP[q]) pAny++; }
  long odd=0, oddhit=0; for(u64 q=1;q<=LIM;q+=2){ odd++; if(hitI[q]||hitP[q]) oddhit++; }
  printf("one-centre words to length %d (%ld words); a_1 = one-centre count check failures: %ld\n",M,words,bad);
  printf("orbit counts <= %llu: odd values hit %ld of %ld; primes hit: I7 %ld, P7 %ld, either %ld of %ld\n",LIM,oddhit,odd,pI,pP,pAny,primes);
  printf("smallest prime orbit counts:"); int k=0; for(u64 q=3;q<=LIM&&k<20;q+=2) if(isprime(q)&&(hitI[q]||hitP[q])){ printf(" %llu",q); k++; } printf("\n");
  return 0; }
