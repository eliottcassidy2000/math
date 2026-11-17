/* Z7-symmetric tournaments on Z7 x {0,1} (two free orbits): orbit k has circulant connection set S_k (S ∪ -S = Z7*, |S| = 3),
   cross arcs: (i,0) -> (j,1) iff j - i in C (C any subset of Z7), else (j,1) -> (i,0).  H is divisible by 7 (free action);
   report the orbit counts H/7 that are prime: each gives a Z7-symmetric realizer of the atom 7p. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 dp[1<<14][14];
static u64 HK(int N,const unsigned*out){ memset(dp,0,sizeof(u64)*(1u<<N)*14); for(int v=0;v<N;v++) dp[1u<<v][v]=1;
  for(unsigned S=1;S<(1u<<N);S++) for(int v=0;v<N;v++){ u64 c=dp[S][v]; if(!c) continue; unsigned o=out[v]&~S; while(o){ int w=__builtin_ctz(o); o&=o-1; dp[S|1u<<w][w]+=c; } }
  u64 h=0; for(int v=0;v<N;v++) h+=dp[(1u<<N)-1][v]; return h; }
static int isprime(u64 q){ if(q<2) return 0; for(u64 d=2;d*d<=q;d++) if(q%d==0) return 0; return 1; }
int main(void){ int halves[8], nh=0; for(int m=0;m<128;m++){ if(m&1) continue; int ok=1, c=0; for(int d=1;d<7;d++){ int a=m>>d&1, b=m>>(7-d)&1; if(a==b) ok=0; c+=a; } if(ok&&c==3) halves[nh++]=m; }
  printf("circulant halves: %d\n",nh); int found=0, total=0; static unsigned char seenp[1<<24]; int nprime=0; u64 minH=~0ULL,maxH=0;
  for(int a=0;a<nh;a++) for(int b=0;b<nh;b++) for(int C=0;C<128;C++){ unsigned out[14]={0};
      for(int i=0;i<7;i++) for(int j=0;j<7;j++){ if(i==j) continue; int d=(j-i+7)%7; if(halves[a]>>d&1) out[i]|=1u<<j; if(halves[b]>>d&1) out[7+i]|=1u<<(7+j); }
      for(int i=0;i<7;i++) for(int j=0;j<7;j++){ int d=(j-i+7)%7; if(C>>d&1) out[i]|=1u<<(7+j); else out[7+j]|=1u<<i; }
      u64 h=HK(14,out); total++; if(h%7){ printf("NOT divisible by 7: %llu\n",h); return 1; }
      if(h<minH) minH=h; if(h>maxH) maxH=h;
      u64 q=h/7; if(isprime(q) && q<(1u<<24) && !seenp[q]){ seenp[q]=1; nprime++; if(nprime<=40) printf("7p atom %llu = 7*%llu (orbit count prime), S0=%d S1=%d C=%d\n",h,q,halves[a],halves[b],C); } }
  printf("configurations %d; H range [%llu, %llu]; all divisible by 7; distinct prime orbit counts: %d\n",total,minH,maxH,nprime); return 0; }
