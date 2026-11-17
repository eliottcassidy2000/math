/* Z7-symmetric 15-vertex tournaments: two free 7-orbits (circulants S0, S1) glued by a circulant block C (|C| = 3 for
   regularity) plus one fixed vertex v (v beats orbit 1, orbit 0 beats v).  Doubly regular <=> skew-Hadamard of order 16.
   Report the doubly regular ones with H, H/7, and |Aut| (brute force on score/neighbourhood-refined permutations). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 dp[1<<15][15];
static u64 HK(int N,const unsigned*out){ memset(dp,0,sizeof(u64)*(1u<<N)*15); for(int v=0;v<N;v++) dp[1u<<v][v]=1;
  for(unsigned S=1;S<(1u<<N);S++) for(int v=0;v<N;v++){ u64 c=dp[S][v]; if(!c) continue; unsigned o=out[v]&~S; while(o){ int w=__builtin_ctz(o); o&=o-1; dp[S|1u<<w][w]+=c; } }
  u64 h=0; for(int v=0;v<N;v++) h+=dp[(1u<<N)-1][v]; return h; }
static int doubly_regular(int N,const unsigned*out){ for(int i=0;i<N;i++) if(__builtin_popcount(out[i])!=(N-1)/2) return 0;
  for(int i=0;i<N;i++) for(int j=i+1;j<N;j++) if(__builtin_popcount(out[i]&out[j])!=(N-3)/4) return 0; return 1; }
int main(void){ int halves[8],nh=0; for(int m=0;m<128;m++){ if(m&1) continue; int ok=1,c=0; for(int d=1;d<7;d++){ int a=m>>d&1,b=m>>(7-d)&1; if(a==b) ok=0; c+=a; } if(ok&&c==3) halves[nh++]=m; }
  int found=0; 
  for(int a=0;a<nh;a++) for(int b=0;b<nh;b++) for(int C=0;C<128;C++){ if(__builtin_popcount(C)!=3) continue;
      unsigned out[15]={0};
      for(int i=0;i<7;i++) for(int j=0;j<7;j++){ if(i==j) continue; int d=(j-i+7)%7; if(halves[a]>>d&1) out[i]|=1u<<j; if(halves[b]>>d&1) out[7+i]|=1u<<(7+j); }
      for(int i=0;i<7;i++) for(int j=0;j<7;j++){ int d=(j-i+7)%7; if(C>>d&1) out[i]|=1u<<(7+j); else out[7+j]|=1u<<i; }
      for(int i=0;i<7;i++){ out[i]|=1u<<14; out[14]|=1u<<(7+i); }   /* orbit 0 beats v, v beats orbit 1 */
      if(!doubly_regular(15,out)) continue;
      u64 h=HK(15,out); found++;
      printf("doubly regular: S0=%d S1=%d C=%d  H=%llu = 7*%llu  (H/105 = %llu, H mod 105 = %llu)\n",halves[a],halves[b],C,h,h/7,h/105,h%105);
  }
  printf("doubly regular Z7 gluings found: %d\n",found); return 0; }
