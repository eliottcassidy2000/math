/* Wave 57 recomputation of known families (credit THM-4876 (6) one centre; Two-Hub Conjecture, wave 5).
   Usage: ./w57_one_centre MAXLEN LIMIT ; ./w57_two_centre MAXLEN LIMIT */
/* Values of H(TT_n + x_w) over all contact words w of length n <= L (Gray-free brute force).
   value(w) = [w1=1] + [wn=0] + #01 + sum_{b<a, w_b=1, w_a=0} 2^(a-b-1). Report odd N <= LIM missed. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc,char**argv){ int L=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10);
  unsigned char *hit=calloc(LIM/2+2,1);
  for(int n=1;n<=L;n++){ u64 W=1ULL<<n;
    for(u64 w=0;w<W;w++){ /* bit i (0-based) = w_{i+1} */
      u64 acc=0,tot=0; int asc=0; int prev=-1;
      for(int a=0;a<n;a++){ int b=(w>>a)&1; if(!b){ tot+=acc; if(tot>LIM) break; } if(prev==0&&b==1) asc++; prev=b; acc=2*acc+b; if(acc>4*LIM) acc=4*LIM; }
      if(tot>LIM) continue;
      u64 v=(u64)((w&1)==1)+(u64)(((w>>(n-1))&1)==0)+asc+tot;
      if(v<=LIM && (v&1)) hit[v/2]=1; }
    if(n>=L-6){ u64 miss=0, firstbig=0; for(u64 N=1;N<=LIM;N+=2) if(!hit[N/2]) miss++;
      printf("L<=%d: odd N <= %llu missed %llu:",n,LIM,miss); int c=0; for(u64 N=1;N<=LIM&&c<40;N+=2) if(!hit[N/2]){ printf(" %llu",N); c++; } printf("\n"); fflush(stdout);} }
  return 0; }
