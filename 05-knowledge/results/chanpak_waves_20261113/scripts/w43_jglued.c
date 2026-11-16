/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
int main(int argc,char**argv){ m=atoi(argv[1]); long hist[64][64]; memset(hist,0,sizeof hist);
  for(unsigned u=0;u<(1u<<m);u++){ ll cA=one_apex(u); if(!cA) continue; for(unsigned w=0;w<(1u<<m);w++){ ll cB=one_apex(w); if(!cB) continue;
    for(int o=0;o<2;o++) if(!both_apex(u,w,o)&&!disjoint_pairs(u,w)) { if(cA<64&&cB<64) hist[cA][cB]++; if(cA>=5||cB>=5) printf("glued with a large part: u=%u w=%u o=%d HA=%lld HB=%lld\n",u,w,o,1+2*cA,1+2*cB); } } }
  for(int a=0;a<64;a++) for(int b=0;b<64;b++) if(hist[a][b]) printf("spine %d: glued tapes with (HA,HB) = (%d,%d): %ld\n",m,1+2*a,1+2*b,hist[a][b]); return 0; }
