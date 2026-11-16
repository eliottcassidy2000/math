/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
int main(int argc,char**argv){ m=atoi(argv[1]); long tot=0, free_=0, onlyAB=0, onlyD=0; ll best=1LL<<60; unsigned bu=0,bw=0; int bo=0;
  for(unsigned u=0;u<(1u<<m);u++){ ll cA=one_apex(u); if(!cA) continue;
    for(unsigned w=0;w<(1u<<m);w++){ ll cB=one_apex(w); if(!cB) continue;
      for(int o=0;o<2;o++){ tot++; ll ab=both_apex(u,w,o), d=disjoint_pairs(u,w);
        if(!ab&&!d) free_++; else if(!d) onlyAB++; else if(!ab) onlyD++;
        /* interaction excess over the glued value: H - (HA + HB - 1) = 2 cAB + 4 dAB */
        ll ex=2*ab+4*d; if(cA>=5&&cB>=5&&ex<best){ best=ex; bu=u; bw=w; bo=o; } } } }
  printf("spine %d: %ld tapes with cA, cB >= 1: glued (cAB=dAB=0) %ld; only cAB>0 %ld; only dAB>0 %ld; least excess 2cAB+4dAB with cA,cB>=5: %lld (u=%u w=%u o=%d, HA=%lld HB=%lld)\n",
         m,tot,free_,onlyAB,onlyD,best,bu,bw,bo,1+2*one_apex(bu),1+2*one_apex(bw)); return 0; }
