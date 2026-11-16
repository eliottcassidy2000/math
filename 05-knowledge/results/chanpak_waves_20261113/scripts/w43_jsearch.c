/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
/* (J) inside two-apex tapes: H = 1 + 2(cA + cB + cAB) + 4 dAB.  (J) gluing <=> cAB = 0 and dAB = 0, so that
   H = HA + HB - 1 with HA = 1 + 2 cA, HB = 1 + 2 cB (each apex's own one-apex count on the same spine).
   Enumerate all (u, w, o) on spines of length m; record which (HA, HB) pairs with HA, HB >= 11 occur glued. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
static int isprime(ll n){ if(n<2) return 0; for(ll d=2;d*d<=n;d++) if(n%d==0) return 0; return 1; }
int main(int argc,char**argv){ m=atoi(argv[1]); int LIM=4000; static char seen[4000][4000]; long glued=0, total=0;
  for(unsigned u=0;u<(1u<<m);u++){ ll cA=one_apex(u); if(1+2*cA<11||1+2*cA>=LIM) continue;
    for(unsigned w=0;w<(1u<<m);w++){ ll cB=one_apex(w); if(1+2*cB<11||1+2*cB>=LIM) continue;
      for(int o=0;o<2;o++){ total++; if(both_apex(u,w,o)==0 && disjoint_pairs(u,w)==0){ glued++; seen[1+2*cA][1+2*cB]=1; } } } }
  long pairs=0, pp=0; int shown=0; for(int a=11;a<LIM;a+=2) for(int b=11;b<LIM;b+=2) if(seen[a][b]){ pairs++; if(isprime(a)&&isprime(b)){ pp++; if(shown<12){ shown++; printf("  glued prime pair p=%d q=%d -> H = p+q-1 = %d\n",a,b,a+b-1);} } }
  printf("spine %d: %ld (u,w,o) with HA,HB in [11,%d); %ld glued (cAB = 0, dAB = 0); %ld distinct (HA,HB) pairs glued, %ld of them prime pairs\n",m,total,LIM,glued,pairs,pp);
  return 0; }
