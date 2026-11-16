/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
/* triangle shield: B beats only spine vertex m-2 (w = all ones except bit m-2); A beats the last two vertices (u bits m-2, m-1 = 0).
   Claim: with A -> B (o = 1), H = 3 * HA exactly, HA = 1 + 2 one_apex(u); with B -> A it is 3 HA + 2 cAB. */
int main(int argc,char**argv){ m=atoi(argv[1]); long ok1=0,n=0; unsigned w=((1u<<m)-1)&~(1u<<(m-2));
  long ok0=0; ll ex0min=1LL<<60, ex0max=0;
  for(unsigned u=0;u<(1u<<(m-2));u++){ /* bits m-2, m-1 of u are 0 */
    ll HA=1+2*one_apex(u); ll h1=formula(u,w,1), h0=formula(u,w,0); n++; if(h1==3*HA) ok1++;
    ll ex=h0-3*HA; if(ex<ex0min) ex0min=ex; if(ex>ex0max) ex0max=ex; }
  printf("spine %d: A->B: H = 3*HA for %ld of %ld A-words; B->A: H - 3*HA ranges over [%lld, %lld] (= 2 cAB)\n",m,ok1,n,ex0min,ex0max); return 0; }
