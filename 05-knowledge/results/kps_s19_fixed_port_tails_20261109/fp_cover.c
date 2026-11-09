/* Fixed-port tails (THM-5085): coverage of I(Omega(T), x) over all cores (T0, u).
 * Core statistics per line (from port_stats.c, optionally prefixed by the order):
 *   [n] code u h2 E2 O2 hu2 h3 E3 O3 hu3
 * State (h, F, D) with hu fixed; letters:
 *   S   : (h, 2F, 0)
 *   I_0 : h += x (F + D)/2 ;            F -> 2F - hu ;          D -> -hu
 *   I_t : h += x 2^(t-1) F (t >= 1) ;   F -> 2^(t+1) F - hu ;   D -> -hu
 * Every letter keeps h nondecreasing, so states with h > LIM are discarded.  A state whose cheapest bridge already
 * overflows is a leaf.
 * Usage: fp_cover STATSFILE X LIM TMAX HAS_ORDER_PREFIX(0/1) [MAXORDER]
 * Prints the misses (v = 1 mod X; X = 2: odd v) and, per core, how many values it alone reaches (to stderr summary).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static unsigned char *hit;
static ll X, LIM; static int TMAX;
/* per-core open-addressing set of states */
static uint64_t *tab; static uint64_t tabsize, tabcount;
static uint64_t mix(uint64_t x){ x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void tab_reset(void){ if (!tab){ tabsize = 1u<<20; tab = calloc(tabsize, 16);} else memset(tab, 0, tabsize*16); tabcount = 0; }
static int tab_insert(uint64_t a, uint64_t b){ /* key (a,b) with a != 0 */
  if (2*(tabcount+1) > tabsize){ /* grow */
    uint64_t old = tabsize, *ot = tab; tabsize *= 2; tab = calloc(tabsize, 16); tabcount = 0;
    for (uint64_t i=0;i<old;i++) if (ot[2*i]) { uint64_t j = mix(ot[2*i]^(ot[2*i+1]*0x9E3779B97F4A7C15ULL)) & (tabsize-1);
      while (tab[2*j]) j = (j+1)&(tabsize-1); tab[2*j]=ot[2*i]; tab[2*j+1]=ot[2*i+1]; tabcount++; }
    free(ot); }
  uint64_t j = mix(a^(b*0x9E3779B97F4A7C15ULL)) & (tabsize-1);
  while (tab[2*j]) { if (tab[2*j]==a && tab[2*j+1]==b) return 0; j = (j+1)&(tabsize-1); }
  tab[2*j]=a; tab[2*j+1]=b; tabcount++; return 1;
}
typedef struct { ll h, F, D; } st;
static st *stack; static size_t scap, slen;
static void push(ll h, ll F, ll D){ if (slen==scap){ scap = scap? 2*scap : 1<<20; stack = realloc(stack, scap*sizeof(st)); } stack[slen].h=h; stack[slen].F=F; stack[slen].D=D; slen++; }
int main(int argc, char **argv){
  FILE *f = fopen(argv[1], "r"); X = atoll(argv[2]); LIM = atoll(argv[3]); TMAX = atoi(argv[4]);
  int pref = atoi(argv[5]); int maxorder = argc > 6 ? atoi(argv[6]) : 99;
  hit = calloc(LIM+1, 1);
  char line[512]; ll ncores = 0, states = 0;
  /* dedupe core statistics */
  ll (*cores)[4] = malloc(sizeof(*cores) * 200000); ll nc = 0;
  while (fgets(line, sizeof line, f)){
    ll a[11]; int k = sscanf(line, "%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld", &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7],&a[8],&a[9],&a[10]);
    ll *r = pref ? a+1 : a; if (pref && a[0] > maxorder) continue; if (k < (pref?11:10)) continue;
    ll h, E, O, hu;
    if (X == 2){ h=r[2]; E=r[3]; O=r[4]; hu=r[5]; } else { h=r[6]; E=r[7]; O=r[8]; hu=r[9]; }
    ll c[4] = { h, E+O, E-O, hu };
    int dup = 0; for (ll i=0;i<nc;i++) if (!memcmp(cores[i], c, sizeof c)) { dup=1; break; }
    if (!dup){ memcpy(cores[nc], c, sizeof c); nc++; }
  }
  for (ll ci=0; ci<nc; ci++){
    ll h0=cores[ci][0], F0=cores[ci][1], D0=cores[ci][2], hu=cores[ci][3];
    if (h0 > LIM) continue; ncores++;
    tab_reset(); slen = 0; push(h0, F0, D0);
    while (slen){
      st s = stack[--slen];
      uint64_t key1 = (uint64_t)s.h*4u + (uint64_t)(s.D==0?0:(s.D==-hu?1:2)) + 1, key2 = (uint64_t)s.F;
      if (!tab_insert(key1, key2)) continue;
      states++; hit[s.h] = 1;
      /* cheapest bridge: I_0 adds x(F+D)/2 >= x(F-hu)/2 ; I_1 adds xF ; stop if all overflow */
      ll off0 = X*(s.F + s.D)/2;          /* I_0 is the cheapest letter (D <= F) */
      if (s.h + off0 > LIM) continue;     /* leaf: every continuation overflows */
      /* S: worthwhile only if a later I_0 (offset x*F after S) can still fit */
      if (s.h + X*s.F <= LIM) push(s.h, 2*s.F, 0);
      for (int t=0; t<=TMAX; t++){
        ll off = (t==0) ? off0 : X*((1LL<<(t-1))*s.F);
        if (s.h + off <= LIM) push(s.h + off, (1LL<<(t+1))*s.F - hu, -hu);
      }
    }
  }
  ll nmiss=0; printf("misses (x=%lld, LIM=%lld, cores=%lld, states=%lld):", X, LIM, ncores, states);
  for (ll v=1; v<=LIM; v++){ if (X==2 ? (v%2==0) : (v%X!=1)) continue; if (!hit[v]){ nmiss++; if (nmiss<=300) printf(" %lld", v);} }
  printf("\ncount %lld\n", nmiss);
  return 0;
}
