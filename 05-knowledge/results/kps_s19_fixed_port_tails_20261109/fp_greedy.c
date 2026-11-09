/* Q4: greedy set cover - few core states whose fixed-port tails reach every target in [LO, LIM]
 * (odd targets for x = 2; targets = 1 mod 3 for x = 3), excluding the permanent gaps.
 * Usage: fp_greedy STATSFILE X LO LIM MAXORDER
 * Output: chosen core states (h0 R0 Delta0 hu) with the number of new targets each covers.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static ll X, LO, LIM;
static uint64_t **sets; static ll nc; static ll (*cores)[4];
static ll W; /* words per bitset over targets */
static ll tid(ll v){ /* target index or -1 */
  if (v < LO || v > LIM) return -1;
  if (X == 2) { if (v % 2 == 0) return -1; return (v - LO) / 2; }
  if (v % X != 1) return -1; return (v - LO) / X;
}
typedef struct { ll h, R, D; } st;
static uint64_t *tab; static uint64_t tabsize, tabcount;
static uint64_t mix(uint64_t x){ x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void tab_reset(void){ if (!tab){ tabsize = 1u<<20; tab = calloc(tabsize, 16);} else memset(tab, 0, tabsize*16); tabcount = 0; }
static int tab_insert(uint64_t a, uint64_t b){
  if (2*(tabcount+1) > tabsize){ uint64_t old = tabsize, *ot = tab; tabsize *= 2; tab = calloc(tabsize, 16); tabcount = 0;
    for (uint64_t i=0;i<old;i++) if (ot[2*i]) { uint64_t j = mix(ot[2*i]^(ot[2*i+1]*0x9E3779B97F4A7C15ULL)) & (tabsize-1);
      while (tab[2*j]) j = (j+1)&(tabsize-1); tab[2*j]=ot[2*i]; tab[2*j+1]=ot[2*i+1]; tabcount++; }
    free(ot); }
  uint64_t j = mix(a^(b*0x9E3779B97F4A7C15ULL)) & (tabsize-1);
  while (tab[2*j]) { if (tab[2*j]==a && tab[2*j+1]==b) return 0; j = (j+1)&(tabsize-1); }
  tab[2*j]=a; tab[2*j+1]=b; tabcount++; return 1; }
static st *stk; static size_t cap, len;
static void push(ll h, ll R, ll D){ if (len == cap){ cap = cap ? 2*cap : 1<<16; stk = realloc(stk, cap*sizeof(st)); } stk[len].h=h; stk[len].R=R; stk[len].D=D; len++; }
int main(int argc, char **argv){
  FILE *f = fopen(argv[1], "r"); X = atoll(argv[2]); LO = atoll(argv[3]); LIM = atoll(argv[4]); int maxo = atoi(argv[5]);
  if (X == 2 && LO % 2 == 0) LO++;
  while (X == 3 && LO % 3 != 1) LO++;
  ll ntar = (X == 2) ? (LIM - LO) / 2 + 1 : (LIM - LO) / X + 1;
  W = (ntar + 63) / 64;
  cores = malloc(sizeof(*cores) * 100000); char line[512];
  while (fgets(line, sizeof line, f)){
    ll a[11]; if (sscanf(line, "%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld", &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7],&a[8],&a[9],&a[10]) < 11) continue;
    if (a[0] > maxo) continue;
    ll *r = a + 1; ll h, E, O, hu;
    if (X == 2){ h=r[2]; E=r[3]; O=r[4]; hu=r[5]; } else { h=r[6]; E=r[7]; O=r[8]; hu=r[9]; }
    ll c[4] = { h, E+O, E-O, hu }; int dup = 0;
    for (ll i=0;i<nc;i++) if (!memcmp(cores[i], c, sizeof c)) { dup = 1; break; }
    if (!dup){ memcpy(cores[nc], c, sizeof c); nc++; }
  }
  sets = malloc(sizeof(uint64_t*) * nc);
  for (ll ci=0; ci<nc; ci++){
    sets[ci] = calloc(W, 8);
    ll h0=cores[ci][0], R0=cores[ci][1], D0=cores[ci][2], hu=cores[ci][3];
    if (h0 > LIM) continue;
    len = 0; push(h0, R0, D0); tab_reset();
    while (len){
      st s = stk[--len];
      if (!tab_insert((uint64_t)s.h*4u + (uint64_t)(s.D==0?0:(s.D==-hu?1:2)) + 1, (uint64_t)s.R)) continue;
      ll t = tid(s.h); if (t >= 0) sets[ci][t>>6] |= 1ULL << (t & 63);
      ll off = X * (s.R + s.D) / 2;
      if (s.h + off > LIM) continue;
      push(s.h + off, 2*s.R - hu, -hu);           /* B */
      if (s.h + X * s.R <= LIM) push(s.h, 2*s.R, 0); /* S, only if a later B can fit */
    }
  }
  uint64_t *unc = malloc(W * 8);
  for (ll i=0;i<W;i++) unc[i] = ~0ULL;
  if (ntar % 64) unc[W-1] = (1ULL << (ntar % 64)) - 1;
  /* remove permanent gaps from the target set */
  ll gaps2[2] = {7, 21}, gaps3[2] = {10, 31};
  for (int g=0; g<2; g++){ ll t = tid(X==2 ? gaps2[g] : gaps3[g]); if (t >= 0) unc[t>>6] &= ~(1ULL << (t&63)); }
  ll remaining = 0; for (ll i=0;i<W;i++) remaining += __builtin_popcountll(unc[i]);
  printf("targets %lld in [%lld, %lld], cores %lld\n", remaining, LO, LIM, nc);
  /* lower bound: cores that are the unique cover of some target are in every cover */
  { ll forced = 0; unsigned char *isf = calloc(nc,1);
    for (ll t=0;t<ntar;t++){ if (!(unc[t>>6]>>(t&63)&1)) continue; ll cnt=0, who=-1;
      for (ll ci=0; ci<nc && cnt<2; ci++) if (sets[ci][t>>6]>>(t&63)&1){ cnt++; who=ci; }
      if (cnt==1 && !isf[who]){ isf[who]=1; forced++; }
      if (cnt==0) printf("target %lld uncovered by every core\n", X==2? LO+2*t : LO+X*t); }
    printf("forced cores (unique cover of some target): %lld\n", forced); }
  /* better lower bound: a packing of targets with pairwise disjoint cover sets (each needs its own core) */
  { ll *cnt = calloc(ntar, sizeof(ll)); ll *ord = malloc(ntar * sizeof(ll)); ll m = 0;
    for (ll t=0;t<ntar;t++){ if (!(unc[t>>6]>>(t&63)&1)) continue;
      for (ll ci=0; ci<nc; ci++) if (sets[ci][t>>6]>>(t&63)&1) cnt[t]++;
      ord[m++] = t; }
    /* sort by cover-set size ascending (simple insertion into buckets) */
    ll maxc = 0; for (ll i=0;i<m;i++) if (cnt[ord[i]]>maxc) maxc=cnt[ord[i]];
    unsigned char *used = calloc(nc,1); ll packed = 0;
    for (ll c=1;c<=maxc;c++) for (ll i=0;i<m;i++){ ll t=ord[i]; if (cnt[t]!=c) continue;
      int ok=1; for (ll ci=0; ci<nc; ci++) if ((sets[ci][t>>6]>>(t&63)&1) && used[ci]) { ok=0; break; }
      if (!ok) continue;
      for (ll ci=0; ci<nc; ci++) if (sets[ci][t>>6]>>(t&63)&1) used[ci]=1;
      packed++; }
    printf("lower bound (targets with pairwise disjoint cover sets): %lld\n", packed); }
  int round = 0;
  while (remaining > 0 && round < 500){
    ll best = -1, bc = 0;
    for (ll ci=0; ci<nc; ci++){ ll c = 0; for (ll i=0;i<W;i++) c += __builtin_popcountll(sets[ci][i] & unc[i]); if (c > bc){ bc = c; best = ci; } }
    if (best < 0) break;
    for (ll i=0;i<W;i++) unc[i] &= ~sets[best][i];
    remaining -= bc; round++;
    printf("%d core (h0=%lld R0=%lld D0=%lld hu=%lld) +%lld left %lld\n", round, cores[best][0], cores[best][1], cores[best][2], cores[best][3], bc, remaining);
  }
  if (remaining > 0){ printf("uncovered:"); ll shown=0; for (ll t=0;t<ntar && shown<50;t++) if (unc[t>>6]>>(t&63)&1){ printf(" %lld", X==2? LO+2*t : LO+X*t); shown++; } printf("\n"); }
  return 0;
}
