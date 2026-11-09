/* aud_fewcycle.c - independent audit enumerator for THM-4940 (audit of opus-2026-11-08-S86).
 * Written from scratch for the audit; it does not reuse the session's code.
 *
 * Enumerates, up to isomorphism, every STRONG tournament with at most M directed odd cycles (alpha_1 <= M),
 * or (mode 1) every STRONG INTERCYCLIC tournament (no two vertex-disjoint cycles) with alpha_1 <= M.
 * Generation: level by level, one-vertex extensions of the previous level.  Completeness rests only on
 *   (a) Moon's deletion lemma: a strong tournament of order n >= 4 has a vertex v with T - v strong;
 *   (b) alpha_1(T - v) <= alpha_1(T) (cycles of T - v are cycles of T); intercyclicity is hereditary.
 * The program does NOT use Moon's cycle-count floor as a stopping rule: it runs until a level is EMPTY, which by
 * (a)+(b) proves that no larger order contributes.
 *
 * Isomorphism: own canonical form = individualisation-refinement with initial colours (score, triangles through v),
 * colour refinement by out-counts, branching on the SMALLEST non-singleton cell, and the lexicographically MAXIMAL
 * leaf code (the session and classgen.c branch on the first cell and keep the minimal code).
 * Odd cycles through the new vertex: depth-first path enumeration with backward-reachability pruning (each explored
 * prefix closes to at least one cycle), capped.  Packing vector: explicit list of odd cycles (vertex masks with
 * multiplicities), then k-sets of pairwise disjoint supports by recursion.  For n <= HCHK the Hamiltonian-path count
 * is computed by subset DP and compared with OCF (H = sum alpha_k 2^k).
 *
 * usage: aud_fewcycle M MODE HCHK > out     MODE 0 = strong, 1 = strong intercyclic
 * line:  n a1 a2 a3 a4 a5 a6 a7 H_ocf H_dp(-1 if not computed) selfconverse code_hex
 * further modes (checks):
 *   aud_fewcycle 0 2 n cls_n.txt   census mode: recompute every class of a classgen list (strong, odd cycles by length,
 *                                  packing vector, OCF against classgen's H, own canonical code, self-converse)
 *   aud_fewcycle 0 3 R inventory   invariance of the canonical form under R random relabellings of every class
 *   aud_fewcycle 0 4 N dir         same-H move graph (H-preserving deletions and single-arc reversals) on the census
 *                                  files dir/census_<n>.txt, n <= N: classes cut off from first creators, H = 9 fibre
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 24
#define CW 5
typedef struct { u64 w[CW]; } code_t;

/* ------------------------------------------------------------------ canonical form */
static int CN;
static uint32_t CADJ[MAXN];
static code_t CBEST; static int CHAVE;

static int code_cmp(const code_t *a, const code_t *b){
  for (int i = 0; i < CW; i++) if (a->w[i] != b->w[i]) return a->w[i] < b->w[i] ? -1 : 1;
  return 0;
}
static code_t make_code(const int *ord){
  code_t c; memset(&c, 0, sizeof c); int b = 0;
  for (int i = 0; i < CN; i++) for (int j = i + 1; j < CN; j++){
    if (CADJ[ord[i]] >> ord[j] & 1) c.w[b >> 6] |= 1ULL << (b & 63);
    b++;
  }
  return c;
}
/* sort helper: vertices by key arrays (lexicographic over len ints) */
static int SIGLEN; static int SIG[MAXN][MAXN + 2];
static int sigcmp(int a, int b){ for (int i = 0; i < SIGLEN; i++) if (SIG[a][i] != SIG[b][i]) return SIG[a][i] < SIG[b][i] ? -1 : 1; return 0; }
/* recolour from SIG: colours = rank of distinct signatures; returns number of colours */
static int recolour(int *col){
  int ord[MAXN]; for (int v = 0; v < CN; v++) ord[v] = v;
  for (int i = 1; i < CN; i++){ int x = ord[i], j = i - 1; while (j >= 0 && sigcmp(ord[j], x) > 0){ ord[j + 1] = ord[j]; j--; } ord[j + 1] = x; }
  int c = 0; col[ord[0]] = 0;
  for (int i = 1; i < CN; i++){ if (sigcmp(ord[i - 1], ord[i]) != 0) c++; col[ord[i]] = c; }
  return c + 1;
}
static int refine(int *col, int k){
  for (;;){
    SIGLEN = k + 1;
    for (int v = 0; v < CN; v++){
      SIG[v][0] = col[v];
      for (int c = 0; c < k; c++) SIG[v][c + 1] = 0;
      uint32_t o = CADJ[v];
      while (o){ int u = __builtin_ctz(o); o &= o - 1; SIG[v][col[u] + 1]++; }
    }
    int k2 = recolour(col);
    if (k2 == k) return k;
    k = k2;
  }
}
static void csearch(int *col, int k){
  k = refine(col, k);
  if (k == CN){
    int ord[MAXN]; for (int v = 0; v < CN; v++) ord[col[v]] = v;
    code_t c = make_code(ord);
    if (!CHAVE || code_cmp(&c, &CBEST) > 0){ CBEST = c; CHAVE = 1; }
    return;
  }
  int sz[MAXN] = {0}; for (int v = 0; v < CN; v++) sz[col[v]]++;
  int t = -1; for (int c = 0; c < k; c++) if (sz[c] > 1 && (t < 0 || sz[c] < sz[t])) t = c;
  for (int v = 0; v < CN; v++) if (col[v] == t){
    int c2[MAXN];
    for (int u = 0; u < CN; u++) c2[u] = col[u] > t ? col[u] + 1 : col[u];
    c2[v] = t + 1;                       /* individualised vertex placed after the rest of its cell */
    csearch(c2, k + 1);
  }
}
static code_t canon(const uint32_t *adj, int n){
  CN = n; memcpy(CADJ, adj, sizeof(uint32_t) * n); CHAVE = 0;
  int col[MAXN];
  /* initial colours: (score, triangles through v) */
  SIGLEN = 2;
  for (int v = 0; v < n; v++){
    int s = __builtin_popcount(adj[v]); int t = 0;
    /* triangles through v: v -> u -> x -> v */
    uint32_t o = adj[v];
    while (o){ int u = __builtin_ctz(o); o &= o - 1; uint32_t w = adj[u]; while (w){ int x = __builtin_ctz(w); w &= w - 1; if (adj[x] >> v & 1) t++; } }
    SIG[v][0] = s; SIG[v][1] = t;
  }
  int k = recolour(col);
  csearch(col, k);
  return CBEST;
}
static void decode(const code_t *c, int n, uint32_t *adj){
  memset(adj, 0, sizeof(uint32_t) * n); int b = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){
    if (c->w[b >> 6] >> (b & 63) & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i;
    b++;
  }
}
/* ------------------------------------------------------------------ cycles */
static int GN; static uint32_t G[MAXN], GIN[MAXN];
static void set_graph(const uint32_t *adj, int n){
  GN = n; memcpy(G, adj, sizeof(uint32_t) * n);
  for (int v = 0; v < n; v++) GIN[v] = 0;
  for (int v = 0; v < n; v++){ uint32_t o = adj[v]; while (o){ int u = __builtin_ctz(o); o &= o - 1; GIN[u] |= 1u << v; } }
}
/* vertices (within 'allowed') that can reach target t using only vertices of 'allowed' */
static uint32_t back_reach(int t, uint32_t allowed){
  uint32_t R = 1u << t, fr = 1u << t;
  while (fr){
    uint32_t nx = 0; uint32_t f = fr;
    while (f){ int w = __builtin_ctz(f); f &= f - 1; nx |= GIN[w]; }
    nx &= allowed & ~R;
    R |= nx; fr = nx;
  }
  return R;
}
/* enumerate cycles through 'start' using only vertices in 'allowed' (start in allowed); odd ones counted,
   optionally recorded; stop once the odd count exceeds cap */
static long long OCNT, OCAP; static int OSTART; static uint32_t OALLOW;
static uint32_t *REC; static long long RECN, RECCAP; static int RECORD;
static void dfs_cyc(int v, uint32_t vis, int len){
  if (OCNT > OCAP) return;
  if ((G[v] >> OSTART & 1) && len >= 3 && (len & 1)){
    OCNT++;
    if (RECORD){ if (RECN == RECCAP){ RECCAP = RECCAP ? 2 * RECCAP : 1024; REC = realloc(REC, RECCAP * sizeof(uint32_t)); } REC[RECN++] = vis; }
    if (OCNT > OCAP) return;
  }
  uint32_t free_ = OALLOW & ~vis;
  uint32_t R = back_reach(OSTART, free_ | (1u << OSTART));
  uint32_t cand = G[v] & free_ & R;
  while (cand){ int u = __builtin_ctz(cand); cand &= cand - 1; dfs_cyc(u, vis | (1u << u), len + 1); if (OCNT > OCAP) return; }
}
static long long odd_through(int z, uint32_t allowed, long long cap, int record){
  OCNT = 0; OCAP = cap; OSTART = z; OALLOW = allowed; RECORD = record;
  dfs_cyc(z, 1u << z, 1);
  return OCNT;
}
/* all odd cycles of the tournament (each once, from its least vertex), recorded as vertex masks */
static long long all_odd_cycles(int n, long long cap){
  RECN = 0; long long tot = 0;
  for (int s = 0; s < n; s++){
    uint32_t allowed = ((n == 32) ? ~0u : ((1u << n) - 1)) & ~((1u << s) - 1);
    long long c = odd_through(s, allowed, cap - tot, 1);
    tot += c; if (tot > cap) return tot;
  }
  return tot;
}
/* packing vector from the recorded cycle masks: group by support, then k-sets of disjoint supports */
static uint32_t SUP[1 << 16]; static u64 MUL[1 << 16]; static int NSUP;
static u64 ALPHA[16];
static void pack_rec(int i, uint32_t used, int k, u64 w){
  ALPHA[k] += w;
  for (int j = i; j < NSUP; j++) if (!(SUP[j] & used)) pack_rec(j + 1, used | SUP[j], k + 1, w * MUL[j]);
}
static void packing_vector(void){
  NSUP = 0;
  for (long long i = 0; i < RECN; i++){
    int f = -1; for (int j = 0; j < NSUP; j++) if (SUP[j] == REC[i]){ f = j; break; }
    if (f < 0){ if (NSUP >= (1 << 16)){ fprintf(stderr, "too many supports\n"); exit(5); } SUP[NSUP] = REC[i]; MUL[NSUP] = 1; NSUP++; }
    else MUL[f]++;
  }
  memset(ALPHA, 0, sizeof ALPHA);
  pack_rec(0, 0, 0, 1);
}
static u64 hpaths(const uint32_t *adj, int n){
  size_t F = (size_t)1 << n; u64 *f = calloc(F * n, sizeof(u64));
  for (int v = 0; v < n; v++) f[((size_t)1 << v) * n + v] = 1;
  for (size_t S = 1; S < F; S++) for (int v = 0; v < n; v++){ u64 x = f[S * n + v]; if (!x) continue;
    uint32_t o = adj[v] & ~(uint32_t)S; while (o){ int u = __builtin_ctz(o); o &= o - 1; f[(S | ((size_t)1 << u)) * n + u] += x; } }
  u64 h = 0; for (int v = 0; v < n; v++) h += f[(F - 1) * n + v]; free(f); return h;
}
/* ------------------------------------------------------------------ hash set */
typedef struct { code_t *key; size_t cap, cnt; } hset;
static u64 hcode(const code_t *c){ u64 h = 0x12345678ULL; for (int i = 0; i < CW; i++){ h ^= c->w[i] + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2); } return h * 0xD6E8FEB86659FD93ULL; }
static int hs_insert(hset *h, const code_t *k){
  if (2 * (h->cnt + 1) > h->cap){                    /* grow */
    hset g; g.cap = h->cap * 2; g.cnt = 0; g.key = malloc(g.cap * sizeof(code_t));
    for (size_t i = 0; i < g.cap; i++) g.key[i].w[CW - 1] = ~0ULL;
    for (size_t i = 0; i < h->cap; i++) if (h->key[i].w[CW - 1] != ~0ULL){
      size_t p = hcode(&h->key[i]) & (g.cap - 1); while (g.key[p].w[CW - 1] != ~0ULL) p = (p + 1) & (g.cap - 1); g.key[p] = h->key[i]; g.cnt++; }
    free(h->key); *h = g;
  }
  size_t p = hcode(k) & (h->cap - 1);
  while (h->key[p].w[CW - 1] != ~0ULL){ if (code_cmp(&h->key[p], k) == 0) return 0; p = (p + 1) & (h->cap - 1); }
  h->key[p] = *k; h->cnt++; return 1;
}
static void hs_init(hset *h){ h->cap = 1 << 12; h->cnt = 0; h->key = malloc(h->cap * sizeof(code_t)); for (size_t i = 0; i < h->cap; i++) h->key[i].w[CW - 1] = ~0ULL; }
static void print_code(const code_t *c){ for (int i = CW - 1; i >= 0; i--) printf("%016llx", c->w[i]); }

/* triangles of a tournament as masks */
static int triangles(const uint32_t *adj, int n, uint32_t *tri){
  int t = 0;
  for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++){
    int ab = adj[a] >> b & 1, bc = adj[b] >> c & 1, ca = adj[c] >> a & 1;
    if (ab == bc && bc == ca) tri[t++] = (1u << a) | (1u << b) | (1u << c);
  }
  return t;
}

static int strongq(const uint32_t *adj, int n){
  uint32_t full = (1u << n) - 1;
  uint32_t R = 1, fr = 1;
  while (fr){ uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= adj[__builtin_ctz(f)]; nx &= full & ~R; R |= nx; fr = nx; }
  if (R != full) return 0;
  uint32_t in[MAXN] = {0}; for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if (adj[v] >> u & 1) in[u] |= 1u << v;
  R = 1; fr = 1;
  while (fr){ uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= in[__builtin_ctz(f)]; nx &= full & ~R; R |= nx; fr = nx; }
  return R == full;
}
/* census mode: read kind-pasteur classgen lines "code H c3 strong SC aut ccode | ..." for order n and recompute
   everything with this file's own code (strong, odd cycles by length, packing vector, OCF against classgen's H,
   own canonical code, self-converse against classgen's SC flag) */
static int census(int n, const char *path){
  FILE *f = fopen(path, "r"); if (!f){ fprintf(stderr, "cannot open %s\n", path); return 1; }
  static char line[1 << 16]; long long cnt = 0;
  while (fgets(line, sizeof line, f)){
    unsigned long long code, H; int c3, st, scf;
    if (sscanf(line, "%llu %llu %d %d %d", &code, &H, &c3, &st, &scf) != 5) continue;
    uint32_t adj[MAXN]; memset(adj, 0, sizeof adj); int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){ if (code >> b & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i; b++; }
    set_graph(adj, n);
    int s = n <= 1 ? 1 : strongq(adj, n);
    if (s != st){ fprintf(stderr, "strong flag disagrees\n"); return 2; }
    all_odd_cycles(n, 1LL << 40);
    long long cl[24] = {0}; for (long long i = 0; i < RECN; i++) cl[__builtin_popcount(REC[i])]++;
    packing_vector();
    u64 Hocf = 0; for (int k = 0; k < 16; k++) Hocf += ALPHA[k] << k;
    if (Hocf != H){ fprintf(stderr, "OCF/H disagreement at order %d code %llu\n", n, code); return 3; }
    code_t c = canon(adj, n);
    uint32_t op[MAXN]; for (int v = 0; v < n; v++) op[v] = ((1u << n) - 1) & ~adj[v] & ~(1u << v);
    code_t cc = canon(op, n); int sc = code_cmp(&cc, &c) == 0;
    if (sc != scf){ fprintf(stderr, "self-converse flag disagrees\n"); return 4; }
    printf("%d %d", n, s); for (int k = 1; k <= 7; k++) printf(" %llu", ALPHA[k]); printf(" %llu %d", Hocf, sc);
    printf(" c3=%lld c5=%lld c7=%lld c9=%lld ", cl[3], cl[5], cl[7], cl[9]); print_code(&c); printf("\n");
    cnt++;
  }
  fclose(f); fprintf(stderr, "census order %d: %lld classes checked (OCF, strong, SC agree with classgen)\n", n, cnt);
  return 0;
}

/* fibre mode: nodes = all classes of orders 1..N (from this program's census files census_<n>.txt);
   edges = H-preserving one-vertex deletions and H-preserving single-arc reversals; report the classes whose
   component contains no first creator (class of least order with that H), per order, and the H = 9 fibre */
typedef struct { code_t c; long long H; int id; int n; } frec;
static frec *FR[MAXN]; static long long FNR[MAXN];
static int fcmp(const void *a, const void *b){ return code_cmp(&((const frec *)a)->c, &((const frec *)b)->c); }
static frec *ffind(int n, const code_t *c){ long long lo = 0, hi = FNR[n]; while (lo < hi){ long long m = (lo + hi) / 2; if (code_cmp(&FR[n][m].c, c) < 0) lo = m + 1; else hi = m; } return (lo < FNR[n] && code_cmp(&FR[n][lo].c, c) == 0) ? &FR[n][lo] : NULL; }
static int *UF; static int uf(int x){ while (UF[x] != x){ UF[x] = UF[UF[x]]; x = UF[x]; } return x; }
static void ujoin(int a, int b){ a = uf(a); b = uf(b); if (a != b) UF[a] = b; }
static int fibres(int N, const char *dir){
  int tot = 0; static char line[4096];
  for (int n = 1; n <= N; n++){
    char fn[2048]; snprintf(fn, sizeof fn, "%s/census_%d.txt", dir, n); FILE *f = fopen(fn, "r"); if (!f){ fprintf(stderr, "missing %s\n", fn); return 1; }
    long long cap = 1024; FR[n] = malloc(cap * sizeof(frec)); FNR[n] = 0;
    while (fgets(line, sizeof line, f)){
      char *tok[32]; int nt = 0; char *p = strtok(line, " \n"); while (p && nt < 32){ tok[nt++] = p; p = strtok(NULL, " \n"); }
      if (nt < 16) continue;
      if (FNR[n] == cap){ cap *= 2; FR[n] = realloc(FR[n], cap * sizeof(frec)); }
      frec *r = &FR[n][FNR[n]++]; r->n = n; r->H = atoll(tok[9]);
      const char *hex = tok[15]; for (int i = 0; i < CW; i++){ char w[17]; memcpy(w, hex + 16 * (CW - 1 - i), 16); w[16] = 0; r->c.w[i] = strtoull(w, 0, 16); }
    }
    fclose(f); qsort(FR[n], FNR[n], sizeof(frec), fcmp);
    for (long long i = 0; i < FNR[n]; i++) FR[n][i].id = tot++;
  }
  UF = malloc(tot * sizeof(int)); for (int i = 0; i < tot; i++) UF[i] = i;
  long long de = 0, fe = 0;
  for (int n = 1; n <= N; n++) for (long long r = 0; r < FNR[n]; r++){
    uint32_t adj[MAXN]; decode(&FR[n][r].c, n, adj); long long H = FR[n][r].H; int id = FR[n][r].id;
    if (n > 1) for (int x = 0; x < n; x++){
      uint32_t d[MAXN]; int k = 0;
      for (int v = 0; v < n; v++) if (v != x){ uint32_t m = 0; int kk = 0; for (int u = 0; u < n; u++) if (u != x){ if (adj[v] >> u & 1) m |= 1u << kk; kk++; } d[k++] = m; }
      code_t c = canon(d, n - 1); frec *q = ffind(n - 1, &c); if (!q){ fprintf(stderr, "lookup failure (deletion)\n"); return 2; }
      if (q->H == H){ ujoin(id, q->id); de++; }
    }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++){
      uint32_t a2[MAXN]; memcpy(a2, adj, sizeof(uint32_t) * n); a2[i] ^= 1u << j; a2[j] ^= 1u << i;
      code_t c = canon(a2, n); frec *q = ffind(n, &c); if (!q){ fprintf(stderr, "lookup failure (flip)\n"); return 2; }
      if (q->H == H){ ujoin(id, q->id); fe++; }
    }
  }
  /* first-creator order of each value */
  #define FHMAX 100000
  static int first[FHMAX]; for (int h = 0; h < FHMAX; h++) first[h] = 0;
  for (int n = 1; n <= N; n++) for (long long r = 0; r < FNR[n]; r++){ long long H = FR[n][r].H; if (H < FHMAX && !first[H]) first[H] = n; }
  char *fc = calloc(tot, 1);
  for (int n = 1; n <= N; n++) for (long long r = 0; r < FNR[n]; r++) if (first[FR[n][r].H] == n) fc[uf(FR[n][r].id)] = 1;
  printf("fibre moves, orders 1..%d: %d classes, H-preserving deletion edges %lld, flip edges %lld\n", N, tot, de, fe);
  for (int n = 1; n <= N; n++){ long long o = 0; for (long long r = 0; r < FNR[n]; r++) if (!fc[uf(FR[n][r].id)]) o++; printf("  order %d: %lld of %lld classes in no component with a first creator of their value\n", n, o, FNR[n]); }
  /* components of the H = 9 fibre */
  int comps[64]; int nc = 0;
  for (int n = 1; n <= N; n++) for (long long r = 0; r < FNR[n]; r++) if (FR[n][r].H == 9){ int c = uf(FR[n][r].id); int seen = 0; for (int i = 0; i < nc; i++) if (comps[i] == c) seen = 1; if (!seen) comps[nc++] = c; }
  printf("  H = 9 fibre: %d components; sizes per order:", nc);
  for (int i = 0; i < nc; i++){ printf(" ["); for (int n = 1; n <= N; n++){ int k = 0; for (long long r = 0; r < FNR[n]; r++) if (FR[n][r].H == 9 && uf(FR[n][r].id) == comps[i]) k++; if (k) printf(" n%d:%d", n, k); } printf(" ]"); }
  printf("\n");
  return 0;
}

int main(int argc, char **argv){
  long long M = argc > 1 ? atoll(argv[1]) : 80;
  int mode = argc > 2 ? atoi(argv[2]) : 0;
  int hchk = argc > 3 ? atoi(argv[3]) : 12;
  if (mode == 2) return census(hchk, argv[4]);
  if (mode == 4) return fibres(hchk, argv[4]);
  if (mode == 3){                     /* invariance test: relabel every class of an inventory file K times */
    FILE *f = fopen(argv[4], "r"); if (!f) return 1; static char line[4096]; long long nc = 0, bad = 0; unsigned long long seed = 88172645463325252ULL;
    while (fgets(line, sizeof line, f)){
      int n; char hex[128]; unsigned long long a[7], Ho; long long Hd; int sc;
      if (sscanf(line, "%d %llu %llu %llu %llu %llu %llu %llu %llu %lld %d %127s", &n, &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &Ho, &Hd, &sc, hex) != 12) continue;
      code_t c; memset(&c, 0, sizeof c);
      for (int i = 0; i < CW; i++){ char w[17]; memcpy(w, hex + 16 * (CW - 1 - i), 16); w[16] = 0; c.w[i] = strtoull(w, 0, 16); }
      uint32_t adj[MAXN]; decode(&c, n, adj);
      code_t c0 = canon(adj, n); if (code_cmp(&c0, &c) != 0){ bad++; continue; }
      for (int rep = 0; rep < hchk; rep++){
        int perm[MAXN]; for (int i = 0; i < n; i++) perm[i] = i;
        for (int i = n - 1; i > 0; i--){ seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; int j = (int)(seed % (unsigned long long)(i + 1)); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        uint32_t b2[MAXN] = {0}; for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if (adj[v] >> u & 1) b2[perm[v]] |= 1u << perm[u];
        code_t c1 = canon(b2, n); if (code_cmp(&c1, &c) != 0){ bad++; break; }
      }
      nc++;
    }
    fclose(f); printf("invariance test: %lld classes, %d random relabellings each, failures %lld\n", nc, hchk, bad); return bad != 0;
  }
  code_t *reps; long long *ra1; size_t nreps;
  uint32_t c3[3] = { 1u << 1, 1u << 2, 1u << 0 };
  reps = malloc(sizeof(code_t)); ra1 = malloc(sizeof(long long)); reps[0] = canon(c3, 3); ra1[0] = 1; nreps = 1;
  long long total = 0;
  for (int n = 3; n < MAXN; n++){
    if (n > 3){
      hset h; hs_init(&h);
      code_t *nr = NULL; long long *na = NULL; size_t nnr = 0, cap = 0;
      long long cand_t3 = 0, cand_dfs = 0;
      for (size_t r = 0; r < nreps; r++){
        uint32_t adj[MAXN]; decode(&reps[r], n - 1, adj); long long base = ra1[r];
        uint32_t fullm = (1u << (n - 1)) - 1;
        uint32_t cover[MAXN]; static uint32_t tri[1 << 14]; int nt = 0;
        if (mode == 1){
          nt = triangles(adj, n - 1, tri);
          for (int p = 0; p < n - 1; p++){ cover[p] = 0; for (int q = 0; q < n - 1; q++){ if (q == p) continue; uint32_t pq = (1u << p) | (1u << q); int ok = 1; for (int i = 0; i < nt; i++) if (!(tri[i] & pq)){ ok = 0; break; } if (ok) cover[p] |= 1u << q; } }
        }
        for (uint32_t P = 1; P < fullm; P++){        /* z beats exactly P; P nonempty and proper <=> strong */
          uint32_t Q = fullm & ~P; long long t3 = 0; int bad = 0;
          for (uint32_t X = P; X; X &= X - 1){ int p = __builtin_ctz(X); uint32_t f = adj[p] & Q; t3 += __builtin_popcount(f); if (mode == 1 && (f & ~cover[p])){ bad = 1; break; } }
          if (bad) continue;
          if (base + t3 > M) continue;
          cand_t3++;
          uint32_t a2[MAXN]; memcpy(a2, adj, sizeof(uint32_t) * (n - 1)); a2[n - 1] = P;
          for (int v = 0; v < n - 1; v++) if (!(P >> v & 1)) a2[v] |= 1u << (n - 1);
          set_graph(a2, n);
          long long thr = odd_through(n - 1, (1u << n) - 1, M - base, 0);
          if (base + thr > M) continue;
          cand_dfs++;
          code_t c = canon(a2, n);
          if (hs_insert(&h, &c)){
            if (nnr == cap){ cap = cap ? 2 * cap : 1024; nr = realloc(nr, cap * sizeof(code_t)); na = realloc(na, cap * sizeof(long long)); }
            nr[nnr] = c; na[nnr] = base + thr; nnr++;
          }
        }
      }
      free(h.key); free(reps); free(ra1); reps = nr; ra1 = na; nreps = nnr;
      fprintf(stderr, "order %d: candidates after t3 %lld, after cycle cap %lld, classes %zu\n", n, cand_t3, cand_dfs, nreps);
    }
    if (nreps == 0){ fprintf(stderr, "order %d EMPTY: generation complete (deletion lemma)\n", n); break; }
    for (size_t r = 0; r < nreps; r++){
      uint32_t adj[MAXN]; decode(&reps[r], n, adj); set_graph(adj, n);
      long long a1 = all_odd_cycles(n, 1LL << 40);
      if (a1 != ra1[r]){ fprintf(stderr, "alpha_1 mismatch at order %d: %lld vs %lld\n", n, a1, ra1[r]); return 3; }
      packing_vector();
      if ((long long)ALPHA[1] != a1){ fprintf(stderr, "support grouping mismatch\n"); return 3; }
      if (mode == 1 && ALPHA[2] != 0){ fprintf(stderr, "intercyclic violated\n"); return 3; }
      u64 Hocf = 0; for (int k = 0; k < 16; k++) Hocf += ALPHA[k] << k;
      long long Hdp = -1;
      if (n <= hchk){ Hdp = (long long)hpaths(adj, n); if ((u64)Hdp != Hocf){ fprintf(stderr, "OCF mismatch at order %d\n", n); return 4; } }
      uint32_t op[MAXN]; for (int v = 0; v < n; v++) op[v] = ((1u << n) - 1) & ~adj[v] & ~(1u << v);
      code_t cc = canon(op, n); int sc = code_cmp(&cc, &reps[r]) == 0;
      printf("%d", n); for (int k = 1; k <= 7; k++) printf(" %llu", ALPHA[k]); printf(" %llu %lld %d ", Hocf, Hdp, sc); print_code(&reps[r]); printf("\n");
    }
    fflush(stdout); total += nreps;
    fprintf(stderr, "order %d: %zu classes written (total %lld)\n", n, nreps, total);
  }
  return 0;
}
