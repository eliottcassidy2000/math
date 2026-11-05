/* aud_rows_brute.c -- independent brute-force audit of THM-4879 (A) (audit S11).
   Everything is counted by enumerating permutations of a vertex set together with cut points, i.e. ordered
   tuples of vertex-disjoint paths covering the set.  No dynamic programming, no code from the audited session.

   State of a rooted digraph (Y, a):  H = #Ham paths, S = #Ham paths starting at a,
     R = #(perm, c) with perm[0..c-1] a path ending at a and perm[c..] a path (possibly empty)
       = sum_{U contains a} E(Y[U], a) H(Y - U),  H(empty) = 1.
   Seam X (.) G (exit a of X, entry u, exit z of G): internal arcs kept, every cross arc G -> X except a -> u.

   Modes:
     exhaustive          all labelled gadgets of order 1..4, all (u,z), x all labelled rooted prefixes of order 1..3
     order5 SEED         all labelled gadgets of order 5, all (u,z), x one random rooted prefix of order 1..3
     random N SEED       N random instances: gadget order 1..5, prefix order 1..4 (total <= 9 vertices)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MV 10

static unsigned long long rng_s = 88172645463325252ULL;
static unsigned long long rnd(void){ rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

/* ---------- state of (Y, a) by permutation brute force ---------- */
static int Y[MV][MV]; static int ny;
static int perm[MV];
static ll bH, bS, bR; static int broot;
static void visit_state(void){
  int n = ny; int pre[MV + 1], suf[MV + 1];
  pre[0] = 1; pre[1] = 1; for (int c = 2; c <= n; c++) pre[c] = pre[c - 1] && Y[perm[c - 2]][perm[c - 1]];
  suf[n] = 1; suf[n - 1] = 1; for (int c = n - 2; c >= 0; c--) suf[c] = suf[c + 1] && Y[perm[c]][perm[c + 1]];
  if (pre[n]) { bH++; if (perm[0] == broot) bS++; }
  for (int c = 1; c <= n; c++) if (perm[c - 1] == broot && pre[c] && suf[c]) bR++;
}
static void gen(int k, int n, void (*cb)(void)){
  if (k == n) { cb(); return; }
  for (int i = k; i < n; i++) { int t = perm[k]; perm[k] = perm[i]; perm[i] = t; gen(k + 1, n, cb); t = perm[k]; perm[k] = perm[i]; perm[i] = t; }
}
static void state(int n, int a, ll out[3]){ ny = n; broot = a; bH = bS = bR = 0; for (int i = 0; i < n; i++) perm[i] = i; gen(0, n, visit_state); out[0] = bH; out[1] = bS; out[2] = bR; }

/* ---------- gadget statistics, all (u,z) at once, by permutation brute force ---------- */
static int G[MV][MV]; static int ng;
static ll h_, E_[MV], s_[MV], P_[MV][MV], Lu_[MV], Lzu_[MV][MV], Rz_[MV], M_[MV][MV], K1_[MV][MV], K2_[MV][MV];
static int seg_ok(int *brk, int i, int j){ /* segment perm[i..j-1] is a path (empty and singletons are paths) */
  for (int t = i; t + 1 < j; t++) if (!brk[t]) return 0; return 1; }
static void visit_gadget(void){
  int g = ng, brk[MV]; for (int t = 0; t + 1 < g; t++) brk[t] = G[perm[t]][perm[t + 1]];
  int full = seg_ok(brk, 0, g);
  if (full) { h_++; E_[perm[g - 1]]++; s_[perm[0]]++; P_[perm[0]][perm[g - 1]]++; }
  /* L_u: (A1 possibly empty, A2 nonempty starting at u) */
  for (int c = 0; c <= g - 1; c++) if (seg_ok(brk, 0, c) && seg_ok(brk, c, g)) Lu_[perm[c]]++;
  /* L_(z,u): A1 nonempty from z, A2 nonempty from u */
  for (int c = 1; c <= g - 1; c++) if (seg_ok(brk, 0, c) && seg_ok(brk, c, g)) Lzu_[perm[0]][perm[c]]++;
  /* R_z(G): P nonempty ending at z, Q possibly empty */
  for (int c = 1; c <= g; c++) if (seg_ok(brk, 0, c) && seg_ok(brk, c, g)) Rz_[perm[c - 1]]++;
  /* M_(z,u): P nonempty ending at z, Q nonempty ending at u */
  for (int c = 1; c <= g - 1; c++) if (seg_ok(brk, 0, c) && seg_ok(brk, c, g)) M_[perm[c - 1]][perm[g - 1]]++;
  /* K1: (A1 possibly empty, A2 nonempty from u to z, Q possibly empty) */
  for (int c1 = 0; c1 < g; c1++) if (seg_ok(brk, 0, c1))
    for (int c2 = c1 + 1; c2 <= g; c2++) if (seg_ok(brk, c1, c2) && seg_ok(brk, c2, g)) K1_[perm[c1]][perm[c2 - 1]]++;
  /* K2: (P nonempty ending at z, Q1 possibly empty, Q3 nonempty starting at u); index [z][u] */
  for (int c1 = 1; c1 <= g - 1; c1++) if (seg_ok(brk, 0, c1))
    for (int c2 = c1; c2 <= g - 1; c2++) if (seg_ok(brk, c1, c2) && seg_ok(brk, c2, g)) K2_[perm[c1 - 1]][perm[c2]]++;
}
static void gadget_stats(int g){
  ng = g; h_ = 0; memset(E_, 0, sizeof E_); memset(s_, 0, sizeof s_); memset(P_, 0, sizeof P_); memset(Lu_, 0, sizeof Lu_);
  memset(Lzu_, 0, sizeof Lzu_); memset(Rz_, 0, sizeof Rz_); memset(M_, 0, sizeof M_); memset(K1_, 0, sizeof K1_); memset(K2_, 0, sizeof K2_);
  for (int i = 0; i < g; i++) perm[i] = i; gen(0, g, visit_gadget);
}
/* the (A) matrix for entry u, exit z */
static void rowmat(int u, int z, ll M[3][3]){
  M[0][0] = h_;       M[0][1] = -E_[u];      M[0][2] = Lu_[u];
  M[1][0] = s_[z];    M[1][1] = -P_[z][u];   M[1][2] = Lzu_[z][u];
  M[2][0] = Rz_[z];   M[2][1] = -M_[z][u];   M[2][2] = K1_[u][z] + K2_[z][u];
}

/* ---------- tournaments ---------- */
static void tour_from_code(int n, unsigned code, int T[MV][MV]){
  int b = 0; for (int i = 0; i < n; i++) { T[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = (code >> b) & 1; b++; T[i][j] = e; T[j][i] = 1 - e; } }
}
static void rand_tour(int n, int T[MV][MV]){ for (int i = 0; i < n; i++) { T[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = rnd() & 1; T[i][j] = e; T[j][i] = 1 - e; } } }

static ll ninst = 0, nbad = 0, chk = 0;
static int Xp[MV][MV];
/* test one instance: prefix (Xp, p, a), gadget G (order g, stats already computed), ports u, z */
static void test(int p, int a, int g, int u, int z){
  ll vin[3], vout[3], M[3][3];
  for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) Y[i][j] = Xp[i][j];
  state(p, a, vin);
  int n = p + g;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i >= p || j >= p) Y[i][j] = 0;
  for (int i = 0; i < g; i++) for (int j = 0; j < g; j++) Y[p + i][p + j] = G[i][j];
  for (int i = 0; i < g; i++) for (int j = 0; j < p; j++) { Y[p + i][j] = 1; Y[j][p + i] = 0; }
  Y[p + u][a] = 0; Y[a][p + u] = 1;
  state(n, p + z, vout);
  rowmat(u, z, M);
  ll pred[3]; for (int r = 0; r < 3; r++) pred[r] = M[r][0] * vin[0] + M[r][1] * vin[1] + M[r][2] * vin[2];
  ninst++; chk += vout[0] * 3 + vout[1] * 5 + vout[2] * 7;
  if (pred[0] != vout[0] || pred[1] != vout[1] || pred[2] != vout[2]) {
    nbad++; if (nbad <= 10) printf("MISMATCH p=%d a=%d g=%d u=%d z=%d pred=(%lld,%lld,%lld) brute=(%lld,%lld,%lld)\n", p, a, g, u, z, pred[0], pred[1], pred[2], vout[0], vout[1], vout[2]);
  }
}

int main(int argc, char **argv){
  if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
  if (!strcmp(argv[1], "exhaustive")) {
    for (int g = 1; g <= 4; g++) {
      int ne = g * (g - 1) / 2; ll before = ninst, badb = nbad;
      for (unsigned gc = 0; gc < (1u << ne); gc++) {
        tour_from_code(g, gc, G); gadget_stats(g);
        for (int u = 0; u < g; u++) for (int z = 0; z < g; z++)
          for (int p = 1; p <= 3; p++) { int pe = p * (p - 1) / 2;
            for (unsigned pc = 0; pc < (1u << pe); pc++) { tour_from_code(p, pc, Xp); for (int a = 0; a < p; a++) test(p, a, g, u, z); } }
      }
      printf("exhaustive: gadget order %d: %lld instances, %lld mismatches\n", g, ninst - before, nbad - badb);
    }
  } else if (!strcmp(argv[1], "order5")) {
    rng_s ^= (unsigned long long)atoll(argv[2]) * 0x9E3779B97F4A7C15ULL;
    ll before = ninst;
    for (unsigned gc = 0; gc < 1024; gc++) {
      tour_from_code(5, gc, G); gadget_stats(5);
      for (int u = 0; u < 5; u++) for (int z = 0; z < 5; z++) { int p = 1 + rnd() % 3; rand_tour(p, Xp); int a = rnd() % p; test(p, a, 5, u, z); }
    }
    printf("order5: all 1024 labelled gadgets x 25 port pairs x 1 random prefix (order 1-3): %lld instances, %lld mismatches\n", ninst - before, nbad);
  } else if (!strcmp(argv[1], "random")) {
    ll N = atoll(argv[2]); rng_s ^= (unsigned long long)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    ll cnt[6][5] = {{0}};
    for (ll t = 0; t < N; t++) {
      int g = 1 + rnd() % 5, p = 1 + rnd() % 4;
      rand_tour(g, G); gadget_stats(g); rand_tour(p, Xp);
      int u = rnd() % g, z = rnd() % g, a = rnd() % p; test(p, a, g, u, z); cnt[g][p]++;
    }
    printf("random: %lld instances (gadget order 1-5 x prefix order 1-4), %lld mismatches; per (g,p):", ninst, nbad);
    for (int g = 1; g <= 5; g++) for (int p = 1; p <= 4; p++) printf(" (%d,%d):%lld", g, p, cnt[g][p]); printf("\n");
  }
  printf("total instances %lld mismatches %lld checksum %lld\n", ninst, nbad, chk);
  return 0;
}
