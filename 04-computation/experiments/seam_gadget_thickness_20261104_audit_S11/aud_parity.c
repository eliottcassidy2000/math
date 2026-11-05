/* aud_parity.c -- audit S11 of THM-4879 (D) and of the proposed addendum (F).
   For every labelled tournament of order 1..NMAX (exhaustive) and random ones of larger order:
     (D)  R_a == S_a (mod 2) for every root a;
     (F1) Pi(W) = number of nonempty directed paths (any vertex set) is odd;
     (F2) for every ordered pair z != u:  L_(z,u) == M_(z,u) == P_(z->u) + P_(u->z) (mod 2);
          L = (A from z, B from u), M = (A ends at z, B ends at u), both nonempty, covering, ordered.
   Also: the number of rows with s_z = P_(z->u) = L_(z,u) = 0 and how many of them have M odd, and the full GF(2)
   relation space among (s_z, e_z, s_u, e_u, P_zu, P_uz, L, Lmix, Mmix, M) (Lmix: A from z, B ends at u;
   Mmix: A ends at z, B from u), homogeneous and affine.
   Counts by subset DP; the DP is validated against permutation brute force on random tournaments of order <= 7.  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MV 13
static int T[MV][MV];
static ll f[1 << MV][MV], g[1 << MV][MV], Hm[1 << MV], pz[1 << MV][MV];
static unsigned long long rng_s = 0x2545F4914F6CDD1DULL;
static unsigned long long rnd(void){ rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

static void dp_all(int n){
  int full = (1 << n) - 1;
  for (int m = 0; m <= full; m++) for (int v = 0; v < n; v++) { f[m][v] = 0; g[m][v] = 0; }
  for (int v = 0; v < n; v++) { f[1 << v][v] = 1; g[1 << v][v] = 1; }
  for (int m = 1; m <= full; m++) for (int v = 0; v < n; v++) {
    ll x = f[m][v], y = g[m][v]; if (!x && !y) continue;
    for (int w = 0; w < n; w++) if (!(m >> w & 1)) { if (x && T[v][w]) f[m | 1 << w][w] += x; if (y && T[w][v]) g[m | 1 << w][w] += y; }
  }
  Hm[0] = 1; for (int m = 1; m <= full; m++) { ll s = 0; for (int v = 0; v < n; v++) s += f[m][v]; Hm[m] = s; }
}
/* Ham paths of the whole set from z, ending at each v: pz[full][v] */
static void dp_from(int n, int z){
  int full = (1 << n) - 1;
  for (int m = 0; m <= full; m++) if (m >> z & 1) for (int v = 0; v < n; v++) pz[m][v] = 0;
  pz[1 << z][z] = 1;
  for (int m = 1; m <= full; m++) if (m >> z & 1) for (int v = 0; v < n; v++) { ll x = pz[m][v]; if (!x) continue;
    for (int w = 0; w < n; w++) if (!(m >> w & 1) && T[v][w]) pz[m | 1 << w][w] += x; }
}
/* permutation brute force for validation */
static int perm[MV], nb; static ll bH, bS[MV], bR[MV], bP[MV][MV];
static void visit(void){
  int n = nb, pre[MV + 1], suf[MV + 1];
  pre[0] = 1; pre[1] = 1; for (int c = 2; c <= n; c++) pre[c] = pre[c - 1] && T[perm[c - 2]][perm[c - 1]];
  suf[n] = 1; suf[n - 1] = 1; for (int c = n - 2; c >= 0; c--) suf[c] = suf[c + 1] && T[perm[c]][perm[c + 1]];
  if (pre[n]) { bH++; bS[perm[0]]++; bP[perm[0]][perm[n - 1]]++; }
  for (int c = 1; c <= n; c++) if (pre[c] && suf[c]) bR[perm[c - 1]]++;
}
static void gen(int k, int n){ if (k == n) { visit(); return; }
  for (int i = k; i < n; i++) { int t = perm[k]; perm[k] = perm[i]; perm[i] = t; gen(k + 1, n); t = perm[k]; perm[k] = perm[i]; perm[i] = t; } }

static ll rowsD = 0, badD = 0, tours = 0, badF1 = 0, rowsF2 = 0, badF2a = 0, badF2b = 0, zeroS = 0, zeroS_Modd = 0;
static int seen[1024];
static void check_tour(int n, int doF2){
  int full = (1 << n) - 1; dp_all(n); tours++;
  /* (F1) */
  ll Pi = 0; for (int m = 1; m <= full; m++) Pi += Hm[m]; if (!(Pi & 1)) badF1++;
  /* (D) */
  for (int a = 0; a < n; a++) { ll R = 0; for (int m = 1; m <= full; m++) if (m >> a & 1) R += f[m][a] * Hm[full ^ m];
    ll S = g[full][a]; rowsD++; if ((R - S) & 1) badD++; }
  if (!doF2 || n < 2) return;
  static ll P[MV][MV];
  for (int z = 0; z < n; z++) { dp_from(n, z); for (int v = 0; v < n; v++) P[z][v] = pz[full][v]; }
  for (int z = 0; z < n; z++) for (int u = 0; u < n; u++) if (z != u) {
    ll L = 0, M = 0, Lm = 0, Mm = 0;
    for (int m = 1; m < full; m++) if ((m >> z & 1) && !(m >> u & 1)) { int c = full ^ m;
      L += g[m][z] * g[c][u]; M += f[m][z] * f[c][u]; Lm += g[m][z] * f[c][u]; Mm += f[m][z] * g[c][u]; }
    ll sz = g[full][z], ez = f[full][z], su = g[full][u], eu = f[full][u];
    rowsF2++;
    if ((L - M) & 1) badF2a++;
    if ((M - P[z][u] - P[u][z]) & 1) badF2b++;
    if (sz == 0 && P[z][u] == 0 && L == 0) { zeroS++; if (M & 1) zeroS_Modd++; }
    int vec = (int)((sz & 1) | (ez & 1) << 1 | (su & 1) << 2 | (eu & 1) << 3 | (P[z][u] & 1) << 4 | (P[u][z] & 1) << 5 | (L & 1) << 6 | (Lm & 1) << 7 | (Mm & 1) << 8 | (M & 1) << 9);
    seen[vec] = 1;
  }
}
/* GF(2) nullspace dimension and basis of the relation space of the observed 10-bit vectors */
static void relations(const char *tag){
  const char *names[10] = {"s_z", "e_z", "s_u", "e_u", "P_zu", "P_uz", "L", "Lmix", "Mmix", "M"};
  /* homogeneous: c with popcount(c & x) even for all observed x */
  printf("%s: homogeneous relations c.x = 0 (mod 2) holding on every row:\n", tag);
  int nh = 0; for (int c = 1; c < 1024; c++) { int ok = 1; for (int x = 0; x < 1024 && ok; x++) if (seen[x] && (__builtin_popcount(c & x) & 1)) ok = 0;
    if (ok) { nh++; printf("   "); int first = 1; for (int i = 0; i < 10; i++) if (c >> i & 1) { printf("%s%s", first ? "" : " + ", names[i]); first = 0; } printf(" == 0\n"); } }
  int dim = 0; while ((1 << dim) < nh + 1) dim++;
  printf("   (%d nonzero relations => relation space of dimension %d)\n", nh, dim);
  int na = 0; for (int c = 1; c < 1024; c++) { int val = -1, ok = 1; for (int x = 0; x < 1024 && ok; x++) if (seen[x]) { int v = __builtin_popcount(c & x) & 1; if (val < 0) val = v; else if (v != val) ok = 0; }
    if (ok && val == 1) na++; }
  printf("   affine relations c.x = 1 holding on every row: %d\n", na);
  int npat = 0; for (int x = 0; x < 1024; x++) npat += seen[x]; printf("   distinct parity patterns observed: %d\n", npat);
}
int main(int argc, char **argv){
  int NMAX = argc > 1 ? atoi(argv[1]) : 7, F2MAX = argc > 2 ? atoi(argv[2]) : 7;
  /* validation of the DP counts against permutation brute force */
  ll vbad = 0, vcnt = 0;
  for (int t = 0; t < 300; t++) { int n = 1 + rnd() % 7;
    for (int i = 0; i < n; i++) { T[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = rnd() & 1; T[i][j] = e; T[j][i] = 1 - e; } }
    dp_all(n); int full = (1 << n) - 1; nb = n; bH = 0; memset(bS, 0, sizeof bS); memset(bR, 0, sizeof bR); memset(bP, 0, sizeof bP);
    for (int i = 0; i < n; i++) perm[i] = i; gen(0, n);
    for (int a = 0; a < n; a++) { ll R = 0; for (int m = 1; m <= full; m++) if (m >> a & 1) R += f[m][a] * Hm[full ^ m];
      vcnt++; if (R != bR[a] || g[full][a] != bS[a] || Hm[full] != bH) vbad++;
      dp_from(n, a); for (int v = 0; v < n; v++) if (pz[full][v] != bP[a][v]) vbad++; }
  }
  printf("DP validation vs permutation brute force: %lld rooted cases, %lld mismatches\n", vcnt, vbad);
  for (int n = 1; n <= NMAX; n++) {
    ll r0 = rowsD, b0 = badD, t0 = tours, f0 = badF1, q0 = rowsF2, a0 = badF2a, c0 = badF2b, z0 = zeroS, zo0 = zeroS_Modd;
    int ne = n * (n - 1) / 2;
    for (unsigned long long code = 0; code < (1ULL << ne); code++) {
      int b = 0; for (int i = 0; i < n; i++) { T[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = (code >> b) & 1; b++; T[i][j] = e; T[j][i] = 1 - e; } }
      check_tour(n, n <= F2MAX);
    }
    printf("order %d (all %lld labelled): (D) rooted %lld viol %lld | (F1) viol %lld | (F2) rows %lld: L!=M %lld, M!=Pzu+Puz %lld | rows s_z=P_zu=L=0: %lld, with M odd: %lld\n",
           n, tours - t0, rowsD - r0, badD - b0, badF1 - f0, rowsF2 - q0, badF2a - a0, badF2b - c0, zeroS - z0, zeroS_Modd - zo0);
    fflush(stdout);
    if (n == 6) { printf("cumulative orders 1-6: (D) rooted %lld; (F2) rows %lld; rows with s_z=P_zu=L=0: %lld (M odd: %lld)\n", rowsD, rowsF2, zeroS, zeroS_Modd); relations("orders 2-6"); }
  }
  if (NMAX >= 7) relations("orders 2-7");
  /* random larger tournaments */
  for (int n = 8; n <= 12; n++) {
    ll r0 = rowsD, b0 = badD, f0 = badF1, q0 = rowsF2, a0 = badF2a, c0 = badF2b; int reps = n <= 10 ? 300 : 40;
    for (int t = 0; t < reps; t++) { for (int i = 0; i < n; i++) { T[i][i] = 0; for (int j = i + 1; j < n; j++) { int e = rnd() & 1; T[i][j] = e; T[j][i] = 1 - e; } }
      check_tour(n, n <= 10); }
    printf("random order %d (%d tournaments): (D) rooted %lld viol %lld | (F1) viol %lld | (F2) rows %lld viol %lld/%lld\n", n, reps, rowsD - r0, badD - b0, badF1 - f0, rowsF2 - q0, badF2a - a0, badF2b - c0);
  }
  printf("TOTAL (D) rooted %lld violations %lld; (F1) tournaments %lld violations %lld; (F2) rows %lld violations %lld + %lld\n", rowsD, badD, tours, badF1, rowsF2, badF2a, badF2b);
  return 0;
}
