// b_thmC.c -- AUDIT check (b): Theorem C of the tilted-ensemble packet. Written from scratch by the auditor.
//
// usage:  gentourng -q m | ./b_thmC m direct        all classes of order m, H(R+v_X) by a DIRECT DP on the (m+1)-vertex
//                                                    tournament for every out-set X (no use of formula (4.1))
//         gentourng -q m | ./b_thmC m fast          values from (4.2) with a_i, mu_ij taken from the two-path-cover matrix,
//                                                    plus one direct DP per class at a pseudo-random X as a control
//         ./b_thmC m labelled                        all 2^C(m,2) labelled tournaments on m vertices, direct mode
//
// Checked for every tournament R:
//   (i)    H(R+v_empty) = H(R+v_V) = H(R)
//   (ii)   (4.1): H(R+v_X) = sum_{i in X} g(i) + sum_{j notin X} f(j) + sum_{j notin X, i in X} M[j][i]      (direct mode, all X)
//   (iii)  (4.2): H(R+v_X) = H(R) + sum_{i in X} a_i - sum_{{i,j} in X} mu_ij  with a_i, mu_ij read off the DIRECT values
//          at singletons and pairs (direct mode, all X);  mu_ij = M[i][j] + M[j][i]
//   (iv)   mu_ij even and >= 2
//   (v)    q(X) = ((H(R+v_X)-H(R))/2 mod 2) equals sum_{i in X} q({i}) + sum_{pairs in X} beta_ij  (quadratic form, polar beta)
//   (vi)   #{q=0} - #{q=1} is 0 or +-2^(m - rank(beta)/2); it is 0 iff q does not vanish on the radical of beta
// Output: rank histogram over classes and over labelled tournaments (weights m!/|Aut|), exact means of the imbalance,
//         and the exact labelled signed sum  sum_R (-1)^((H(R)-1)/2) (#{q=0}-#{q=1})  (= 2^C(m+1,2) (P(H=1)-P(H=3)) at order m+1).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXM 10
static int m, full;
static int adj[MAXM][MAXM], inm[MAXM], outm[MAXM];
static i64 E[1 << MAXM][MAXM], B[1 << MAXM][MAXM];
static i64 Dv[1 << MAXM], Dw[1 << MAXM][MAXM];
static i64 Mx[MAXM][MAXM], fj[MAXM], gi[MAXM], a[MAXM], mu[MAXM][MAXM], val[1 << MAXM];
static unsigned char qv[1 << MAXM];

static void base_dp(void) {
  for (int S = 1; S <= full; S++) for (int v = 0; v < m; v++) {
    if (!(S >> v & 1)) { E[S][v] = B[S][v] = 0; continue; }
    int R = S ^ (1 << v);
    if (!R) { E[S][v] = B[S][v] = 1; continue; }
    i64 e = 0, b = 0;
    for (int u = 0; u < m; u++) if (R >> u & 1) { if (adj[u][v]) e += E[R][u]; else b += B[R][u]; }
    E[S][v] = e; B[S][v] = b;      /* E: paths on S ending at v ; B: paths on S starting at v */
  }
}
/* direct count of Hamiltonian paths of R + v_X (v beats exactly X): DP over the states that contain v */
static i64 direct_ext(int X) {
  Dv[0] = 1;
  for (int S = 1; S <= full; S++) {
    i64 t = 0; int uu = S & ~X;                      /* u -> v iff u not in X */
    while (uu) { int u = __builtin_ctz(uu); uu &= uu - 1; t += E[S][u]; }
    Dv[S] = t;
    for (int w = 0; w < m; w++) {
      if (!(S >> w & 1)) { Dw[S][w] = 0; continue; }
      int R = S ^ (1 << w); i64 s = 0;
      if (X >> w & 1) s += Dv[R];                    /* v -> w */
      int pp = R & inm[w];
      while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; s += Dw[R][u]; }
      Dw[S][w] = s;
    }
  }
  i64 tot = Dv[full]; for (int w = 0; w < m; w++) tot += Dw[full][w];
  return tot;
}
static int perm_[MAXM], used_, sc[MAXM];
static i64 autcount(int k) {
  if (k == m) return 1; i64 c = 0;
  for (int w = 0; w < m; w++) {
    if (used_ >> w & 1) continue; if (sc[w] != sc[k]) continue; int ok = 1;
    for (int j = 0; j < k && ok; j++) if (adj[j][k] != adj[perm_[j]][w]) ok = 0;
    if (ok) { perm_[k] = w; used_ |= 1 << w; c += autcount(k + 1); used_ &= ~(1 << w); }
  }
  return c;
}
int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: b_thmC m direct|fast|labelled\n"); return 2; }
  m = atoi(argv[1]); full = (1 << m) - 1;
  int direct = !strcmp(argv[2], "direct") || !strcmp(argv[2], "labelled"), labelled = !strcmp(argv[2], "labelled");
  int npairs = m * (m - 1) / 2; i64 mfact = 1; for (int i = 2; i <= m; i++) mfact *= i;
  u64 nt = 0, fail_ends = 0, fail41 = 0, fail42 = 0, fail_cover = 0, fail_even = 0, fail_ge2 = 0, fail_quad = 0, fail_gauss = 0, fail_rad = 0, fail_control = 0;
  i64 minmu = 1LL << 60, maxmu = 0;
  u64 cls[MAXM + 2][2]; i64 lab[MAXM + 2][2]; memset(cls, 0, sizeof cls); memset(lab, 0, sizeof lab);
  i64 signed_lab = 0, signed_plain_lab = 0, labtotal = 0;
  char line[256]; u64 code = 0, ncodes = labelled ? (1ULL << npairs) : 0;
  for (;;) {
    if (labelled) { if (code >= ncodes) break; int k = 0;
      for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { int b = code >> k & 1; k++; adj[i][j] = b; adj[j][i] = !b; }
      code++; }
    else { if (!fgets(line, sizeof line, stdin)) break; if (line[0] != '0' && line[0] != '1') continue; int k = 0;
      for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { int b = (line[k++] == '1'); adj[i][j] = b; adj[j][i] = !b; } }
    for (int i = 0; i < m; i++) { adj[i][i] = 0; inm[i] = outm[i] = 0; }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (adj[i][j]) { outm[i] |= 1 << j; inm[j] |= 1 << i; }
    nt++;
    base_dp();
    i64 H = 0; for (int j = 0; j < m; j++) { fj[j] = E[full][j]; gi[j] = B[full][j]; H += fj[j]; }
    for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue;
      i64 s = 0; for (int S = 1; S < full; S++) if ((S >> j & 1) && !(S >> i & 1)) s += E[S][j] * B[full ^ S][i];
      Mx[j][i] = s; }
    if (direct) {
      for (int X = 0; X <= full; X++) val[X] = direct_ext(X);
      for (int i = 0; i < m; i++) a[i] = val[1 << i] - H;
      for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) mu[i][j] = (i == j) ? 0 : a[i] + a[j] - (val[(1 << i) | (1 << j)] - H);
      for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (i != j && mu[i][j] != Mx[i][j] + Mx[j][i]) fail_cover++;
      for (int X = 0; X <= full; X++) {
        i64 v41 = 0, v42 = H;
        for (int i = 0; i < m; i++) { if (X >> i & 1) { v41 += gi[i]; v42 += a[i]; for (int j = 0; j < i; j++) if (X >> j & 1) v42 -= mu[i][j]; }
                                      else { v41 += fj[i]; for (int t = 0; t < m; t++) if (X >> t & 1) v41 += Mx[i][t]; } }
        if (v41 != val[X]) fail41++; if (v42 != val[X]) fail42++;
      }
    } else {
      for (int i = 0; i < m; i++) { a[i] = gi[i] - fj[i]; for (int j = 0; j < m; j++) if (j != i) a[i] += Mx[j][i]; }
      for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) mu[i][j] = (i == j) ? 0 : Mx[i][j] + Mx[j][i];
      val[0] = H;
      for (int X = 1; X <= full; X++) { int t = 31 - __builtin_clz(X); int Y = X ^ (1 << t); i64 v = val[Y] + a[t];
        int yy = Y; while (yy) { int i = __builtin_ctz(yy); yy &= yy - 1; v -= mu[i][t]; } val[X] = v; }
      int Xc = (int)((nt * 2654435761ULL >> 7) & (u64)full);
      if (direct_ext(Xc) != val[Xc]) fail_control++;
    }
    if (val[0] != H || val[full] != H) fail_ends++;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (mu[i][j] & 1) fail_even++; if (mu[i][j] < 2) fail_ge2++;
      if (mu[i][j] < minmu) minmu = mu[i][j]; if (mu[i][j] > maxmu) maxmu = mu[i][j]; }
    int brow[MAXM]; for (int i = 0; i < m; i++) { brow[i] = 0; for (int j = 0; j < m; j++) if (j != i && ((mu[i][j] / 2) & 1)) brow[i] |= 1 << j; }
    i64 diff = 0;
    for (int X = 0; X <= full; X++) { i64 d = val[X] - H; if (d & 1) { fail_quad++; } qv[X] = (unsigned char)((d / 2) & 1); diff += qv[X] ? -1 : 1; }
    for (int X = 0; X <= full; X++) { int pr = 0; for (int i = 0; i < m; i++) if (X >> i & 1) { pr ^= qv[1 << i]; pr ^= __builtin_parity(brow[i] & X & ((1 << i) - 1)); }
      if (pr != qv[X]) fail_quad++; }
    int rr[MAXM], rk = 0; for (int i = 0; i < m; i++) rr[i] = brow[i];
    for (int c = 0; c < m; c++) { int p = -1; for (int i = rk; i < m; i++) if (rr[i] >> c & 1) { p = i; break; } if (p < 0) continue;
      int t = rr[p]; rr[p] = rr[rk]; rr[rk] = t; for (int i = 0; i < m; i++) if (i != rk && (rr[i] >> c & 1)) rr[i] ^= rr[rk]; rk++; }
    i64 ad = diff < 0 ? -diff : diff;
    if ((rk & 1) || !(ad == 0 || ad == (1LL << (m - rk / 2)))) fail_gauss++;
    int qnonzero_on_rad = 0, radsize = 0;
    for (int Z = 0; Z <= full; Z++) { int inrad = 1; for (int i = 0; i < m && inrad; i++) if (__builtin_parity(brow[i] & Z)) inrad = 0;
      if (inrad) { radsize++; if (qv[Z]) qnonzero_on_rad = 1; } }
    if (radsize != (1 << (m - rk))) fail_rad++;
    if ((ad == 0) != (qnonzero_on_rad != 0)) fail_rad++;
    i64 w;
    if (labelled) w = 1; else { for (int i = 0; i < m; i++) sc[i] = __builtin_popcount(outm[i]); used_ = 0; i64 au = autcount(0); if (mfact % au) { fprintf(stderr, "aut does not divide m!\n"); return 1; } w = mfact / au; }
    cls[rk][ad != 0]++; lab[rk][ad != 0] += w; labtotal += w;
    signed_lab += w * diff * (((H - 1) / 2) & 1 ? -1 : 1); signed_plain_lab += w * diff;
  }
  printf("order %d (%s): tournaments %llu, labelled weight total %lld (2^C(m,2) = %llu)\n", m, argv[2], nt, labtotal, 1ULL << npairs);
  printf("  failures: ends %llu, (4.1) %llu, (4.2) %llu, mu=M+M^T %llu, mu even %llu, mu>=2 %llu, quadratic form %llu, Gauss sum %llu, radical criterion %llu, direct control %llu\n",
         fail_ends, fail41, fail42, fail_cover, fail_even, fail_ge2, fail_quad, fail_gauss, fail_rad, fail_control);
  printf("  min mu %lld, max mu %lld\n", minmu, maxmu);
  double mc = 0, ml = 0;
  for (int r = 0; r <= m; r++) if (cls[r][0] || cls[r][1]) {
    printf("  rank %d: classes zero-bias %llu, bias 2^-%d %llu | labelled zero-bias %lld, bias 2^-%d %lld\n", r, cls[r][0], r / 2, cls[r][1], lab[r][0], r / 2, lab[r][1]);
    mc += (double)cls[r][1] / (double)(1 << (r / 2)); ml += (double)lab[r][1] / (double)(1 << (r / 2)); }
  printf("  mean |imbalance|: classes %.6f ; labelled %.6f ; labelled mean of 2^(-rank/2) = ", mc / nt, ml / labtotal);
  { double s = 0; for (int r = 0; r <= m; r++) s += (double)(lab[r][0] + lab[r][1]) / (double)(1 << (r / 2)); printf("%.6f\n", s / labtotal); }
  printf("  labelled signed sum  sum_R (-1)^((H(R)-1)/2) diff_R = %lld  => P(H=1 mod 4) - P(H=3 mod 4) at order %d = %lld / 2^%d = %.10f\n",
         signed_lab, m + 1, signed_lab, npairs + m, (double)signed_lab / (double)(1ULL << (npairs + m)));
  printf("  labelled unsigned-by-H sum  sum_R diff_R = %lld (mean %.6f)\n", signed_plain_lab, (double)signed_plain_lab / (double)(1ULL << (npairs + m)));
  return 0;
}
