/* a09_gadgets.c (S20 audit): independent one-arc gadget matrices by brute force.
   For a gadget (G, u, z), M is defined by (h,e,R)(N, z) = M (h,e,R)(T, v), N = T => G with the arc u -> v.
   We build N for five bases (T, v) (single vertex; TT_2 at its source; 3-cycle; TT_3 at its source; the strong
   4-tournament), count (h, e_z, R_z) of N by a subset DP, solve the 3x3 system on the first three bases (Cramer,
   integrality checked) and verify linearity on the other two.
   usage: a09_gadgets labelled GMAX            -> all labelled gadgets of order 1..GMAX
          a09_gadgets classes GMAX CLSDIR      -> class representatives cls_g.txt, g = 1..GMAX, all (u, z)
   output lines: g code u z m00 m01 m02 m10 m11 m12 m20 m21 m22 ; stderr: counts, failures */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MAXN 12
static int A[MAXN][MAXN];
static ll f[1 << MAXN][MAXN], s[1 << MAXN][MAXN], Hs[1 << MAXN];
/* triple (H, e_z, R_z) for every z of tournament A on n vertices */
static void stats(int n, ll *H, ll *e, ll *R){
  int N = 1 << n;
  Hs[0] = 1;
  for (int S = 1; S < N; S++) {
    ll tot = 0;
    for (int v = 0; v < n; v++) {
      if (!(S >> v & 1)) continue;
      int P = S & ~(1 << v); ll a = 0, b = 0;
      if (!P) { a = b = 1; }
      else for (int w = 0; w < n; w++) if (P >> w & 1) { if (A[w][v]) a += f[P][w]; if (A[v][w]) b += s[P][w]; }
      f[S][v] = a; s[S][v] = b; tot += a;
    }
    Hs[S] = tot;
  }
  *H = Hs[N - 1];
  for (int z = 0; z < n; z++) {
    e[z] = f[N - 1][z]; ll r = 0;
    for (int S = 0; S < N; S++) if (S >> z & 1) r += s[S][z] * Hs[(N - 1) & ~S];
    R[z] = r;
  }
}
static int nb = 5, bn[5], bv[5], bA[5][4][4];
static ll btrip[5][3];
static void setbases(void){
  memset(bA, 0, sizeof bA);
  bn[0] = 1; bv[0] = 0;
  bn[1] = 2; bv[1] = 0; bA[1][0][1] = 1;
  bn[2] = 3; bv[2] = 0; bA[2][0][1] = 1; bA[2][1][2] = 1; bA[2][2][0] = 1;
  bn[3] = 3; bv[3] = 0; bA[3][0][1] = bA[3][0][2] = bA[3][1][2] = 1;
  /* strong 4-tournament: 0->1->2->3->0, 0->2, 1->3 */
  bn[4] = 4; bv[4] = 2; bA[4][0][1] = bA[4][1][2] = bA[4][2][3] = bA[4][3][0] = bA[4][0][2] = bA[4][1][3] = 1;
  for (int b = 0; b < nb; b++) {
    memset(A, 0, sizeof A);
    for (int i = 0; i < bn[b]; i++) for (int j = 0; j < bn[b]; j++) A[i][j] = bA[b][i][j];
    ll H, e[MAXN], R[MAXN]; stats(bn[b], &H, e, R);
    btrip[b][0] = H; btrip[b][1] = e[bv[b]]; btrip[b][2] = R[bv[b]];
    fprintf(stderr, "base %d: (h,e,R) = (%lld,%lld,%lld)\n", b, H, btrip[b][0] ? e[bv[b]] : 0, R[bv[b]]);
  }
}
static ll det3(ll a[3][3]){
  return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
static long nfail = 0, nout = 0;
static int G[MAXN][MAXN];
static void do_gadget(int g, unsigned long long code){
  /* out[b][z][k] for every entry u */
  static ll out[MAXN][5][MAXN][3];
  for (int u = 0; u < g; u++) for (int b = 0; b < nb; b++) {
    int t = bn[b], n = t + g;
    memset(A, 0, sizeof A);
    for (int i = 0; i < t; i++) for (int j = 0; j < t; j++) A[i][j] = bA[b][i][j];
    for (int i = 0; i < g; i++) for (int j = 0; j < g; j++) A[t + i][t + j] = G[i][j];
    for (int i = 0; i < t; i++) for (int k = 0; k < g; k++) A[i][t + k] = 1;
    A[bv[b]][t + u] = 0; A[t + u][bv[b]] = 1;
    ll H, e[MAXN], R[MAXN]; stats(n, &H, e, R);
    for (int z = 0; z < g; z++) { out[u][b][z][0] = H; out[u][b][z][1] = e[t + z]; out[u][b][z][2] = R[t + z]; }
  }
  ll X[3][3]; for (int b = 0; b < 3; b++) for (int k = 0; k < 3; k++) X[b][k] = btrip[b][k];
  ll D = det3(X);
  for (int u = 0; u < g; u++) for (int z = 0; z < g; z++) {
    ll M[3][3]; int ok = 1;
    for (int row = 0; row < 3; row++) {
      for (int col = 0; col < 3; col++) {
        ll Y[3][3]; memcpy(Y, X, sizeof Y);
        for (int b = 0; b < 3; b++) Y[b][col] = out[u][b][z][row];
        ll num = det3(Y);
        if (num % D) ok = 0;
        M[row][col] = num / D;
      }
      for (int b = 3; b < nb; b++) {
        ll pred = M[row][0] * btrip[b][0] + M[row][1] * btrip[b][1] + M[row][2] * btrip[b][2];
        if (pred != out[u][b][z][row]) ok = 0;
      }
    }
    if (!ok) { nfail++; fprintf(stderr, "FAIL g=%d code=%llu u=%d z=%d\n", g, code, u, z); }
    printf("%d %llu %d %d", g, code, u, z);
    for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) printf(" %lld", M[r][c]);
    printf("\n"); nout++;
  }
}
static void decode(int g, unsigned long long code){
  memset(G, 0, sizeof G); int k = 0;
  for (int i = 0; i < g; i++) for (int j = i + 1; j < g; j++) { if (code >> k & 1) G[i][j] = 1; else G[j][i] = 1; k++; }
}
int main(int argc, char **argv){
  setbases();
  int GMAX = atoi(argv[2]);
  if (!strcmp(argv[1], "labelled")) {
    for (int g = 1; g <= GMAX; g++) {
      int np = g * (g - 1) / 2;
      for (unsigned long long code = 0; code < (1ULL << np); code++) { decode(g, code); do_gadget(g, code); }
    }
  } else {
    for (int g = 1; g <= GMAX; g++) {
      char fn[1024]; snprintf(fn, sizeof fn, "%s/cls_%d.txt", argv[3], g);
      FILE *fp = fopen(fn, "r"); if (!fp) { perror(fn); return 1; }
      char line[4096];
      while (fgets(line, sizeof line, fp)) { unsigned long long code; if (sscanf(line, "%llu", &code) != 1) continue; decode(g, code); do_gadget(g, code); }
      fclose(fp);
    }
  }
  fprintf(stderr, "gadgets written %ld, failures %ld\n", nout, nfail);
  return 0;
}
