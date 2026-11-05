/* aud_tt.c -- audit S11 of THM-4879 (B) and (C): two-port transitive gadgets TT_k, k <= 9, every (u, z).
   Three independent routes to the transfer matrix:
     (i)   prefix fit: v(X (.) TT_k, z) for X = singleton (1,1,1), TT_2 at its sink (1,0,2), TT_3 at its sink (1,0,4),
           counted by subset DP on the seamed tournament (validated against permutation brute force for k <= 6),
           then M = Vout Vin^{-1}; a fourth prefix, the directed triangle (3,1,3), is used as a hold-out check;
     (ii)  the (A) row statistics counted by permutation brute force on TT_k (ordered path tuples via cut points);
     (iii) the closed form stated in THM-4879 (B).
   Also prints det for every case (C).                                                              */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MV 13
static int Y[MV][MV];
/* ---------- subset DP state of (Y, a) on n vertices ---------- */
static ll dp[1 << MV][MV], ds[1 << MV][MV], Hm[1 << MV];
static void state_dp(int n, int a, ll out[3]){
  int full = (1 << n) - 1;
  for (int m = 0; m <= full; m++) for (int v = 0; v < n; v++) { dp[m][v] = 0; ds[m][v] = 0; }
  for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
  ds[1 << a][a] = 1;
  for (int m = 1; m <= full; m++) for (int v = 0; v < n; v++) {
    ll x = dp[m][v], y = ds[m][v]; if (!x && !y) continue;
    for (int w = 0; w < n; w++) if (!(m >> w & 1) && Y[v][w]) { dp[m | 1 << w][w] += x; ds[m | 1 << w][w] += y; }
  }
  Hm[0] = 1; for (int m = 1; m <= full; m++) { ll s = 0; for (int v = 0; v < n; v++) s += dp[m][v]; Hm[m] = s; }
  ll R = 0; for (int m = 1; m <= full; m++) if (m >> a & 1) R += dp[m][a] * Hm[full ^ m];
  ll S = 0; for (int v = 0; v < n; v++) S += ds[full][v];
  out[0] = Hm[full]; out[1] = S; out[2] = R;
}
/* ---------- permutation brute force state (for validation) ---------- */
static int perm[MV]; static ll bH, bS, bR; static int broot, ny;
static void visit_state(void){
  int n = ny, pre[MV + 1], suf[MV + 1];
  pre[0] = 1; pre[1] = 1; for (int c = 2; c <= n; c++) pre[c] = pre[c - 1] && Y[perm[c - 2]][perm[c - 1]];
  suf[n] = 1; suf[n - 1] = 1; for (int c = n - 2; c >= 0; c--) suf[c] = suf[c + 1] && Y[perm[c]][perm[c + 1]];
  if (pre[n]) { bH++; if (perm[0] == broot) bS++; }
  for (int c = 1; c <= n; c++) if (perm[c - 1] == broot && pre[c] && suf[c]) bR++;
}
static void gen(int k, int n, void (*cb)(void)){
  if (k == n) { cb(); return; }
  for (int i = k; i < n; i++) { int t = perm[k]; perm[k] = perm[i]; perm[i] = t; gen(k + 1, n, cb); t = perm[k]; perm[k] = perm[i]; perm[i] = t; }
}
static void state_brute(int n, int a, ll out[3]){ ny = n; broot = a; bH = bS = bR = 0; for (int i = 0; i < n; i++) perm[i] = i; gen(0, n, visit_state); out[0] = bH; out[1] = bS; out[2] = bR; }
/* ---------- gadget statistics by permutation brute force (all (u,z) in one pass) ---------- */
static int G[MV][MV]; static int ng;
static ll h_, E_[MV], s_[MV], P_[MV][MV], Lu_[MV], Lzu_[MV][MV], Rz_[MV], M_[MV][MV], K1_[MV][MV], K2_[MV][MV];
static void visit_gadget(void){
  int g = ng, brk[MV], okc[MV + 1][MV + 1];
  for (int t = 0; t + 1 < g; t++) brk[t] = G[perm[t]][perm[t + 1]];
  for (int i = 0; i <= g; i++) { okc[i][i] = 1; for (int j = i + 1; j <= g; j++) okc[i][j] = (j == i + 1) ? 1 : (okc[i][j - 1] && brk[j - 2]); }
  if (okc[0][g]) { h_++; E_[perm[g - 1]]++; s_[perm[0]]++; P_[perm[0]][perm[g - 1]]++; }
  for (int c = 0; c <= g - 1; c++) if (okc[0][c] && okc[c][g]) Lu_[perm[c]]++;
  for (int c = 1; c <= g - 1; c++) if (okc[0][c] && okc[c][g]) Lzu_[perm[0]][perm[c]]++;
  for (int c = 1; c <= g; c++) if (okc[0][c] && okc[c][g]) Rz_[perm[c - 1]]++;
  for (int c = 1; c <= g - 1; c++) if (okc[0][c] && okc[c][g]) M_[perm[c - 1]][perm[g - 1]]++;
  for (int c1 = 0; c1 < g; c1++) if (okc[0][c1]) for (int c2 = c1 + 1; c2 <= g; c2++) if (okc[c1][c2] && okc[c2][g]) K1_[perm[c1]][perm[c2 - 1]]++;
  for (int c1 = 1; c1 <= g - 1; c1++) if (okc[0][c1]) for (int c2 = c1; c2 <= g - 1; c2++) if (okc[c1][c2] && okc[c2][g]) K2_[perm[c1 - 1]][perm[c2]]++;
}
static void gadget_stats(int g){
  ng = g; h_ = 0; memset(E_, 0, sizeof E_); memset(s_, 0, sizeof s_); memset(P_, 0, sizeof P_); memset(Lu_, 0, sizeof Lu_);
  memset(Lzu_, 0, sizeof Lzu_); memset(Rz_, 0, sizeof Rz_); memset(M_, 0, sizeof M_); memset(K1_, 0, sizeof K1_); memset(K2_, 0, sizeof K2_);
  for (int i = 0; i < g; i++) perm[i] = i; gen(0, g, visit_gadget);
}
/* ---------- closed form (THM-4879 (B)) ---------- */
static ll pw(ll b, int e){ ll r = 1; while (e-- > 0) r *= b; return r; }
static void closed(int k, int u, int z, ll M[3][3]){
  M[0][0] = 1; M[0][1] = -(u == k - 1); M[0][2] = pw(2, k - 1 - u);
  if (k == 1) { M[1][0] = 1; M[1][1] = -1; M[1][2] = 0; }
  else if (z == 0) { if (u == 0) { M[1][0] = 1; M[1][1] = 0; M[1][2] = 0; } else { M[1][0] = M[0][0]; M[1][1] = M[0][1]; M[1][2] = M[0][2]; } }
  else { if (u == 0) { M[1][0] = 0; M[1][1] = 0; M[1][2] = pw(2, k - 1 - z); } else { M[1][0] = 0; M[1][1] = 0; M[1][2] = 0; } }
  ll mu = 0; if (z == k - 1 && k - 1 > u) mu = pw(2, u); else if (u == k - 1 && k - 1 > z) mu = pw(2, z);
  ll c = (z <= u) ? pw(2, k - 1 - u + z) : pw(2, k - z + u) * pw(3, z - u - 1);
  M[2][0] = pw(2, z); M[2][1] = -mu; M[2][2] = c;
}
static ll det3(ll M[3][3]){ return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]); }
/* seam builder: prefix kind p (0 singleton, 1 TT2 exit sink, 2 TT3 exit sink, 3 directed triangle exit vertex 0) */
static int pre_n[4] = {1, 2, 3, 3};
static int build(int p, int k, int u, int z, int brute, ll out[3]){
  int pn = pre_n[p], n = pn + k, a;
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Y[i][j] = 0;
  if (p == 1) { Y[0][1] = 1; a = 1; }
  else if (p == 2) { Y[0][1] = Y[0][2] = Y[1][2] = 1; a = 2; }
  else if (p == 3) { Y[0][1] = Y[1][2] = Y[2][0] = 1; a = 0; }
  else a = 0;
  for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) Y[pn + i][pn + j] = 1;   /* TT_k: i -> j iff i < j */
  for (int i = 0; i < k; i++) for (int j = 0; j < pn; j++) { Y[pn + i][j] = 1; Y[j][pn + i] = 0; }
  Y[pn + u][a] = 0; Y[a][pn + u] = 1;
  if (brute) state_brute(n, pn + z, out); else state_dp(n, pn + z, out);
  return n;
}
int main(int argc, char **argv){
  int K = argc > 1 ? atoi(argv[1]) : 9;
  /* validate the DP against permutation brute force on seamed TT graphs with k <= 6 and all four prefixes */
  ll nval = 0, badval = 0;
  for (int k = 1; k <= 6; k++) for (int u = 0; u < k; u++) for (int z = 0; z < k; z++) for (int p = 0; p < 4; p++) {
    ll a[3], b[3]; build(p, k, u, z, 0, a); build(p, k, u, z, 1, b); nval++; if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2]) badval++; }
  printf("DP vs permutation brute force on seamed TT graphs (k <= 6, 4 prefixes): %lld cases, %lld mismatches\n", nval, badval);
  /* Vin columns: (1,1,1), (1,0,2), (1,0,4); integer adjugate */
  ll V[3][3] = {{1, 1, 1}, {1, 0, 0}, {1, 2, 4}}, adjV[3][3]; ll dV = det3(V);
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { int r1 = (j + 1) % 3, r2 = (j + 2) % 3, c1 = (i + 1) % 3, c2 = (i + 2) % 3; adjV[i][j] = V[r1][c1] * V[r2][c2] - V[r1][c2] * V[r2][c1]; }
  int cases = 0, bad_i = 0, bad_ii = 0, bad_hold = 0, bad_C = 0;
  for (int k = 1; k <= K; k++) {
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) G[i][j] = (i < j);
    gadget_stats(k);
    for (int u = 0; u < k; u++) for (int z = 0; z < k; z++) {
      ll Vo[3][3], Mi[3][3], Mii[3][3], Mc[3][3];
      for (int p = 0; p < 3; p++) { ll o[3]; build(p, k, u, z, 0, o); for (int r = 0; r < 3; r++) Vo[r][p] = o[r]; }
      for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) { ll s = 0; for (int t = 0; t < 3; t++) s += Vo[r][t] * adjV[t][c]; if (s % dV) { printf("non-integer fit k=%d u=%d z=%d\n", k, u, z); return 1; } Mi[r][c] = s / dV; }
      Mii[0][0] = h_; Mii[0][1] = -E_[u]; Mii[0][2] = Lu_[u];
      Mii[1][0] = s_[z]; Mii[1][1] = -P_[z][u]; Mii[1][2] = Lzu_[z][u];
      Mii[2][0] = Rz_[z]; Mii[2][1] = -M_[z][u]; Mii[2][2] = K1_[u][z] + K2_[z][u];
      closed(k, u, z, Mc);
      /* hold-out: directed triangle prefix (3,1,3) */
      ll ho[3]; build(3, k, u, z, 0, ho); ll pr[3]; for (int r = 0; r < 3; r++) pr[r] = Mc[r][0] * 3 + Mc[r][1] * 1 + Mc[r][2] * 3;
      int di = memcmp(Mi, Mc, sizeof Mi) != 0, dii = memcmp(Mii, Mc, sizeof Mii) != 0, dh = (pr[0] != ho[0] || pr[1] != ho[1] || pr[2] != ho[2]);
      ll d = det3(Mc); int predicted_unimod = (k == 1) || (u == 0 && z == k - 1);
      int dC = predicted_unimod ? (d != 1) : (d != 0);
      cases++; bad_i += di; bad_ii += dii; bad_hold += dh; bad_C += dC;
      if (di || dii || dh || dC) printf("DIFF k=%d u=%d z=%d (fit %d, brute-(A) %d, holdout %d, rank %d)\n", k, u, z, di, dii, dh, dC);
      if (k <= 4 || (k == 9 && (u == 0 || z == 8 || u == z))) {
        printf("k=%d u=%d z=%d  M=[%lld %lld %lld; %lld %lld %lld; %lld %lld %lld] det=%lld\n", k, u, z, Mc[0][0], Mc[0][1], Mc[0][2], Mc[1][0], Mc[1][1], Mc[1][2], Mc[2][0], Mc[2][1], Mc[2][2], d);
      }
    }
  }
  printf("TT_k, k <= %d: %d (k,u,z) cases; closed form vs prefix fit: %d diffs; vs brute-force (A) statistics: %d diffs; "
         "hold-out C3 prefix: %d diffs; rank rule (C) (det 1 iff k=1 or (u,z)=(0,k-1), else det 0): %d violations\n", K, cases, bad_i, bad_ii, bad_hold, bad_C);
  return 0;
}
