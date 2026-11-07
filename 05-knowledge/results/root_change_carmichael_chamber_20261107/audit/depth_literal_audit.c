/* depth_literal_audit.c -- auditor's own code (audit of kind-pasteur-2026-11-07-S14).
 * Part 1: exact root-change depth with NO cap on the number of letters J.
 *   depth(N) = min #J(w) over all words w over {O,Z,J} with F(t.w) = N (F = D-E+U strictly increases along every
 *   letter, so the DFS pruned at F > X is exhaustive). Also one witness word of minimal depth (shortest, then
 *   lexicographically first in DFS order O < Z < J).
 * Part 2: literal tournaments. Builder (Session200 carmichael-fibers verify.py, re-implemented): start with vertices
 *   0 -> 1 and root 0; a letter adds a vertex n that is beaten by every old vertex, except that for Z and J it beats
 *   the current root; J then makes n the root. A word w closed by J gives the C = 1 candidate.
 *   Hamiltonian paths H by subset DP; Hamiltonian cycles C by subset DP from vertex 0.
 *   (a) every word w of length <= LMAX: H(graph(w)) = D(t.w), H(graph(wJ)) = F(t.w), C(graph(wJ)) = 1, and
 *       d(t.w) = H(graph(w) - root);
 *   (b) the targets: minimal depth, witness, literal (H, C).
 * usage: depth_literal_audit XDEPTH LMAX target...
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
typedef struct { ll D, E, U, d; } st;
static inline st mv(st s, char c){
  st r;
  if (c == 'O'){ r.D = s.D; r.E = 0; r.U = 2*s.U; r.d = s.d; }
  else if (c == 'Z'){ r.D = s.D - s.E + s.U; r.E = s.d; r.U = 2*s.U - s.d; r.d = s.d; }
  else { r.D = s.D - s.E + s.U; r.E = s.D - s.E; r.U = s.D + s.U; r.d = s.D; }
  return r;
}
static inline ll FF(st s){ return s.D - s.E + s.U; }

static ll XD; static unsigned char *dep; static char **wit; static int *witlen;
static char wb[4096]; static ll nodes = 0, incr_fail = 0;
static void dfs(st s, int j, int len){
  ll F = FF(s);
  if (F > XD) return;
  nodes++;
  ll i = F >> 1;
  if (j < dep[i] || (j == dep[i] && len < witlen[i])){
    dep[i] = (unsigned char)j; witlen[i] = len;
    free(wit[i]); wit[i] = malloc(len + 1); memcpy(wit[i], wb, len); wit[i][len] = 0;
  }
  if (len >= 4000){ fprintf(stderr, "too deep\n"); exit(2); }
  const char *al = "OZJ";
  for (int a = 0; a < 3; a++){
    st q = mv(s, al[a]);
    if (FF(q) <= F) incr_fail++;
    if (j == 254 && al[a] == 'J') continue;
    wb[len] = al[a]; dfs(q, j + (al[a] == 'J'), len + 1);
  }
}

/* literal graphs: adjacency as bitmasks out[v] (bit u set iff v -> u) */
static int build(const char *w, uint32_t *out){
  int n = 2, root = 0; out[0] = 2; out[1] = 0;
  for (const char *p = w; *p; p++){
    int v = n++; out[v] = 0;
    for (int u = 0; u < v; u++){
      if (*p != 'O' && u == root) out[v] |= 1u << u; else out[u] |= 1u << v;
    }
    if (*p == 'J') root = v;
  }
  return n | (root << 8);
}
static ll *dp = NULL; static size_t dpcap = 0;
static void need_dp(int n){ size_t sz = ((size_t)1 << n) * n; if (sz > dpcap){ free(dp); dp = malloc(sz * sizeof(ll)); dpcap = sz; } }
static ll ham_paths(const uint32_t *out, int n, uint32_t mask){  /* paths covering exactly the vertex set 'mask' */
  need_dp(n);
  size_t full = (size_t)1 << n;
  memset(dp, 0, full * n * sizeof(ll));
  for (int v = 0; v < n; v++) if (mask >> v & 1) dp[((size_t)1 << v) * n + v] = 1;
  ll tot = 0;
  for (size_t S = 1; S < full; S++){
    if ((S & ~(size_t)mask)) continue;
    for (int v = 0; v < n; v++){ ll x = dp[S * n + v]; if (!x) continue;
      if (S == mask){ tot += x; continue; }
      uint32_t o = out[v] & mask & ~(uint32_t)S;
      while (o){ int u = __builtin_ctz(o); o &= o - 1; dp[(S | ((size_t)1 << u)) * n + u] += x; } }
  }
  return tot;
}
static ll ham_cycles(const uint32_t *out, int n){
  need_dp(n);
  size_t full = (size_t)1 << n;
  memset(dp, 0, full * n * sizeof(ll));
  dp[1 * n + 0] = 1;
  ll tot = 0;
  for (size_t S = 1; S < full; S += 2){
    for (int v = 0; v < n; v++){ ll x = dp[S * n + v]; if (!x) continue;
      if (S == full - 1){ if (v != 0 && (out[v] & 1)) tot += x; continue; }
      uint32_t o = out[v] & ~(uint32_t)S;
      while (o){ int u = __builtin_ctz(o); o &= o - 1; dp[(S | ((size_t)1 << u)) * n + u] += x; } }
  }
  return tot;
}
static st state_of(const char *w){ st s = {1, 0, 2, 1}; for (const char *p = w; *p; p++) s = mv(s, *p); return s; }

static ll lw_n = 0, lw_fail = 0; static int LMAX;
static void literal_words(char *w, int len){
  if (len > 0 || 1){
    st s = state_of(w);
    uint32_t out[40]; int r = build(w, out); int n = r & 255, root = r >> 8;
    ll H = ham_paths(out, n, (1u << n) - 1);
    ll hd = ham_paths(out, n, ((1u << n) - 1) & ~(1u << root));
    char wj[64]; strcpy(wj, w); strcat(wj, "J");
    uint32_t out2[40]; int r2 = build(wj, out2); int n2 = r2 & 255;
    ll H2 = ham_paths(out2, n2, (1u << n2) - 1), C2 = ham_cycles(out2, n2);
    lw_n++;
    if (H != s.D || hd != s.d || H2 != FF(s) || C2 != 1){ lw_fail++; if (lw_fail < 5) fprintf(stderr, "literal mismatch %s: H %lld D %lld | Hroot %lld d %lld | H(wJ) %lld F %lld C %lld\n", w, H, s.D, hd, s.d, H2, FF(s), C2); }
  }
  if (len == LMAX) return;
  const char *al = "OZJ";
  for (int a = 0; a < 3; a++){ w[len] = al[a]; w[len + 1] = 0; literal_words(w, len + 1); w[len] = 0; }
}

int main(int argc, char **argv){
  if (argc < 3){ fprintf(stderr, "usage: depth_literal_audit XDEPTH LMAX target...\n"); return 1; }
  XD = atoll(argv[1]); LMAX = atoi(argv[2]);
  ll H = XD / 2 + 2;
  dep = malloc(H); memset(dep, 255, H); wit = calloc(H, sizeof(char*)); witlen = malloc(H * sizeof(int));
  for (ll i = 0; i < H; i++) witlen[i] = 1 << 30;
  st t = {1, 0, 2, 1};
  dfs(t, 0, 0);
  dep[0] = 0;   /* N = 1: the empty counted word */
  printf("Part 1: all words (any number of J) with F(t.w) <= %lld: %lld prefixes visited; F-increase failures %lld\n", XD, nodes, incr_fail);
  ll cnt[256] = {0}; ll none = 0;
  printf("odd N <= %lld with no word at all:", XD);
  for (ll N = 1; N <= XD; N += 2){ unsigned char e = dep[N >> 1]; if (e == 255){ none++; printf(" %lld", N); } else cnt[e]++; }
  printf("\n  count %lld\n", none);
  printf("depth histogram (number of odd N <= %lld with exact minimal #J = k, N = 1 counted at k = 0):", XD);
  for (int k = 0; k < 255; k++){ if (cnt[k]) printf(" k=%d:%lld", k, cnt[k]); }
  printf("\n");
  printf("values of depth >= 4:"); for (ll N = 3; N <= XD; N += 2){ unsigned char e = dep[N >> 1]; if (e != 255 && e >= 4) printf(" %lld(d%d)", N, e); } printf("\n");
  /* Part 2a */
  char w[64] = {0};
  literal_words(w, 0);
  printf("Part 2a: literal builder vs state machine on all %lld words of length 0..%d: H(w) = D, H(w - root) = d, H(wJ) = F, C(wJ) = 1; failures %lld\n", lw_n, LMAX, lw_fail);
  /* Part 2b */
  printf("Part 2b: target, exact minimal #J (no cap), my witness w (closed by J), vertices, literal (H, C); then the packet's witness\n");
  for (int a = 3; a < argc; a++){
    char *tok = argv[a]; char *col = strchr(tok, ':');
    ll N = atoll(tok);
    if (N > XD){ printf("%lld beyond XDEPTH\n", N); continue; }
    unsigned char e = dep[N >> 1];
    if (e == 255){ printf("%lld: no word at all\n", N); continue; }
    char ww[256]; snprintf(ww, sizeof ww, "%sJ", wit[N >> 1]);
    uint32_t out[40]; int r = build(ww, out); int n = r & 255;
    ll Hh = ham_paths(out, n, (1u << n) - 1), C = ham_cycles(out, n);
    printf("%lld depth %d my witness %s (%d vertices) literal (H, C) = (%lld, %lld)", N, e, wit[N >> 1], n, Hh, C);
    if (col){ char pw[256]; snprintf(pw, sizeof pw, "%sJ", col + 1);
      st ps = state_of(col + 1); int nj = 0; for (char *p = col + 1; *p; p++) nj += (*p == 'J');
      uint32_t o2[40]; int r2 = build(pw, o2); int n2 = r2 & 255;
      ll H2 = ham_paths(o2, n2, (1u << n2) - 1), C2 = ham_cycles(o2, n2);
      printf(" | packet witness %s: #J %d, F %lld, %d vertices, literal (H, C) = (%lld, %lld)", col + 1, nj, FF(ps), n2, H2, C2); }
    printf("\n");
  }
  return 0;
}
