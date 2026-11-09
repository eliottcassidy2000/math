/* tight_audit.c - independent audit of the tight-line counts of THM-4951 (4) (written from scratch by the auditor).
 *
 * A line {t, t + 1} is tight iff nu(t) + nu(t + 1) = k = floor(n/3).  By THM-4951 (1) (equality case of
 * nu >= P(A), nu(t + 1) >= P(B)), the triple word x_a = [tile (a, a + 2) reversed] of the end with x_0 = 1 is then a
 * tight colouring and nu(t) = P(A), nu(t + 1) = P(B).  Conversely a line whose word is tight and with nu(t) <= P(A),
 * nu(t + 1) <= P(B) is tight.  So the auditor enumerates, for every TIGHT colouring x of [0, n-3] with x_0 = 1, the tilings with that triple word, vertex by vertex,
 * pruning with the LOOKAHEAD bound
 *     nu(T[0..v]) + P(A cap [v+1, n-3]) <= P(A)     and the same for u = t + 1 and B,
 * which is valid because the windows {a, a+1, a+2} with a >= v + 1 avoid [0, v].
 * This design differs from the packet's tight_lines.c (which searches all words at once with the bound
 * nu(t) + nu(u) <= k on prefixes and no lookahead).
 * Output: per word (tight words only appear), counts black/blue; totals by nu-pair and colour.
 * usage: tight_audit n [SLICE NSLICES]   (slices deal the words round-robin)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 16
typedef long long ll;

static int n, L, K;
static uint32_t ot[MAXN], ou[MAXN];
static int bit[MAXN][MAXN];
static int word[MAXN], PA, PB, sufA[MAXN + 2], sufB[MAXN + 2];
static ll cntw[2];
static ll pairc[8][8][2];

static int P_of(const int *x, int from, int to, int colour) { /* greedy: max 3-separated subset of {a in [from,to]: x_a = colour} */
  int c = 0, last = -100;
  for (int a = from; a <= to; a++) if (x[a] == colour && a - last >= 3) { c++; last = a; }
  return c;
}

/* nu on vertices {0..V-1}: triangles grouped by highest vertex; branch on the highest free vertex */
static uint32_t tri[MAXN][256]; static int ntri[MAXN];
static int bestp;
static void pack(uint32_t freev, int cur, int cap) {
  if (cur > bestp) bestp = cur;
  if (bestp >= cap || !freev) return;
  if (cur + __builtin_popcount(freev) / 3 <= bestp) return;
  int h = 31 - __builtin_clz(freev);
  for (int k = 0; k < ntri[h]; k++) if ((tri[h][k] & freev) == tri[h][k]) pack(freev & ~tri[h][k], cur + 1, cap);
  pack(freev & ~(1u << h), cur, cap);
}
static int nu(const uint32_t *out, int V, int cap) {
  uint32_t in[MAXN];
  for (int v = 0; v < V; v++) in[v] = 0;
  for (int v = 0; v < V; v++) { uint32_t o = out[v] & ((1u << V) - 1); while (o) { int u = __builtin_ctz(o); o &= o - 1; in[u] |= 1u << v; } }
  for (int h = 0; h < V; h++) {
    ntri[h] = 0;
    uint32_t low = (1u << h) - 1, A = out[h] & low;
    while (A) { int a = __builtin_ctz(A); A &= A - 1; uint32_t B = out[a] & in[h] & low;
      while (B) { int b = __builtin_ctz(B); B &= B - 1; tri[h][ntri[h]++] = (1u << h) | (1u << a) | (1u << b); } }
  }
  bestp = 0; pack((1u << V) - 1, 0, cap);
  return bestp;
}

static void dfs(int v) {
  if (v == n) {
    int a = nu(ot, n, PA + 1), b = nu(ou, n, PB + 1);
    if (a != PA || b != PB) return;           /* cannot happen after the pruning at v = n - 1; kept as a guard */
    int blue = 1;
    for (int j = 2; j < n && blue; j++) for (int i = 0; i + 2 <= j; i++) if (bit[i][j] != bit[n - 1 - j][n - 1 - i]) { blue = 0; break; }
    cntw[blue]++;
    int lo = a < b ? a : b, hi = a < b ? b : a;
    pairc[lo][hi][blue]++;
    return;
  }
  int free_tiles = v - 2;                     /* tiles (i, v), i = 0..v-3 free; (v-2, v) fixed by the word */
  for (uint32_t mask = 0; mask < (1u << free_tiles); mask++) {
    ot[v] = ou[v] = 0;
    for (int i = 0; i < v; i++) { ot[i] &= ~(1u << v); ou[i] &= ~(1u << v); }
    ot[v - 1] |= 1u << v; ou[v - 1] |= 1u << v;
    for (int i = 0; i <= v - 2; i++) {
      int r = (i == v - 2) ? word[v - 2] : (int)((mask >> i) & 1u);
      bit[i][v] = r;
      if (r) { ot[v] |= 1u << i; ou[i] |= 1u << v; } else { ot[i] |= 1u << v; ou[v] |= 1u << i; }
    }
    int la = sufA[v + 1], lb = sufB[v + 1];
    int a = nu(ot, v + 1, PA - la + 1);
    if (a + la > PA) continue;
    int b = nu(ou, v + 1, PB - lb + 1);
    if (b + lb > PB) continue;
    dfs(v + 1);
  }
  ot[v] = ou[v] = 0;
  for (int i = 0; i < v; i++) { ot[i] &= ~(1u << v); ou[i] &= ~(1u << v); }
}

int main(int argc, char **argv) {
  n = atoi(argv[1]);
  int slice = argc > 3 ? atoi(argv[2]) : 0, ns = argc > 3 ? atoi(argv[3]) : 1;
  L = n - 2; K = n / 3;
  ll tot[2] = {0, 0}; int ntight = 0, nrealised = 0, npal = 0, npal_blue = 0;
  int idx = 0;
  for (uint32_t w = 0; w < (1u << L); w++) {
    if (!(w & 1u)) continue;                  /* x_0 = 1 */
    for (int a = 0; a < L; a++) word[a] = (w >> a) & 1u;
    PA = P_of(word, 0, L - 1, 1); PB = P_of(word, 0, L - 1, 0);
    int tight = (PA + PB == K);
    int pal = 1;
    for (int a = 0; a < L; a++) if (word[a] != word[L - 1 - a]) pal = 0;
    if (!tight) continue;                     /* THM-4951 (1): no tight line has a non-tight word (cross-checked at n <= 9 against the full census of audit_lines.c) */
    if ((idx++ % ns) != slice) continue;
    for (int v = 0; v <= L + 2; v++) { sufA[v] = (v <= L - 1) ? P_of(word, v, L - 1, 1) : 0; sufB[v] = (v <= L - 1) ? P_of(word, v, L - 1, 0) : 0; }
    /* vertices 0, 1, 2 with tile (0, 2) reversed (x_0 = 1) */
    memset(ot, 0, sizeof ot); memset(ou, 0, sizeof ou); memset(bit, 0, sizeof bit);
    ot[0] = ou[0] = 1u << 1; ot[1] = ou[1] = 1u << 2;
    ot[2] |= 1u << 0; ou[0] |= 1u << 2; bit[0][2] = 1;
    cntw[0] = cntw[1] = 0;
    dfs(3);
    ntight++;
    if (pal) npal++;
    if (cntw[0] + cntw[1]) nrealised++;
    if (pal && cntw[1]) npal_blue++;
    tot[0] += cntw[0]; tot[1] += cntw[1];
    printf("word ");
    for (int a = 0; a < L; a++) putchar('0' + word[a]);
    printf(" P(A)=%d P(B)=%d %s black %lld blue %lld\n", PA, PB, pal ? "palindrome" : "", cntw[0], cntw[1]);
    fflush(stdout);
  }
  printf("n=%d slice %d/%d: tight colourings (x_0 = 1) %d, realised %d, palindromic %d (realised by a blue line %d); "
         "tight lines %lld (blue %lld)\n",
         n, slice, ns, ntight, nrealised, npal, npal_blue, tot[0] + tot[1], tot[1]);
  for (int i = 0; i < 8; i++) for (int j = i; j < 8; j++) for (int c = 0; c < 2; c++)
    if (pairc[i][j][c]) printf("  pair (%d,%d) %s %lld\n", i, j, c ? "blue" : "black", pairc[i][j][c]);
  return 0;
}
