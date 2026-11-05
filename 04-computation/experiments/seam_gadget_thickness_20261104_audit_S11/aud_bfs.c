/* aud_bfs.c -- audit S11: re-implementation (own code) of the two BFS procedures described in HYP-9471 and the note.
   seam  X MATFILE... : 3D BFS of SEAM states (H,S,R) under the DISTINCT matrices read from the files (deduplicated),
                        seeds (1,1,1),(1,0,2),(3,1,3),(1,0,4),(1,0,8),(1,0,16); a state is kept iff H <= X and R <= X;
                        an odd H <= X counts as attained if any image of a kept state has that H (whatever its R).
                        Prints: #kept states, the missed odd H <= X, the minimum slope R/H over kept S = 0 states,
                        the two band fills (fraction of lattice points H odd, R even, 2H <= R <= 2.75H realised as
                        kept S = 0 states) for H in [801,999] and [1601,1799], and their ratio.
   plane X SEEDFILE MATFILE... : pure plane BFS of (H,R) under the distinct nonsingular 2x2 restrictions
                        [[M00,M02],[M20,M22]] of the S-row-zero matrices; seeds from SEEDFILE ("H R" per line);
                        same box and bands; also slope statistics of the states with H in [801,1799] and the average
                        of D(t) = sum (a + b t)^-2 over the maps, at the slopes t = R/H of those states.
   Completeness of the box truncation: for every gadget H' = hH - E_u S + L_u R >= R (since E_u <= h, S <= H,
   L_u >= 1) and H' >= H (since R >= S); so every ancestor of a state with H, R <= X has H, R <= X.            */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long long ll; typedef unsigned long long ull;
static ll (*MM)[9]; static int nm;
static int cmp9(const void *a, const void *b){ const ll *x = a, *y = b; for (int i = 0; i < 9; i++) { if (x[i] < y[i]) return -1; if (x[i] > y[i]) return 1; } return 0; }
static int cmp4(const void *a, const void *b){ const ll *x = a, *y = b; for (int i = 0; i < 4; i++) { if (x[i] < y[i]) return -1; if (x[i] > y[i]) return 1; } return 0; }
static void load(int nf, char **fs){
  int cap = 600000; MM = malloc(sizeof(ll) * 9 * cap); nm = 0;
  for (int i = 0; i < nf; i++) { FILE *f = fopen(fs[i], "r"); if (!f) { fprintf(stderr, "cannot open %s\n", fs[i]); exit(1); }
    ll m[9]; while (fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld", &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]) == 9) { memcpy(MM[nm++], m, sizeof m); } fclose(f); }
  qsort(MM, nm, sizeof(ll) * 9, cmp9); int k = 0; for (int i = 0; i < nm; i++) if (k == 0 || cmp9(MM[i], MM[k - 1])) { memcpy(MM[k], MM[i], sizeof(ll) * 9); k++; } nm = k;
}
/* packed key: 21 bits each */
static ull *tab; static size_t cap;
static inline ull key(ll h, ll s, ll r){ return ((ull)h << 42) | ((ull)s << 21) | (ull)r; }
static inline size_t hh(ull k){ k ^= k >> 31; k *= 0x9E3779B97F4A7C15ULL; k ^= k >> 29; return (size_t)k & (cap - 1); }
static int ins(ull k){ size_t i = hh(k); while (tab[i]) { if (tab[i] == k) return 0; i = (i + 1) & (cap - 1); } tab[i] = k; return 1; }
static int has(ull k){ size_t i = hh(k); while (tab[i]) { if (tab[i] == k) return 1; i = (i + 1) & (cap - 1); } return 0; }
int main(int argc, char **argv){
  if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
  ll X = atoll(argv[2]);
  if (!strcmp(argv[1], "seam")) {
    load(argc - 3, argv + 3);
    cap = 1ULL << 26; tab = calloc(cap, sizeof(ull)); ull *q = malloc(sizeof(ull) * (cap / 2)); size_t qh = 0, qt = 0, cnt = 0;
    char *hit = calloc(X + 1, 1);
    ll seeds[6][3] = {{1, 1, 1}, {1, 0, 2}, {3, 1, 3}, {1, 0, 4}, {1, 0, 8}, {1, 0, 16}};
    for (int i = 0; i < 6; i++) { ull k = key(seeds[i][0], seeds[i][1], seeds[i][2]); if (seeds[i][0] <= X && seeds[i][2] <= X && ins(k)) { q[qt++] = k; cnt++; hit[seeds[i][0]] = 1; } }
    while (qh < qt) { ull k = q[qh++]; ll H = k >> 42, S = (k >> 21) & 0x1FFFFF, R = k & 0x1FFFFF;
      for (int i = 0; i < nm; i++) { ll *A = MM[i]; ll h = A[0] * H + A[1] * S + A[2] * R; if (h > X) continue;
        if (h <= 0) { fprintf(stderr, "nonpositive H\n"); return 1; }
        hit[h] = 1; ll r = A[6] * H + A[7] * S + A[8] * R; if (r > X) continue; ll s = A[3] * H + A[4] * S + A[5] * R;
        if (s < 0 || r < 0 || s > h) { fprintf(stderr, "invalid state\n"); return 1; }
        ull kk = key(h, s, r); if (ins(kk)) { if (qt >= cap / 2) { fprintf(stderr, "queue full\n"); return 1; } q[qt++] = kk; cnt++; } } }
    printf("seam BFS X=%lld distinct matrices=%d kept states=%zu\n", X, nm, cnt);
    ll miss = 0; printf("missed odd H <= X:"); for (ll n = 1; n <= X; n += 2) if (!hit[n]) { miss++; if (miss <= 40) printf(" %lld", n); } printf("  (count %lld)\n", miss);
    double minsl = 1e9; ll np = 0; for (size_t i = 0; i < qt; i++) { ull k = q[i]; ll H = k >> 42, S = (k >> 21) & 0x1FFFFF, R = k & 0x1FFFFF; if (S == 0) { np++; double t = (double)R / H; if (t < minsl) minsl = t; } }
    double minsl_band = 1e9; for (size_t i = 0; i < qt; i++) { ull k = q[i]; ll H = k >> 42, S = (k >> 21) & 0x1FFFFF, R = k & 0x1FFFFF; if (S == 0 && H >= 801) { double t = (double)R / H; if (t < minsl_band) minsl_band = t; } }
    printf("kept S=0 states %lld, min slope R/H %.4f (over those with H >= 801: %.4f)\n", np, minsl, minsl_band);
    if (X >= 4948) { double fl[2]; ll lo[2] = {801, 1601}, hi[2] = {999, 1799};
      for (int b = 0; b < 2; b++) { ll tot = 0, have = 0; for (ll H = lo[b]; H <= hi[b]; H += 2) { ll r1 = (ll)floor(2.75 * H); for (ll R = 2 * H; R <= r1; R += 2) { tot++; have += has(key(H, 0, R)); } } fl[b] = (double)have / tot;
        printf("band H in [%lld,%lld], R in [2H,2.75H]: lattice points %lld, realised %lld, fill %.4f\n", lo[b], hi[b], tot, have, fl[b]); }
      printf("ratio fill2/fill1 = %.4f\n", fl[1] / fl[0]); }
    return 0;
  }
  if (!strcmp(argv[1], "plane")) {
    load(argc - 4, argv + 4);
    ll (*A)[4] = malloc(sizeof(ll) * 4 * nm); int na = 0;
    for (int i = 0; i < nm; i++) { ll *m = MM[i]; if (m[3] || m[4] || m[5]) continue; if (m[0] * m[8] - m[2] * m[6] == 0) continue; A[na][0] = m[0]; A[na][1] = m[2]; A[na][2] = m[6]; A[na][3] = m[8]; na++; }
    qsort(A, na, sizeof(ll) * 4, cmp4); int k = 0; for (int i = 0; i < na; i++) if (k == 0 || cmp4(A[i], A[k - 1])) { memcpy(A[k], A[i], sizeof(ll) * 4); k++; } na = k;
    size_t W = (size_t)(X + 1) * (size_t)(X + 1); unsigned char *seen = calloc(W / 8 + 1, 1);
#define IDX(h, r) ((size_t)(h) * (size_t)(X + 1) + (size_t)(r))
#define GET(id) ((seen[(id) >> 3] >> ((id) & 7)) & 1)
#define SET(id) (seen[(id) >> 3] |= (unsigned char)(1 << ((id) & 7)))
    ll *qq = malloc(sizeof(ll) * 2 * 40000000); size_t qs = 0, qe = 0;
    FILE *sf = fopen(argv[3], "r"); ll h0, r0; int nseed = 0;
    while (fscanf(sf, "%lld %lld", &h0, &r0) == 2) { nseed++; if (h0 <= X && r0 <= X && !GET(IDX(h0, r0))) { SET(IDX(h0, r0)); qq[2 * qe] = h0; qq[2 * qe + 1] = r0; qe++; } } fclose(sf);
    while (qs < qe) { ll h = qq[2 * qs], r = qq[2 * qs + 1]; qs++;
      for (int i = 0; i < na; i++) { ll h2 = A[i][0] * h + A[i][1] * r; if (h2 > X) continue; ll r2 = A[i][2] * h + A[i][3] * r; if (r2 > X) continue;
        if (r2 < 0 || h2 <= 0 || (r2 & 1) || !(h2 & 1)) { fprintf(stderr, "parity/sign violation\n"); return 1; }
        size_t id = IDX(h2, r2); if (GET(id)) continue; SET(id); qq[2 * qe] = h2; qq[2 * qe + 1] = r2; qe++; } }
    printf("plane BFS X=%lld seeds=%d maps=%d states=%zu\n", X, nseed, na, qe);
    double fl[2]; ll lo[2] = {801, 1601}, hi[2] = {999, 1799};
    for (int b = 0; b < 2; b++) { ll tot = 0, have = 0; for (ll H = lo[b]; H <= hi[b]; H += 2) { ll r1 = (ll)floor(2.75 * H); for (ll R = 2 * H; R <= r1; R += 2) { tot++; have += GET(IDX(H, R)); } } fl[b] = (double)have / tot;
      printf("band H in [%lld,%lld]: lattice points %lld realised %lld fill %.4f\n", lo[b], hi[b], tot, have, fl[b]); }
    printf("ratio = %.4f\n", fl[1] / fl[0]);
    /* slope statistics over states with 801 <= H <= 1799 */
    ll ns = 0; ll histo[64] = {0}; double *sl = malloc(sizeof(double) * qe); double Dsum = 0, Dmin = 1e9; ll nbelow = 0;
    for (size_t i = 0; i < qe; i++) { ll h = qq[2 * i], r = qq[2 * i + 1]; if (h < 801 || h > 1799) continue; double t = (double)r / h; sl[ns++] = t; int fb = (int)floor(t); if (fb < 64) histo[fb]++;
      double D = 0; for (int j = 0; j < na; j++) { double g = A[j][0] + A[j][1] * t; D += 1.0 / (g * g); } Dsum += D; if (D < Dmin) Dmin = D; if (D < 1) nbelow++; }
    int cmpd(const void *a, const void *b){ double x = *(const double *)a, y = *(const double *)b; return x < y ? -1 : x > y; }
    qsort(sl, ns, sizeof(double), cmpd);
    printf("states with H in [801,1799]: %lld; slope min %.4f median %.4f p90 %.4f max %.4f; histogram floor(t):", ns, ns ? sl[0] : 0, ns ? sl[ns / 2] : 0, ns ? sl[(ll)(0.9 * ns)] : 0, ns ? sl[ns - 1] : 0);
    for (int b = 0; b < 64; b++) if (histo[b]) printf(" %d:%lld", b, histo[b]); printf("\n");
    printf("average over those states of D(t) for this alphabet: %.4f (states with D(t) < 1: %lld)\n", ns ? Dsum / ns : 0, nbelow);
    const char *dump = getenv("AUD_DUMP");
    if (dump) { FILE *o = fopen(dump, "w"); for (size_t i = 0; i < qe; i++) { ll h = qq[2 * i], r = qq[2 * i + 1]; if (h >= 801 && h <= 1799) fprintf(o, "%lld %lld\n", h, r); } fclose(o);
      FILE *o2 = fopen("aud_plane_maps_dump.txt", "w"); for (int j = 0; j < na; j++) fprintf(o2, "%lld %lld %lld %lld\n", A[j][0], A[j][1], A[j][2], A[j][3]); fclose(o2); }
    return 0;
  }
  return 1;
}
