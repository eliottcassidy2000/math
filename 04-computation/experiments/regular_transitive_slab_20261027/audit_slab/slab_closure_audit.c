/* slab_closure_audit.c -- independent re-derivation of the THM-4819 (B) slab census at odd orders n <= 11 by forward
 * closure from TT under R-descending native moves (the characterisation (A3)(iii), which slab_audit.c checks against
 * full class-graph BFS distances at orders 5 and 7).  Shares no code with the producer; in particular the canonical
 * form is the author's own individualisation-refinement (no nauty):
 *   colour by score, refine to the coarsest equitable colouring (signature = own colour + out-neighbour counts per
 *   colour, new colours = ranks of sorted signatures), and if not discrete individualise every vertex of the first
 *   non-singleton colour in turn; the code is the least upper-triangle bit string over all leaves.
 * Steps: read every class (gentourng -q n); keep the slab classes (L + R = c); closure from TT level by level using every
 * Hamilton path from a vertex of score > m to a vertex of score < m (exactly the R-descending moves); every image must be
 * a slab class one level up.  Prints per level: slab / reached / unreached split (source or sink, strong with lambda = 1,
 * other).  Build: cc -O2 -o slab_closure_audit slab_closure_audit.c   Run: gentourng -q 9 | ./slab_closure_audit 9
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, m, c;
typedef struct { uint16_t o[16]; } Tour;

/* ---------- canonical form ---------- */
static uint64_t best_code;
static void refine(const Tour *T, int *col) {
    int k_old = -1;
    for (;;) {
        int k = 0; for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
        if (k == k_old) return;
        k_old = k;
        /* signature: col, then counts of out-neighbours in each colour */
        int sig[16][17];
        for (int v = 0; v < n; v++) {
            sig[v][0] = col[v];
            for (int j = 0; j < k; j++) sig[v][1 + j] = 0;
            for (int w = 0; w < n; w++) if (T->o[v] >> w & 1) sig[v][1 + col[w]]++;
        }
        /* rank distinct signatures */
        int order[16]; for (int v = 0; v < n; v++) order[v] = v;
        for (int i = 1; i < n; i++) for (int j = i; j > 0; j--) {
            int a = order[j - 1], b = order[j], cmp = 0;
            for (int t = 0; t <= k && !cmp; t++) cmp = (sig[a][t] > sig[b][t]) - (sig[a][t] < sig[b][t]);
            if (cmp > 0) { order[j - 1] = b; order[j] = a; } else break;
        }
        int newcol[16], r = 0;
        for (int i = 0; i < n; i++) {
            if (i > 0) { int a = order[i - 1], b = order[i], same = 1; for (int t = 0; t <= k; t++) if (sig[a][t] != sig[b][t]) { same = 0; break; } if (!same) r++; }
            newcol[order[i]] = r;
        }
        memcpy(col, newcol, sizeof(int) * n);
    }
}
static void search(const Tour *T, int *col) {
    refine(T, col);
    int k = 0; for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
    if (k == n) {
        int at[16]; for (int v = 0; v < n; v++) at[col[v]] = v;
        uint64_t code = 0; int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (T->o[at[i]] >> at[j] & 1) code |= 1ull << b; b++; }
        if (code < best_code) best_code = code;
        return;
    }
    /* first colour with >= 2 vertices */
    int cnt[16] = {0}; for (int v = 0; v < n; v++) cnt[col[v]]++;
    int target = 0; while (cnt[target] < 2) target++;
    for (int w = 0; w < n; w++) {
        if (col[w] != target) continue;
        int c2[16];
        for (int v = 0; v < n; v++) c2[v] = 2 * col[v] + 1;
        c2[w] = 2 * col[w];
        /* compress to ranks */
        int used[40] = {0}; for (int v = 0; v < n; v++) used[c2[v]] = 1;
        int rk[40], r = 0; for (int x = 0; x < 40; x++) { rk[x] = r; if (used[x]) r++; }
        for (int v = 0; v < n; v++) c2[v] = rk[c2[v]];
        search(T, c2);
    }
}
static uint64_t canon(const Tour *T) {
    int col[16]; for (int v = 0; v < n; v++) col[v] = __builtin_popcount(T->o[v]);
    /* compress score colours to ranks */
    int used[17] = {0}; for (int v = 0; v < n; v++) used[col[v]] = 1;
    int rk[17], r = 0; for (int x = 0; x < 17; x++) { rk[x] = r; if (used[x]) r++; }
    for (int v = 0; v < n; v++) col[v] = rk[col[v]];
    best_code = ~0ull; search(T, col); return best_code;
}

/* ---------- hash map code -> id ---------- */
#define HBITS 22
static uint64_t hkey[1u << HBITS]; static int hval[1u << HBITS]; static char hused[1u << HBITS];
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return x; }
static int hput(uint64_t k, int v) { uint32_t h = (uint32_t)(mix(k) & ((1u << HBITS) - 1)); while (hused[h]) { if (hkey[h] == k) return 0; h = (h + 1) & ((1u << HBITS) - 1); } hused[h] = 1; hkey[h] = k; hval[h] = v; return 1; }
static int hget(uint64_t k) { uint32_t h = (uint32_t)(mix(k) & ((1u << HBITS) - 1)); while (hused[h]) { if (hkey[h] == k) return hval[h]; h = (h + 1) & ((1u << HBITS) - 1); } return -1; }

static Tour *reps; static int *lev, *reached, nslab;
static int Lof(const int *s0) { int s[16]; memcpy(s, s0, sizeof s); for (int i = 1; i < n; i++) for (int j = i; j > 0 && s[j] < s[j-1]; j--) { int t = s[j]; s[j] = s[j-1]; s[j-1] = t; } int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static int Rof(const int *s) { int r = 0; for (int v = 0; v < n; v++) r += abs(s[v] - m); return r / 2; }
static Tour native(const Tour *T, const int *Q) {
    Tour R; uint16_t full = (uint16_t)((1u << n) - 1);
    for (int u = 0; u < n; u++) R.o[u] = (uint16_t)(~T->o[u] & full & ~(1u << u));
    for (int i = 0; i + 1 < n; i++) { int a = Q[i], b = Q[i + 1]; R.o[a] |= (uint16_t)(1u << b); R.o[b] &= (uint16_t)~(1u << a); }
    return R;
}
static int lambda_of(const Tour *T) {
    int best = 99;
    for (uint32_t X = 1; X + 1 < (1u << n); X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T->o[x] & ~X & ((1u << n) - 1)); if (f < best) best = f; }
    return best;
}
static int *nextf, nnext, curL; static long long moves;
static int Q[16], sc_[16];
static void paths(const Tour *T, int depth, uint32_t used) {
    if (depth == n) {
        if (sc_[Q[n - 1]] >= m) return;
        moves++;
        Tour R = native(T, Q); int t = hget(canon(&R));
        if (t < 0) { fprintf(stderr, "image outside the slab (or canonical-form bug)\n"); exit(3); }
        if (lev[t] != curL + 1) { fprintf(stderr, "image not one level up\n"); exit(3); }
        if (!reached[t]) { reached[t] = 1; nextf[nnext++] = t; }
        return;
    }
    int last = Q[depth - 1];
    for (int w = 0; w < n; w++) if (!(used >> w & 1) && (T->o[last] >> w & 1)) { Q[depth] = w; paths(T, depth + 1, used | 1u << w); }
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); m = (n - 1) / 2; c = (n * n - 1) / 8;
    int cap = 1 << 21; reps = malloc(sizeof(Tour) * cap); lev = malloc(sizeof(int) * cap); reached = calloc(cap, sizeof(int)); nextf = malloc(sizeof(int) * cap);
    char line[256]; long all = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        all++;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.o[i] |= (uint16_t)(1u << j); else T.o[j] |= (uint16_t)(1u << i); k++; }
        int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T.o[v]);
        int L = Lof(s), R = Rof(s);
        if (L + R < c) { fprintf(stderr, "slab inequality violated\n"); return 5; }
        if (L + R != c) continue;
        if (!hput(canon(&T), nslab)) { fprintf(stderr, "two gentourng classes share a code: canonical-form bug\n"); return 4; }
        reps[nslab] = T; lev[nslab] = L; nslab++;
    }
    int *front = malloc(sizeof(int) * cap), nf = 0;
    for (int i = 0; i < nslab; i++) if (lev[i] == 0) { reached[i] = 1; front[nf++] = i; }
    int lvcount[64] = {0};
    for (curL = 0; curL < c; curL++) {
        nnext = 0;
        for (int f = 0; f < nf; f++) {
            const Tour *T = &reps[front[f]];
            for (int v = 0; v < n; v++) sc_[v] = __builtin_popcount(T->o[v]);
            for (int a = 0; a < n; a++) if (sc_[a] > m) { Q[0] = a; paths(T, 1, 1u << a); }
        }
        memcpy(front, nextf, sizeof(int) * nnext); nf = nnext; lvcount[curL + 1] = nnext;
    }
    printf("n=%d: classes %ld; slab classes %d; R-descending moves explored %lld\n", n, all, nslab, moves);
    long tot_reached = 0;
    for (int L = 0; L <= c; L++) {
        int tot = 0, re = 0, ss = 0, l1 = 0, ot = 0;
        for (int i = 0; i < nslab; i++) if (lev[i] == L) {
            tot++; if (reached[i]) { re++; continue; }
            int src = 0; for (int v = 0; v < n; v++) { int s = __builtin_popcount(reps[i].o[v]); if (s == 0 || s == n - 1) src = 1; }
            if (src) ss++; else if (lambda_of(&reps[i]) == 1) l1++; else ot++;
        }
        tot_reached += re;
        printf("   (%d,%d): %d / %d (%d, %d, %d)\n", L, c - L, tot, re, ss, l1, ot);
    }
    printf("  geodesic classes in total: %ld\n", tot_reached);
    return 0;
}
