/* slab_audit.c -- independent audit of THM-4819 (A) and (B) at small odd orders (n = 5, 7; n <= 8 feasible).
 * Shares no code with the producer.  Canonical forms: lexicographically least upper-triangle code over all
 * relabellings that list the vertices in non-decreasing score order (an isomorphism-invariant set of labellings).
 *
 * Steps: read all classes from gentourng (-q n); build the FULL undirected class move graph (every Hamilton path of every
 * class, native move = reverse every arc not on the path); BFS distances dTT (from the transitive class) and dREG (from
 * the regular classes).  Checks:
 *   (A1)  L + R >= c for every class, and = c iff every sorted score s_(i) lies between i and m;
 *   (B8)  dREG = R for every class (inherited, re-checked);  dTT >= L for every class (THM-4810, re-checked);
 *   (A2)  dTT(regular) = c for the regular classes and every class on a TT-REG geodesic (dTT + dREG = c) is in the slab;
 *   (A3)  inside the slab: {on a geodesic} = {dTT = L} = {reachable from TT by R-descending moves};
 *   (A4)  every class with a source or sink has R >= m and no L-descending move, and none with L >= 1 is on a geodesic.
 * Prints the (B) table: per level L, slab classes / geodesic classes / failing (source or sink, strong with lambda = 1,
 * other), lambda = min over nonempty proper X of #arcs X -> V\X.
 * Build: cc -O2 -o slab_audit slab_audit.c     Run: gentourng -q 7 | ./slab_audit 7
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXC 8000
static int n, m, c;
typedef struct { uint8_t o[16]; } Tour;  /* n <= 8: out-neighbour bitmask per vertex */

static Tour cls[MAXC]; static uint32_t code_of[MAXC]; static int ncls;
static int sc[MAXC][16];
static int adjc[MAXC][700], nadj[MAXC];

static uint32_t canon(const Tour *T) {
    int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(T->o[v]);
    /* order vertices by score; enumerate permutations inside equal-score cells */
    int ord[16], k = 0;
    for (int val = 0; val < n; val++) for (int v = 0; v < n; v++) if (s[v] == val) ord[k++] = v;
    int cellstart[16];
    for (int i = 0; i < n; i++) cellstart[i] = (i > 0 && s[ord[i]] == s[ord[i - 1]]) ? cellstart[i - 1] : i;
    uint32_t best = 0xffffffffu;
    /* iterative generation: Heap's algorithm per cell is fiddly; use recursion */
    int cur[16]; int usedv[16] = {0};
    /* recursive lambda substitute */
    int depth = 0; int choice[16];
    (void)depth; (void)choice;
    /* simple recursion via explicit function */
    {
        int stackpos[17]; int pos = 0; stackpos[0] = -1;
        /* candidates for position pos: vertices in the cell of position pos not used */
        while (pos >= 0) {
            if (pos == n) {
                uint32_t code = 0; int b = 0;
                for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (T->o[cur[i]] >> cur[j] & 1) code |= 1u << b; b++; }
                if (code < best) best = code;
                pos--; if (pos >= 0) usedv[cur[pos]] = 0; continue;
            }
            int cs = cellstart[pos], ce = cs; while (ce < n && cellstart[ce] == cs) ce++;
            int next = -1;
            for (int idx = stackpos[pos] + 1; idx < ce - cs; idx++) { int v = ord[cs + idx]; if (!usedv[v]) { next = idx; break; } }
            if (next < 0) { stackpos[pos] = -1; pos--; if (pos >= 0) usedv[cur[pos]] = 0; continue; }
            stackpos[pos] = next; cur[pos] = ord[cs + next]; usedv[cur[pos]] = 1; pos++; stackpos[pos] = -1;
        }
    }
    return best;
}
static int find_class(uint32_t code) {
    int lo = 0, hi = ncls - 1;
    while (lo <= hi) { int mid = (lo + hi) / 2; if (code_of[mid] == code) return mid; if (code_of[mid] < code) lo = mid + 1; else hi = mid - 1; }
    return -1;
}
static int cmpu(const void *a, const void *b) { uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b; return x < y ? -1 : x > y; }

static Tour native(const Tour *T, const int *Q) {
    Tour R; uint8_t full = (uint8_t)((1u << n) - 1);
    for (int u = 0; u < n; u++) R.o[u] = (uint8_t)(~T->o[u] & full & ~(1u << u));
    for (int i = 0; i + 1 < n; i++) { int a = Q[i], b = Q[i + 1]; R.o[a] |= (uint8_t)(1u << b); R.o[b] &= (uint8_t)~(1u << a); }
    return R;
}
static int Lof(const int *s0) { int s[16]; memcpy(s, s0, sizeof s); for (int i = 1; i < n; i++) for (int j = i; j > 0 && s[j] < s[j-1]; j--) { int t = s[j]; s[j] = s[j-1]; s[j-1] = t; } int l = 0; for (int i = 0; i < n; i++) l += abs(s[i] - i); return l / 2; }
static int Rof(const int *s) { int r = 0; for (int v = 0; v < n; v++) r += abs(s[v] - m); return r / 2; }
static int between_ok(const int *s0) { int s[16]; memcpy(s, s0, sizeof s); for (int i = 1; i < n; i++) for (int j = i; j > 0 && s[j] < s[j-1]; j--) { int t = s[j]; s[j] = s[j-1]; s[j-1] = t; }
    for (int i = 0; i < n; i++) { int lo = i < m ? i : m, hi = i < m ? m : i; if (s[i] < lo || s[i] > hi) return 0; } return 1; }
static int lambda_of(const Tour *T) {
    int best = 99;
    for (uint32_t X = 1; X + 1 < (1u << n); X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T->o[x] & ~X & ((1u << n) - 1)); if (f < best) best = f; }
    return best;
}

static int curQ[16], cid;
static void add_nb(int t) { for (int i = 0; i < nadj[cid]; i++) if (adjc[cid][i] == t) return; adjc[cid][nadj[cid]++] = t; }
static long long nmoves;
static int Lc[MAXC], Rc[MAXC];
static long long ldesc_from_srcsink;
static void paths_rec(const Tour *T, int depth, uint32_t used) {
    if (depth == n) {
        Tour R = native(T, curQ); int t = find_class(canon(&R));
        if (t < 0) { fprintf(stderr, "unknown class\n"); exit(3); }
        nmoves++; add_nb(t);
        return;
    }
    int last = curQ[depth - 1];
    for (int w = 0; w < n; w++) if (!(used >> w & 1) && (T->o[last] >> w & 1)) { curQ[depth] = w; paths_rec(T, depth + 1, used | 1u << w); }
}
static int dist[2][MAXC];
static void bfs(int which, int (*src)(int)) {
    static int q[MAXC]; int h = 0, tl = 0;
    for (int i = 0; i < ncls; i++) { dist[which][i] = -1; if (src(i)) { dist[which][i] = 0; q[tl++] = i; } }
    while (h < tl) { int u = q[h++]; for (int k = 0; k < nadj[u]; k++) { int t = adjc[u][k]; if (dist[which][t] < 0) { dist[which][t] = dist[which][u] + 1; q[tl++] = t; } } }
}
static int is_TT(int i) { return Lc[i] == 0; }
static int is_REG(int i) { return Rc[i] == 0; }

int main(int argc, char **argv) {
    n = atoi(argv[1]); m = (n - 1) / 2; c = (n * n - 1) / 8;
    char line[256]; Tour tmp[MAXC]; int nt = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.o[i] |= (uint8_t)(1u << j); else T.o[j] |= (uint8_t)(1u << i); k++; }
        tmp[nt++] = T;
    }
    uint32_t codes[MAXC]; for (int i = 0; i < nt; i++) codes[i] = canon(&tmp[i]);
    /* sort classes by code; check distinctness (gentourng gives one per class) */
    int idx[MAXC]; for (int i = 0; i < nt; i++) idx[i] = i;
    for (int i = 1; i < nt; i++) for (int j = i; j > 0 && codes[idx[j]] < codes[idx[j-1]]; j--) { int t = idx[j]; idx[j] = idx[j-1]; idx[j-1] = t; }
    ncls = nt; for (int i = 0; i < nt; i++) { cls[i] = tmp[idx[i]]; code_of[i] = codes[idx[i]]; }
    for (int i = 1; i < ncls; i++) if (code_of[i] == code_of[i-1]) { fprintf(stderr, "duplicate class from gentourng?\n"); return 4; }
    (void)cmpu;
    for (int i = 0; i < ncls; i++) { for (int v = 0; v < n; v++) sc[i][v] = __builtin_popcount(cls[i].o[v]); Lc[i] = Lof(sc[i]); Rc[i] = Rof(sc[i]); }
    for (cid = 0; cid < ncls; cid++) {
        for (int s = 0; s < n; s++) { curQ[0] = s; paths_rec(&cls[cid], 1, 1u << s); }
    }
    bfs(0, is_TT); bfs(1, is_REG);
    int a1 = 1, eqchar = 1, b8 = 1, lowL = 1, a2 = 1, a3 = 1, a4 = 1, connected = 1;
    for (int i = 0; i < ncls; i++) {
        if (dist[0][i] < 0 || dist[1][i] < 0) connected = 0;
        if (Lc[i] + Rc[i] < c) a1 = 0;
        if ((Lc[i] + Rc[i] == c) != between_ok(sc[i])) eqchar = 0;
        if (dist[1][i] != Rc[i]) b8 = 0;
        if (dist[0][i] >= 0 && dist[0][i] < Lc[i]) lowL = 0;
        if (Rc[i] == 0 && dist[0][i] != c) a2 = 0;
        if (dist[0][i] + dist[1][i] == c && Lc[i] + Rc[i] != c) a2 = 0;
    }
    /* R-descending closure from TT */
    static int reach[MAXC]; static int q[MAXC]; int h = 0, tl = 0;
    for (int i = 0; i < ncls; i++) { reach[i] = 0; if (Lc[i] == 0) { reach[i] = 1; q[tl++] = i; } }
    while (h < tl) { int u = q[h++]; for (int k = 0; k < nadj[u]; k++) { int t = adjc[u][k]; if (Rc[t] == Rc[u] - 1 && !reach[t]) { reach[t] = 1; q[tl++] = t; } } }
    int geo_total = 0, slab_total = 0;
    for (int i = 0; i < ncls; i++) {
        int slab = (Lc[i] + Rc[i] == c);
        int g1 = (dist[0][i] + dist[1][i] == c), g2 = slab && dist[0][i] == Lc[i], g3 = reach[i];
        if (slab) slab_total++;
        if (g1) geo_total++;
        if (slab && !(g1 == g2 && g2 == g3)) a3 = 0;
        if (!slab && (g1 || g3)) a3 = 0;
        int ss = 0; for (int v = 0; v < n; v++) if (sc[i][v] == 0 || sc[i][v] == n - 1) ss = 1;
        if (ss) {
            if (Rc[i] < m) a4 = 0;
            for (int k = 0; k < nadj[i]; k++) if (Lc[adjc[i][k]] < Lc[i]) a4 = 0;
            if (Lc[i] >= 1 && g1) a4 = 0;
        }
    }
    printf("n=%d: classes %d; slab classes %d; on a TT-REG geodesic %d; move graph connected (all classes reach TT and REG): %d; native moves enumerated %lld\n",
           n, ncls, slab_total, geo_total, connected, nmoves);
    printf("  (A1) L+R >= c for all classes: %d; equality iff s_(i) between i and m: %d\n", a1, eqchar);
    printf("  (B8) dREG = R for all classes: %d;  dTT >= L for all classes: %d\n", b8, lowL);
    printf("  (A2) regular classes at dTT = c and geodesic classes inside the slab: %d\n", a2);
    printf("  (A3) slab: on-geodesic = (dTT = L) = R-descending reachable from TT: %d\n", a3);
    printf("  (A4) source/sink classes: R >= m, no L-descending move, never on a geodesic when L >= 1: %d\n", a4);
    printf("  level (L,R): slab / geodesic (failing: source or sink, strong lambda = 1, other)\n");
    for (int Lv = 0; Lv <= c; Lv++) {
        int tot = 0, geo = 0, ss = 0, l1 = 0, ot = 0, Rmin_fail = 99;
        for (int i = 0; i < ncls; i++) {
            if (Lc[i] + Rc[i] != c || Lc[i] != Lv) continue;
            tot++;
            if (dist[0][i] + dist[1][i] == c) { geo++; continue; }
            int s_ = 0; for (int v = 0; v < n; v++) if (sc[i][v] == 0 || sc[i][v] == n - 1) s_ = 1;
            if (s_) ss++; else if (lambda_of(&cls[i]) == 1) l1++; else ot++;
            if (Rc[i] < Rmin_fail) Rmin_fail = Rc[i];
        }
        printf("   (%d,%d): %d / %d (%d, %d, %d)\n", Lv, c - Lv, tot, geo, ss, l1, ot);
    }
    int ok_obs = 1;
    for (int i = 0; i < ncls; i++) if (Lc[i] + Rc[i] == c && Rc[i] <= m - 1 && dist[0][i] + dist[1][i] != c) ok_obs = 0;
    printf("  every slab class with R <= m-1 is on a geodesic: %d\n", ok_obs);
    /* details of the 'other' failures */
    for (int i = 0; i < ncls; i++) {
        if (Lc[i] + Rc[i] != c || dist[0][i] + dist[1][i] == c) continue;
        int s_ = 0; for (int v = 0; v < n; v++) if (sc[i][v] == 0 || sc[i][v] == n - 1) s_ = 1;
        if (s_ || lambda_of(&cls[i]) == 1) continue;
        int ss2[16]; memcpy(ss2, sc[i], sizeof ss2);
        printf("   'other' failure: level (%d,%d), dTT %d, lambda %d, scores", Lc[i], Rc[i], dist[0][i], lambda_of(&cls[i]));
        for (int v = 0; v < n; v++) printf(" %d", ss2[v]);
        printf("\n");
    }
    return 0;
}
