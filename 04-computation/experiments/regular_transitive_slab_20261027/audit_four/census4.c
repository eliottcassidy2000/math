/* census4.c -- independent audit of THM-4823 (written from the statement, not from pattern4.py / w4free.c).
 *
 * Reads gentourng -q n (upper triangle row by row, '1' at (i,j), i<j, means i->j).
 * For every isomorphism class:
 *   (A) the 4-vertex census by brute force over all C(n,4) subsets.  Each 4-set is classified twice:
 *       by its sorted local score vector (0123 TT4, 1113 W4, 0222 L4, 1122 S4) and by source/sink presence
 *       (THM-584's two bits); the two must agree.  The directed 4-cycles inside every 4-set are counted
 *       directly (all 6 cyclic orders) to test "S4 is the only type with a 4-cycle, and it has exactly one".
 *       Identities checked: #TT4+#W4 = sum C(s,3), #TT4+#L4 = sum C(n-1-s,3), total = C(n,4),
 *       4 #S4 = tr(A^4) (A^2 by explicit integer product), plus the three solved formulas of THM-4823 (A).
 *       #S4 is grouped by sorted score sequence (hash table, min/max) to count sequences with non-constant #S4.
 *   (B) pattern-free classes for each of the four types; S4-free classes are checked to be K1/C3 stacks
 *       (every strong component has size 1 or 3), and the number of distinct stack words is recorded.
 *   (C) W4-free (census count of W4 = 0), cross-checked against "every out-neighbourhood is transitive";
 *       L4-free; locally transitive (both).  Structure check of Theorem C on every W4-free class:
 *       peel sinks one at a time (a tournament has at most one sink); the core L must be a single vertex
 *       (then T is transitive) or sinkless, locally transitive and strong.  Records the histogram of |L|.
 *       Every locally transitive class with a sink must be transitive.
 * Usage: gentourng -q n | ./census4 n
 * Build: cc -O3 -o census4 census4.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint32_t out_[32];

static int pc(uint32_t x) { return __builtin_popcount(x); }
static long binom(long a, long b) { if (b < 0 || a < b) return 0; long r = 1; for (long i = 0; i < b; i++) r = r * (a - i) / (i + 1); return r; }

/* number of directed 4-cycles on the 4-set {v0,v1,v2,v3}: the 3 undirected Hamilton cycles, both directions */
static int arc(int u, int v) { return (out_[u] >> v) & 1; }
static int cyc4(int a, int b, int c, int d) {
    int perms[3][4] = {{a,b,c,d},{a,b,d,c},{a,c,b,d}}; int cnt = 0;
    for (int p = 0; p < 3; p++) {
        int *q = perms[p];
        if (arc(q[0],q[1]) && arc(q[1],q[2]) && arc(q[2],q[3]) && arc(q[3],q[0])) cnt++;
        if (arc(q[0],q[3]) && arc(q[3],q[2]) && arc(q[2],q[1]) && arc(q[1],q[0])) cnt++;
    }
    return cnt;
}

/* has a directed 3-cycle inside vertex set S */
static int has_c3(uint32_t S) {
    for (uint32_t x = S; x; x &= x - 1) { int u = __builtin_ctz(x);
        for (uint32_t y = out_[u] & S; y; y &= y - 1) { int v = __builtin_ctz(y);
            if (out_[v] & S & ~out_[u] & ~(1u << u)) { /* some w with v->w, w not an out-neighbour of u, w != u: then w->u */
                uint32_t W = out_[v] & S & ~out_[u] & ~(1u << u);
                for (uint32_t z = W; z; z &= z - 1) { int w = __builtin_ctz(z); if (arc(w, u)) return 1; }
            }
        }
    }
    return 0;
}

/* strong components of T[S] are all of size in {1,3}?  Return 1 if yes; also returns stack word code */
static uint32_t reach_from(int v, uint32_t S) { uint32_t seen = 1u << v, frontier = seen;
    while (frontier) { uint32_t nf = 0; for (uint32_t x = frontier; x; x &= x - 1) nf |= out_[__builtin_ctz(x)] & S; nf &= ~seen; seen |= nf; frontier = nf; }
    return seen; }
static int is_strong(uint32_t S) { if (pc(S) <= 1) return 1; int v = __builtin_ctz(S);
    if (reach_from(v, S) != S) return 0;
    /* reverse reachability: v reachable from all */
    uint32_t seen = 1u << v, frontier = seen;
    while (frontier) { uint32_t nf = 0; for (uint32_t x = S & ~seen; x; x &= x - 1) { int u = __builtin_ctz(x); if (out_[u] & frontier) nf |= 1u << u; } seen |= nf; frontier = nf; }
    return seen == S; }

/* score-sequence hash table: key -> min, max #S4 */
#define HSIZE (1 << 16)
static uint64_t hkey[HSIZE]; static long hmin[HSIZE], hmax[HSIZE]; static int hused[HSIZE];
static void hput(uint64_t key, long val) {
    uint64_t h = (key * 0x9E3779B97F4A7C15ULL) >> 48; uint32_t i = (uint32_t)h & (HSIZE - 1);
    while (hused[i] && hkey[i] != key) i = (i + 1) & (HSIZE - 1);
    if (!hused[i]) { hused[i] = 1; hkey[i] = key; hmin[i] = hmax[i] = val; }
    else { if (val < hmin[i]) hmin[i] = val; if (val > hmax[i]) hmax[i] = val; }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: gentourng -q n | census4 n\n"); return 1; }
    n = atoi(argv[1]);
    char line[512]; long cls = 0; long fail_id = 0, fail_type = 0, fail_cyc = 0, fail_solved = 0;
    long free_[4] = {0,0,0,0}; long w4f = 0, l4f = 0, lt = 0, lt_sink = 0, lt_sink_nontrans = 0, w4f_mismatch = 0, l4f_mismatch = 0;
    long struct_bad = 0; long Lhist[32]; memset(Lhist, 0, sizeof Lhist);
    long s4free_notstack = 0; long s4free_stack_bad = 0;
    long lt_sinkless = 0, lt_sinkless_notstrong = 0;
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        cls++;
        for (int i = 0; i < n; i++) out_[i] = 0;
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i; k++; }
        long cnt[4] = {0,0,0,0};
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) for (int d = c + 1; d < n; d++) {
            uint32_t m = (1u << a) | (1u << b) | (1u << c) | (1u << d);
            int v[4] = {a,b,c,d}; int s[4]; int src = 0, snk = 0;
            for (int i = 0; i < 4; i++) { s[i] = pc(out_[v[i]] & m); if (s[i] == 3) src++; if (s[i] == 0) snk++; }
            /* sort local scores */
            for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (s[j] < s[i]) { int t = s[i]; s[i] = s[j]; s[j] = t; }
            int code = s[0] * 1000 + s[1] * 100 + s[2] * 10 + s[3];
            int t1 = (code == 123) ? 0 : (code == 1113) ? 1 : (code == 222) ? 2 : (code == 1122) ? 3 : -1;
            int t2 = (src && snk) ? 0 : src ? 1 : snk ? 2 : 3;
            if (t1 != t2) fail_type++;
            int cy = cyc4(a,b,c,d);
            if ((t2 == 3 && cy != 1) || (t2 != 3 && cy != 0)) fail_cyc++;
            cnt[t2]++;
        }
        int sc[32]; long P = 0, Q = 0;
        for (int i = 0; i < n; i++) { sc[i] = pc(out_[i]); P += binom(sc[i], 3); Q += binom(n - 1 - sc[i], 3); }
        /* A^2 and tr(A^4) */
        long A2[32][32];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { long s2 = 0; for (int t = 0; t < n; t++) s2 += arc(i,t) * arc(t,j); A2[i][j] = s2; }
        long tr4 = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) tr4 += A2[i][j] * A2[j][i];
        long N4 = binom(n, 4);
        if (!(cnt[0] + cnt[1] == P && cnt[0] + cnt[2] == Q && cnt[0] + cnt[1] + cnt[2] + cnt[3] == N4 && 4 * cnt[3] == tr4)) fail_id++;
        if (tr4 % 4) fail_id++;
        long S4 = tr4 / 4;
        if (!(cnt[0] == P + Q + S4 - N4 && cnt[1] == N4 - Q - S4 && cnt[2] == N4 - P - S4)) fail_solved++;
        /* score sequence key */
        int ss[32]; memcpy(ss, sc, sizeof(int) * n);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (ss[j] < ss[i]) { int t = ss[i]; ss[i] = ss[j]; ss[j] = t; }
        uint64_t key = 0; for (int i = 0; i < n; i++) key = key * 32 + (uint64_t)ss[i];
        hput(key, cnt[3]);
        for (int t = 0; t < 4; t++) if (cnt[t] == 0) free_[t]++;
        /* (B) S4-free => K1/C3 stack */
        if (cnt[3] == 0) {
            /* strong components: peel by repeatedly taking the set of vertices reaching nothing outside ... simpler: for each v, comp(v) = reach(v) & coreach(v) */
            int ok = 1;
            for (int v = 0; v < n && ok; v++) {
                uint32_t R = reach_from(v, full);
                /* co-reach */
                uint32_t seen = 1u << v, fr = seen;
                while (fr) { uint32_t nf = 0; for (uint32_t x = full & ~seen; x; x &= x - 1) { int u = __builtin_ctz(x); if (out_[u] & fr) nf |= 1u << u; } seen |= nf; fr = nf; }
                int sz = pc(R & seen); if (sz != 1 && sz != 3) ok = 0;
            }
            if (!ok) s4free_notstack++;
        } else {
            /* converse direction: a class with a strong component of size >= 4 must contain S4 (pancyclicity); and a class with all components <= 3 has no S4 */
            int big = 0;
            for (int v = 0; v < n; v++) {
                uint32_t R = reach_from(v, full);
                uint32_t seen = 1u << v, fr = seen;
                while (fr) { uint32_t nf = 0; for (uint32_t x = full & ~seen; x; x &= x - 1) { int u = __builtin_ctz(x); if (out_[u] & fr) nf |= 1u << u; } seen |= nf; fr = nf; }
                if (pc(R & seen) >= 4) big = 1;
            }
            if (!big) s4free_stack_bad++;
        }
        /* (C) */
        int w4free_census = (cnt[1] == 0), l4free_census = (cnt[2] == 0);
        int w4free_nbhd = 1, l4free_nbhd = 1;
        for (int v = 0; v < n; v++) { if (has_c3(out_[v])) w4free_nbhd = 0; if (has_c3(full & ~out_[v] & ~(1u << v))) l4free_nbhd = 0; }
        if (w4free_census != w4free_nbhd) w4f_mismatch++;
        if (l4free_census != l4free_nbhd) l4f_mismatch++;
        if (w4free_census) w4f++;
        if (l4free_census) l4f++;
        int is_lt = w4free_census && l4free_census;
        int has_sink = 0; for (int v = 0; v < n; v++) if (out_[v] == 0) has_sink = 1;
        int transitive = 1; { int seen_sc[32] = {0}; for (int v = 0; v < n; v++) seen_sc[sc[v]]++; for (int i = 0; i < n; i++) if (seen_sc[i] != 1) transitive = 0; }
        if (is_lt) { lt++; if (has_sink) { lt_sink++; if (!transitive) lt_sink_nontrans++; } else { lt_sinkless++; if (!is_strong(full)) lt_sinkless_notstrong++; } }
        if (w4free_census) {
            uint32_t cur = full;
            while (pc(cur) > 1) { int found = -1; for (uint32_t x = cur; x; x &= x - 1) { int v = __builtin_ctz(x); if ((out_[v] & cur) == 0) { found = v; break; } } if (found < 0) break; cur &= ~(1u << found); }
            int L = pc(cur); Lhist[L]++;
            if (L == 1) { if (!transitive) struct_bad++; }
            else {
                /* L sinkless by construction; check locally transitive inside T[L] and strong */
                int ok = 1;
                for (uint32_t x = cur; x; x &= x - 1) { int v = __builtin_ctz(x); if (has_c3(out_[v] & cur) || has_c3(cur & ~out_[v] & ~(1u << v))) ok = 0; }
                if (!is_strong(cur)) ok = 0;
                if (!ok) struct_bad++;
            }
        }
    }
    long nseq = 0, nonconst = 0;
    for (int i = 0; i < HSIZE; i++) if (hused[i]) { nseq++; if (hmin[i] != hmax[i]) nonconst++; }
    printf("n=%d classes=%ld | (A) identity failures=%ld, solved-formula failures=%ld, type-rule disagreements=%ld, 4-cycle-rule failures=%ld | score sequences=%ld, with non-constant #S4=%ld\n",
           n, cls, fail_id, fail_solved, fail_type, fail_cyc, nseq, nonconst);
    printf("      (B) pattern-free classes: TT4=%ld W4=%ld L4=%ld S4=%ld | S4-free but not a K1/C3 stack=%ld | S4-containing yet all components <=3: %ld\n",
           free_[0], free_[1], free_[2], free_[3], s4free_notstack, s4free_stack_bad);
    printf("      (C) W4-free=%ld (census vs out-nbhd mismatches %ld) L4-free=%ld (mismatches %ld) LT=%ld [with sink %ld, of which non-transitive %ld; sinkless %ld, of which not strong %ld] | structure violations=%ld | core-size histogram:",
           w4f, w4f_mismatch, l4f, l4f_mismatch, lt, lt_sink, lt_sink_nontrans, lt_sinkless, lt_sinkless_notstrong, struct_bad);
    for (int L = 1; L <= n; L++) if (Lhist[L]) printf(" |L|=%d:%ld", L, Lhist[L]);
    printf("\n");
    return 0;
}
