/* Audit E (2026-10-14), THM-4655 statements 3-4 / HYP-9279: independent recomputation of the balance delta(T) of
   uniform Hamiltonian paths of strong tournaments.
   Method (different from tbal.c's forward/backward subset DP): explicit depth-first enumeration of every Hamiltonian
   path, with exact integer tallies bef[x][y] = #paths with x before y.  Input: DIGRAPH6 lines (gentourng -q -z -c n),
   decoded here from the digraph6 specification (not the ascii triangle format used by tbal.c).
   Also: strong connectivity (Warshall), automorphism group order by backtracking (orbit-counting check
   sum_T n!/|Aut T| = number of LABELLED strong tournaments s(n) from the recurrence
   2^C(n,2) = sum_k C(n,k) s(k) 2^C(n-k,2)), minimum feedback arc set (subset DP), and the largest transitive
   subtournament (so "transitive plus k" with k = n - that size) for every minimiser.
   Modes:  ./strong_dfs_balance d6 < digraph6-lines         (isomorphism classes from gentourng)
           ./strong_dfs_balance labelled n                   (all 2^C(n,2) labelled tournaments, n <= 7)        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static int n;
static unsigned out[16];          /* out-neighbour bitmasks */
static int path[16];
static u64 bef[16][16];
static u64 H;

static void dfs(int depth, unsigned used) {
    if (depth == n) {
        H++;
        for (int i = 0; i < n; i++) { int x = path[i]; for (int j = i + 1; j < n; j++) bef[x][path[j]]++; }
        return;
    }
    unsigned cand = out[path[depth - 1]] & ~used;
    while (cand) { int v = __builtin_ctz(cand); cand &= cand - 1; path[depth] = v; dfs(depth + 1, used | 1u << v); }
}

static int strong(void) {
    unsigned R[16];
    for (int i = 0; i < n; i++) R[i] = out[i] | 1u << i;
    for (int k = 0; k < n; k++) for (int i = 0; i < n; i++) if (R[i] >> k & 1) R[i] |= R[k];
    unsigned full = (1u << n) - 1;
    for (int i = 0; i < n; i++) if (R[i] != full) return 0;
    return 1;
}

/* automorphisms by backtracking */
static int img[16], usedimg; static u64 autc; static int score[16];
static void autbt(int i) {
    if (i == n) { autc++; return; }
    for (int c = 0; c < n; c++) if (!(usedimg >> c & 1) && score[c] == score[i]) {
        int ok = 1;
        for (int j = 0; j < i && ok; j++) {
            int a = out[i] >> j & 1, b = out[c] >> img[j] & 1;   /* i->j must map to c->img[j] */
            if (a != b) ok = 0;
        }
        if (ok) { img[i] = c; usedimg |= 1u << c; autbt(i + 1); usedimg &= ~(1u << c); }
    }
}
static u64 autorder(void) {
    for (int i = 0; i < n; i++) score[i] = __builtin_popcount(out[i]);
    autc = 0; usedimg = 0; autbt(0); return autc;
}

static int minfas(void) {
    int N = 1 << n; static int f[1 << 12];
    f[0] = 0;
    for (int S = 1; S < N; S++) {
        int best = 1 << 30;
        for (int v = 0; v < n; v++) if (S >> v & 1) {
            int R = S ^ (1 << v);
            int c = f[R] + __builtin_popcount(out[v] & R);   /* v placed last in S: arcs v -> earlier are backward */
            if (c < best) best = c;
        }
        f[S] = best;
    }
    return f[N - 1];
}
static int maxtrans(void) {
    int best = 0;
    for (int S = 1; S < (1 << n); S++) {
        int k = __builtin_popcount(S); if (k <= best) continue;
        unsigned seen = 0; int ok = 1;
        for (int v = 0; v < n && ok; v++) if (S >> v & 1) {
            int d = __builtin_popcount(out[v] & S);
            if (seen >> d & 1) ok = 0; seen |= 1u << d;      /* transitive iff internal scores are distinct */
        }
        if (ok) best = k;
    }
    return best;
}

static u64 binom(int a, int b) { u64 r = 1; for (int i = 1; i <= b; i++) r = r * (a - b + i) / i; return r; }

/* statistics */
static long long total = 0, nonstrong = 0, below = 0, balarc = 0;
static u64 bestn = 1, bestd = 1, arcn = 1, arcd = 1;
static char minis[64][128]; static int nmin = 0, nmin_total = 0;
static u128 labelled_sum = 0; static u64 sumH = 0, maxH = 0;
static long long nonbalz = 0;
static long long vertexwise_ok = 0; static int vertexwise_fail_shown = 0;
static u64 vbest_num = 1, vbest_den = 1; static char vbest_tag[128];

static void process(const char *tag) {
    if (!strong()) { nonstrong++; return; }
    total++;
    memset(bef, 0, sizeof bef); H = 0;
    for (int s = 0; s < n; s++) { path[0] = s; dfs(1, 1u << s); }
    sumH += H; if (H > maxH) maxH = H;
    u64 bm = 0; int hasarc = 0; u64 amin = H;
    for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) {
        if (bef[x][y] + bef[y][x] != H) { printf("TALLY ERROR\n"); exit(1); }
        u64 c = bef[x][y], mn = c < H - c ? c : H - c; if (mn > bm) bm = mn;
        int t = (out[x] >> y & 1) ? x : y, h = t == x ? y : x;
        u64 ct = bef[t][h];
        if (3 * ct >= H && 3 * ct <= 2 * H) hasarc = 1;
        if (ct < amin) amin = ct;
    }
    if (3 * bm < H) below++;
    if (hasarc) balarc++;
    /* HYP-9279(B): every vertex v has a partner w with 1/3 <= P(v before w) <= 2/3.
       vertex balance b(v) = max_w min(P, 1-P); tournament vertex balance = min_v b(v); track the minimum over T. */
    {
        u64 vmin_num = 1, vmin_den = 1; int allv = 1;
        for (int v = 0; v < n; v++) {
            u64 bv = 0;
            for (int w = 0; w < n; w++) if (w != v) { u64 c = bef[v][w], mn = c < H - c ? c : H - c; if (mn > bv) bv = mn; }
            if (3 * bv < H) allv = 0;
            if ((u128)bv * vmin_den < (u128)vmin_num * H) { vmin_num = bv; vmin_den = H; }
        }
        if (allv) vertexwise_ok++; else if (vertexwise_fail_shown++ < 3) printf("VERTEX WITHOUT BALANCED PARTNER: %s\n", tag);
        if ((u128)vmin_num * vbest_den < (u128)vbest_num * vmin_den) { vbest_num = vmin_num; vbest_den = vmin_den;
            strncpy(vbest_tag, tag, 127); vbest_tag[127] = 0; }
    }
    if ((u128)amin * arcd < (u128)arcn * H) { arcn = amin; arcd = H; }
    if ((u128)bm * bestd < (u128)bestn * H) { bestn = bm; bestd = H; nmin = 0; nmin_total = 0; }
    if ((u128)bm * bestd == (u128)bestn * H) {
        nmin_total++;
        if (nmin < 64) { strncpy(minis[nmin], tag, 127); minis[nmin][127] = 0; nmin++; }
    }
    if (tag[0] == '&') { u64 a = autorder(); u64 f = 1; for (int i = 2; i <= n; i++) f *= i; labelled_sum += f / a;
        if (f % a) printf("AUT ORDER DOES NOT DIVIDE n!\n"); }
}

static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    if (argc >= 3 && !strcmp(argv[1], "labelled")) {
        n = atoi(argv[2]);
        int E = n * (n - 1) / 2;
        for (u64 code = 0; code < (1ULL << E); code++) {
            memset(out, 0, sizeof out); int k = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code >> k & 1) out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
            char tag[64]; snprintf(tag, 64, "code=%llu", code);
            process(tag);
        }
    } else {
        char line[512];
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '&') continue;
            line[strcspn(line, "\r\n")] = 0;
            n = line[1] - 63;
            memset(out, 0, sizeof out);
            int bad = 0;
            for (int k = 0; k < n * n; k++) {
                int byte = line[2 + k / 6] - 63, bit = byte >> (5 - k % 6) & 1;
                int i = k / n, j = k % n;
                if (bit) { if (i == j) bad = 1; out[i] |= 1u << j; }
            }
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (((out[i] >> j) & 1) == ((out[j] >> i) & 1)) bad = 1;
            if (bad) { printf("NOT A TOURNAMENT: %s\n", line); continue; }
            process(line);
        }
    }
    u64 g = gcd(bestn, bestd), ga = gcd(arcn, arcd);
    printf("n=%d: strong tournaments %lld (non-strong skipped %lld); delta < 1/3: %lld; with a balanced arc: %lld\n",
           n, total, nonstrong, below, balarc);
    printf("  sum H = %llu, max H = %llu\n", sumH, maxH);
    printf("  minimum delta = %llu/%llu = %.6f, attained by %d tournaments; min over arcs of p(tail first) = %llu/%llu = %.6f\n",
           bestn / g, bestd / g, (double)bestn / bestd, nmin_total, arcn / ga, arcd / ga, (double)arcn / arcd);
    {   u64 gv = gcd(vbest_num, vbest_den);
        printf("  HYP-9279(B): tournaments in which EVERY vertex has a balanced partner: %lld of %lld; "
               "min over T of min over v of max_w min(P,1-P) = %llu/%llu = %.6f (at %s)\n",
               vertexwise_ok, total, vbest_num / gv, vbest_den / gv, (double)vbest_num / vbest_den, vbest_tag); }
    if (labelled_sum) {
        /* labelled strong count by recurrence */
        static u128 s[16]; s[1] = 1;
        for (int m = 2; m <= n; m++) {
            u128 t = (u128)1 << (m * (m - 1) / 2);
            for (int k = 1; k < m; k++) t -= (u128)binom(m, k) * s[k] * ((u128)1 << ((m - k) * (m - k - 1) / 2));
            s[m] = t;
        }
        printf("  orbit-count check: sum n!/|Aut| = %llu ; labelled strong tournaments s(%d) = %llu ; %s\n",
               (u64)labelled_sum, n, (u64)s[n], labelled_sum == s[n] ? "EQUAL (gentourng list complete and irredundant)" : "MISMATCH");
    }
    /* analyse minimisers */
    for (int t = 0; t < nmin && t < 20; t++) {
        const char *L = minis[t];
        if (L[0] == '&') {
            n = L[1] - 63; memset(out, 0, sizeof out);
            for (int k = 0; k < n * n; k++) { int byte = L[2 + k / 6] - 63; if (byte >> (5 - k % 6) & 1) out[k / n] |= 1u << (k % n); }
        } else {
            u64 code = strtoull(L + 5, 0, 10); memset(out, 0, sizeof out); int k = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code >> k & 1) out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
        }
        int sc[16]; for (int i = 0; i < n; i++) sc[i] = __builtin_popcount(out[i]);
        /* sort scores */
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (sc[j] < sc[i]) { int q = sc[i]; sc[i] = sc[j]; sc[j] = q; }
        memset(bef, 0, sizeof bef); H = 0; for (int s0 = 0; s0 < n; s0++) { path[0] = s0; dfs(1, 1u << s0); }
        printf("  minimiser %s: H=%llu scores (", L, H);
        for (int i = 0; i < n; i++) printf("%d%s", sc[i], i + 1 < n ? "," : ")");
        printf(" minFAS=%d transitive-plus-%d\n", minfas(), n - maxtrans());
    }
    return 0;
}
