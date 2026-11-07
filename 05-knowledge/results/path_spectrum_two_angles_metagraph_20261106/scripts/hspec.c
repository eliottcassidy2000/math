/* hspec.c (klein-2026-11-06) -- level-resolved Hamiltonian-path census of tournaments.
   usage:  hspec n [hc] < gentourng-output
   Reads gentourng's default format (upper triangle row by row, one tournament per line, '1' at position (i,j), i<j,
   meaning i -> j).  For every tournament it computes
     h   = number of Hamiltonian paths (subset DP, pull form: f[S][v] = sum over u in S-v with u->v of f[S-v][u]),
     c3  = number of cyclic triangles = C(n,3) - sum_v C(s_v, 2)  (s_v = out-degree),
     str = 1 if strongly connected (bitmask closure from vertex 0 both ways),
     hc  = number of Hamiltonian cycles (only with the "hc" argument; cycles through vertex 0 = paths from 0 ending at
           an in-neighbour of 0, second DP over subsets containing 0).
   Output (stdout):
     "T h c3 str count"      one line per realised (h, c3, strong) triple
     "W h witness"           the first tournament read with value h (gentourng line)
     "S h witness"           the first strong tournament read with value h
     "C h hc count"          (hc mode) one line per realised (h, hc) pair
     "# ..." summary lines.
   Values of h are assumed below HMAX = 2^18 (checked; the program aborts otherwise). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef unsigned int u32;
typedef u32 cell;
#define HMAX (1u << 18)
#define C3MAX 64
#define HCMAX 4096
static int n; static u64 adj[32], radj[32];
static cell *f;
static u64 count_paths(void) {
    u64 N = 1ULL << n;
    for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
    for (u64 S = 1; S < N; S++) {
        if (!(S & (S - 1))) continue;                 /* singletons already set */
        u64 vs = S; cell *dst = f + S * n;
        while (vs) {
            int v = __builtin_ctzll(vs); vs &= vs - 1;
            u64 R = S & ~(1ULL << v); u64 in = radj[v] & R; u64 c = 0;
            const cell *row = f + R * n;
            while (in) { int u = __builtin_ctzll(in); in &= in - 1; c += row[u]; }
            dst[v] = c;
        }
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += f[(N - 1) * n + v]; return h;
}
static u64 count_cycles(void) {           /* Hamiltonian cycles = paths from vertex 0 covering all, ending at an in-neighbour of 0 */
    u64 N = 1ULL << n; u64 h = 0;
    /* reuse the table: f[S][v] for S containing 0, paths starting at 0 */
    for (u64 S = 1; S < N; S += 2) {      /* S contains vertex 0 */
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue;
            if (v == 0) { f[S * n + 0] = (S == 1); continue; }
            u64 R = S & ~(1ULL << v); u64 in = radj[v] & R; u64 c = 0; const cell *row = f + R * n;
            while (in) { int u = __builtin_ctzll(in); in &= in - 1; c += row[u]; }
            f[S * n + v] = c;
        }
    }
    u64 in0 = radj[0]; while (in0) { int u = __builtin_ctzll(in0); in0 &= in0 - 1; h += f[(N - 1) * n + u]; }
    return h;
}
static int is_strong(void) {
    u64 all = (1ULL << n) - 1, seen = 1, frontier = 1;
    while (frontier) { u64 nf = 0; u64 x = frontier; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= adj[v]; } nf &= ~seen; seen |= nf; frontier = nf; }
    if (seen != all) return 0;
    seen = 1; frontier = 1;
    while (frontier) { u64 nf = 0; u64 x = frontier; while (x) { int v = __builtin_ctzll(x); x &= x - 1; nf |= radj[v]; } nf &= ~seen; seen |= nf; frontier = nf; }
    return seen == all;
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: hspec n [hc] < gentourng output\n"); return 1; }
    n = atoi(argv[1]); int want_hc = (argc > 2 && !strcmp(argv[2], "hc"));
    if (argc > 4 && !strcmp(argv[2], "witness")) {   /* hspec n witness h k : print the first tournament with value h at level k */
        u64 wh = strtoull(argv[3], 0, 10); int wk = atoi(argv[4]); u64 N = 1ULL << n; f = malloc(sizeof(cell) * N * n);
        char line[1024]; int L0 = n * (n - 1) / 2;
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0; if (L != L0) continue;
            memset(adj, 0, sizeof adj); memset(radj, 0, sizeof radj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') { adj[i] |= 1ULL << j; radj[j] |= 1ULL << i; } else { adj[j] |= 1ULL << i; radj[i] |= 1ULL << j; } }
            int c3 = n * (n - 1) * (n - 2) / 6; for (int v = 0; v < n; v++) { int s = __builtin_popcountll(adj[v]); c3 -= s * (s - 1) / 2; }
            if (c3 != wk) continue;
            if (count_paths() == wh) { printf("order %d h %llu level %d strong %d: %s  scores", n, wh, wk, is_strong(), line); for (int v = 0; v < n; v++) printf(" %d", __builtin_popcountll(adj[v])); printf("\n"); return 0; }
        }
        printf("order %d h %llu level %d: no tournament\n", n, wh, wk); return 1;
    }
    if (n < 2 || n > 20) { fprintf(stderr, "n out of range\n"); return 1; }
    u64 N = 1ULL << n; f = malloc(sizeof(cell) * N * n); if (!f) { fprintf(stderr, "alloc failed\n"); return 2; }
    u32 *cnt = calloc((size_t)HMAX * C3MAX * 2, sizeof(u32)); if (!cnt) { fprintf(stderr, "alloc failed\n"); return 2; }
    u32 *chc = want_hc ? calloc((size_t)HMAX * HCMAX, sizeof(u32)) : NULL; if (want_hc && !chc) { fprintf(stderr, "alloc failed\n"); return 2; }
    char **wit = calloc(HMAX, sizeof(char *)), **wits = calloc(HMAX, sizeof(char *));
    char line[1024]; u64 total = 0, nstrong = 0, maxh = 0, maxhc = 0; int L0 = n * (n - 1) / 2;
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (L != L0) continue;
        memset(adj, 0, sizeof adj); memset(radj, 0, sizeof radj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (line[idx++] == '1') { adj[i] |= 1ULL << j; radj[j] |= 1ULL << i; } else { adj[j] |= 1ULL << i; radj[i] |= 1ULL << j; }
        }
        u64 h = count_paths(); total++;
        if (h >= HMAX) { fprintf(stderr, "h = %llu exceeds HMAX\n", h); return 3; }
        int c3 = n * (n - 1) * (n - 2) / 6; for (int v = 0; v < n; v++) { int s = __builtin_popcountll(adj[v]); c3 -= s * (s - 1) / 2; }
        int st = is_strong(); nstrong += st; if (h > maxh) maxh = h;
        cnt[((size_t)h * C3MAX + c3) * 2 + st]++;
        if (!wit[h]) wit[h] = strdup(line);
        if (st && !wits[h]) wits[h] = strdup(line);
        if (want_hc) { u64 hc = st ? count_cycles() : 0; if (hc >= HCMAX) { fprintf(stderr, "hc = %llu exceeds HCMAX\n", hc); return 3; } chc[(size_t)h * HCMAX + hc]++; if (hc > maxhc) maxhc = hc; }
    }
    printf("# order %d: %llu tournaments read, %llu strong, max h = %llu\n", n, total, nstrong, maxh);
    if (want_hc) printf("# max hc = %llu\n", maxhc);
    u64 distinct = 0, distinct_s = 0;
    for (u64 h = 0; h < HMAX; h++) {
        int any = 0, anys = 0;
        for (int c3 = 0; c3 < C3MAX; c3++) for (int st = 0; st < 2; st++) {
            u32 c = cnt[((size_t)h * C3MAX + c3) * 2 + st]; if (!c) continue; any = 1; if (st) anys = 1;
            printf("T %llu %d %d %u\n", h, c3, st, c);
        }
        distinct += any; distinct_s += anys;
    }
    for (u64 h = 0; h < HMAX; h++) if (wit[h]) printf("W %llu %s\n", h, wit[h]);
    for (u64 h = 0; h < HMAX; h++) if (wits[h]) printf("S %llu %s\n", h, wits[h]);
    if (want_hc) for (u64 h = 0; h < HMAX; h++) for (u64 hc = 0; hc < HCMAX; hc++) { u32 c = chc[(size_t)h * HCMAX + hc]; if (c) printf("C %llu %llu %u\n", h, hc, c); }
    printf("# distinct h = %llu, distinct strong h = %llu\n", distinct, distinct_s);
    return 0;
}
