/* two_apex.c (klein-2026-11-15) -- the Hamiltonian-path spectrum of two-apex tapes.
   usage: two_apex_small m c [HMAX]  -- only tapes where apex B beats at most c spine vertices (popcount b <= c)
   Spine: transitive 0 < 1 < ... < m-1 (i -> j for i < j).  Apex A beats spine vertex i iff bit i of a is set, else
   i -> A; apex B likewise with b; the arc between A and B is A -> B (orientation o = 0) or B -> A (o = 1).  All
   2^(2m+1) tapes are counted by the subset DP; the output lists every attained value H < HMAX (default 2^22) with the
   number of (a, b, o) triples attaining it, one line "V h count minpop orient" per value (minpop = least over realisations of min(popcount a, popcount b): 0 means a one-apex tape; orient = bitmask of orientations seen), plus a summary.  The one-apex family is the
   sub-family b = 0 with o = 0 (B a sink). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned int u32;
static int n; static u64 adj[32], radj[32]; static u32 *f;
static u64 count_paths(void) {
    u64 Nn = 1ULL << n;
    for (int v = 0; v < n; v++) f[(1ULL << v) * n + v] = 1;
    for (u64 S = 1; S < Nn; S++) {
        if (!(S & (S - 1))) continue;
        u64 vs = S; u32 *dst = f + S * n;
        while (vs) { int v = __builtin_ctzll(vs); vs &= vs - 1; u64 R = S & ~(1ULL << v); u64 in = radj[v] & R; u64 c = 0; const u32 *row = f + R * n;
            while (in) { int u = __builtin_ctzll(in); in &= in - 1; c += row[u]; } dst[v] = (u32)c; }
    }
    u64 hh = 0; for (int v = 0; v < n; v++) hh += f[(Nn - 1) * n + v]; return hh;
}
int main(int argc, char **argv) {
    int m = atoi(argv[1]); int cmax = atoi(argv[2]); u64 HMAX = argc > 3 ? strtoull(argv[3], 0, 10) : (1ULL << 22);
    n = m + 2; f = malloc(sizeof(u32) * (1ULL << n) * n);
    u32 *cnt = calloc(HMAX, sizeof(u32)); unsigned char *minpop = malloc(HMAX); memset(minpop, 255, HMAX); unsigned char *orient = calloc(HMAX, 1); u64 maxh = 0, total = 0;
    int Av = m, Bv = m + 1;
    for (u64 a = 0; a < (1ULL << m); a++) for (u64 b = 0; b < (1ULL << m); b++) { if (__builtin_popcountll(b) > cmax) continue; for (int o = 0; o < 2; o++) {
        memset(adj, 0, sizeof adj); memset(radj, 0, sizeof radj);
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { adj[i] |= 1ULL << j; radj[j] |= 1ULL << i; }
        for (int i = 0; i < m; i++) {
            if (a >> i & 1) { adj[Av] |= 1ULL << i; radj[i] |= 1ULL << Av; } else { adj[i] |= 1ULL << Av; radj[Av] |= 1ULL << i; }
            if (b >> i & 1) { adj[Bv] |= 1ULL << i; radj[i] |= 1ULL << Bv; } else { adj[i] |= 1ULL << Bv; radj[Bv] |= 1ULL << i; }
        }
        if (o == 0) { adj[Av] |= 1ULL << Bv; radj[Bv] |= 1ULL << Av; } else { adj[Bv] |= 1ULL << Av; radj[Av] |= 1ULL << Bv; }
        u64 hh = count_paths(); total++;
        if (hh < HMAX) { cnt[hh]++; int pb = __builtin_popcountll(b), pa = __builtin_popcountll(a); int pm = pb < pa ? pb : pa; if (pm < minpop[hh]) minpop[hh] = pm; orient[hh] |= 1 << o; } if (hh > maxh) maxh = hh;
    } }
    u64 distinct = 0; for (u64 v = 0; v < HMAX; v++) if (cnt[v]) distinct++;
    printf("# two-apex tapes with popcount(b) <= %d, spine %d: %llu tapes, max H = %llu, distinct values below %llu: %llu\n", cmax, m, total, maxh, HMAX, distinct);
    for (u64 v = 0; v < HMAX; v++) if (cnt[v]) printf("V %llu %u %d %d\n", v, cnt[v], minpop[v], orient[v]);
    return 0;
}
