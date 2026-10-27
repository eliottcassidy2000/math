// scan12.c -- independent audit re-implementation of the THM-4812 (a) order-12 invariant scans.
// Differences from the repository's tc_invariant_scan.cpp: sigma is a different conjugate (cycles placed on the LAST
// points, each cycle run backwards), enumeration in plain binary order with the tournament rebuilt from scratch,
// strong connectivity by Landau's score criterion, Held-Karp in pull form, |Aut| by an orbit-stabilizer backtracking
// counter (no nauty), and |Aut| is computed for EVERY strong invariant tournament (no candidate filter); the
// repository's candidate filter (H = t*a, t in {7,21}, a odd, a >= amin, ord(sigma) | a) is replayed only to compare
// counters.
// Build: cc -O3 -march=native -o scan12 scan12.c
// Run:   ./scan12 n amin c1 c2 ...   [slice r m]
#include "../common/tourlib.h"
#include <time.h>

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: scan12 n amin c1 c2 ... [slice r m]\n"); return 1; }
    int n = atoi(argv[1]); long long amin = atoll(argv[2]);
    int cyc[16], nc = 0; long long sr = 0, sm = 1;
    for (int i = 3; i < argc; i++) { if (!strcmp(argv[i], "slice")) { sr = atoll(argv[i + 1]); sm = atoll(argv[i + 2]); break; } cyc[nc++] = atoi(argv[i]); }
    int sigma[32]; for (int i = 0; i < n; i++) sigma[i] = i;
    int tot = 0; for (int i = 0; i < nc; i++) tot += cyc[i];
    if (tot > n) { fprintf(stderr, "too many points\n"); return 1; }
    // place cycles on the last tot points: start at n - tot; each cycle on [a, a+c) maps a+i -> a + (i-1 mod c)
    int pos = n - tot; long long ord = 1;
    for (int i = 0; i < nc; i++) {
        int c = cyc[i];
        for (int j = 0; j < c; j++) sigma[pos + j] = pos + (j - 1 + c) % c;
        pos += c;
        long long a = ord, b = c; while (b) { long long t = a % b; a = b; b = t; } ord = ord / a * c;
    }
    // check sigma has the right cycle type
    // pair orbits under <sigma>
    static int orbit_of[32][32]; for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) orbit_of[u][v] = -1;
    int K = 0; static int ou[600][64], ov[600][64], olen[600];
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        if (orbit_of[u][v] >= 0) continue;
        int a = u, b = v, L = 0;
        while (orbit_of[a][b] < 0) {
            if (orbit_of[b][a] >= 0) { fprintf(stderr, "orbit contains a pair and its reverse\n"); return 1; }
            orbit_of[a][b] = K; ou[K][L] = a; ov[K][L] = b; L++;
            int na = sigma[a], nb = sigma[b]; a = na; b = nb;
        }
        if (!(a == u && b == v)) { fprintf(stderr, "orbit did not close\n"); return 1; }
        // mark reversed pairs as belonging to the same unordered orbit (for the scan over unordered pairs)
        for (int j = 0; j < L; j++) orbit_of[ov[K][j]][ou[K][j]] = K;
        olen[K] = L; K++;
    }
    printf("scan12 n=%d type", n); for (int i = 0; i < nc; i++) printf(" (%d)", cyc[i]);
    printf(" ord=%lld pair-orbits=%d invariant=2^%d amin=%lld slice=%lld/%lld\n", ord, K, K, amin, sr, sm);
    fflush(stdout);
    uint64_t *f = (uint64_t *)malloc(sizeof(uint64_t) * ((size_t)1 << n) * n);
    long long strongc = 0, strongbfs_mismatch = 0, candidates = 0, hits = 0, divfail = 0, nontrivAut = 0;
    uint64_t minH = ~0ull; long long minTc = -1;
    long long autHistAll[4096]; memset(autHistAll, 0, sizeof autHistAll); long long autBig = 0;
    long long candHist[4096]; memset(candHist, 0, sizeof candHist);
    long long tcSmall[64]; memset(tcSmall, 0, sizeof tcSmall);
    clock_t t0 = clock();
    for (uint64_t code = 0; code < (1ull << K); code++) {
        if ((long long)(code % (uint64_t)sm) != sr) continue;
        uint32_t out[32]; memset(out, 0, sizeof out);
        for (int k = 0; k < K; k++) {
            int rev = (code >> k) & 1ull;
            for (int j = 0; j < olen[k]; j++) { int a = ou[k][j], b = ov[k][j]; if (rev) out[b] |= 1u << a; else out[a] |= 1u << b; }
        }
        int st = strong_landau(n, out);
        if ((code & 1023) == 0) { if (st != strong_bfs(n, out)) strongbfs_mismatch++; }
        if (!st) continue;
        strongc++;
        uint32_t in[32]; make_in(n, out, in);
        uint64_t H = ham_pull(n, in, f);
        if (H < minH) minH = H;
        long long A = aut_order(n, out);
        if (A < 4096) autHistAll[A]++; else autBig++;
        if (A > 1) nontrivAut++;
        if (H % (uint64_t)A) { divfail++; continue; }
        long long tc = (long long)(H / (uint64_t)A);
        if (minTc < 0 || tc < minTc) minTc = tc;
        if (tc < 64) tcSmall[tc]++;
        if (tc == 7 || tc == 21) { hits++; printf("HIT tc=%lld H=%llu |Aut|=%lld code=%llx\n", tc, (unsigned long long)H, A, (unsigned long long)code); fflush(stdout); }
        // replay the repository's candidate filter, to compare its counters
        for (int t = 7; t <= 21; t += 14) {
            if (H % t) continue; uint64_t a = H / t;
            if (!(a & 1) || (long long)a < amin || a % (uint64_t)ord) continue;
            candidates++; if (A < 4096) candHist[A]++;
        }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("strong=%lld minH=%llu minTc=%lld hits(tc in {7,21})=%lld H%%Aut failures=%lld strong-test mismatches(sampled)=%lld nontrivial-Aut=%lld\n",
           strongc, (unsigned long long)minH, minTc, hits, divfail, strongbfs_mismatch, nontrivAut);
    printf("replayed repository filter: candidates=%lld, |Aut| histogram:", candidates);
    for (int a = 1; a < 4096; a++) if (candHist[a]) printf(" %d:%lld", a, candHist[a]);
    printf("\n|Aut| histogram of all strong invariant:");
    for (int a = 1; a < 4096; a++) if (autHistAll[a]) printf(" %d:%lld", a, autHistAll[a]);
    if (autBig) printf(" >=4096:%lld", autBig);
    printf("\nsmall tc values (<64) among strong invariant:");
    for (int t = 1; t < 64; t++) if (tcSmall[t]) printf(" %d:%lld", t, tcSmall[t]);
    printf("\ntime %.1f s\n", secs);
    return 0;
}
