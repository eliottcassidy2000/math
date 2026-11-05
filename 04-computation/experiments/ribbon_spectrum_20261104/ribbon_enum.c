/* ribbon_enum.c -- exhaustive DFS over span-w ribbon words with H <= N, recording the Hamiltonian-path
   counts of the STRONG words.  Inputs: the profile transfer system (export_matrices.py) and the
   condensation automaton (export_automaton.py; verified against brute-force reachability).
   Pruning: H(prefix) <= H(word) (insertion lemma), so prefixes with H > N are cut; prefixes whose
   condensation state is dead (a strong component before the window that can never be merged) are cut;
   depth is bounded by the Busch floor of strong tournaments (passed as maxdepth).
   usage: ribbon_enum <matrix file> <automaton file> <w> <N> <maxdepth> <out prefix>                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#define MAXD 128
#define MAXL 8
#define MAXC 256
typedef unsigned __int128 u128;
static int d, nl, nf, w, maxdepth, nc, na;
static int fin[MAXD], lmask[MAXL], nnz[MAXL];
static int *Mj[MAXL], *Mi[MAXL];
static uint64_t *Mm[MAXL];
static int trans[MAXC][MAXL];
static uint8_t accept[MAXC];
static uint64_t N;
static uint8_t *attained;
static u128 *witness;
static uint64_t nodes = 0, maxentry = 0, leaves = 0;
static long depth_hist[80], alive_hist[80];

static inline uint64_t Hval(const uint64_t *v) { uint64_t s = 0; for (int k = 0; k < nf; k++) s += v[fin[k]]; return s; }

static void dfs(const uint64_t *v, int depth, int cstate, u128 packed) {
    nodes++;
    alive_hist[depth]++;
    if ((nodes & ((1ull << 26) - 1)) == 0) { fprintf(stderr, "  nodes %llu depth %d\n", (unsigned long long)nodes, depth); }
    if (depth >= maxdepth) return;
    int q = depth + 1;
    for (int li = 0; li < nl; li++) {
        int m = lmask[li];
        if (q >= 2 && (m >> (q - 1))) continue;      /* spans must be < q */
        if (q == 1 && m) continue;
        int c2 = trans[cstate][li];
        if (c2 < 0) continue;                        /* dead for strong connectivity */
        uint64_t nv[MAXD];
        memset(nv, 0, sizeof(uint64_t) * d);
        const int *J = Mj[li], *I = Mi[li]; const uint64_t *MM = Mm[li];
        for (int k = 0; k < nnz[li]; k++) nv[J[k]] += MM[k] * v[I[k]];
        uint64_t H = Hval(nv);
        if (H > N) continue;
        for (int k = 0; k < d; k++) if (nv[k] > maxentry) maxentry = nv[k];
        u128 np = packed | ((u128)m << (w * depth));
        if (accept[c2]) {
            leaves++;
            if (!(H & 1)) { fprintf(stderr, "EVEN H=%llu ?!\n", (unsigned long long)H); exit(2); }
            uint64_t idx = H >> 1;
            if (!((attained[idx >> 3] >> (idx & 7)) & 1)) { attained[idx >> 3] |= (uint8_t)(1u << (idx & 7)); witness[idx] = np; depth_hist[q]++; }
        }
        dfs(nv, depth + 1, c2, np);
    }
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s matrixfile automatonfile w N maxdepth outprefix\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("matrix file"); return 1; }
    w = atoi(argv[3]); N = strtoull(argv[4], 0, 10); maxdepth = atoi(argv[5]);
    if (w * maxdepth > 120) { fprintf(stderr, "maxdepth too large for packing\n"); return 1; }
    if (fscanf(f, "%d %d %d", &d, &nl, &nf) != 3) return 1;
    for (int k = 0; k < nf; k++) if (fscanf(f, "%d", &fin[k]) != 1) return 1;
    for (int li = 0; li < nl; li++) {
        if (fscanf(f, "%d %d", &lmask[li], &nnz[li]) != 2) return 1;
        Mj[li] = malloc(sizeof(int) * nnz[li]); Mi[li] = malloc(sizeof(int) * nnz[li]); Mm[li] = malloc(sizeof(uint64_t) * nnz[li]);
        for (int k = 0; k < nnz[li]; k++) { unsigned long long m; if (fscanf(f, "%d %d %llu", &Mj[li][k], &Mi[li][k], &m) != 3) return 1; Mm[li][k] = m; }
    }
    fclose(f);
    f = fopen(argv[2], "r"); if (!f) { perror("automaton file"); return 1; }
    if (fscanf(f, "%d %d", &nc, &na) != 2) return 1;
    for (int k = 0; k < na; k++) { int s; if (fscanf(f, "%d", &s) != 1) return 1; accept[s] = 1; }
    for (int s = 0; s < nc; s++) for (int li = 0; li < nl; li++) if (fscanf(f, "%d", &trans[s][li]) != 1) return 1;
    fclose(f);
    uint64_t nodd = N / 2 + 1;
    attained = calloc(nodd / 8 + 1, 1); witness = calloc(nodd, sizeof(u128));
    uint64_t v0[MAXD]; memset(v0, 0, sizeof v0); v0[0] = 1;
    clock_t t0 = clock();
    dfs(v0, 0, 0, 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    uint64_t count = 0; for (uint64_t i = 0; i < nodd; i++) if ((attained[i >> 3] >> (i & 7)) & 1) count++;
    char name[512]; snprintf(name, sizeof name, "%s_misses.txt", argv[6]);
    FILE *o = fopen(name, "w");
    fprintf(o, "# strong span-%d ribbon words, H <= %llu, maxdepth %d: nodes %llu, strong leaves %llu, attained odd values %llu of %llu, max state entry %llu, %.1f s\n",
            w, (unsigned long long)N, maxdepth, (unsigned long long)nodes, (unsigned long long)leaves, (unsigned long long)count, (unsigned long long)nodd, (unsigned long long)maxentry, secs);
    fprintf(o, "# new values first attained at length n:");
    for (int q = 0; q < 80; q++) if (depth_hist[q]) fprintf(o, " %d:%ld", q, depth_hist[q]);
    fprintf(o, "\n# alive prefixes per depth:");
    for (int q = 0; q < 80; q++) if (alive_hist[q]) fprintf(o, " %d:%ld", q, alive_hist[q]);
    fprintf(o, "\n# odd misses:\n");
    uint64_t nm = 0;
    for (uint64_t i = 0; i < nodd; i++) if (!((attained[i >> 3] >> (i & 7)) & 1)) { fprintf(o, "%llu\n", (unsigned long long)(2 * i + 1)); nm++; }
    fclose(o);
    snprintf(name, sizeof name, "%s_witnesses.txt", argv[6]);
    o = fopen(name, "w");
    fprintf(o, "# H  word (letters as span masks, vertex 1 first; bit s-1 = back arc of span s); strong words only; H <= 20000\n");
    for (uint64_t i = 0; i < nodd && 2 * i + 1 <= 20000; i++) if (witness[i]) {
        fprintf(o, "%llu ", (unsigned long long)(2 * i + 1));
        u128 p = witness[i]; int q = 0;
        while (p) { fprintf(o, "%s%d", q ? "." : "", (int)(p & ((1u << w) - 1))); p >>= w; q++; }
        fprintf(o, "\n");
    }
    fclose(o);
    printf("done: nodes %llu, strong leaves %llu, attained %llu/%llu, misses %llu, max entry %llu, %.1f s\n", (unsigned long long)nodes, (unsigned long long)leaves, (unsigned long long)count, (unsigned long long)nodd, (unsigned long long)nm, (unsigned long long)maxentry, secs);
    return 0;
}
