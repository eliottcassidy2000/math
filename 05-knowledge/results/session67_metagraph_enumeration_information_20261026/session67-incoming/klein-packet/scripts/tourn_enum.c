/* tourn_enum.c -- exact per-class invariants for n-tournaments read from gentourng output.
 *
 * Reads gentourng's default ascii format (upper triangle row by row; char '1' at pair (i,j), i<j, means i->j).
 * For every class: scores, c3 (cyclic triples), strong components (Landau cuts), the modular
 * decomposition tree (linear nodes = strong components; prime nodes = maximal proper modules of a
 * strong tournament, quotient canonicalised with nauty), |Aut| (nauty), self-converse flag,
 * optional fas (Slater index) and H (Hamiltonian paths) by subset DP, and the canonical form of the
 * 3-cycle hypergraph (nauty on the vertex/edge incidence graph, 2 colours).
 *
 * Also collects, for every distinct prime quotient Q, the cycle-type histogram of Aut(Q) and of
 * Aut+-(Q) (automorphisms and anti-automorphisms), for Polya substitution downstream.
 *
 * Build: cc -O2 -I/opt/homebrew/include/nauty tourn_enum.c -L/opt/homebrew/lib -lnauty -o tourn_enum
 * Use:   gentourng n | ./tourn_enum n [--fas] [--hpaths] [--c3canon] [--table FILE] [--primes FILE]
 */
#define MAXN 64
#include "nauty.h"
#include "nautinv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int N;
static uint32_t out[32];
static int do_fas = 0, do_h = 0, do_c3 = 0;
static FILE *ftab = NULL, *fprimes = NULL;

static int popc(uint32_t x) { return __builtin_popcount(x); }

/* ---------------- nauty canonical form of a tournament on k vertices (masks q[]) ---------------- */
static int gen_count;
static int gens[256][MAXN];
static void autoproc(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n) {
    (void)orbits; (void)numorbits; (void)stabvertex;
    if (gen_count < 256) { memcpy(gens[gen_count], perm, n * sizeof(int)); gen_count++; }
    (void)count;
}

/* returns packed canonical upper triangle; fills grp (group order), lab (canonical labelling) */
static uint64_t canon_tourn(const uint32_t *q, int k, double *grp, int *lab_out, int want_gens) {
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;
    gen_count = 0;
    if (want_gens) options.userautomproc = autoproc;
    EMPTYGRAPH(g, 1, k);
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) if (q[i] >> j & 1) g[i] |= bit[j];
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, k, cg);
    uint64_t key = 0; int bitpos = 0;
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { if (cg[i] & bit[j]) key |= (1ULL << bitpos); bitpos++; }
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10.0; *grp = s; }
    if (lab_out) memcpy(lab_out, lab, k * sizeof(int));
    return key;
}

/* ---------------- group enumeration from generators (k<=16, perms packed 4 bits/entry) ---------------- */
#define HSZ (1<<16)
static uint64_t hkeys[HSZ]; static unsigned char hused[HSZ];
static uint64_t pack_perm(const int *p, int k) { uint64_t v = 0; for (int i = 0; i < k; i++) v |= (uint64_t)p[i] << (4*i); return v; }
static int hins(uint64_t v) { uint32_t h = (uint32_t)((v * 0x9E3779B97F4A7C15ULL) >> 48); while (hused[h]) { if (hkeys[h] == v) return 0; h = (h + 1) & (HSZ-1); } hused[h] = 1; hkeys[h] = v; return 1; }
static uint64_t grp_elems[65536]; static int grp_n;
static void enumerate_group(int k) {
    memset(hused, 0, sizeof hused); grp_n = 0;
    int id[MAXN]; for (int i = 0; i < k; i++) id[i] = i;
    grp_elems[grp_n++] = pack_perm(id, k); hins(grp_elems[0]);
    for (int head = 0; head < grp_n; head++) {
        int p[MAXN]; for (int i = 0; i < k; i++) p[i] = (grp_elems[head] >> (4*i)) & 15;
        for (int gi = 0; gi < gen_count; gi++) {
            int r[MAXN]; for (int i = 0; i < k; i++) r[i] = gens[gi][p[i]];
            uint64_t v = pack_perm(r, k); if (hins(v)) { if (grp_n < 65536) grp_elems[grp_n++] = v; }
        }
    }
}
static void cycle_type(const int *p, int k, int *ct) { /* ct[len] = number of cycles of that length */
    int seen[MAXN] = {0}; for (int i = 1; i <= k; i++) ct[i] = 0;
    for (int i = 0; i < k; i++) if (!seen[i]) { int len = 0, j = i; while (!seen[j]) { seen[j] = 1; j = p[j]; len++; } ct[len]++; }
}

/* ---------------- prime quotient registry ---------------- */
typedef struct { int k; uint64_t key; long count; } prime_rec;
#define PCAP (1<<23)
static prime_rec *primes; static int nprimes = 0; static int *pslot; static int noreg = 0, compact = 0;
static int find_prime(int k, uint64_t key) { uint64_t h = (key * 0x9E3779B97F4A7C15ULL + (uint64_t)k * 0xC2B2AE3D27D4EB4FULL) >> 41; for (;;) { h &= PCAP - 1; int s = pslot[h]; if (s < 0) return -1; if (primes[s].k == k && primes[s].key == key) return s; h++; } }
static int add_prime(int k, uint64_t key) { if (nprimes >= PCAP/2) { fprintf(stderr, "registry full\n"); exit(5); } uint64_t h = (key * 0x9E3779B97F4A7C15ULL + (uint64_t)k * 0xC2B2AE3D27D4EB4FULL) >> 41; for (;;) { h &= PCAP - 1; if (pslot[h] < 0) { pslot[h] = nprimes; primes[nprimes].k = k; primes[nprimes].key = key; primes[nprimes].count = 0; return nprimes++; } h++; } }

static void emit_prime(int k, uint64_t key, const uint32_t *q) {
    /* called once per new prime quotient: write cycle-type histograms of Aut and Aut+- */
    double grp; int lab[MAXN];
    uint64_t kk = canon_tourn(q, k, &grp, lab, 1);
    (void)kk;
    enumerate_group(k);
    if (!fprimes) return;
    fprintf(fprimes, "PRIME k=%d key=%llu aut=%d", k, (unsigned long long)key, grp_n);
    /* histogram of cycle types: print each element's cycle type as partition */
    fprintf(fprimes, " AUT:");
    for (int e = 0; e < grp_n; e++) { int p[MAXN]; for (int i = 0; i < k; i++) p[i] = (grp_elems[e] >> (4*i)) & 15; int ct[MAXN+1]; cycle_type(p, k, ct); fprintf(fprimes, " "); for (int l = 1; l <= k; l++) if (ct[l]) fprintf(fprimes, "%d^%d,", l, ct[l]); }
    /* anti-automorphisms: isomorphism Q -> Q* if self-converse */
    uint32_t qc[MAXN]; for (int i = 0; i < k; i++) { qc[i] = 0; for (int j = 0; j < k; j++) if (i != j && !(q[i] >> j & 1)) qc[i] |= 1u << j; }
    int labc[MAXN]; double grpc; uint64_t keyc = canon_tourn(qc, k, &grpc, labc, 0);
    fprintf(fprimes, " SC=%d", keyc == key);
    if (keyc == key) {
        int psi[MAXN]; for (int i = 0; i < k; i++) psi[lab[i]] = labc[i];
        fprintf(fprimes, " ANTI:");
        for (int e = 0; e < grp_n; e++) { int p[MAXN]; for (int i = 0; i < k; i++) p[i] = (grp_elems[e] >> (4*i)) & 15; int r[MAXN]; for (int i = 0; i < k; i++) r[i] = psi[p[i]]; int ct[MAXN+1]; cycle_type(r, k, ct); fprintf(fprimes, " "); for (int l = 1; l <= k; l++) if (ct[l]) fprintf(fprimes, "%d^%d,", l, ct[l]); }
    }
    /* also the quotient itself as an arc list for reconstruction */
    fprintf(fprimes, " ARCS:"); for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) if (q[i] >> j & 1) fprintf(fprimes, "%d>%d,", i, j);
    fprintf(fprimes, "\n");
}

/* ---------------- strong components of T|V, bottom first ---------------- */
static int scc_split(uint32_t V, uint32_t *comps) {
    int verts[32], m = 0; for (int v = 0; v < N; v++) if (V >> v & 1) verts[m++] = v;
    int sc[32]; for (int i = 0; i < m; i++) sc[i] = popc(out[verts[i]] & V);
    /* insertion sort by score */
    for (int i = 1; i < m; i++) { int v = verts[i], s = sc[i], j = i - 1; while (j >= 0 && sc[j] > s) { verts[j+1] = verts[j]; sc[j+1] = sc[j]; j--; } verts[j+1] = v; sc[j+1] = s; }
    int nc = 0, cum = 0; uint32_t cur = 0;
    for (int i = 0; i < m; i++) { cum += sc[i]; cur |= 1u << verts[i]; if (cum == (i+1)*i/2) { comps[nc++] = cur; cur = 0; } }
    return nc;
}

/* ---------------- modular decomposition ---------------- */
static uint32_t closure(uint32_t V, uint32_t S) {
    int changed = 1;
    while (changed) { changed = 0; uint32_t rest = V & ~S; while (rest) { int z = __builtin_ctz(rest); rest &= rest - 1; uint32_t a = out[z] & S; if (a != 0 && a != S) { S |= 1u << z; changed = 1; } } }
    return S;
}
static char sig[4096]; static int siglen;
static int n_linear, n_prime; static double prod_kfact; static int lin_k_hist[33];
static void sigcat(const char *s) { int l = strlen(s); if (siglen + l < 4000) { memcpy(sig + siglen, s, l); siglen += l; sig[siglen] = 0; } }
static double fact(int k) { double f = 1; for (int i = 2; i <= k; i++) f *= i; return f; }

static void decompose(uint32_t V) {
    int m = popc(V);
    if (m == 1) { sigcat("."); return; }
    uint32_t comps[32]; int nc = scc_split(V, comps);
    char buf[64];
    if (nc >= 2) {
        n_linear++; prod_kfact *= fact(nc); lin_k_hist[nc]++;
        sprintf(buf, "L%d[", nc); sigcat(buf);
        for (int c = 0; c < nc; c++) { decompose(comps[c]); if (c + 1 < nc) sigcat(","); }
        sigcat("]");
        return;
    }
    /* strong, m>=3: maximal proper modules */
    uint32_t M[32]; int verts[32], mm = 0; for (int v = 0; v < N; v++) if (V >> v & 1) verts[mm++] = v;
    for (int i = 0; i < mm; i++) { int x = verts[i]; M[i] = 1u << x; for (int j = 0; j < mm; j++) if (j != i) { uint32_t C = closure(V, (1u << x) | (1u << verts[j])); if (C != V) M[i] |= C; } }
    /* blocks = distinct M[i]; check partition */
    uint32_t blocks[32]; int nb = 0;
    for (int i = 0; i < mm; i++) { int found = 0; for (int b = 0; b < nb; b++) { if (blocks[b] == M[i]) { found = 1; break; } if (blocks[b] & M[i]) { fprintf(stderr, "PARTITION VIOLATION\n"); exit(2); } } if (!found) blocks[nb++] = M[i]; }
    if (nb < 3) { fprintf(stderr, "quotient with %d blocks on strong set\n", nb); exit(3); }
    /* quotient on representatives (lowest vertex of each block) */
    uint32_t q[32]; int rep[32]; for (int b = 0; b < nb; b++) rep[b] = __builtin_ctz(blocks[b]);
    for (int a = 0; a < nb; a++) { q[a] = 0; for (int b = 0; b < nb; b++) if (a != b && (out[rep[a]] >> rep[b] & 1)) q[a] |= 1u << b; }
    uint64_t key = canon_tourn(q, nb, NULL, NULL, 0);
    if (!noreg) { int pi = find_prime(nb, key); if (pi < 0) { pi = add_prime(nb, key); emit_prime(nb, key, q); } primes[pi].count++; }
    n_prime++;
    sprintf(buf, "P%d:%llu[", nb, (unsigned long long)key); sigcat(buf);
    for (int b = 0; b < nb; b++) { decompose(blocks[b]); if (b + 1 < nb) sigcat(","); }
    sigcat("]");
}

/* ---------------- fas and H by subset DP ---------------- */
static int fas_dp(void) {
    int S, full = (1 << N) - 1; static int f[1 << 12];
    f[0] = 0;
    for (S = 1; S <= full; S++) { int best = 1 << 30; uint32_t s = S; while (s) { int v = __builtin_ctz(s); s &= s - 1; int c = f[S & ~(1 << v)] + popc(out[v] & (S & ~(1 << v))); if (c < best) best = c; } f[S] = best; }
    return f[full];
}
static double h_dp(void) {
    int full = (1 << N) - 1; static double h[1 << 12][12];
    for (int S = 1; S <= full; S++) for (int v = 0; v < N; v++) h[S][v] = 0;
    for (int v = 0; v < N; v++) h[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) { uint32_t s = S; while (s) { int v = __builtin_ctz(s); s &= s - 1; double hv = h[S][v]; if (hv == 0) continue; uint32_t nx = out[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; h[S | (1 << w)][w] += hv; } } }
    double tot = 0; for (int v = 0; v < N; v++) tot += h[full][v]; return tot;
}

/* ---------------- 3-cycle hypergraph canonical hash ---------------- */
static uint64_t c3_canon_hash(int *c3out) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options); statsblk stats; options.getcanon = TRUE; options.defaultptn = FALSE;
    int e = 0; int tot;
    /* count cyclic triples */
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) for (int k = j + 1; k < N; k++) {
        int ij = out[i] >> j & 1, jk = out[j] >> k & 1, ki = out[k] >> i & 1;
        if ((ij && jk && ki) || (!ij && !jk && !ki)) e++; }
    tot = N + e; if (c3out) *c3out = e;
    EMPTYGRAPH(g, 1, tot);
    int idx = N;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) for (int k = j + 1; k < N; k++) {
        int ij = out[i] >> j & 1, jk = out[j] >> k & 1, ki = out[k] >> i & 1;
        if ((ij && jk && ki) || (!ij && !jk && !ki)) { g[i] |= bit[idx]; g[idx] |= bit[i]; g[j] |= bit[idx]; g[idx] |= bit[j]; g[k] |= bit[idx]; g[idx] |= bit[k]; idx++; } }
    for (int i = 0; i < tot; i++) { lab[i] = i; ptn[i] = 1; } ptn[N - 1] = 0; ptn[tot - 1] = 0;
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, tot, cg);
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < tot; i++) { h ^= (uint64_t)cg[i]; h *= 1099511628211ULL; h ^= h >> 29; }
    return h ^ ((uint64_t)tot << 56);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    N = atoi(argv[1]);
    for (int a = 2; a < argc; a++) { if (!strcmp(argv[a], "--fas")) do_fas = 1; else if (!strcmp(argv[a], "--hpaths")) do_h = 1; else if (!strcmp(argv[a], "--c3canon")) do_c3 = 1; else if (!strcmp(argv[a], "--table")) ftab = fopen(argv[++a], "w"); else if (!strcmp(argv[a], "--primes")) fprimes = fopen(argv[++a], "w"); else if (!strcmp(argv[a], "--noreg")) noreg = 1; else if (!strcmp(argv[a], "--compact")) compact = 1; }
    primes = malloc(sizeof(prime_rec) * (PCAP/2)); pslot = malloc(sizeof(int) * PCAP); for (int i = 0; i < PCAP; i++) pslot[i] = -1;
    char line[1024]; long idx = 0;
    long cnt_sc = 0, cnt_strong = 0, cnt_linroot = 0, cnt_primeroot = 0, cnt_prime_itself = 0;
    static long c3hist[256]; static long fashist[128]; static long auth[4096]; static long nchist[33];
    if (ftab && !compact) fprintf(ftab, "idx\tkey\taut\tsc\tscores\tc3\tfas\tH\tncomp\tcompsizes\tlin\tprime\tlogN\tsig\tc3hash\n");
    if (ftab && compact) fprintf(ftab, "aut\tsc\tc3\tncomp\tlin\tprime\tlogN\tisprime\tc3hash\n");
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (L != N * (N - 1) / 2) continue;
        for (int i = 0; i < N; i++) out[i] = 0;
        int p = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (line[p++] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
        idx++;
        int sc[32]; int ssum2 = 0; for (int i = 0; i < N; i++) { sc[i] = popc(out[i]); ssum2 += sc[i] * (sc[i] - 1) / 2; }
        int c3 = N * (N - 1) * (N - 2) / 6 - ssum2;
        double grp; int lab[MAXN]; uint64_t key = canon_tourn(out, N, &grp, lab, 0);
        uint32_t conv[32]; for (int i = 0; i < N; i++) { conv[i] = 0; for (int j = 0; j < N; j++) if (i != j && !(out[i] >> j & 1)) conv[i] |= 1u << j; }
        uint64_t keyc = canon_tourn(conv, N, NULL, NULL, 0); int isc = (keyc == key);
        uint32_t comps[32]; int nc = scc_split((1u << N) - 1, comps);
        siglen = 0; sig[0] = 0; n_linear = 0; n_prime = 0; prod_kfact = 1;
        decompose((1u << N) - 1);
        int fas = do_fas ? fas_dp() : -1; double H = do_h ? h_dp() : -1;
        int c3b = c3; uint64_t c3h = do_c3 ? c3_canon_hash(&c3b) : 0;
        if (do_c3 && c3b != c3) { fprintf(stderr, "c3 mismatch\n"); return 4; }
        cnt_sc += isc; if (nc == 1) cnt_strong++; if (nc >= 2) cnt_linroot++; else if (N >= 3) { cnt_primeroot++; if (n_prime == 1 && n_linear == 0 && sig[0] == 'P' && strchr(sig, '[') && !strchr(sig, 'L') && strstr(sig, "[.") ) { /* prime itself iff all children leaves */ int allleaf = 1; for (char *s = strchr(sig, '[') + 1; *s && *s != ']'; s++) if (*s != '.' && *s != ',') { allleaf = 0; break; } if (allleaf) cnt_prime_itself++; } }
        c3hist[c3]++; if (fas >= 0) fashist[fas]++; { int g = (int)(grp + 0.5); if (g < 4096) auth[g]++; } nchist[nc]++;
        int isprime = 0; if (N >= 3 && nc == 1 && n_prime == 1 && n_linear == 0) { isprime = 1; for (char *s = strchr(sig, '[') + 1; *s && *s != ']'; s++) if (*s != '.' && *s != ',') { isprime = 0; break; } }
        if (ftab && compact) { double l = 0; char *s = sig; while (*s) { if (*s == 'L') { int k = atoi(s + 1); for (int i = 2; i <= k; i++) l += __builtin_log2(i); } s++; } fprintf(ftab, "%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%d\t%016llx\n", (int)(grp + 0.5), isc, c3, nc, n_linear, n_prime, l + n_prime, isprime, (unsigned long long)c3h); }
        if (ftab && !compact) {
            /* sorted scores */
            int s2[32]; memcpy(s2, sc, sizeof(int) * N); for (int i = 1; i < N; i++) { int v = s2[i], j = i - 1; while (j >= 0 && s2[j] > v) { s2[j+1] = s2[j]; j--; } s2[j+1] = v; }
            fprintf(ftab, "%ld\t%llu\t%d\t%d\t", idx, (unsigned long long)key, (int)(grp + 0.5), isc);
            for (int i = 0; i < N; i++) fprintf(ftab, "%d%s", s2[i], i + 1 < N ? "," : "");
            fprintf(ftab, "\t%d\t%d\t%.0f\t%d\t", c3, fas, H, nc);
            for (int c = 0; c < nc; c++) fprintf(ftab, "%d%s", popc(comps[c]), c + 1 < nc ? "," : "");
            double logN = n_prime;
            /* log2 of class size N(T) = prod k! * 2^prime */
            { double l = 0; char *s = sig; while (*s) { if (*s == 'L') { int k = atoi(s + 1); for (int i = 2; i <= k; i++) l += __builtin_log2(i); } s++; } logN = l + n_prime; }
            fprintf(ftab, "\t%d\t%d\t%.6f\t%s\t%016llx\n", n_linear, n_prime, logN, sig, (unsigned long long)c3h);
        }
    }
    printf("n=%d classes=%ld strong=%ld selfconverse=%ld linear_root=%ld prime_root=%ld prime_tournaments=%ld distinct_prime_quotients=%d\n", N, idx, cnt_strong, cnt_sc, cnt_linroot, cnt_primeroot, cnt_prime_itself, nprimes);
    printf("c3_levels:"); for (int k = 0; k < 256; k++) if (c3hist[k]) printf(" %d:%ld", k, c3hist[k]); printf("\n");
    if (do_fas) { printf("fas_layers:"); for (int k = 0; k < 128; k++) if (fashist[k]) printf(" %d:%ld", k, fashist[k]); printf("\n"); }
    printf("aut_hist:"); for (int g = 0; g < 4096; g++) if (auth[g]) printf(" %d:%ld", g, auth[g]); printf("\n");
    printf("ncomp_hist:"); for (int k = 1; k <= 32; k++) if (nchist[k]) printf(" %d:%ld", k, nchist[k]); printf("\n");
    printf("prime_quotients_by_k:"); { long byk[33] = {0}; for (int i = 0; i < nprimes; i++) byk[primes[i].k]++; for (int k = 3; k <= N; k++) if (byk[k]) printf(" %d:%ld", k, byk[k]); } printf("\n");
    if (fprimes) { fprintf(fprimes, "#USAGE"); for (int i = 0; i < nprimes; i++) fprintf(fprimes, " %d:%llu:%ld", primes[i].k, (unsigned long long)primes[i].key, primes[i].count); fprintf(fprimes, "\n"); fclose(fprimes); }
    if (ftab) fclose(ftab);
    return 0;
}
