/* Strong-tailed near-join grammar, recursive.  Start: one vertex.  Moves on T (order n):
   PAD: T => w (new sink);  TIP(t0, C): (T => C) + t with t beating only t0 in T, C a strong core
   (single vertex, C3, strong 4-tournament, the six strong 5-tournaments).
   Child count via H' = H(T)H(C) + 2(r_t0 F_e(C) + s_t0 F_o(C)); children kept up to order NMAX are expanded,
   deduplicated by (n, H, sorted multiset of (r_v, s_v)).  Usage: grammar NMAX LIMIT */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
#define MAXN 13
static ull Fs[1 << MAXN], dpE[(1 << MAXN) * MAXN], dpS[(1 << MAXN) * MAXN];
static void regs(int n, const unsigned *out, ull *H, ull *r, ull *s) {
    int N = 1 << n;
    memset(dpE, 0, sizeof(ull) * N * n); memset(dpS, 0, sizeof(ull) * N * n);
    unsigned in[MAXN] = {0};
    for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if (out[v] >> u & 1) in[u] |= 1u << v;
    for (int v = 0; v < n; v++) { dpE[(1u << v) * n + v] = 1; dpS[(1u << v) * n + v] = 1; }
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
        ull c = dpE[S * n + v]; if (c) { unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpE[(S | 1u << u) * n + u] += c; } }
        ull d = dpS[S * n + v]; if (d) { unsigned t = in[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpS[(S | 1u << u) * n + u] += d; } } }
    Fs[0] = 1;
    for (int S = 1; S < N; S++) { ull a = 0; for (int v = 0; v < n; v++) a += dpE[S * n + v]; Fs[S] = a; }
    *H = Fs[N - 1];
    for (int v = 0; v < n; v++) {
        ull B = 0;
        for (int S = 1; S < N; S++) if (S >> v & 1) B += dpS[S * n + v] * Fs[(N - 1) & ~S];
        ull he = dpE[(N - 1) * n + v];
        r[v] = (B - he) / 2; s[v] = (B + he) / 2;
    }
}
typedef struct { int m; unsigned out[6]; ull H, Fe, Fo; } Core;
static Core cores[16]; static int ncores;
static void core_F(Core *c) {
    int m = c->m; ull r[MAXN], s[MAXN]; ull H; regs(m, c->out, &H, r, s); c->H = H;
    /* F_e, F_o: sums of H(C - Q) over paths Q (empty in F_e) by parity of |Q| = sum over start vertices of (r_v, s_v) */
    ull Fe = H, Fo = 0; for (int v = 0; v < m; v++) { Fe += r[v]; Fo += s[v]; } c->Fe = Fe; c->Fo = Fo;
}
static int is_strong(int m, const unsigned *out) {
    unsigned reach = 1, fr = 1; while (fr) { unsigned nf = 0; for (int v = 0; v < m; v++) if (fr >> v & 1) nf |= out[v]; nf &= ~reach; reach |= nf; fr = nf; }
    if (reach != (1u << m) - 1) return 0;
    unsigned in[6] = {0}; for (int v = 0; v < m; v++) for (int u = 0; u < m; u++) if (out[v] >> u & 1) in[u] |= 1u << v;
    reach = 1; fr = 1; while (fr) { unsigned nf = 0; for (int v = 0; v < m; v++) if (fr >> v & 1) nf |= in[v]; nf &= ~reach; reach |= nf; fr = nf; }
    return reach == (1u << m) - 1;
}
static unsigned canon(int m, const unsigned *out) {      /* min bit-key over all permutations (m <= 5) */
    int p[6]; for (int i = 0; i < m; i++) p[i] = i; unsigned best = ~0u;
    while (1) { unsigned key = 0; int b = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++, b++) if (out[p[i]] >> p[j] & 1) key |= 1u << b;
        if (key < best) best = key;
        int i = m - 2; while (i >= 0 && p[i] > p[i + 1]) i--; if (i < 0) break;
        int j = m - 1; while (p[j] < p[i]) j--; int t = p[i]; p[i] = p[j]; p[j] = t;
        for (int a = i + 1, z = m - 1; a < z; a++, z--) { t = p[a]; p[a] = p[z]; p[z] = t; } }
    return best;
}
static void build_cores(void) {
    cores[0].m = 1; cores[0].out[0] = 0; core_F(&cores[0]); ncores = 1;
    for (int m = 3; m <= 5; m++) { unsigned seen[64]; int ns = 0; int np = m * (m - 1) / 2;
        for (unsigned bits = 0; bits < (1u << np); bits++) { unsigned out[6] = {0}; int b = 0;
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++, b++) { if (bits >> b & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
            if (!is_strong(m, out)) continue; unsigned c = canon(m, out); int dup = 0;
            for (int k = 0; k < ns; k++) if (seen[k] == c) dup = 1; if (dup) continue; seen[ns++] = c;
            cores[ncores].m = m; memcpy(cores[ncores].out, out, sizeof out); core_F(&cores[ncores]); ncores++; } }
}
typedef struct { int n; unsigned out[MAXN]; } Tour;
static Tour *cur, *nxt; static int ncur, nnxt, CAP = 400000;
#define HS (1 << 22)
static ull *hkeys;
static int insert_key(ull k) { ull i = k & (HS - 1); while (hkeys[i]) { if (hkeys[i] == k) return 0; i = (i + 1) & (HS - 1); } hkeys[i] = k; return 1; }
static unsigned char *reached; static ull LIMIT;
int main(int argc, char **argv) {
    int NMAX = atoi(argv[1]); LIMIT = strtoull(argv[2], 0, 10);
    reached = calloc(LIMIT + 1, 1); hkeys = calloc(HS, sizeof(ull));
    cur = malloc(sizeof(Tour) * CAP); nxt = malloc(sizeof(Tour) * CAP);
    build_cores();
    printf("cores: %d (H, Fe, Fo):", ncores); for (int k = 0; k < ncores; k++) printf(" (%llu,%llu,%llu)", cores[k].H, cores[k].Fe, cores[k].Fo); printf("\n");
    cur[0].n = 1; cur[0].out[0] = 0; ncur = 1; reached[1] = 1;
    if (argc > 3) { FILE *f = fopen(argv[3], "r"); char line[512];
        while (fgets(line, sizeof line, f)) { char *q = strstr(line, "n="); if (!q) continue;
            int n = atoi(q + 2); char *pp = strstr(line, "parents"); if (!pp) continue; pp += 7;
            int par[MAXN]; for (int j = 1; j < n; j++) par[j] = (int)strtol(pp, &pp, 10);
            Tour T; T.n = n; for (int i = 0; i < n; i++) T.out[i] = 0;
            for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (par[j] == i) T.out[j] |= 1u << i; else T.out[i] |= 1u << j; }
            ull H, r[MAXN], s[MAXN]; regs(n, T.out, &H, r, s); if (H <= LIMIT) reached[H] = 1; cur[ncur++] = T; }
        fclose(f); printf("seeds loaded: %d\n", ncur - 1); }
    for (int level = 0; ncur > 0; level++) {
        nnxt = 0;
        for (int a = 0; a < ncur; a++) {
            Tour *T = &cur[a]; int n = T->n; ull H, r[MAXN], s[MAXN]; regs(n, T->out, &H, r, s);
            for (int k = -1; k < ncores; k++) for (int t0 = 0; t0 < n; t0++) {
                if (k < 0 && t0 > 0) break;                       /* PAD: one child */
                int m = k < 0 ? 1 : cores[k].m; int n2 = n + m + (k < 0 ? 0 : 1);
                ull Hc = k < 0 ? H : H * cores[k].H + 2 * (r[t0] * cores[k].Fe + s[t0] * cores[k].Fo);
                if (Hc <= LIMIT) reached[Hc] = 1;
                if (n2 > NMAX || Hc > LIMIT) continue;
                Tour C2; C2.n = n2; memcpy(C2.out, T->out, sizeof(unsigned) * n);
                for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) C2.out[i] |= 1u << (n + j);
                for (int j = 0; j < m; j++) C2.out[n + j] = k < 0 ? 0 : (cores[k].out[j] << n);
                if (k >= 0) { int t = n + m; C2.out[t] = 1u << t0;
                    for (int v = 0; v < n + m; v++) if (v != t0) C2.out[v] |= 1u << t; }
                ull H2, r2[MAXN], s2[MAXN]; regs(n2, C2.out, &H2, r2, s2);
                if (H2 != Hc && k >= 0) { printf("FORMULA MISMATCH %llu %llu\n", H2, Hc); return 1; }
                /* dedupe key */
                ull pr[MAXN]; for (int v = 0; v < n2; v++) pr[v] = r2[v] * 1000003ULL + s2[v];
                for (int x = 1; x < n2; x++) { ull y = pr[x]; int z = x - 1; while (z >= 0 && pr[z] > y) { pr[z + 1] = pr[z]; z--; } pr[z + 1] = y; }
                ull key = 1469598103934665603ULL ^ (ull)n2 * 1099511628211ULL ^ H2 * 7919ULL;
                for (int v = 0; v < n2; v++) key = (key ^ pr[v]) * 1099511628211ULL;
                if (!key) key = 1;
                if (insert_key(key) && nnxt < CAP) nxt[nnxt++] = C2;
            }
        }
        Tour *tmp = cur; cur = nxt; nxt = tmp; ncur = nnxt;
        int miss = 0; ull bound = 0;
        for (ull h = 1; h <= LIMIT; h += 2) if (!reached[h] && h != 7 && h != 21) { if (!bound) bound = h; miss++; }
        printf("after level %d: frontier %d; first odd miss (besides 7, 21): %llu; misses <= %llu: %d\n", level, ncur, bound, LIMIT, miss);
        fflush(stdout);
    }
    printf("odd misses <= %llu:", LIMIT); for (ull h = 1; h <= LIMIT; h += 2) if (!reached[h]) printf(" %llu", h); printf("\n");
}
