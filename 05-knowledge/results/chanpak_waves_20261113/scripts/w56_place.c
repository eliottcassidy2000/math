/* Place an arbitrary odd target p (e.g. a prime) as p = h(A) + 2^k B(A,t0): A from a pool of one-parent trees
   (n <= NB), t0 any vertex, k >= 1; the tournament is (A => M_k) + t with t beating only t0 (a near-join).
   Verification by direct DP when the order is <= 22.  Usage: place NB p1 p2 ...  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
#define MAXN 23
static ull Fs[1 << 12], dpE[(1 << 12) * 12], dpS[(1 << 12) * 12];
static void regs(int n, const unsigned *out, ull *H, ull *B) {
    int N = 1 << n; unsigned in[12] = {0};
    for (int v = 0; v < n; v++) for (int u = 0; u < n; u++) if (out[v] >> u & 1) in[u] |= 1u << v;
    memset(dpE, 0, sizeof(ull) * N * n); memset(dpS, 0, sizeof(ull) * N * n);
    for (int v = 0; v < n; v++) { dpE[(1u << v) * n + v] = 1; dpS[(1u << v) * n + v] = 1; }
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) {
        ull c = dpE[S * n + v]; if (c) { unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpE[(S | 1u << u) * n + u] += c; } }
        ull d = dpS[S * n + v]; if (d) { unsigned t = in[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dpS[(S | 1u << u) * n + u] += d; } } }
    Fs[0] = 1; for (int S = 1; S < N; S++) { ull a = 0; for (int v = 0; v < n; v++) a += dpE[S * n + v]; Fs[S] = a; }
    *H = Fs[N - 1];
    for (int v = 0; v < n; v++) { ull b = 0; for (int S = 1; S < N; S++) if (S >> v & 1) b += dpS[S * n + v] * Fs[(N - 1) & ~S]; B[v] = b; }
}
typedef struct { ull H, B; int n, t0; signed char par[12]; } Ent;
static Ent *pool; static int npool, cap = 2000000;
#define HS (1 << 22)
static int *ht;
static ull hk(ull H, ull B) { ull k = H * 0x9E3779B97F4A7C15ULL ^ (B + 0x632BE59BD9B4E019ULL) * 0xBF58476D1CE4E5B9ULL; return k ^ (k >> 31); }
static int find(ull H, ull B) { ull i = hk(H, B) & (HS - 1); while (ht[i] >= 0) { if (pool[ht[i]].H == H && pool[ht[i]].B == B) return ht[i]; i = (i + 1) & (HS - 1); } return -1; }
static void add(Ent e) { ull i = hk(e.H, e.B) & (HS - 1); while (ht[i] >= 0) { if (pool[ht[i]].H == e.H && pool[ht[i]].B == e.B) return; i = (i + 1) & (HS - 1); }
    if (npool < cap) { pool[npool] = e; ht[i] = npool++; } }
static int p[12]; static unsigned out[12];
static void rec(int j, int n) {
    if (j == n) { for (int i = 0; i < n; i++) out[i] = 0;
        for (int b = 1; b < n; b++) for (int i = 0; i < b; i++) { if (p[b] == i) out[b] |= 1u << i; else out[i] |= 1u << b; }
        ull H, B[12]; regs(n, out, &H, B);
        for (int v = 0; v < n; v++) { Ent e; e.H = H; e.B = B[v]; e.n = n; e.t0 = v; for (int i = 0; i < n; i++) e.par[i] = p[i]; add(e); } return; }
    for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
static ull Hfull(int n, const unsigned *o) {    /* DP for verification, n <= 22 */
    int N = 1 << n; ull *dp = calloc((size_t)N * n, sizeof(ull));
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { ull c = dp[(size_t)S * n + v]; if (!c) continue;
        unsigned t = o[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; free(dp); return s; }
int main(int argc, char **argv) {
    int NB = atoi(argv[1]); pool = malloc(sizeof(Ent) * cap); ht = malloc(sizeof(int) * HS); memset(ht, -1, sizeof(int) * HS);
    p[0] = -1; for (int n = 1; n <= NB; n++) rec(1, n);
    printf("pool: %d distinct (H, B) rooted pairs from one-parent trees with <= %d vertices\n", npool, NB);
    /* distinct H values */
    for (int a = 2; a < argc; a++) {
        ull target = strtoull(argv[a], 0, 10); int done = 0;
        static const ull HC[8] = {1, 3, 5, 9, 11, 15, 13, 15}, PC[8] = {1, 12, 32, 84, 100, 108, 104, 120};
        for (int c = 0; c < 8 && !done; c++) for (int k = (c == 0 ? 1 : 0); k < 62 && !done; k++) { ull m = PC[c] << k; if (m > target) break;
            for (int e = 0; e < npool && !done; e++) { ull h = pool[e].H; if (h * HC[c] >= target || (target - h * HC[c]) % m) continue;
                int f = find(h, (target - h * HC[c]) / m); if (f < 0) continue;
                if (c > 0) { printf("p=%llu = h*%llu + B*%llu*2^%d with h=%llu, B=%llu (strong core of count %llu, then a chain of %d)\n", target, HC[c], PC[c], k, h, (target - h * HC[c]) / m, HC[c], k); done = 1; break; }
                Ent *E = &pool[f]; int n = E->n, N2 = n + k + 1;
                printf("p=%llu = h + 2^k B with h=%llu, B=%llu, k=%d: base tree parents", target, h, E->B, k);
                for (int j = 1; j < n; j++) printf(" %d", E->par[j]); printf(", root %d; order %d", E->t0, N2);
                if (N2 <= 22) { unsigned o[MAXN] = {0};
                    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (E->par[j] == i) o[j] |= 1u << i; else o[i] |= 1u << j; }
                    for (int i = 0; i < n; i++) for (int j = 0; j < k; j++) o[i] |= 1u << (n + j);
                    for (int j = 0; j < k; j++) for (int j2 = j + 1; j2 < k; j2++) o[n + j] |= 1u << (n + j2);
                    int t = n + k; o[t] = 1u << E->t0; for (int v = 0; v < n + k; v++) if (v != E->t0) o[v] |= 1u << t;
                    ull Hv = Hfull(N2, o); printf("; direct DP count %llu %s", Hv, Hv == target ? "OK" : "MISMATCH"); }
                printf("\n"); done = 1; } }
        if (!done) printf("p=%llu: no placement with this pool\n", target); }
}
