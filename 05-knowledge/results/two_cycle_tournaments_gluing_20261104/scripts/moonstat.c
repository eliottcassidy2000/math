/* moonstat.c (klein-2026-11-04): deletion statistics of two-cycle tournaments. stdin: "hc digraph6" lines; for each line with hc = 2,
   for every vertex v with Y - v strong and hc(Y - v) <= 2 ("deletable"), report hc(Y - v), |O| (O = out-neighbourhood of v),
   the number of cyclic intervals of O along one Hamiltonian cycle of Y - v, k1 = #out-degree-1 vertices of Y - v inside O,
   head1/tail1 = how many intervals of O start / end (in cycle order) with such a vertex.  Output per class:
   "n od1 ndel : hcp,|O|,nint,k1,head1,tail1 ..." */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 *dp;
static int parse_d6(const char *s, u64 *a) {
    if (s[0] != '&') return -1; int n = s[1] - 63; int bits = n * n; const char *p = s + 2;
    for (int i = 0; i < n; i++) a[i] = 0;
    for (int k = 0; k < bits; k++) { int ch = p[k / 6] - 63; int b = (ch >> (5 - k % 6)) & 1; if (b) a[k / n] |= 1ULL << (k % n); }
    return n;
}
static int strong(int m, u64 *a) {
    u64 all = (1ULL << m) - 1, seen = 1, fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < m; v++) if (fr >> v & 1) nf |= a[v]; nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != all) return 0;
    u64 r[32]; for (int v = 0; v < m; v++) r[v] = 0;
    for (int u = 0; u < m; u++) for (int v = 0; v < m; v++) if (a[u] >> v & 1) r[v] |= 1ULL << u;
    seen = 1; fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < m; v++) if (fr >> v & 1) nf |= r[v]; nf &= ~seen; seen |= nf; fr = nf; }
    return seen == all;
}
/* hc of tournament a on m vertices, reconstruct one cycle into cyc[] */
static u64 hcycle(int m, u64 *a, int *cyc) {
    u64 N = 1ULL << m; memset(dp, 0, sizeof(u64) * N * m); dp[1 * m + 0] = 1;
    for (u64 S = 1; S < N; S += 2) for (int v = 0; v < m; v++) { u64 c = dp[S * m + v]; if (!c) continue; u64 nx = a[v] & ~S; while (nx) { int u = __builtin_ctzll(nx); nx &= nx - 1; dp[(S | (1ULL << u)) * m + u] += c; } }
    u64 hc = 0; int endv = -1; for (int v = 1; v < m; v++) if (a[v] & 1ULL) { hc += dp[(N - 1) * m + v]; if (dp[(N - 1) * m + v]) endv = v; }
    if (hc) { u64 S = N - 1; int cur = endv; cyc[m - 1] = cur;
        for (int pos = m - 1; pos > 0; pos--) { u64 Sp = S & ~(1ULL << cur); int found = -1;
            for (int i = 0; i < m; i++) if ((Sp >> i & 1) && dp[Sp * m + i] && (a[i] >> cur & 1)) { found = i; break; }
            cur = found; S = Sp; cyc[pos - 1] = cur; } }
    return hc;
}
int main(void) {
    char line[512]; u64 adj[32], sub[32]; int cyc[32];
    dp = malloc(sizeof(u64) * (1ULL << 20) * 20);
    while (fgets(line, sizeof line, stdin)) {
        int h; char d6[256]; if (sscanf(line, "%d %255s", &h, d6) != 2 || h != 2) continue;
        int n = parse_d6(d6, adj); int m = n - 1;
        int od1 = 0, id1 = 0, ncut = 0; for (int v = 0; v < n; v++) if (__builtin_popcountll(adj[v]) == 1) od1++;
        for (int w = 0; w < n; w++) { int d = 0; for (int u = 0; u < n; u++) if (adj[u] >> w & 1) d++; if (d == 1) id1++; }
        for (int u = 0; u < n; u++) for (int w = 0; w < n; w++) if (adj[u] >> w & 1) {
            /* is u -> w a cut arc: remove it (not reverse) and test strong connectivity of the digraph */
            u64 a2[32]; for (int i = 0; i < n; i++) a2[i] = adj[i]; a2[u] &= ~(1ULL << w);
            u64 all = (1ULL << n) - 1, seen = 1, fr = 1; while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= a2[v]; nf &= ~seen; seen |= nf; fr = nf; }
            int ok = (seen == all); if (ok) { u64 r[32]; for (int v = 0; v < n; v++) r[v] = 0; for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if (a2[x] >> y & 1) r[y] |= 1ULL << x;
                seen = 1; fr = 1; while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= r[v]; nf &= ~seen; seen |= nf; fr = nf; } ok = (seen == all); }
            if (!ok) ncut++; }
        printf("%d %d %d %d", n, od1, id1, ncut); int ndel = 0; char buf[4096] = "";
        for (int v = 0; v < n; v++) {
            int map[32], k = 0, inv[32]; for (int u = 0; u < n; u++) if (u != v) { inv[u] = k; map[k++] = u; }
            for (int i = 0; i < m; i++) { sub[i] = 0; for (int j = 0; j < m; j++) if (adj[map[i]] >> map[j] & 1) sub[i] |= 1ULL << j; }
            if (!strong(m, sub)) continue;
            u64 hp = hcycle(m, sub, cyc); if (hp < 1 || hp > 2) continue;
            ndel++;
            u64 O = 0; for (int i = 0; i < m; i++) if (adj[v] >> map[i] & 1) O |= 1ULL << i;
            int Osz = __builtin_popcountll(O), nsw = 0, k1 = 0, head1 = 0, tail1 = 0;
            for (int i = 0; i < m; i++) { int a_ = (O >> cyc[i]) & 1, b_ = (O >> cyc[(i + 1) % m]) & 1; if (a_ != b_) nsw++;
                if (a_ && __builtin_popcountll(sub[cyc[i]]) == 1) { k1++; if (!b_) tail1++; if (!((O >> cyc[(i + m - 1) % m]) & 1)) head1++; } }
            int se = 0, ie = 0; for (int i = 0; i < m; i++) { if (__builtin_popcountll(sub[i]) == 1) se++; int d = 0; for (int j = 0; j < m; j++) if (sub[j] >> i & 1) d++; if (d == 1) ie++; }
            char t[80]; snprintf(t, sizeof t, " %llu,%d,%d,%d,%d,%d,%d,%d,%d", hp, Osz, nsw / 2, k1, head1, tail1, se, ie, v); strcat(buf, t);
        }
        printf(" %d :%s\n", ndel, buf);
    }
    return 0;
}
