/* AUDIT program (independent of the author's code).
 * Theorem Q of the note: exhaustive check on ALL labelled tournaments on 2m+1 vertices that are
 * invariant under the fixed anti-involution psi = (1 2)(3 4)...(2m-1 2m), psi(0) = 0.
 *
 * Everything is computed from the definitions:
 *   - invariant cycles are found by enumerating every directed cycle (DFS rooted at the least vertex)
 *     and testing psi(C) = C arc by arc (the half-path form of the note is NOT used);
 *   - H by subset dynamic programming;
 *   - the fold (Z, D_1) from its definition, k-vertex paths of Z ending in D_1 by DFS.
 * Checks per tournament:
 *   (i)   F_k = #{k-vertex paths of Z ending in D_1} (mod 2), k = 1..m          [Theorem Q (2)]
 *   (ii)  c_{2k+1} = F_k (mod 2)                                                  [pairing of non-invariant cycles]
 *   (iii) H = 1 + 2 c_odd (mod 4) and H = 1 + 2 #{paths of Z ending in D_1} (mod 4) [Theorem Q (3)]
 *   (iv)  Z is a tournament; every invariant odd cycle passes through vertex 0
 *   (v)   D_1 empty => all F_k even; D_1 full => F_k = C(m,k) mod 2            [Corollaries (a), (b)]
 *   (vi)  balance: for fixed arcs outside the pairs, exactly half of the 2^m choices of the arcs
 *         inside the pairs give H = 1 (mod 4)                                    [Corollary (c)]
 * usage: aud_thmq m     (m = 1..4)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 11
static int n, m;
static int adj[MAXN][MAXN];
static int psi[MAXN];
static long long cyc[MAXN + 1];   /* cycles by length */
static long long inv[MAXN + 1];   /* invariant cycles by length */
static long long inv_missing0;    /* invariant odd cycles avoiding vertex 0 */
static int path[MAXN], onpath[MAXN], root;

static void test_cycle(int len) {
    /* path[0..len-1] is a directed cycle (path[len-1] -> path[0]) */
    cyc[len]++;
    int succ[MAXN];
    int i;
    for (i = 0; i < n; i++) succ[i] = -1;
    for (i = 0; i < len; i++) succ[path[i]] = path[(i + 1) % len];
    /* psi(C) = C iff for every arc u -> w of C, psi(w) -> psi(u) is an arc of C */
    for (i = 0; i < len; i++) {
        int u = path[i], w = succ[u];
        if (succ[psi[w]] != psi[u]) return;
    }
    inv[len]++;
    if (len & 1) {
        int has0 = 0;
        for (i = 0; i < len; i++) if (path[i] == 0) has0 = 1;
        if (!has0) inv_missing0++;
    }
}

static void dfs(int len) {
    int u = path[len - 1], w;
    if (len >= 3 && adj[u][root]) test_cycle(len);
    for (w = root + 1; w < n; w++) {
        if (!onpath[w] && adj[u][w]) {
            onpath[w] = 1; path[len] = w;
            dfs(len + 1);
            onpath[w] = 0;
        }
    }
}

static unsigned long long hampaths(void) {
    static unsigned long long dp[1 << MAXN][MAXN];
    int full = (1 << n) - 1, s, v, w;
    for (s = 0; s <= full; s++) for (v = 0; v < n; v++) dp[s][v] = 0;
    for (v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (s = 1; s <= full; s++)
        for (v = 0; v < n; v++) {
            unsigned long long x = dp[s][v];
            if (!x) continue;
            for (w = 0; w < n; w++)
                if (!(s >> w & 1) && adj[v][w]) dp[s | 1 << w][w] += x;
        }
    unsigned long long h = 0;
    for (v = 0; v < n; v++) h += dp[full][v];
    return h;
}

/* fold */
static int Z[MAXN][MAXN], inD1[MAXN];
static long long zpaths[MAXN + 1];  /* k-vertex paths of Z ending in D_1 */
static int zon[MAXN];
static void zdfs(int last, int len) {
    int j;
    if (inD1[last]) zpaths[len]++;
    for (j = 0; j < m; j++)
        if (!zon[j] && Z[last][j]) { zon[j] = 1; zdfs(j, len + 1); zon[j] = 0; }
}

int main(int argc, char **argv) {
    m = argc > 1 ? atoi(argv[1]) : 2;
    if (m < 1 || m > 4) { fprintf(stderr, "m in 1..4\n"); return 1; }
    n = 2 * m + 1;
    int i, j, k;
    psi[0] = 0;
    for (i = 0; i < m; i++) { psi[1 + 2 * i] = 2 + 2 * i; psi[2 + 2 * i] = 1 + 2 * i; }
    int nbits_outer = m + m * (m - 1);   /* arcs at vertex 0, and two bits per pair of pairs */
    long long total = 0, viol_Q2 = 0, viol_pair = 0, viol_H_codd = 0, viol_H_fold = 0, viol_Ztour = 0;
    long long viol_a = 0, viol_b = 0, viol_bal = 0, not_antiinv = 0, missing0 = 0;
    long long cntH1 = 0, cntH3 = 0;
    long long binom[8][8];
    for (i = 0; i < 8; i++) { binom[i][0] = 1; for (j = 1; j < 8; j++) binom[i][j] = (i == 0) ? 0 : binom[i-1][j-1] + binom[i-1][j]; }
    unsigned long long outer;
    for (outer = 0; outer < (1ULL << nbits_outer); outer++) {
        int h1_count = 0;
        unsigned int inner;
        for (inner = 0; inner < (1u << m); inner++) {
            memset(adj, 0, sizeof adj);
            int bit = 0;
            for (i = 0; i < m; i++) {
                int a = 1 + 2 * i, b = 2 + 2 * i;
                if (outer >> bit & 1) { adj[0][b] = 1; adj[a][0] = 1; } else { adj[0][a] = 1; adj[b][0] = 1; }
                bit++;
                if (inner >> i & 1) adj[b][a] = 1; else adj[a][b] = 1;
            }
            for (i = 0; i < m; i++) for (j = i + 1; j < m; j++) {
                int a = 1 + 2 * i, b = 1 + 2 * j;           /* arc {a,b}; its image is {psi b, psi a} */
                if (outer >> bit & 1) { adj[a][b] = 1; adj[psi[b]][psi[a]] = 1; }
                else                  { adj[b][a] = 1; adj[psi[a]][psi[b]] = 1; }
                bit++;
                int c = 2 + 2 * j;                           /* arc {a,c}; its image is {psi c, psi a} */
                if (outer >> bit & 1) { adj[a][c] = 1; adj[psi[c]][psi[a]] = 1; }
                else                  { adj[c][a] = 1; adj[psi[a]][psi[c]] = 1; }
                bit++;
            }
            /* sanity: tournament and anti-involution */
            int ok = 1;
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
                if (i != j && adj[i][j] + adj[j][i] != 1) ok = 0;
                if (i != j && adj[i][j] != adj[psi[j]][psi[i]]) ok = 0;
            }
            if (!ok) { not_antiinv++; continue; }
            total++;
            /* cycles from the definition */
            for (i = 0; i <= n; i++) cyc[i] = inv[i] = 0;
            inv_missing0 = 0;
            for (root = 0; root < n; root++) {
                memset(onpath, 0, sizeof onpath);
                path[0] = root; onpath[root] = 1;
                dfs(1);
            }
            missing0 += inv_missing0;
            unsigned long long H = hampaths();
            long long codd = 0, F = 0;
            for (k = 1; k <= m; k++) { codd += cyc[2 * k + 1]; F += inv[2 * k + 1]; }
            /* fold */
            int plus[MAXN], minus[MAXN];
            for (i = 0; i < m; i++) {
                int a = 1 + 2 * i, b = 2 + 2 * i;
                if (adj[0][a]) { plus[i] = a; minus[i] = b; } else { plus[i] = b; minus[i] = a; }
                inD1[i] = adj[plus[i]][minus[i]];
            }
            for (i = 0; i < m; i++) for (j = 0; j < m; j++) {
                Z[i][j] = 0;
                if (i != j) Z[i][j] = (adj[plus[i]][plus[j]] + adj[plus[i]][minus[j]]) == 1;
            }
            for (i = 0; i < m; i++) for (j = i + 1; j < m; j++) if (Z[i][j] + Z[j][i] != 1) viol_Ztour++;
            for (k = 0; k <= m; k++) zpaths[k] = 0;
            for (i = 0; i < m; i++) { memset(zon, 0, sizeof zon); zon[i] = 1; zdfs(i, 1); }
            long long zp = 0; int nd1 = 0;
            for (k = 1; k <= m; k++) zp += zpaths[k];
            for (i = 0; i < m; i++) nd1 += inD1[i];
            for (k = 1; k <= m; k++) {
                if ((inv[2 * k + 1] ^ zpaths[k]) & 1) viol_Q2++;
                if ((cyc[2 * k + 1] ^ inv[2 * k + 1]) & 1) viol_pair++;
                if (nd1 == 0 && (inv[2 * k + 1] & 1)) viol_a++;
                if (nd1 == m && ((inv[2 * k + 1] ^ binom[m][k]) & 1)) viol_b++;
            }
            if ((H & 3) != ((1 + 2 * codd) & 3)) viol_H_codd++;
            if ((H & 3) != (unsigned long long)((1 + 2 * zp) & 3)) viol_H_fold++;
            if ((H & 3) == 1) { cntH1++; h1_count++; } else cntH3++;
        }
        if (h1_count != (1 << (m - 1))) viol_bal++;
    }
    printf("order %d (m = %d): %lld labelled tournaments invariant under the fixed anti-involution (not anti-invariant by construction error: %lld)\n",
           n, m, total, not_antiinv);
    printf("  Z not a tournament: %lld\n", viol_Ztour);
    printf("  invariant odd cycles avoiding the fixed vertex: %lld\n", missing0);
    printf("  violations of F_k = #{k-vertex paths of Z ending in D_1} (mod 2), any k: %lld\n", viol_Q2);
    printf("  violations of c_(2k+1) = F_k (mod 2): %lld\n", viol_pair);
    printf("  violations of H = 1 + 2 c_odd (mod 4): %lld\n", viol_H_codd);
    printf("  violations of H = 1 + 2 #{paths of Z ending in D_1} (mod 4): %lld\n", viol_H_fold);
    printf("  violations of corollary (a) [D_1 empty => F_k even]: %lld; of (b) [D_1 full => F_k = C(m,k)]: %lld\n", viol_a, viol_b);
    printf("  violations of corollary (c) [exactly half of the D_1 choices give H = 1 mod 4]: %lld of %llu outer configurations\n",
           viol_bal, 1ULL << nbits_outer);
    printf("  H = 1 (mod 4): %lld; H = 3 (mod 4): %lld\n", cntH1, cntH3);
    return 0;
}
