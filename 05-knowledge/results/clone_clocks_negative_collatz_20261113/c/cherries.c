/* cherries.c -- targeted search for the C2-line two-clock inside the two-cherry family.
   Skeleton: vertex u; in-set I (|I| = p) and out-set O (|O| = q), n = p + q + 1.  All I -> u -> O, all I -> O except
   exactly four reversed arcs forming two in-cherries s1,s2 -> t1 and s3,s4 -> t2 (t1 != t2 in I, s's distinct in O).
   T[I] and T[O] are transitive with at most B reversed arcs each (B = 0, 1 or 2).  For every skeleton whose u is a
   ternary clock of digit 1, every other vertex v is tested; C2 = (sigma, sig2) = (1, 0).
   usage: cherries n B part nparts   ; prints CLOCK1 counts, GEN lines for digit-1 skeletons and C2LINE hits.      */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int adj[24][24];
static ll *dp; static size_t dpcap = 0;
static void ensure(size_t need) { if (need > dpcap) { free(dp); dpcap = need; dp = malloc(dpcap * sizeof(ll)); } }
static int table2(int u, int v, ll t[32][32]) {      /* v < 0 : single-module table in t[i][0] */
    int inu = 0, inv = 0; for (int x = 0; x < n; x++) { if (adj[x][u]) inu++; if (v >= 0 && adj[x][v]) inv++; }
    int CU, CV;
    if (v < 0) { CU = inu + 1; CV = 0; }
    else if (adj[v][u]) { CV = inv + 1; CU = inu + inv + 1; } else { CU = inu + 1; CV = inv + inu + 1; }
    if (CU > 30) CU = 30; if (CV > 30) CV = 30;
    int others[24], m = 0; for (int x = 0; x < n; x++) if (x != u && x != v) others[m++] = x;
    int pos[24]; for (int i = 0; i < m; i++) pos[others[i]] = i;
    size_t S = (size_t)1 << m; size_t W = (size_t)(CU + 1) * (CV + 1) * n;
    ensure(S * W); memset(dp, 0, S * W * sizeof(ll));
    #define I(mask,cu,cv,last) ((size_t)(mask) * W + ((size_t)(cu) * (CV + 1) + (cv)) * n + (last))
    dp[I(0,1,0,u)] = 1; if (v >= 0) dp[I(0,0,1,v)] = 1; for (int i = 0; i < m; i++) dp[I((size_t)1 << i,0,0,others[i])] = 1;
    memset(t, 0, sizeof(ll) * 32 * 32);
    for (size_t mask = 0; mask < S; mask++) for (int cu = 0; cu <= CU; cu++) for (int cv = 0; cv <= CV; cv++) for (int last = 0; last < n; last++) {
        ll val = dp[I(mask,cu,cv,last)]; if (!val) continue;
        if (mask == S - 1) t[cu][cv] += val;
        for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
            if (w == u) { if (last == u || cu == CU) continue; dp[I(mask,cu+1,cv,u)] += val; }
            else if (w == v) { if (last == v || cv == CV) continue; dp[I(mask,cu,cv+1,v)] += val; }
            else { size_t b = (size_t)1 << pos[w]; if (mask & b) continue; dp[I(mask|b,cu,cv,w)] += val; } }
    }
    return 1;
}
static long nsk = 0, nclock1 = 0, ngen = 0, nc2 = 0;
static void test_skeleton(int u, int *Iv, int p, int *Ov, int q, const char *tag) {
    ll t[32][32]; nsk++;
    table2(u, -1, t);
    for (int i = 4; i < 32; i++) if (t[i][0]) return;
    ll A = t[3][0]; if (A <= 0 || t[2][0] != 3*A || t[1][0] != 3*A - 1) return;   /* digit-1 clock */
    nclock1++;
    for (int v = 0; v < n; v++) { if (v == u) continue;
        table2(u, v, t);
        int bad = 0; for (int i = 1; i < 32; i++) for (int ip = 1; ip < 32; ip++) if ((i >= 4 || ip >= 3) && t[i][ip]) bad = 1;
        if (bad || t[3][2] <= 0 || t[2][1] != 3*t[3][1] || t[2][2] != 3*t[3][2]) continue;
        ll sig = 3*t[3][1] - t[1][1], sig2 = 3*t[3][2] - t[1][2];
        ngen++;
        printf("GEN n=%d %s v=%d A=%lld sigma=%lld al=%lld sig2=%lld rows=", n, tag, v, t[3][1], sig, t[3][2], sig2);
        for (int a = 0; a < n; a++) { int r = 0; for (int b = 0; b < n; b++) if (adj[a][b]) r |= 1 << b; printf("%d%s", r, a + 1 < n ? "," : "\n"); }
        if (sig == 1 && sig2 == 0) { nc2++; printf("C2LINE\n"); fflush(stdout); }
    }
}
/* enumerate near-transitive orders on a set: list of back-arc sets with <= B arcs among k vertices */
static int nback(int k, int B, int out[][2][2], int *cnt) {
    int pairs[64][2], np = 0; for (int hi = 1; hi < k; hi++) for (int lo = 0; lo < hi; lo++) { pairs[np][0] = hi; pairs[np][1] = lo; np++; }
    int c = 0;
    out[c][0][0] = -1; out[c][1][0] = -1; c++;
    if (B >= 1) for (int a = 0; a < np; a++) { out[c][0][0] = pairs[a][0]; out[c][0][1] = pairs[a][1]; out[c][1][0] = -1; c++; }
    if (B >= 2) for (int a = 0; a < np; a++) for (int b = a+1; b < np; b++) { out[c][0][0] = pairs[a][0]; out[c][0][1] = pairs[a][1]; out[c][1][0] = pairs[b][0]; out[c][1][1] = pairs[b][1]; c++; }
    *cnt = c; return c;
}
static int bI[4000][2][2], bO[4000][2][2];
int main(int argc, char **argv) {
    n = atoi(argv[1]); int B = atoi(argv[2]); int part = argc > 3 ? atoi(argv[3]) : 0, nparts = argc > 4 ? atoi(argv[4]) : 1;
    long counter = 0;
    for (int p = 2; p <= n - 5; p++) { int q = n - 1 - p; if (q < 4) continue;
        int cI, cO; nback(p, B, bI, &cI); nback(q, B, bO, &cO);
        /* vertex labels: I = 0..p-1 (transitive order), u = p, O = p+1..n-1 */
        int u = p;
        for (int t1 = 0; t1 < p; t1++) for (int t2 = t1 + 1; t2 < p; t2++)
        for (int s1 = 0; s1 < q; s1++) for (int s2 = s1 + 1; s2 < q; s2++)
        for (int s3 = 0; s3 < q; s3++) for (int s4 = s3 + 1; s4 < q; s4++) {
            if (s3 == s1 || s3 == s2 || s4 == s1 || s4 == s2) continue;
            for (int iI = 0; iI < cI; iI++) for (int iO = 0; iO < cO; iO++) {
                if (counter++ % nparts != part) continue;
                memset(adj, 0, sizeof adj);
                for (int a = 0; a < p; a++) for (int b = a + 1; b < p; b++) adj[a][b] = 1;
                for (int k = 0; k < 2; k++) if (bI[iI][k][0] >= 0) { int hi = bI[iI][k][0], lo = bI[iI][k][1]; adj[lo][hi] = 0; adj[hi][lo] = 1; }
                for (int a = 0; a < q; a++) for (int b = a + 1; b < q; b++) adj[p+1+a][p+1+b] = 1;
                for (int k = 0; k < 2; k++) if (bO[iO][k][0] >= 0) { int hi = p+1+bO[iO][k][0], lo = p+1+bO[iO][k][1]; adj[lo][hi] = 0; adj[hi][lo] = 1; }
                for (int a = 0; a < p; a++) { adj[a][u] = 1; for (int b = 0; b < q; b++) adj[a][p+1+b] = 1; }
                for (int b = 0; b < q; b++) adj[u][p+1+b] = 1;
                int S[4] = {s1, s2, s3, s4}, T[4] = {t1, t1, t2, t2};
                for (int k = 0; k < 4; k++) { int o = p + 1 + S[k], i = T[k]; adj[i][o] = 0; adj[o][i] = 1; }
                char tag[64]; snprintf(tag, sizeof tag, "p=%d u=%d", p, u);
                test_skeleton(u, 0, p, 0, q, tag);
            }
        }
    }
    fprintf(stderr, "n=%d B=%d part %d/%d: skeletons %ld digit1 %ld gen %ld c2 %ld\n", n, B, part, nparts, nsk, nclock1, ngen, nc2);
}
