/* backscan.c -- near-transitive skeleton scan: vertices 0..n-1 in transitive order (i->j for i<j) with a set of
   B back arcs (hi->lo) reversed, |B| <= maxb.  For each u with matching number 2: ternary clock test
   (a_4=0, a_2=3a_3), digit sigma = 3a_3 - a_1.  For sigma = 1 hits, every v is tested for the TWO-CLOCK condition.
   Also reports general two-clock tables (X_1, X_2 both ternary-shaped with zero 2^j terms), printing digits.
   usage: backscan n maxb [part nparts]                                                                     */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int adj[16][16];
static ll *dp1, *dp2;
static void visits1(int u, ll *a) {
    int others[16], m = 0; for (int x = 0; x < n; x++) if (x != u) others[m++] = x;
    int pos[16]; for (int i = 0; i < m; i++) pos[others[i]] = i;
    size_t S = 1u << m;
    #define I1(mask,c,last) ((((size_t)(mask))*6+(c))*n+(last))
    memset(dp1, 0, sizeof(ll) * S * 6 * n);
    dp1[I1(0,1,u)] = 1; for (int i = 0; i < m; i++) dp1[I1(1u<<i,0,others[i])] = 1;
    for (int c = 0; c <= 5; c++) a[c] = 0;
    for (size_t mask = 0; mask < S; mask++) for (int c = 0; c <= 5; c++) for (int last = 0; last < n; last++) {
        ll val = dp1[I1(mask,c,last)]; if (!val) continue;
        if (mask == S-1) a[c] += val;
        for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
            if (w == u) { if (last == u || c == 5) continue; dp1[I1(mask,c+1,u)] += val; }
            else { size_t b = 1u << pos[w]; if (mask & b) continue; dp1[I1(mask|b,c,w)] += val; } }
    }
}
static void visits2(int u, int v, ll t[6][5]) {
    int others[16], m = 0; for (int x = 0; x < n; x++) if (x != u && x != v) others[m++] = x;
    int pos[16]; for (int i = 0; i < m; i++) pos[others[i]] = i;
    size_t S = 1u << m;
    #define I2(mask,cu,cv,last) (((((size_t)(mask))*6+(cu))*5+(cv))*n+(last))
    memset(dp2, 0, sizeof(ll) * S * 6 * 5 * n);
    dp2[I2(0,1,0,u)] = 1; dp2[I2(0,0,1,v)] = 1; for (int i = 0; i < m; i++) dp2[I2(1u<<i,0,0,others[i])] = 1;
    memset(t, 0, sizeof(ll)*30);
    for (size_t mask = 0; mask < S; mask++) for (int cu = 0; cu <= 5; cu++) for (int cv = 0; cv <= 4; cv++) for (int last = 0; last < n; last++) {
        ll val = dp2[I2(mask,cu,cv,last)]; if (!val) continue;
        if (mask == S-1) t[cu][cv] += val;
        for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
            if (w == u) { if (last == u || cu == 5) continue; dp2[I2(mask,cu+1,cv,u)] += val; }
            else if (w == v) { if (last == v || cv == 4) continue; dp2[I2(mask,cu,cv+1,v)] += val; }
            else { size_t b = 1u << pos[w]; if (mask & b) continue; dp2[I2(mask|b,cu,cv,w)] += val; } }
    }
}
static int matchnum(int u) { /* max matching O(u)->I(u), capped at 3 (augmenting paths) */
    int O[16], I[16], no = 0, ni = 0;
    for (int x = 0; x < n; x++) if (x != u) { if (adj[u][x]) O[no++] = x; else I[ni++] = x; }
    int matchI[16]; for (int i = 0; i < ni; i++) matchI[i] = -1; int res = 0;
    for (int o = 0; o < no; o++) {
        int seen[16] = {0};
        /* DFS augment */
        int stack_o[32], sp = 0; int parentI[16]; for (int i = 0; i < ni; i++) parentI[i] = -2;
        /* simple recursive via function pointer not allowed in C easily; use iterative BFS over alternating paths */
        int fromO[16]; for (int i = 0; i < 16; i++) fromO[i] = -1;
        int queue[64], qh = 0, qt = 0; queue[qt++] = o; int found = -1;
        int prevI[16]; for (int i = 0; i < ni; i++) prevI[i] = -1;
        int visitedO[16] = {0}; visitedO[o] = 1;
        while (qh < qt && found < 0) { int oo = queue[qh++];
            for (int i = 0; i < ni; i++) if (!seen[i] && adj[O[oo]][I[i]]) { seen[i] = 1; prevI[i] = oo;
                if (matchI[i] < 0) { found = i; break; }
                int o2 = matchI[i]; if (!visitedO[o2]) { visitedO[o2] = 1; queue[qt++] = o2; } } }
        if (found >= 0) { res++; int i = found; while (1) { int oo = prevI[i]; int prevmatch = -1;
                for (int k = 0; k < ni; k++) if (matchI[k] == oo) prevmatch = k;
                matchI[i] = oo; if (prevmatch < 0) break; i = prevmatch; } }
        (void)stack_o; (void)sp; (void)parentI; (void)fromO;
    }
    return res;
}
static long hits_t = 0, hits_two = 0, hits_gen = 0;
static void test(int nb, int *bk) {
    for (int u = 0; u < n; u++) {
        if (matchnum(u) != 2) continue;
        ll a[6]; visits1(u, a);
        if (a[4] || a[5] || a[3] <= 0 || a[2] != 3*a[3]) continue;
        ll A = a[3], sig = 3*a[3] - a[1];
        hits_t++;
        printf("TERN n=%d u=%d A=%lld sigma=%lld back=", n, u, A, sig); for (int i = 0; i < nb; i++) printf("%d>%d,", bk[2*i], bk[2*i+1]); printf("\n");
        for (int v = 0; v < n; v++) { if (v == u) continue;
            if (matchnum(v) != 1) continue;
            ll t[6][5]; visits2(u, v, t);
            int bad = 0; for (int i = 1; i <= 5; i++) for (int ip = 1; ip <= 4; ip++) if ((i >= 4 || ip >= 3) && t[i][ip]) bad = 1;
            if (bad || t[3][2] <= 0) continue;
            ll al = t[3][2];
            if (t[2][2] != 3*al) continue;           /* X_2 has no 2^j term */
            ll sig2 = 3*al - t[1][2];                /* X_2 = al 3^j - sig2 */
            hits_gen++;
            printf("GEN n=%d u=%d v=%d A=%lld sigma=%lld al=%lld sig2=%lld  H(j,k)=3^j(%lld*2^k+%lld) - %lld*2^k - (%lld-2*%lld) back=",
                   n, u, v, A, sig, al, sig2, al, A-2*al, sig2, sig, sig2);
            for (int i = 0; i < nb; i++) printf("%d>%d,", bk[2*i], bk[2*i+1]); printf("\n");
            if (sig == 1 && sig2 == 0) { hits_two++; printf("TWO!! n=%d u=%d v=%d A=%lld al=%lld beta=%lld\n", n, u, v, A, al, A - 2*al); }
        }
    }
    fflush(stdout);
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); int maxb = atoi(argv[2]); int part = argc > 3 ? atoi(argv[3]) : 0, np = argc > 4 ? atoi(argv[4]) : 1;
    dp1 = malloc(sizeof(ll) * ((size_t)1 << (n-1)) * 6 * n);
    dp2 = malloc(sizeof(ll) * ((size_t)1 << (n-2)) * 6 * 5 * n);
    int pairs[256][2], np2 = 0; for (int hi = 1; hi < n; hi++) for (int lo = 0; lo < hi; lo++) { pairs[np2][0] = hi; pairs[np2][1] = lo; np2++; }
    long counter = 0;
    for (int nb = 1; nb <= maxb; nb++) {
        int idx[16]; for (int i = 0; i < nb; i++) idx[i] = i;
        while (1) {
            if (counter++ % np == part) {
                memset(adj, 0, sizeof adj);
                for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) adj[i][j] = 1;
                int bk[32];
                for (int i = 0; i < nb; i++) { int hi = pairs[idx[i]][0], lo = pairs[idx[i]][1]; adj[lo][hi] = 0; adj[hi][lo] = 1; bk[2*i] = hi; bk[2*i+1] = lo; }
                test(nb, bk);
            }
            int k = nb - 1; while (k >= 0 && idx[k] == np2 - nb + k) k--;
            if (k < 0) break; idx[k]++; for (int i = k+1; i < nb; i++) idx[i] = idx[i-1] + 1;
        }
    }
    fprintf(stderr, "n=%d maxb=%d part %d/%d: tern=%ld gen=%ld two=%ld\n", n, maxb, part, np, hits_t, hits_gen, hits_two);
}
