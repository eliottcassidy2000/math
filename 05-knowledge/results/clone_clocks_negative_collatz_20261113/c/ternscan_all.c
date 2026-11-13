/* ternscan.c -- exhaustive scan (gentourng input) for transitive-clone TERNARY CLOCKS at a vertex u:
   visit numbers a_1..a_3 of u with a_4 = 0 (matching number 2) and a_2 = 3 a_3, so that
        H(S[u <- TT_m]) = A 3^m + C,   A = a_3,  C = a_1 - 3 a_3        (all m >= 1)
   i.e. an exact C2-type recursion H -> 3H + 2 sigma with digit sigma = -C (fixed point -sigma).
   For every hit with C = -1 (sigma = 1, the C2 line through -1) also test every second vertex v for the
   TWO-CLOCK condition a[.][1] = (3A-1,3A,A), a[.][2] = (3al,3al,al), no u-4 or v-3 visits:
        H(S[u<-TT_j, v<-TT_k]) = 3^j (al 2^k + A - 2 al) - 1.
   Output lines:  TERN n code u A C   |  TWO n code u v A al beta  |  NEAR (v-twice sector off by small amount) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int adj[16][16];
static ll *dp1;
/* single-module visit numbers with cap 4 */
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
static ll *dp2;
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
static int matching2(int u) { /* is the max matching from O(u) to I(u) of size exactly 2 ? (>=2 and <3) */
    int O[16], I[16], no = 0, ni = 0;
    for (int x = 0; x < n; x++) if (x != u) { if (adj[u][x]) O[no++] = x; else I[ni++] = x; }
    /* brute force: size>=2 exists? size>=3 exists? */
    int has2 = 0, has3 = 0;
    for (int a = 0; a < no; a++) for (int b = a+1; b < no; b++)
      for (int p = 0; p < ni; p++) for (int q = 0; q < ni; q++) { if (p == q) continue;
        if (adj[O[a]][I[p]] && adj[O[b]][I[q]]) { has2 = 1;
            for (int c = b+1; c < no && !has3; c++) for (int r = 0; r < ni; r++) if (r != p && r != q && adj[O[c]][I[r]]) { has3 = 1; break; }
            if (!has3) { /* also c<a or between */
                for (int c = 0; c < no && !has3; c++) { if (c == a || c == b) continue; for (int r = 0; r < ni; r++) if (r != p && r != q && adj[O[c]][I[r]]) { has3 = 1; break; } }
            }
        } }
    return has2 && !has3;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]);
    dp1 = malloc(sizeof(ll) * ((size_t)1 << (n-1)) * 6 * n);
    dp2 = malloc(sizeof(ll) * ((size_t)1 << (n-2)) * 6 * 5 * n);
    char line[512]; long cnt = 0, tern = 0, two = 0;
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && line[L-1] == '\n') line[--L] = 0;
        if (L != n*(n-1)/2) continue; cnt++;
        memset(adj, 0, sizeof adj); int idx = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[idx++] == '1') adj[i][j] = 1; else adj[j][i] = 1; }
        for (int u = 0; u < n; u++) {
            if (!matching2(u)) continue;
            ll a[6]; visits1(u, a);
            if (a[4] || a[5] || a[3] <= 0) continue;
            if (a[2] != 3*a[3]) continue;
            ll A = a[3], C = a[1] - 3*a[3];
            tern++; printf("TERN %d %s u=%d A=%lld C=%lld\n", n, line, u, A, C);
            /* all digits */
            for (int v = 0; v < n; v++) { if (v == u) continue;
                ll t[6][5]; visits2(u, v, t);
                int bad = 0; for (int i = 1; i <= 5; i++) for (int ip = 1; ip <= 4; ip++) if ((i >= 4 || ip >= 3) && t[i][ip]) bad = 1;
                if (bad) continue;
                ll al = t[3][2]; if (al <= 0) continue;
                if (t[2][2] == 3*al) { printf("GEN %d %s u=%d v=%d A=%lld sigma=%lld al=%lld sig2=%lld\n", n, line, u, v, A, -C, al, 3*al - t[1][2]); if (C == -1 && t[1][2] == 3*al) { two++; printf("TWO\n"); } }
                
            }
        }
        fflush(stdout);
    }
    fprintf(stderr, "n=%d tournaments=%ld tern=%ld two=%ld\n", n, cnt, tern, two);
}
