/* extscan.c -- extend the Session-197 ternary skeleton S60 (8 vertices, u = 3) by e = 1 or 2 new vertices with every
   adjacency pattern; test every vertex v (new or old, v != u) for the exact two-clock condition with u:
   a[.][1] = (3A - sigma, 3A, A), a[.][2] = (3al - sig2, 3al, al) (ternary-shaped), report (sigma, sig2); TWO = (1, 0). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int adj[16][16];
static ll *dp2;
static void visits2(int u, int v, ll t[8][8], int CU, int CV) {
    int others[16], m = 0; for (int x = 0; x < n; x++) if (x != u && x != v) others[m++] = x;
    int pos[16]; for (int i = 0; i < m; i++) pos[others[i]] = i;
    size_t S = 1u << m;
    #define I2(mask,cu,cv,last) (((((size_t)(mask))*8+(cu))*8+(cv))*n+(last))
    memset(dp2, 0, sizeof(ll) * S * 64 * n);
    dp2[I2(0,1,0,u)] = 1; dp2[I2(0,0,1,v)] = 1; for (int i = 0; i < m; i++) dp2[I2(1u<<i,0,0,others[i])] = 1;
    memset(t, 0, sizeof(ll)*64);
    for (size_t mask = 0; mask < S; mask++) for (int cu = 0; cu < 8; cu++) for (int cv = 0; cv < 8; cv++) for (int last = 0; last < n; last++) {
        ll val = dp2[I2(mask,cu,cv,last)]; if (!val) continue;
        if (mask == S-1) t[cu][cv] += val;
        for (int w = 0; w < n; w++) { if (!adj[last][w]) continue;
            if (w == u) { if (last == u || cu == CU) continue; dp2[I2(mask,cu+1,cv,u)] += val; }
            else if (w == v) { if (last == v || cv == CV) continue; dp2[I2(mask,cu,cv+1,v)] += val; }
            else { size_t b = 1u << pos[w]; if (mask & b) continue; dp2[I2(mask|b,cu,cv,w)] += val; } }
    }
}
int main(int argc, char **argv) {
    int e = atoi(argv[1]);
    int base[8] = {2+4+8+64+128, 4+8+16+32+64+128, 8+16+32, 16+32+64+128, 1+32+64+128, 1+64+128, 4+128, 4};
    /* S60 rows: 0->{1,2,3,6,7},1->{2..7},2->{3,4,5},3->{4,5,6,7},4->{0,5,6,7},5->{0,6,7},6->{2,7},7->{2} */
    int u = 3; n = 8 + e;
    dp2 = malloc(sizeof(ll) * ((size_t)1 << (n-2)) * 64 * n);
    long pats = 1L << (8*e + (e==2 ? 1 : 0)); long two = 0, gen = 0;
    for (long p = 0; p < pats; p++) {
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) if (base[i] >> j & 1) adj[i][j] = 1;
        for (int k = 0; k < e; k++) for (int j = 0; j < 8; j++) { int b = (p >> (8*k + j)) & 1; if (b) adj[8+k][j] = 1; else adj[j][8+k] = 1; }
        if (e == 2) { if ((p >> 16) & 1) adj[8][9] = 1; else adj[9][8] = 1; }
        for (int v = 0; v < n; v++) { if (v == u) continue;
            ll t[8][8]; visits2(u, v, t, 7, 7);
            int bad = 0; for (int i = 1; i < 8; i++) for (int ip = 1; ip < 8; ip++) if ((i >= 4 || ip >= 3) && t[i][ip]) bad = 1;
            if (bad) continue;
            ll A = t[3][1], al = t[3][2]; if (A <= 0 || al <= 0) continue;
            if (t[2][1] != 3*A || t[2][2] != 3*al) continue;
            ll sig = 3*A - t[1][1], sig2 = 3*al - t[1][2]; gen++;
            printf("GEN e=%d pat=%ld v=%d A=%lld al=%lld sigma=%lld sig2=%lld c11=%lld c12=%lld\n", e, p, v, A, al, sig, sig2, 2*sig2 - sig, -sig2);
            if (sig == 1 && sig2 == 0) { two++; printf("TWO!!\n"); }
        }
    }
    fprintf(stderr, "e=%d patterns=%ld gen=%ld two=%ld\n", e, pats, gen, two);
}
