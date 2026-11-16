/* allsearch2.c <n> <k>: all rooted tournaments of order n (gentourng on stdin), extensions by k = 1, 2 or 3 new vertices
   (arbitrary contacts to old vertices; among new vertices every orientation).  Per rooted base: profile (H,E,O,g); whether
   some extension has count 3H+2 (has2) or 3H+4 (has4); whether the three digits share a common (X',g') (weak = complete
   after a trailing sink) or a common (E',O',g').  Prints ANY-COMMON lines and a table by (D==0, H). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[16];
static i64 Hm[1<<12], rs[1<<12]; static i64 dp[1<<12][12], dpr[1<<12][12];
static void profile(int root, i64 *H, i64 *E, i64 *O, i64 *g) {
    int N = 1 << n;
    for (int m = 0; m < N; m++) { Hm[m] = 0; for (int l = 0; l < n; l++) dp[m][l] = 0, dpr[m][l] = 0; }
    for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
    dpr[1<<root][root] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) {
        i64 a = dp[m][l], b = dpr[m][l]; if (!a && !b) continue;
        u64 nx = adj[l] & ~(u64)m;
        while (nx) { int v = __builtin_ctzll(nx); nx &= nx - 1; dp[m|(1<<v)][v] += a; dpr[m|(1<<v)][v] += b; }
    }
    for (int m = 1; m < N; m++) { i64 s = 0, r = 0; for (int l = 0; l < n; l++) s += dp[m][l], r += dpr[m][l]; Hm[m] = s; rs[m] = r; }
    Hm[0] = 1; *H = Hm[N-1]; *g = Hm[(N-1) & ~(1<<root)];
    i64 e = 0, o = 0;
    for (int m = 1; m < N; m++) if ((m >> root) & 1) { i64 w = rs[m] * Hm[(N-1) & ~m]; if (__builtin_popcount(m) % 2 == 0) e += w; else o += w; }
    *E = e; *O = o;
}
typedef struct { int d; i64 X, g, E, O; } Hit;
int main(int argc, char **argv) {
    int nn = atoi(argv[1]), K = atoi(argv[2]); char line[256]; long tcount = 0;
    long tab[2][64][4] = {{{0}}}; /* [D==0][H][bases, has2, has4, common] */
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0;
        if (!len) continue; tcount++;
        int base[16] = {0}; int k = 0;
        for (int i = 0; i < nn; i++) for (int j = i+1; j < nn; j++) { if (line[k++] == '1') base[i] |= 1<<j; else base[j] |= 1<<i; }
        for (int root = 0; root < nn; root++) {
            n = nn; for (int i = 0; i < nn; i++) adj[i] = base[i];
            i64 H0,E0,O0,g0; profile(root,&H0,&E0,&O0,&g0);
            int dz = (E0 - O0 == 0); int hb = H0 < 63 ? (int)H0 : 63;
            tab[dz][hb][0]++;
            int N0 = 1 << nn; Hit *hits = malloc(sizeof(Hit) * 4000000); int cnt = 0;
            long total = 1; for (int t = 0; t < K; t++) total *= N0; long orient = 1 << (K*(K-1)/2);
            for (long code = 0; code < total; code++) for (long oc = 0; oc < orient; oc++) {
                n = nn + K; long c = code; int cm[3];
                for (int t = 0; t < K; t++) { cm[t] = c % N0; c /= N0; }
                for (int i = 0; i < nn; i++) { adj[i] = base[i]; for (int t = 0; t < K; t++) if (!((cm[t]>>i)&1)) adj[i] |= 1ULL<<(nn+t); }
                for (int t = 0; t < K; t++) adj[nn+t] = cm[t];
                int bit = 0; for (int t = 0; t < K; t++) for (int u = t+1; u < K; u++) { if ((oc>>bit)&1) adj[nn+u] |= 1ULL<<(nn+t); else adj[nn+t] |= 1ULL<<(nn+u); bit++; }
                i64 H,E,O,g; profile(root,&H,&E,&O,&g); i64 d = H - 3*H0;
                if (d == 0 || d == 2 || d == 4) { hits[cnt].d = d; hits[cnt].X = E+O; hits[cnt].g = g; hits[cnt].E = E; hits[cnt].O = O; cnt++; }
            }
            int has2 = 0, has4 = 0; for (int i = 0; i < cnt; i++) { if (hits[i].d == 2) has2 = 1; if (hits[i].d == 4) has4 = 1; }
            tab[dz][hb][1] += has2; tab[dz][hb][2] += has4;
            int common = 0, commonS = 0;
            for (int i = 0; i < cnt && !common; i++) if (hits[i].d == 0) { int f2 = 0, f4 = 0; for (int j = 0; j < cnt; j++) if (hits[j].X == hits[i].X && hits[j].g == hits[i].g) { if (hits[j].d == 2) f2 = 1; if (hits[j].d == 4) f4 = 1; } if (f2 && f4) common = 1; }
            for (int i = 0; i < cnt && !commonS; i++) if (hits[i].d == 0) { int f2 = 0, f4 = 0; for (int j = 0; j < cnt; j++) if (hits[j].E == hits[i].E && hits[j].O == hits[i].O && hits[j].g == hits[i].g) { if (hits[j].d == 2) f2 = 1; if (hits[j].d == 4) f4 = 1; } if (f2 && f4) commonS = 1; }
            tab[dz][hb][3] += common;
            if (common) printf("ANY-COMMON %s root=%d (H,E,O,g)=(%lld,%lld,%lld,%lld) D=%lld strictcommon=%d\n", line, root, H0,E0,O0,g0, E0-O0, commonS);
            free(hits);
        }
    }
    printf("order %d, %d added vertices, %ld tournaments. Table [D==0?][H]: bases, with a 3H+2 extension, with a 3H+4 extension, with a common (X',g') triple\n", nn, K, tcount);
    for (int dz = 0; dz < 2; dz++) for (int h = 0; h < 64; h++) if (tab[dz][h][0]) printf("  D%s H=%d: %ld %ld %ld %ld\n", dz ? "=0 " : "!=0", h, tab[dz][h][0], tab[dz][h][1], tab[dz][h][2], tab[dz][h][3]);
    long s[2][4] = {{0}}; for (int dz = 0; dz < 2; dz++) for (int h = 0; h < 64; h++) for (int c = 0; c < 4; c++) s[dz][c] += tab[dz][h][c];
    printf("  totals D!=0: %ld bases, has2 %ld, has4 %ld, common %ld; D=0: %ld bases, has2 %ld, has4 %ld, common %ld\n", s[0][0],s[0][1],s[0][2],s[0][3], s[1][0],s[1][1],s[1][2],s[1][3]);
    return 0;
}
