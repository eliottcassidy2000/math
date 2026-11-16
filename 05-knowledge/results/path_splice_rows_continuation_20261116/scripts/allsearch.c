/* allsearch.c <n>: read gentourng -q output for order n on stdin; for every tournament and every root, compute the
   profile (H,E,O,g); tabulate the law D + g by (strong, root terminal); enumerate all one- and two-vertex extensions
   (fixed root) and report sources that have extensions of counts 3H, 3H+2, 3H+4 sharing a continuation:
   weak = (X', g') = (3X, 3g); strict = (E',O',g') = (3E,3O,3g).  Also report pairs {3H+2, 3H+4} sharing (3X,3g). */
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
static u64 reachfrom(u64 start, u64 *A) { u64 r = start, p = 0; while (r != p) { p = r; for (int v = 0; v < n; v++) if ((r>>v)&1) r |= A[v]; } return r; }
static int strongT(void) { u64 radj[16]; for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
    u64 full = (1ULL<<n)-1; return reachfrom(1, adj) == full && reachfrom(1, radj) == full; }
static int terminal_root(int root) { /* root in the terminal strong component <=> every vertex reaches root */
    u64 radj[16]; for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
    return reachfrom(1ULL<<root, radj) == ((1ULL<<n)-1); }
int main(int argc, char **argv) {
    int nn = atoi(argv[1]); char line[256]; long tcount = 0;
    long law[2][2][3] = {{{0}}}; /* [strong][terminal][D+g <0, =0, >0] */
    long lawDg[2][2] = {{0}};    /* D == -g counts */
    long srcWeak = 0, srcStrict = 0, pairWeak = 0, pairStrict = 0, srcWeakStrong = 0;
    static i64 keyX[4096*2*2+10], keyG[4096*2*2+10]; /* not used */
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0;
        if (!len) continue; tcount++;
        int base[16] = {0}; int k = 0;
        for (int i = 0; i < nn; i++) for (int j = i+1; j < nn; j++) { if (line[k++] == '1') base[i] |= 1<<j; else base[j] |= 1<<i; }
        for (int root = 0; root < nn; root++) {
            n = nn; for (int i = 0; i < nn; i++) adj[i] = base[i];
            i64 H0,E0,O0,g0; profile(root,&H0,&E0,&O0,&g0);
            int st = strongT(), tr = terminal_root(root);
            i64 Dg = (E0 - O0) + g0; law[st][tr][Dg < 0 ? 0 : (Dg == 0 ? 1 : 2)]++; if (Dg == 0) lawDg[st][tr]++;
            /* extensions: collect (d, X', g', E', O') lists */
            int N0 = 1 << nn; int cap = 2*N0*N0 + N0; 
            i64 *lx = malloc(sizeof(i64)*cap), *lg = malloc(sizeof(i64)*cap), *le = malloc(sizeof(i64)*cap), *lo = malloc(sizeof(i64)*cap); int *ld = malloc(sizeof(int)*cap); int *lstrong = malloc(sizeof(int)*cap); int cnt = 0;
            int lawonly = (argc > 2 && !strcmp(argv[2], "law")); int twomax = lawonly ? 0 : ((argc > 2 && !strcmp(argv[2], "one")) ? 1 : 2); for (int two = 0; two < twomax; two++) for (int cx = 0; cx < N0; cx++) for (int cy = 0; cy < (two ? N0 : 1); cy++) for (int o = 0; o < (two ? 2 : 1); o++) {
                n = nn + 1 + two;
                for (int i = 0; i < nn; i++) { adj[i] = base[i]; if (!((cx>>i)&1)) adj[i] |= 1ULL<<nn; if (two && !((cy>>i)&1)) adj[i] |= 1ULL<<(nn+1); }
                adj[nn] = cx; if (two) { adj[nn+1] = cy; if (o == 0) adj[nn] |= 1ULL<<(nn+1); else adj[nn+1] |= 1ULL<<nn; }
                i64 H,E,O,g; profile(root,&H,&E,&O,&g); i64 d = H - 3*H0;
                if (d == 0 || d == 2 || d == 4) { ld[cnt] = d; lx[cnt] = E+O; lg[cnt] = g; le[cnt] = E; lo[cnt] = O; lstrong[cnt] = strongT(); cnt++; }
            }
            int has[3] = {0,0,0}, hasS[3] = {0,0,0}, hasSt[3] = {0,0,0};
            for (int i = 0; i < cnt; i++) { if (lx[i] == 3*(E0+O0) && lg[i] == 3*g0) { has[ld[i]/2] = 1; if (lstrong[i]) hasSt[ld[i]/2] = 1; } if (le[i] == 3*E0 && lo[i] == 3*O0 && lg[i] == 3*g0) hasS[ld[i]/2] = 1; }
            /* also any common continuation (not necessarily tripled) shared by all three digits: */
            int anyCommon = 0;
            for (int i = 0; i < cnt && !anyCommon; i++) if (ld[i] == 0) { int f2 = 0, f4 = 0; for (int j = 0; j < cnt; j++) if (lx[j] == lx[i] && lg[j] == lg[i]) { if (ld[j] == 2) f2 = 1; if (ld[j] == 4) f4 = 1; } if (f2 && f4) anyCommon = 1; }
            if (has[0] && has[1] && has[2]) { srcWeak++; if (st) srcWeakStrong++; printf("WEAK-TRIPLE %s root=%d strong=%d (H,E,O,g)=(%lld,%lld,%lld,%lld) strongchildren=%d%d%d\n", line, root, st, H0,E0,O0,g0, hasSt[0],hasSt[1],hasSt[2]); }
            if (hasS[0] && hasS[1] && hasS[2]) { srcStrict++; printf("STRICT-TRIPLE %s root=%d strong=%d (H,E,O,g)=(%lld,%lld,%lld,%lld)\n", line, root, st, H0,E0,O0,g0); }
            if (has[1] && has[2]) pairWeak++; if (hasS[1] && hasS[2]) pairStrict++;
            if (anyCommon) printf("ANY-COMMON %s root=%d strong=%d (H,E,O,g)=(%lld,%lld,%lld,%lld)\n", line, root, st, H0,E0,O0,g0);
            free(lx); free(lg); free(le); free(lo); free(ld); free(lstrong);
        }
    }
    printf("order %d: %ld tournaments; law D+g: [strong][terminal root] -> (neg, zero, pos):", nn, tcount);
    for (int s = 0; s < 2; s++) for (int t = 0; t < 2; t++) printf(" s%dt%d=(%ld,%ld,%ld)", s, t, law[s][t][0], law[s][t][1], law[s][t][2]);
    printf("\nsources with the three digits 3H,3H+2,3H+4 under common tripled (X,g): %ld (%ld strong); under tripled (E,O,g): %ld; pairs {+2,+4} weak: %ld strict: %ld\n", srcWeak, srcWeakStrong, srcStrict, pairWeak, pairStrict);
    return 0;
}
