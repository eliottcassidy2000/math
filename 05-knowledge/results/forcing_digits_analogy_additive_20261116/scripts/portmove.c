/* portmove.c <n>: P107 complement.  For every tournament class of order n (gentourng on stdin) and every root x that is
   READY (defect b = s - r - g = 0, i.e. O - E - g = 0, wave 45 §6), with charge a = 2r + 1 - H (r = E): enumerate all
   extensions by one or two new vertices with arbitrary contacts (every orientation among the new vertices) and every new
   root y (old or new); test the binary port move H' = 2H + 1 with a' = a and b' = 0 (P107), and the weaker variants.
   Prints each hit and the counts. */
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
static void seedmode(int argc, char **argv) {
    int nn = atoi(argv[2]), x = atoi(argv[3]); int base[16]; for (int i = 0; i < nn; i++) base[i] = atoi(argv[4+i]);
    n = nn; for (int i = 0; i < nn; i++) adj[i] = base[i];
    i64 H0,E0,O0,g0; profile(x,&H0,&E0,&O0,&g0); i64 a0 = 2*E0+1-H0, b0 = O0-E0-g0;
    printf("seed n=%d root=%d (H,E,O,g)=(%lld,%lld,%lld,%lld) charge a=%lld defect b=%lld; target H'=2H+1+a=%lld\n", nn, x, H0,E0,O0,g0, a0, b0, 2*H0+1+a0);
    long hits = 0, cnt2h1 = 0; int N0 = 1 << nn; int maxtwo = nn <= 9 ? 2 : 1;
    for (int two = 0; two < maxtwo; two++) for (int cx = 0; cx < N0; cx++) for (int cy = 0; cy < (two ? N0 : 1); cy++) for (int o = 0; o < (two ? 2 : 1); o++) {
        n = nn + 1 + two;
        for (int i = 0; i < nn; i++) { adj[i] = base[i]; if (!((cx>>i)&1)) adj[i] |= 1ULL<<nn; if (two && !((cy>>i)&1)) adj[i] |= 1ULL<<(nn+1); }
        adj[nn] = cx; if (two) { adj[nn+1] = cy; if (o == 0) adj[nn] |= 1ULL<<(nn+1); else adj[nn+1] |= 1ULL<<nn; }
        i64 H,E,O,g; profile(x,&H,&E,&O,&g); if (H != 2*H0 + 1 + a0) continue; cnt2h1++;
        for (int y = 0; y < n; y++) { profile(y,&H,&E,&O,&g); if (2*E+1-H == a0 && O-E-g == 0) { hits++; if (hits <= 6) printf("  HIT masks cx=%d cy=%d o=%d newroot=%d (H',E',O',g')=(%lld,%lld,%lld,%lld)\n", cx, cy, o, y, H,E,O,g); } }
    }
    printf("seed done: extensions with H'=2H+1+a: %ld; root-moving ready same-charge hits: %ld\n", cnt2h1, hits);
}
int main(int argc, char **argv) {
    if (!strcmp(argv[1], "seed")) { seedmode(argc, argv); return 0; }
    int nn = atoi(argv[1]); int maxk = (argc > 2 && !strcmp(argv[2], "one")) ? 1 : 2; char line[256]; long classes = 0, readyroots = 0, hitsP107 = 0, hits2H1ready = 0, hits2H1charge = 0, hits2H1 = 0;
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
        int base[16] = {0}; int k = 0;
        for (int i = 0; i < nn; i++) for (int j = i+1; j < nn; j++) { if (line[k++] == '1') base[i] |= 1<<j; else base[j] |= 1<<i; }
        for (int x = 0; x < nn; x++) {
            n = nn; for (int i = 0; i < nn; i++) adj[i] = base[i];
            i64 H0,E0,O0,g0; profile(x,&H0,&E0,&O0,&g0);
            i64 b0 = O0 - E0 - g0, a0 = 2*E0 + 1 - H0;
            if (b0 != 0) continue; readyroots++;
            int N0 = 1 << nn;
            for (int two = 0; two < maxk; two++) for (int cx = 0; cx < N0; cx++) for (int cy = 0; cy < (two ? N0 : 1); cy++) for (int o = 0; o < (two ? 2 : 1); o++) {
                n = nn + 1 + two;
                for (int i = 0; i < nn; i++) { adj[i] = base[i]; if (!((cx>>i)&1)) adj[i] |= 1ULL<<nn; if (two && !((cy>>i)&1)) adj[i] |= 1ULL<<(nn+1); }
                adj[nn] = cx; if (two) { adj[nn+1] = cy; if (o == 0) adj[nn] |= 1ULL<<(nn+1); else adj[nn+1] |= 1ULL<<nn; }
                /* count H' once (root-independent) */
                i64 H,E,O,g; profile(x,&H,&E,&O,&g);
                if (H != 2*H0 + 1 + a0) continue;   /* r -> 2r + 1 at the same charge: H' = 2H + 1 + a */
                hits2H1++;
                for (int y = 0; y < n; y++) {
                    profile(y,&H,&E,&O,&g); i64 a = 2*E + 1 - H, b = O - E - g;
                    if (b == 0) hits2H1ready++;
                    if (a == a0) hits2H1charge++;
                    if (b == 0 && a == a0) { hitsP107++; printf("P107-HIT charge=%lld base=%s root=%d (H,E,O,g)=(%lld,%lld,%lld,%lld) a=%lld -> masks cx=%d cy=%d o=%d newroot=%d (H',E',O',g')=(%lld,%lld,%lld,%lld)\n", a0, line, x, H0,E0,O0,g0, a0, cx, cy, o, y, H,E,O,g); }
                }
            }
        }
    }
    printf("order %d: %ld classes, %ld ready roots; extensions (<= 2 vertices) with H' = 2H + 1 + a: %ld (root-independent count); of the (extension, new root) pairs: ready %ld, same charge %ld, both (P107) %ld\n", nn, classes, readyroots, hits2H1, hits2H1ready, hits2H1charge, hitsP107);
    return 0;
}
