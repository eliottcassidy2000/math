/* profile.c: rooted profiles (H, E, O, g) of tournaments on n <= 12 vertices (adjacency as out-masks, vertex 0 = low bit),
   and the fresh-digit continuation search.
   Usage:
     profile check                      -- verify the literal ternary-deposit source and children (fresh_unit_prime_gluing_20261116)
     profile search1 <n> <root> m0..m_{n-1}   -- all one-vertex extensions x (x beats mask c of old vertices, 2^n choices):
                                            print those with H' = 3H + d, d in {0,2,4}, grouped by continuation
     profile search2 <n> <root> m0..m_{n-1}   -- all two-vertex extensions (x -> y, masks cx, cy): same report
   Profile: H = Hamiltonian paths; g = H(T - root); E (O) = sum over root-start paths P with an even (odd) number of
   vertices of H(T - V(P)), H(empty) = 1 (the ordered two-path covers of the bridge protocol). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef long long i64;
static int n; static u64 adj[16];
static i64 Hm[1<<13];           /* H(T[mask]) */
static i64 dp[1<<13][13];       /* paths covering mask ending at last */
static i64 rs[1<<13];           /* root-start Hamiltonian paths of T[mask] (mask contains root) */
static i64 dpr[1<<13][13];
static void profile(int root, i64 *H, i64 *E, i64 *O, i64 *g) {
    int N = 1 << n;
    for (int m = 0; m < N; m++) { Hm[m] = 0; for (int l = 0; l < n; l++) dp[m][l] = 0, dpr[m][l] = 0; }
    for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
    dpr[1<<root][root] = 1;
    for (int m = 1; m < N; m++) {
        for (int l = 0; l < n; l++) {
            i64 a = dp[m][l], b = dpr[m][l];
            if (!a && !b) continue;
            u64 nx = adj[l] & ~(u64)m;
            while (nx) { int v = __builtin_ctzll(nx); nx &= nx - 1; dp[m | (1<<v)][v] += a; dpr[m | (1<<v)][v] += b; }
        }
    }
    for (int m = 1; m < N; m++) { i64 s = 0, r = 0; for (int l = 0; l < n; l++) s += dp[m][l], r += dpr[m][l]; Hm[m] = s; rs[m] = r; }
    Hm[0] = 1;
    *H = Hm[N-1]; *g = Hm[(N-1) & ~(1<<root)];
    i64 e = 0, o = 0;
    for (int m = 1; m < N; m++) if ((m >> root) & 1) {
        i64 w = rs[m] * Hm[(N-1) & ~m];
        if (__builtin_popcount(m) % 2 == 0) e += w; else o += w;
    }
    *E = e; *O = o;
}
static int strong(void) { /* strongly connected? */
    int N = n; u64 reach = 1, prev = 0;
    while (reach != prev) { prev = reach; for (int v = 0; v < N; v++) if ((reach >> v) & 1) reach |= adj[v]; }
    if (reach != ((1ULL<<N)-1)) return 0;
    u64 radj[16]; for (int v = 0; v < N; v++) { radj[v] = 0; for (int u = 0; u < N; u++) if ((adj[u] >> v) & 1) radj[v] |= 1ULL<<u; }
    reach = 1; prev = 0;
    while (reach != prev) { prev = reach; for (int v = 0; v < N; v++) if ((reach >> v) & 1) reach |= radj[v]; }
    return reach == ((1ULL<<N)-1);
}
static void setbase(int nn, int *m) { n = nn; for (int i = 0; i < n; i++) adj[i] = (u64)m[i]; }
int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "check")) {
        int src[8] = {118,124,248,113,224,192,128,11};
        int c0[9] = {374,380,504,369,224,448,128,11,208};
        int c2[9] = {118,372,504,115,224,192,128,11,249};
        int c4[9] = {118,124,504,369,224,192,128,11,243};
        i64 H,E,O,g;
        setbase(8, src); profile(7,&H,&E,&O,&g); printf("source  (H,E,O,g) = (%lld,%lld,%lld,%lld) expected (163,84,89,5) strong=%d\n",H,E,O,g,strong());
        setbase(9, c0); profile(7,&H,&E,&O,&g); printf("child 0 (H,E,O,g) = (%lld,%lld,%lld,%lld) expected (489,252,267,15) strong=%d\n",H,E,O,g,strong());
        setbase(9, c2); profile(7,&H,&E,&O,&g); printf("child 2 (H,E,O,g) = (%lld,%lld,%lld,%lld) expected (491,252,267,15) strong=%d\n",H,E,O,g,strong());
        setbase(9, c4); profile(7,&H,&E,&O,&g); printf("child 4 (H,E,O,g) = (%lld,%lld,%lld,%lld) expected (493,252,267,15) strong=%d\n",H,E,O,g,strong());
        return 0;
    }
    int nn = atoi(argv[2]), root = atoi(argv[3]); int base[16];
    for (int i = 0; i < nn; i++) base[i] = atoi(argv[4+i]);
    i64 H0,E0,O0,g0; setbase(nn, base); profile(root,&H0,&E0,&O0,&g0);
    printf("source n=%d root=%d (H,E,O,g)=(%lld,%lld,%lld,%lld) X=%lld D=%lld strong=%d\n", nn, root, H0,E0,O0,g0, E0+O0, E0-O0, strong());
    int two = !strcmp(argv[1], "search2");
    int N0 = 1 << nn;
    /* report lines: d, strong flag, masks, (H',E',O',g') */
    long found[3] = {0,0,0}, foundstrong[3] = {0,0,0};
    for (int cx = 0; cx < N0; cx++) {
        int cyMax = two ? N0 : 1;
        for (int cy = 0; cy < cyMax; cy++) for (int o = 0; o < (two ? 2 : 1); o++) {
            /* new vertex x = nn, y = nn+1; x beats old vertices in cx, y beats old vertices in cy; o: x->y (0) or y->x (1) */
            n = nn + (two ? 2 : 1);
            for (int i = 0; i < nn; i++) {
                adj[i] = (u64)base[i];
                if (!((cx >> i) & 1)) adj[i] |= 1ULL << nn;
                if (two && !((cy >> i) & 1)) adj[i] |= 1ULL << (nn+1);
            }
            adj[nn] = (u64)cx; if (two) { adj[nn+1] = (u64)cy; if (o == 0) adj[nn] |= 1ULL << (nn+1); else adj[nn+1] |= 1ULL << nn; }
            i64 H,E,O,g; profile(root,&H,&E,&O,&g);
            i64 d = H - 3*H0;
            if (d == 0 || d == 2 || d == 4) {
                int st = strong(); found[d/2]++; if (st) foundstrong[d/2]++;
                printf("d=%lld strong=%d masks=%d%s%s%s (H,E,O,g)=(%lld,%lld,%lld,%lld) X=%lld D=%lld tripled(E,O,g)=%d tripled(X,g)=%d\n",
                       d, st, cx, two ? "," : "", two ? (char[]){0} : "", "", H,E,O,g, E+O, E-O,
                       (E==3*E0 && O==3*O0 && g==3*g0), (E+O==3*(E0+O0) && g==3*g0));
                if (two) printf("   cy=%d orient=%d\n", cy, o);
            }
        }
    }
    printf("counts: d=0: %ld (%ld strong), d=2: %ld (%ld strong), d=4: %ld (%ld strong)\n", found[0],foundstrong[0],found[1],foundstrong[1],found[2],foundstrong[2]);
    return 0;
}
