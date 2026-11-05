/*
 * aud_bank.c  (S78 audit, independent code)
 * Rooted-triple bank from ALL labelled tournaments of order n = 1..NMAX, rooted at vertex 0
 * (every rooted tournament is isomorphic to one rooted at vertex 0).
 * Triple: (h, e_0, R_0) with R_0 = sum_{X subset V-0} h(T[X]) * s_0(T[V-X]), h(empty) = 1.
 * Pull-style DP: end[S][v] = sum_{u in S-v, u->v} end[S-v][u]; start[S][v] = sum_{u in S-v, v->u} start[S-v][u].
 * Prints distinct triples "h e R" (sorted by h, e, R) to stdout; summary to stderr.
 * Usage: aud_bank NMAX > triples.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
typedef struct { ll h, e, R; } Tr;
static int cmp(const void *a, const void *b) {
    const Tr *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->e != y->e) return x->e < y->e ? -1 : 1;
    if (x->R != y->R) return x->R < y->R ? -1 : 1;
    return 0;
}
int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 7;
    size_t cap = 1 << 20, cnt = 0; Tr *all = malloc(sizeof(Tr) * cap);
    for (int n = 1; n <= NMAX; n++) {
        int na = n * (n - 1) / 2, full = (1 << n) - 1;
        static ll end[1 << 8][8], start[1 << 8][8], hS[1 << 8];
        size_t before = cnt;
        for (long code = 0; code < (1L << na); code++) {
            int adj[8][8]; int k = 0;           /* adj[i][j] = 1 iff i -> j */
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if ((code >> k) & 1) adj[i][j] = 1; else adj[j][i] = 1; }
            hS[0] = 1;
            for (int S = 1; S <= full; S++) {
                ll t = 0;
                for (int v = 0; v < n; v++) {
                    if (!((S >> v) & 1)) { end[S][v] = 0; start[S][v] = 0; continue; }
                    int P = S & ~(1 << v);
                    if (P == 0) { end[S][v] = 1; start[S][v] = 1; t += 1; continue; }
                    ll ee = 0, ss = 0;
                    for (int u = 0; u < n; u++) if ((P >> u) & 1) {
                        if (adj[u][v]) ee += end[P][u];     /* ... u -> v */
                        if (adj[v][u]) ss += start[P][u];   /* v -> u ... */
                    }
                    end[S][v] = ee; start[S][v] = ss; t += ee;
                }
                hS[S] = t;
            }
            ll h = hS[full], e = end[full][0], R = 0;
            for (int X = 0; X <= full; X++) if (!(X & 1)) R += hS[X] * start[full & ~X][0];
            if (cnt == cap) { cap *= 2; all = realloc(all, sizeof(Tr) * cap); }
            all[cnt].h = h; all[cnt].e = e; all[cnt].R = R; cnt++;
            /* compact periodically */
            if (cnt == cap) {
                qsort(all, cnt, sizeof(Tr), cmp);
                size_t m = 0; for (size_t i = 0; i < cnt; i++) if (m == 0 || cmp(&all[i], &all[m - 1]) != 0) all[m++] = all[i];
                cnt = m;
            }
        }
        qsort(all, cnt, sizeof(Tr), cmp);
        size_t m = 0; for (size_t i = 0; i < cnt; i++) if (m == 0 || cmp(&all[i], &all[m - 1]) != 0) all[m++] = all[i];
        cnt = m;
        fprintf(stderr, "after order %d: %zu distinct triples (labelled tournaments %ld)\n", n, cnt, 1L << na);
        (void)before;
    }
    for (size_t i = 0; i < cnt; i++) printf("%lld %lld %lld\n", all[i].h, all[i].e, all[i].R);
    return 0;
}
