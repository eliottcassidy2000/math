/* lbfs7.c -- BFS over all labelled 7-tournaments (2^21 states) under native moves, from a given start.
 * State bit for pair (i<j): 1 if i->j.  Reports distance to the nearest transitive tournament and the
 * distance distribution.  Usage: ./lbfs7 paley | ./lbfs7 tt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define N 7
#define NP 21
static int pidx[N][N];
static uint8_t *dist; static uint32_t *queue;
static int path[N]; static uint32_t curadj[N]; static uint32_t curmask; static long qtail;
static void push(uint32_t m, int d) { if (dist[m] == 255) { dist[m] = d; queue[qtail++] = m; } }
static void visit(int depth, uint32_t used, int d) {
    if (depth == N) {
        uint32_t m = (~curmask) & ((1u << NP) - 1);
        for (int i = 0; i + 1 < N; i++) { int a = path[i], b = path[i+1]; int bit = pidx[a][b]; if (curmask >> bit & 1) m |= 1u << bit; else m &= ~(1u << bit); }
        push(m, d + 1); return;
    }
    int last = path[depth-1];
    for (int v = 0; v < N; v++) if (!(used >> v & 1) && (curadj[last] >> v & 1)) { path[depth] = v; visit(depth + 1, used | (1u << v), d); }
}
int main(int argc, char **argv) {
    int p = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { pidx[i][j] = p; pidx[j][i] = p; p++; }
    dist = malloc(1u << NP); memset(dist, 255, 1u << NP); queue = malloc(sizeof(uint32_t) << NP);
    uint32_t start = 0;
    if (!strcmp(argv[1], "paley")) { for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { int r = (j - i) % 7; if (r == 1 || r == 2 || r == 4) start |= 1u << pidx[i][j]; } }
    else { for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) start |= 1u << pidx[i][j]; }
    qtail = 0; push(start, 0); long qhead = 0;
    while (qhead < qtail) {
        uint32_t m = queue[qhead++]; int d = dist[m]; curmask = m;
        for (int i = 0; i < N; i++) curadj[i] = 0;
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (m >> pidx[i][j] & 1) curadj[i] |= 1u << j; else curadj[j] |= 1u << i; }
        for (int v = 0; v < N; v++) { path[0] = v; visit(1, 1u << v, d); }
    }
    long hist[64] = {0}; int best_tt = 99, best_paley = 99;
    for (uint32_t m = 0; m < (1u << NP); m++) {
        if (dist[m] == 255) { hist[63]++; continue; } hist[dist[m]]++;
        int s[N] = {0}; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (m >> pidx[i][j] & 1) s[i]++; else s[j]++; }
        int seen = 0, tt = 1; for (int i = 0; i < N; i++) { if (seen >> s[i] & 1) { tt = 0; break; } seen |= 1 << s[i]; }
        if (tt && dist[m] < best_tt) best_tt = dist[m];
        int reg = 1; for (int i = 0; i < N; i++) if (s[i] != 3) reg = 0;
        if (reg) { /* Paley iff every vertex's out-neighbourhood induces a 3-cycle (doubly regular) */
            int dr = 1; for (int i = 0; i < N && dr; i++) { uint32_t out = 0; for (int j = 0; j < N; j++) if (j != i && (m >> pidx[i][j] & 1) == (i < j)) out |= 1u << j; /* placeholder */ }
        }
    }
    printf("start=%s reachable distance histogram:", argv[1]); for (int d = 0; d < 40; d++) if (hist[d]) printf(" %d:%ld", d, hist[d]); printf(" unreachable:%ld\n", hist[63]);
    printf("distance to nearest transitive tournament = %d\n", best_tt);
    uint32_t ttstd = 0, ttop = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) ttstd |= 1u << pidx[i][j];
    printf("labelled distance to TT (i<j order) = %d, to its converse = %d\n", dist[ttstd], dist[ttop]);
    long cnt6 = 0, cnt7 = 0; for (uint32_t m = 0; m < (1u << NP); m++) { int s[N] = {0}; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (m >> pidx[i][j] & 1) s[i]++; else s[j]++; } int seen = 0, tt = 1; for (int i = 0; i < N; i++) { if (seen >> s[i] & 1) { tt = 0; break; } seen |= 1 << s[i]; } if (tt) { if (dist[m] == 6) cnt6++; else if (dist[m] == 7) cnt7++; } }
    printf("transitive labellings at distance 6: %ld, at distance 7: %ld (of 5040)\n", cnt6, cnt7);
    return 0;
}
