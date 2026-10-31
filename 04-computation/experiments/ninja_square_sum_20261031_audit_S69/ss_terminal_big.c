/* S69 audit (independent referee, 2026-10-31): randomized search for terminal-form square-sum cycles at large N.
   Written from scratch; does not reuse terminal_cycles_large_20261031.c.

   Terminal form: a Hamiltonian path of Q_N that starts 1,3 and ends at T = 8 (N even) or T = 15 (N odd); closing it to 1
   (8+1 = 9, 15+1 = 16) gives a Hamiltonian cycle through 8-1-3 or 15-1-3.

   Search: DFS from 1,3 with
     - avail[v] = number of neighbours of v in S = U + {cur} (U = unvisited), kept incrementally: a move cur -> w only
       removes cur from S, so only the neighbours of the old cur change;
     - dead state if some unvisited v != T has avail[v] <= 1, or avail[T] = 0;
     - forced move if an unvisited non-target neighbour w of cur has avail[w] = 2 (two such: dead);
     - candidate order: increasing avail (Warnsdorff), ties broken by a per-attempt pseudo-random key;
     - node budget per attempt, restarts with fresh seeds and doubling budget.
   Every found path is printed in full; correctness is NOT trusted here: verify_terminal_big.py re-checks every path.

   usage: ss_terminal_big LO HI [BUDGET0] [MAXATTEMPTS]   (prints one line per N:  N=.. T=.. attempts=.. nodes=.. path=..)
   Build: gcc -O2 -Wl,--stack,268435456 -o ss_terminal_big.exe ss_terminal_big.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXV 8192
static int N, T;
static int *adjl[MAXV], deg[MAXV];
static int avail[MAXV];
static unsigned char inU[MAXV];
static int path[MAXV];
static unsigned int key[MAXV];
static long long nodes, budget;
static unsigned long long rs;

static unsigned int rnd(void) {   /* xorshift64* */
    rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27;
    return (unsigned int)((rs * 2685821657736338717ULL) >> 32);
}

static void build(int n) {
    N = n;
    for (int v = 1; v <= n; v++) {
        int tmp[256], c = 0;
        for (long long r = 1; r * r <= 2LL * n; r++) {
            long long w = r * r - v;
            if (w >= 1 && w <= n && w != v) tmp[c++] = (int)w;
        }
        free(adjl[v]);
        adjl[v] = malloc(sizeof(int) * (c ? c : 1));
        memcpy(adjl[v], tmp, sizeof(int) * c);
        deg[v] = c;
    }
}

static int dfs(int len) {
    if (++nodes > budget) return -1;
    int cur = path[len - 1];
    if (len == N) return cur == T;
    if (cur == T) return 0;
    int cand[256], nc = 0, forced = 0, nforced = 0;
    for (int k = 0; k < deg[cur]; k++) {
        int w = adjl[cur][k];
        if (!inU[w]) continue;
        if (w != T && avail[w] == 2) { forced = w; nforced++; }
    }
    if (nforced > 1) return 0;
    if (forced) cand[nc++] = forced;
    else {
        for (int k = 0; k < deg[cur]; k++) {
            int w = adjl[cur][k];
            if (!inU[w]) continue;
            if (w == T && len != N - 1) continue;
            cand[nc++] = w;
        }
        /* insertion sort by (avail, key) */
        for (int i = 1; i < nc; i++) {
            int x = cand[i], j = i;
            while (j > 0 && (avail[cand[j - 1]] > avail[x] || (avail[cand[j - 1]] == avail[x] && key[cand[j - 1]] > key[x]))) {
                cand[j] = cand[j - 1]; j--;
            }
            cand[j] = x;
        }
    }
    for (int i = 0; i < nc; i++) {
        int w = cand[i];
        inU[w] = 0;
        int bad = 0;
        for (int k = 0; k < deg[cur]; k++) {
            int v = adjl[cur][k];
            avail[v]--;
            if (inU[v]) {
                if (v != T && avail[v] <= 1) bad = 1;
                if (v == T && avail[v] == 0) bad = 1;
            }
        }
        int r = 0;
        if (!bad) { path[len] = w; r = dfs(len + 1); }
        for (int k = 0; k < deg[cur]; k++) avail[adjl[cur][k]]++;
        inU[w] = 1;
        if (r) return r;
    }
    return 0;
}

static int attempt(int n, unsigned long long seed, long long b) {
    rs = seed * 0x9E3779B97F4A7C15ULL + 12345;
    if (!rs) rs = 1;
    for (int v = 1; v <= n; v++) key[v] = rnd();
    memset(inU, 0, sizeof(unsigned char) * (n + 2));
    for (int v = 1; v <= n; v++) inU[v] = 1;
    inU[1] = inU[3] = 0;
    path[0] = 1; path[1] = 3;
    for (int v = 1; v <= n; v++) {
        int c = 0;
        for (int k = 0; k < deg[v]; k++) { int w = adjl[v][k]; if (inU[w] || w == 3) c++; }
        avail[v] = c;
    }
    for (int v = 1; v <= n; v++) if (inU[v]) {
        if (v != T && avail[v] <= 1) return 0;
        if (v == T && avail[v] == 0) return 0;
    }
    nodes = 0; budget = b;
    return dfs(2);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s LO HI [BUDGET0] [MAXATTEMPTS]\n", argv[0]); return 1; }
    int lo = atoi(argv[1]), hi = atoi(argv[2]);
    long long b0 = argc > 3 ? atoll(argv[3]) : 200000;
    int maxa = argc > 4 ? atoi(argv[4]) : 40;
    if (hi >= MAXV) { fprintf(stderr, "HI too large\n"); return 1; }
    for (int n = lo; n <= hi; n++) {
        build(n);
        T = (n % 2 == 0) ? 8 : 15;
        int res = 0, a; long long b = b0, total = 0;
        for (a = 1; a <= maxa; a++) {
            res = attempt(n, (unsigned long long)n * 1000003ULL + (unsigned long long)a, b);
            total += nodes < b ? nodes : b;
            if (res == 1 || res == 0) break;      /* 0 = complete search without success (no restart helps) */
            if (b < 1000000000LL) b *= 2;
        }
        if (res == 1) {
            printf("N=%d T=%d attempts=%d nodes=%lld path=", n, T, a, nodes);
            for (int i = 0; i < n; i++) printf("%d%s", path[i], i + 1 < n ? "," : "\n");
        } else if (res == 0) printf("N=%d T=%d EXHAUSTED (complete search, attempt %d)\n", n, T, a);
        else printf("N=%d T=%d UNDECIDED after %d attempts (%lld nodes)\n", n, T, maxa, total);
        fflush(stdout);
    }
    return 0;
}
