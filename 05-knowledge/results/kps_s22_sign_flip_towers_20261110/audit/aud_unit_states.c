/* aud_unit_states.c - all pairs (B, O), B a LABELLED tournament of order n = argv[1] (<= 7), O any vertex set
 * (empty included), with offset r = +1 or -1: prints the state as "r g M E" with g = H(B) + r (the format of
 * runs/e2_states_order_le7.txt).  Output is to be piped through sort -u.  Auditor's own subset recursion. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef long long i64;
#define MAXN 7
static i64 endat[1 << MAXN][MAXN], startat[1 << MAXN][MAXN], tot[1 << MAXN];
int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int nb = n * (n - 1) / 2, full = (1 << n) - 1;
    u64 ncodes = 1ull << nb;
    for (u64 code = 0; code < ncodes; code++) {
        unsigned out[MAXN] = {0}, in[MAXN] = {0};
        int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (code >> b & 1) { out[i] |= 1u << j; in[j] |= 1u << i; } else { out[j] |= 1u << i; in[i] |= 1u << j; }
            b++;
        }
        tot[0] = 1;
        for (int S = 1; S <= full; S++) {
            i64 sum = 0;
            for (int v = 0; v < n; v++) { endat[S][v] = 0; startat[S][v] = 0; }
            if ((S & (S - 1)) == 0) { int v = __builtin_ctz(S); endat[S][v] = 1; startat[S][v] = 1; sum = 1; }
            else {
                unsigned t = S;
                while (t) {
                    int v = __builtin_ctz(t); t &= t - 1;
                    int R = S ^ (1 << v);
                    i64 c = 0; unsigned xs = R & in[v];
                    while (xs) { int x = __builtin_ctz(xs); xs &= xs - 1; c += endat[R][x]; }
                    endat[S][v] = c; sum += c;
                    c = 0; unsigned ys = R & out[v];
                    while (ys) { int y = __builtin_ctz(ys); ys &= ys - 1; c += startat[R][y]; }
                    startat[S][v] = c;
                }
            }
            tot[S] = sum;
        }
        i64 H = tot[full];
        i64 Q[MAXN][MAXN], U[MAXN];
        memset(Q, 0, sizeof Q); memset(U, 0, sizeof U);
        for (int S = 0; S < full; S++) {
            int R = full ^ S;
            for (int o = 0; o < n; o++) if (R >> o & 1) {
                i64 s = startat[R][o]; if (!s) continue;
                U[o] += tot[S] * s;
                if (S) for (int e = 0; e < n; e++) if ((S >> e & 1) && endat[S][e]) Q[e][o] += endat[S][e] * s;
            }
        }
        for (int O = 0; O <= full; O++) {
            i64 W = 0, M = 0, E = 0;
            for (int e = 0; e < n; e++) if (O >> e & 1) {
                M += U[e]; E += endat[full][e];
                for (int f = 0; f < n; f++) if (f != e && (O >> f & 1)) W += Q[e][f];
            }
            i64 r = H - W;
            if (r == 1 || r == -1) printf("%lld %lld %lld %lld\n", r, H + r, M, E);
        }
    }
    return 0;
}
