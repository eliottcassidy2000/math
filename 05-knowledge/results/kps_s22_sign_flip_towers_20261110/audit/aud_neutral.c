/* aud_neutral.c - auditor's independent check of THM-5130 (G) / HYP-9662 (D):
 * over ALL LABELLED tournaments X of order n (n = argv[1], <= 7) and every nonempty vertex set P,
 *     kappa_X(P) = Pi_2(X) - 2 W_X(P) - 2 W3_X(P)
 * with W_X(P)  = sum over ordered pairs e != f in P of #(D1, D2) nonempty disjoint covering X, D1 ending at e, D2 starting at f,
 *      W3_X(P) = same with a third nonempty path D3.
 * Reports, for the gadgets with kappa = 0: |P|, whether W_X(P) = 2 H(X) (i.e. r_X = -H(X)), H(X), and for H(X) = 1
 * the distance of the two vertices in the transitive order.  No class lists are used (labelled exhaustion).
 */
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
    u64 neutral = 0, bad_size = 0, bad_r = 0, trans = 0, bad_dist = 0, total = 0;
    u64 trans_d1 = 0, trans_d2 = 0, trans_pairs_other_zero = 0;
    /* histogram of (H) among neutral gadgets: small table */
    static u64 histH[4096];
    for (u64 code = 0; code < ncodes; code++) {
        unsigned out[MAXN] = {0}, in[MAXN] = {0};
        int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (code >> b & 1) { out[i] |= 1u << j; in[j] |= 1u << i; }
            else { out[j] |= 1u << i; in[i] |= 1u << j; }
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
        i64 H = tot[full], Pi = 0;
        for (int S = 0; S <= full; S++) Pi += tot[S] * tot[full ^ S];
        /* K2[e][f], K3[e][f] */
        i64 K2[MAXN][MAXN], K3[MAXN][MAXN];
        memset(K2, 0, sizeof K2); memset(K3, 0, sizeof K3);
        for (int S = 1; S <= full; S++) {
            int rest = full ^ S;
            for (int Y = rest; Y; Y = (Y - 1) & rest) {
                int Z = rest ^ Y;
                i64 tz = tot[Z];
                unsigned se = S;
                while (se) {
                    int e = __builtin_ctz(se); se &= se - 1;
                    i64 ce = endat[S][e]; if (!ce) continue;
                    unsigned sf = Y;
                    while (sf) {
                        int f = __builtin_ctz(sf); sf &= sf - 1;
                        i64 cf = startat[Y][f]; if (!cf) continue;
                        if (Z == 0) K2[e][f] += ce * cf; else K3[e][f] += ce * cf * tz;
                    }
                }
            }
        }
        for (int P = 1; P <= full; P++) {
            i64 w = 0, w3 = 0;
            for (int e = 0; e < n; e++) if (P >> e & 1) for (int f = 0; f < n; f++) if (f != e && (P >> f & 1)) { w += K2[e][f]; w3 += K3[e][f]; }
            total++;
            if (Pi - 2 * w - 2 * w3 == 0) {
                neutral++;
                int sz = __builtin_popcount(P);
                if (sz != 2) bad_size++;
                if (w != 2 * H) bad_r++;
                if (H < 4096) histH[H]++;
                if (H == 1) {
                    trans++;
                    if (sz == 2) {
                        int e = __builtin_ctz(P), f = __builtin_ctz(P & (P - 1));
                        int d = abs(__builtin_popcount(out[e]) - __builtin_popcount(out[f]));
                        if (d == 1) trans_d1++; else if (d == 2) trans_d2++; else bad_dist++;
                    }
                }
            } else if (H == 1 && __builtin_popcount(P) == 2) {
                int e = __builtin_ctz(P), f = __builtin_ctz(P & (P - 1));
                int d = abs(__builtin_popcount(out[e]) - __builtin_popcount(out[f]));
                if (d <= 2) trans_pairs_other_zero++;   /* a transitive pair at distance 1 or 2 that is NOT neutral */
            }
        }
    }
    printf("n=%d labelled tournaments=%llu gadgets (nonempty P)=%llu\n", n, ncodes, total);
    printf("  neutral (kappa=0): %llu ; with |P| != 2: %llu ; with W != 2H (r_X != -H): %llu\n", neutral, bad_size, bad_r);
    printf("  neutral with H(X)=1: %llu (distance 1: %llu, distance 2: %llu, other distance: %llu) ; transitive pairs at distance <= 2 that are not neutral: %llu\n",
           trans, trans_d1, trans_d2, bad_dist, trans_pairs_other_zero);
    printf("  H(X) values among neutral gadgets:");
    for (int h = 0; h < 4096; h++) if (histH[h]) printf(" %d:%llu", h, histH[h]);
    printf("\n");
    return 0;
}
