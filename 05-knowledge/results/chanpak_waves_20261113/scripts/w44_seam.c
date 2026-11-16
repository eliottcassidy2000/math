/* Wave 44: the least-resistance seam.  T = (triangle x->y->z->x) placed FIRST over the two-apex window tape T'
   (every triangle vertex beats every spine vertex), except: apex A beats the triangle vertices in Sa, apex B beats
   those in Sb (the other triangle vertices beat A, B).  Exact count by subset DP.
   Compared with: 3 * H(T'_{0u}) where T'_{0u} prepends one bottom spine site b with A -> b and b -> B
   (computed with the k-apex transfer).  Hub seam theorem (README wave 44 §2): Sa = {x}, Sb = {} gives
   H = 3 H(0u) + 2 [u_0 = 0] H(T' - A).  For one apex H(T' - A) = 1; for this window tape H(T' - A) = H_B = 9 for all
   four settings (w44_facts.py), so the predicted difference is 18 [u_0 = 0]; the printed "one-apex prediction" is
   2 [u_0 = 0] and the receipts show the difference is 0 or 18 accordingly.
   usage: seam h Sa Sb      (Sa, Sb bit masks over x=1, y=2, z=4) */
#include <stdio.h>
#include <stdlib.h>
#include "w44_kapex.h"
static const int et[4] = {45, 51, 21, 45}, es[4] = {68, 98, 8, 8}, eo[4] = {0, 1, 0, 1};
static int h;
static ull dpH(int n, unsigned *out) {
    static ull dp[1 << 19][19]; int N = 1 << n;
    for (int S = 0; S < N; S++) memset(dp[S], 0, sizeof(ull) * n);
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { ull c = dp[S][v]; if (!c) continue;
        unsigned t = out[v] & ~S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[S | 1u << u][u] += c; } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[N - 1][v]; return s;
}
static void bits(unsigned head, int c, int i, int *u, int *w) {
    if (i < h) { *u = head >> i & 1; *w = 1; } else { int k = i - h; *u = et[c] >> k & 1; *w = !(es[c] >> k & 1); }
}
static ull seamH(unsigned head, int c, int Sa, int Sb) {
    int m = h + 7, A = m, B = m + 1, X = m + 2, n = m + 5; unsigned out[24] = {0};
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) out[i] |= 1u << j;
    for (int i = 0; i < m; i++) { int u, w; bits(head, c, i, &u, &w);
        if (u) out[i] |= 1u << A; else out[A] |= 1u << i;
        if (w) out[i] |= 1u << B; else out[B] |= 1u << i; }
    if (eo[c]) out[A] |= 1u << B; else out[B] |= 1u << A;
    for (int t = 0; t < 3; t++) { int v = X + t; out[v] |= 1u << (X + (t + 1) % 3);
        for (int i = 0; i < m; i++) out[v] |= 1u << i;
        if (Sa >> t & 1) out[A] |= 1u << v; else out[v] |= 1u << A;
        if (Sb >> t & 1) out[B] |= 1u << v; else out[v] |= 1u << B; }
    return dpH(n, out);
}
static ull tapeH(unsigned head, int c, int prepend) {   /* two-apex tape, optionally with a bottom site b: A -> b, b -> B */
    int m = h + 7, off = prepend ? 1 : 0; M_ = m + off; K_ = 2;
    if (prepend) beat_[0] = 0 | (1 << 1);
    for (int i = 0; i < m; i++) { int u, w; bits(head, c, i, &u, &w); beat_[i + off] = u | (w << 1); }
    arc_[0][1] = eo[c]; arc_[1][0] = !eo[c];
    return kapex_count();
}
int main(int argc, char **argv) {
    h = atoi(argv[1]); int Sa = atoi(argv[2]), Sb = atoi(argv[3]);
    kapex_init_perms(2);
    int nh = 1 << h; long win6 = 0, form = 0; long long dmin = 1LL << 60, dmax = -(1LL << 60);
    for (int hd = 0; hd < nh; hd++) {
        ull s0 = seamH(hd, 0, Sa, Sb); int ok6 = 1;
        for (int c = 0; c < 4; c++) {
            ull s = c ? seamH(hd, c, Sa, Sb) : s0;
            if (s != s0 + 6 * c) ok6 = 0;
            ull r = tapeH(hd, c, 1);
            long long d = (long long)s - 3 * (long long)r;
            if (d < dmin) dmin = d; if (d > dmax) dmax = d;
            long long pred = (hd & 1) ? 0 : 2;
            if (d == pred) form++;
            if (hd < 4 && c < 2) printf("  head %2d c %d: seamed %llu, 3*H(0u) %llu, difference %lld (one-apex prediction %lld)\n", hd, c, s, 3 * r, d, pred);
        }
        win6 += ok6;
    }
    printf("h=%d Sa=%d Sb=%d: seamed window steps exactly 6 for %ld of %d heads; seamed - 3 H(0u) in [%lld, %lld]; equals 2[u_0=0] in %ld of %d cases\n",
           h, Sa, Sb, win6, nh, dmin, dmax, form, 4 * nh);
    return 0;
}
