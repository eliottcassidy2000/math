/* seams9.c: contact seams and arc-flip seams from bases of orders 3..9 (gentourng -q streams).  Per class: h, e_x, s_x,
   the pair-cover matrix M(y, x), B_x = s_x + sum_y M(y, x); values H(T + v_X) for all X by the quadratic form, built
   incrementally by the contact digit Delta_X(x) = (B_x - e_x) - sum_{y in X} (M(y, x) + M(x, y)); a contact seam is
   |Delta| = 2 with centre val[X] + Delta/2.  Arc-flip seams: Delta = M(b, a) - M(a, b).  Reports coverage of even centres
   below LIM, the misses below 2000 and their residues mod 8, the contact-digit histogram, the arc-flip coverage, and the
   +2 window (3H + 2, 3H + 4) for the odd counts H <= 1000 realised as contact values. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
#define LIM 10000
#define VLIM (1 << 20)
static int n; static u64 adj[12], radj[12]; static i64 Hm[1<<10], dp[1<<10][10], dpr[1<<10][10], M[10][10], val[1<<10];
static unsigned char cen[LIM + 2], fcen[LIM + 2], vals[VLIM + 2]; static long dig[4100]; static long seams = 0, flips = 0, fseams = 0;
static void dps(u64 *A, i64 d[][10]) { int N = 1 << n; for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) d[m][l] = 0; for (int v = 0; v < n; v++) d[1<<v][v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = d[m][l]; if (!a) continue; u64 nx = A[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; d[m|(1<<v)][v] += a; } } }
int main(void) {
    char line[256];
    for (n = 3; n <= 9; n++) {
        char cmd[64]; snprintf(cmd, 64, "gentourng -q %d", n); FILE *g = popen(cmd, "r"); int N = 1 << n, full = N - 1;
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue;
            for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
            for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
            dps(adj, dp); dps(radj, dpr);
            for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
            i64 h = Hm[full]; i64 e[10], B[10];
            for (int y = 0; y < n; y++) for (int x = 0; x < n; x++) M[y][x] = 0;
            for (int m = 1; m < full; m++) { int c = full & ~m; for (int y = 0; y < n; y++) { i64 a = dp[m][y]; if (!a) continue; for (int x = 0; x < n; x++) { i64 b = dpr[c][x]; if (b) M[y][x] += a * b; } } }
            for (int x = 0; x < n; x++) { e[x] = dp[full][x]; B[x] = dpr[full][x]; for (int y = 0; y < n; y++) B[x] += M[y][x]; }
            val[0] = h;
            for (int X = 0; X < N; X++) {
                if (val[X] < VLIM && val[X] > 0) vals[val[X]] = 1;
                for (int x = 0; x < n; x++) {
                    if ((X >> x) & 1) continue;
                    i64 D = B[x] - e[x]; for (int y = 0; y < n; y++) if ((X >> y) & 1) D -= M[y][x] + M[x][y];
                    int X2 = X | (1 << x); if (!(X2 & (X2 - 1) & ~X) || 1) val[X2] = val[X] + D;   /* consistent by symmetry */
                    i64 ad = D < 0 ? -D : D; if (ad / 2 < 4100) dig[ad / 2]++;
                    if (ad == 2) { seams++; i64 c = val[X] + D / 2; if (c > 0 && c <= LIM) cen[c] = 1; }
                }
            }
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if ((adj[a] >> b) & 1) { flips++; i64 D = M[b][a] - M[a][b]; if (D == 2 || D == -2) { fseams++; i64 c = h + D / 2; if (c > 0 && c <= LIM) fcen[c] = 1; } }
        }
        pclose(g); fprintf(stderr, "order %d done\n", n);
    }
    long cov = 0, miss = 0; int res[8] = {0}; printf("contact seams from bases of orders 3-9: %ld seam pairs; even centres in [4, %d] missed:", seams, LIM);
    for (int c = 4; c <= LIM; c += 2) { if (cen[c]) cov++; else { miss++; res[c % 8]++; if (c <= 2000) printf(" %d", c); } }
    printf("\n  covered %ld of %d; missed below %d: %ld (printed those <= 2000); residues mod 8 of the misses: 0:%d 2:%d 4:%d 6:%d\n", cov, LIM/2 - 1, LIM, miss, res[0], res[2], res[4], res[6]);
    long tot = 0; for (int d = 0; d < 4100; d++) tot += dig[d]; int gap = 0; while (gap < 4100 && dig[gap]) gap++;
    printf("  contact digits |Delta|/2: every value in [0, %d] occurs; share digit 0: %.4f, digit 1: %.4f, digit 2: %.4f\n", gap - 1, (double)dig[0]/tot, (double)dig[1]/tot, (double)dig[2]/tot);
    long fcov = 0; printf("arc-flip seams from orders 3-9: %ld flips, %ld seams; even centres in [4, 2900] missed:", flips, fseams);
    for (int c = 4; c <= 2900; c += 2) { if (fcen[c]) fcov++; else printf(" %d", c); } printf("\n  covered %ld of 1449\n", fcov);
    long nh = 0, nw = 0; printf("the +2 window: odd counts H <= 1000 realised as contact values whose centre 3(H+1) is NOT a contact seam centre:");
    for (int H = 1; H <= 1000; H += 2) if (vals[H]) { nh++; if (3 * H + 3 <= LIM && cen[3 * H + 3]) nw++; else printf(" %d", H); }
    printf("\n  windows realised: %ld of %ld\n", nw, nh);
    return 0;
}
