/* h_spectrum_n9_cover_boxeph_S153.c  (HYP-8225)
 *
 * Complete n=9 h-spectrum by double augmentation cover:
 * 456 canonical n=7 reps -> 58,368 n=8 candidates (cover all 6880 classes)
 * -> x 256 patterns = 14,942,208 n=9 tournaments covering ALL 191,536
 * n=9 classes.  For each: h (Hamiltonian path count, subset DP) and
 * strongness.  Outputs: attained-value bitmap, strong bitmap, min strong h,
 * max h.  Tests: 7/21 still absent; f(9) = 75 (Moon-Busch); near-top holes.
 *
 * build: gcc -O2 -o /tmp/n9cover h_spectrum_n9_cover_boxeph_S153.c
 * boxeph-2026-07-20-S153.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define N7 7
#define N8 8
#define N9 9
#define HMAX 65536

static uint8_t attained[HMAX], strong_att[HMAX];

static inline uint32_t hpaths9(const int adj[N9]) {
    static uint32_t dp[1 << N9][N9];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N9; v++) dp[1 << v][v] = 1;
    for (int S = 1; S < (1 << N9); S++) {
        for (int v = 0; v < N9; v++) {
            uint32_t c = dp[S][v];
            if (!c || !((S >> v) & 1)) continue;
            int m = adj[v] & ~S;
            while (m) {
                int w = __builtin_ctz(m);
                dp[S | (1 << w)][w] += c;
                m &= m - 1;
            }
        }
    }
    uint32_t tot = 0;
    for (int v = 0; v < N9; v++) tot += dp[(1 << N9) - 1][v];
    return tot;
}

static inline int is_strong9(const int adj[N9]) {
    int full = (1 << N9) - 1;
    for (int s = 0; s < N9; s++) {
        int seen = 1 << s, stack[N9], sp = 0;
        stack[sp++] = s;
        while (sp) {
            int u = stack[--sp];
            int m = adj[u] & ~seen;
            while (m) {
                int w = __builtin_ctz(m);
                seen |= 1 << w;
                stack[sp++] = w;
                m &= m - 1;
            }
        }
        if (seen != full) return 0;
    }
    return 1;
}

int main(void) {
    FILE *f = fopen("05-knowledge/results/n7_class_reps_boxeph_S152.txt", "r");
    if (!f) { fprintf(stderr, "reps file missing\n"); return 1; }
    long reps[456]; int nreps = 0;
    while (nreps < 456 && fscanf(f, "%ld", &reps[nreps]) == 1) nreps++;
    fclose(f);
    printf("n=7 reps: %d\n", nreps);

    int pi = 0, pairs[21][2];
    for (int i = 0; i < N7; i++)
        for (int j = i + 1; j < N7; j++) { pairs[pi][0] = i; pairs[pi][1] = j; pi++; }

    long done = 0;
    for (int r = 0; r < nreps; r++) {
        int adj7[N7]; memset(adj7, 0, sizeof(adj7));
        for (int k = 0; k < 21; k++) {
            int i = pairs[k][0], j = pairs[k][1];
            if ((reps[r] >> k) & 1) adj7[i] |= 1 << j; else adj7[j] |= 1 << i;
        }
        for (int p8 = 0; p8 < (1 << N7); p8++) {
            int adj8[N8];
            for (int v = 0; v < N7; v++)
                adj8[v] = adj7[v] | ((!((p8 >> v) & 1)) ? (1 << N7) : 0);
            adj8[N7] = p8;
            for (int p9 = 0; p9 < (1 << N8); p9++) {
                int adj9[N9];
                for (int v = 0; v < N8; v++)
                    adj9[v] = adj8[v] | ((!((p9 >> v) & 1)) ? (1 << N8) : 0);
                adj9[N8] = p9;
                uint32_t h = hpaths9(adj9);
                if (h < HMAX) {
                    attained[h] = 1;
                    if (!strong_att[h] && is_strong9(adj9)) strong_att[h] = 1;
                }
                done++;
            }
        }
        if ((r + 1) % 50 == 0) {
            fprintf(stderr, "  rep %d/456 (%ld tournaments)\n", r + 1, done);
            fflush(stderr);
        }
    }
    printf("tournaments processed: %ld\n", done);
    int mx = 0; for (int v = 0; v < HMAX; v++) if (attained[v]) mx = v;
    int msf = -1; for (int v = 0; v < HMAX; v++) if (strong_att[v]) { msf = v; break; }
    int mss = 0; for (int v = 0; v < HMAX; v++) if (strong_att[v]) mss = v;
    printf("max h at n=9: %d\n", mx);
    printf("MIN strong h(9): %d   (Moon-Busch predicts 75)\n", msf);
    printf("max strong h: %d\n", mss);
    printf("h=7:%d h=21:%d h=63:%d h=105:%d h=189:%d\n",
           attained[7], attained[21], attained[63], attained[105], attained[189]);
    printf("missing odd values up to max:\n");
    int cnt = 0;
    for (int v = 1; v <= mx; v += 2)
        if (!attained[v]) { printf("%d ", v); if (++cnt % 20 == 0) printf("\n"); }
    printf("\n(total missing odds: %d)\n", cnt);
    printf("strong holes between floor and max-strong:\n");
    cnt = 0;
    for (int v = msf; v <= mss; v += 2)
        if (!strong_att[v]) { printf("%d ", v); if (++cnt % 20 == 0) printf("\n"); }
    printf("\n(total strong holes: %d)\nDONE\n", cnt);
    return 0;
}
