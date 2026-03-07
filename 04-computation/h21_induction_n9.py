"""
h21_induction_n9.py — Test induction hypothesis for H=21 at n=9

Key hypothesis: H=21 impossible at n=9 because:
1. If source/sink exists: reduces to n=8 (impossible by prior)
2. If no source/sink: alpha_1 + 2*i_2 != 10

Uses C for speed. Random sampling.

Author: opus-2026-03-07-S41
"""

import os, sys, tempfile, subprocess

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define N 9
#define NBITS (N*(N-1)/2)

static inline int arc(unsigned long long T, int i, int j) {
    if (i < j) {
        int pos = i*(2*N-i-1)/2 + (j-i-1);
        return (T >> pos) & 1;
    } else {
        int pos = j*(2*N-j-1)/2 + (i-j-1);
        return 1 - ((T >> pos) & 1);
    }
}

int find_3cycles(unsigned long long T, int *masks) {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = arc(T,a,b), bc = arc(T,b,c), ca = arc(T,c,a);
                if ((ab && bc && ca) || ((1-ab) && (1-bc) && (1-ca))) {
                    masks[count++] = (1<<a)|(1<<b)|(1<<c);
                }
            }
    return count;
}

/* For n=9, 5-cycles are expensive. Use fast canonical counting. */
int find_5cycles(unsigned long long T, int *masks) {
    int count = 0;
    int vv[5], p[5];

    for (int a = 0; a < N; a++)
    for (int b = a+1; b < N; b++)
    for (int c = b+1; c < N; c++)
    for (int d = c+1; d < N; d++)
    for (int e = d+1; e < N; e++) {
        vv[0]=a; vv[1]=b; vv[2]=c; vv[3]=d; vv[4]=e;
        int mask = (1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e);

        /* Count directed 5-cycles with canonical form starting at min vertex (=a) */
        int nseen = 0;
        p[0] = a;
        for (int i1 = 1; i1 < 5; i1++) {
            p[1] = vv[i1];
            for (int i2 = 1; i2 < 5; i2++) {
                if (i2==i1) continue;
                p[2] = vv[i2];
                for (int i3 = 1; i3 < 5; i3++) {
                    if (i3==i1||i3==i2) continue;
                    p[3] = vv[i3];
                    int i4;
                    for (i4 = 1; i4 < 5; i4++) {
                        if (i4!=i1 && i4!=i2 && i4!=i3) break;
                    }
                    p[4] = vv[i4];

                    int ok = 1;
                    for (int q = 0; q < 5; q++) {
                        if (!arc(T, p[q], p[(q+1)%5])) { ok=0; break; }
                    }
                    if (ok) nseen++;
                }
            }
        }
        for (int k = 0; k < nseen; k++)
            masks[count++] = mask;
    }
    return count;
}

/* Held-Karp for H(T) */
long long hamiltonian_paths(unsigned long long T) {
    /* dp[mask][v] = number of paths ending at v visiting exactly mask */
    static long long dp[1<<N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1<<v][v] = 1;
    for (int mask = 1; mask < (1<<N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1<<v))) continue;
            if (dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (arc(T, v, u))
                    dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    int full = (1<<N)-1;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

unsigned long long rand64() {
    unsigned long long r = 0;
    for (int i = 0; i < 4; i++)
        r = (r << 16) ^ (unsigned long long)rand();
    return r;
}

int main(int argc, char **argv) {
    long long nsamples = 2000000;
    if (argc > 1) nsamples = atoll(argv[1]);

    srand(time(NULL) ^ getpid());

    int c3_masks[200], c5_masks[1000];
    int all_masks[2000];

    long long h21_count = 0;
    long long no_src_sink_low_alpha = 0;
    int alpha_i2_hist[20][50]; /* [alpha1][i2] counts */
    memset(alpha_i2_hist, 0, sizeof(alpha_i2_hist));
    long long total_low_alpha = 0;

    for (long long s = 0; s < nsamples; s++) {
        unsigned long long T = rand64() & ((1ULL<<NBITS)-1);

        /* Check scores for source/sink */
        int scores[N];
        for (int i = 0; i < N; i++) {
            scores[i] = 0;
            for (int j = 0; j < N; j++) {
                if (i != j && arc(T, i, j)) scores[i]++;
            }
        }
        int has_src_sink = 0;
        for (int i = 0; i < N; i++) {
            if (scores[i] == 0 || scores[i] == N-1) { has_src_sink = 1; break; }
        }

        int t3 = find_3cycles(T, c3_masks);
        if (t3 > 10) continue;

        int t5 = find_5cycles(T, c5_masks);
        int alpha1 = t3 + t5;

        /* Skip 7-cycles and 9-cycles for speed — they only add to alpha1 */
        /* If alpha1 already > 10, skip */
        if (alpha1 > 12) continue;

        /* For alpha1 <= 10 range, compute H directly via Held-Karp */
        long long H = hamiltonian_paths(T);

        if (H == 21) {
            h21_count++;
            printf("*** H=21 FOUND at n=9! T=%llu scores=", T);
            for (int i = 0; i < N; i++) printf("%d ", scores[i]);
            printf("t3=%d t5=%d has_src_sink=%d\n", t3, t5, has_src_sink);
            fflush(stdout);
        }

        /* Track alpha1 with H in range for understanding */
        if (alpha1 <= 12) {
            /* Actually compute i_2 for low alpha1 */
            int ncycles = 0;
            for (int i = 0; i < t3; i++) all_masks[ncycles++] = c3_masks[i];
            for (int i = 0; i < t5; i++) all_masks[ncycles++] = c5_masks[i];

            int i2 = 0;
            for (int i = 0; i < ncycles; i++)
                for (int j = i+1; j < ncycles; j++)
                    if ((all_masks[i] & all_masks[j]) == 0)
                        i2++;

            /* H from OCF (approximate, missing 7 and 9-cycles): */
            /* H_approx = 1 + 2*alpha1 + 4*i2 + ... */
            /* The actual H is from Held-Karp */

            if (alpha1 <= 12 && i2 < 50) {
                if (alpha1 < 20) alpha_i2_hist[alpha1][i2]++;
            }
            total_low_alpha++;

            if (!has_src_sink && alpha1 <= 10) {
                no_src_sink_low_alpha++;
                printf("NO_SRC_SINK: a1=%d i2=%d H=%lld t3=%d t5=%d scores=",
                       alpha1, i2, H, t3, t5);
                for (int i = 0; i < N; i++) printf("%d ", scores[i]);
                printf("\n");
                fflush(stdout);
            }
        }

        if ((s+1) % 200000 == 0) {
            fprintf(stderr, "Progress: %lld/%lld. h21=%lld low_alpha=%lld no_src_sink_low=%lld\n",
                    s+1, nsamples, h21_count, total_low_alpha, no_src_sink_low_alpha);
        }
    }

    printf("\n=== RESULTS (n=9, %lld samples) ===\n", nsamples);
    printf("H=21 found: %lld\n", h21_count);
    printf("Total with alpha1<=12 (3+5 only): %lld\n", total_low_alpha);
    printf("No source/sink with alpha1<=10: %lld\n", no_src_sink_low_alpha);

    printf("\nalpha1 vs i2 distribution (3+5 cycles only):\n");
    for (int a = 0; a < 13; a++) {
        int any = 0;
        for (int i = 0; i < 50; i++) if (alpha_i2_hist[a][i]) any = 1;
        if (!any) continue;
        printf("  alpha1=%d:", a);
        for (int i = 0; i < 50; i++)
            if (alpha_i2_hist[a][i])
                printf(" i2=%d:%d", i, alpha_i2_hist[a][i]);
        printf("\n");
    }

    return 0;
}
"""

def main():
    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, "h21_n9.c")
    exe = os.path.join(tmpdir, "h21_n9")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running 2M random n=9 tournaments...")
    sys.stdout.flush()

    proc = subprocess.run([exe, "2000000"], capture_output=True, text=True, timeout=600)
    # Print results
    lines = proc.stdout.strip().split('\n')
    for l in lines:
        print(l)
    if proc.stderr:
        stderr_lines = proc.stderr.strip().split('\n')
        for l in stderr_lines[-3:]:
            print(l, file=sys.stderr)

if __name__ == "__main__":
    main()
