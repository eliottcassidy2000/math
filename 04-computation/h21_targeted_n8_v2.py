"""
h21_targeted_n8_v2.py — Fixed 5-cycle counting for H=21 search at n=8

Bug fix: use min-vertex start for canonical directed 5-cycle deduplication.

Author: opus-2026-03-07-S41
"""

import os, sys, tempfile, subprocess

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define N 8
#define NBITS (N*(N-1)/2)

static inline int arc(unsigned int T, int i, int j) {
    if (i < j) {
        int pos = i*(2*N-i-1)/2 + (j-i-1);
        return (T >> pos) & 1;
    } else {
        int pos = j*(2*N-j-1)/2 + (i-j-1);
        return 1 - ((T >> pos) & 1);
    }
}

int find_3cycles(unsigned int T, int *masks) {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = arc(T,a,b), bc = arc(T,b,c), ca = arc(T,c,a);
                int ac = 1-ca, cb = 1-bc, ba = 1-ab;
                if ((ab && bc && ca) || (ac && cb && ba)) {
                    masks[count++] = (1<<a)|(1<<b)|(1<<c);
                }
            }
    return count;
}

/* Count directed 5-cycles properly.
   For each 5-vertex subset, enumerate all directed 5-cycles.
   Canonicalize: start at the minimum vertex (unique among 5 rotations).
   Use a hash set to track seen canonical forms.
*/
int find_5cycles(unsigned int T, int *masks) {
    int count = 0;
    int vv[5], p[5];

    for (int a = 0; a < N; a++)
    for (int b = a+1; b < N; b++)
    for (int c = b+1; c < N; c++)
    for (int d = c+1; d < N; d++)
    for (int e = d+1; e < N; e++) {
        vv[0]=a; vv[1]=b; vv[2]=c; vv[3]=d; vv[4]=e;
        int mask = (1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e);

        /* Track canonical forms for this vertex set.
           A canonical form is (min_vertex, v1, v2, v3, v4) where
           min_vertex is the start. Store as a packed int. */
        int seen[12]; /* max 12 directed 5-cycles on 5 vertices */
        int nseen = 0;

        for (int i0 = 0; i0 < 5; i0++) {
            p[0] = vv[i0];
            /* Only consider rotations starting at the minimum vertex */
            if (p[0] != a) continue; /* a is min since a < b < c < d < e */

            for (int i1 = 0; i1 < 5; i1++) {
                if (i1==i0) continue;
                p[1] = vv[i1];
                for (int i2 = 0; i2 < 5; i2++) {
                    if (i2==i0||i2==i1) continue;
                    p[2] = vv[i2];
                    for (int i3 = 0; i3 < 5; i3++) {
                        if (i3==i0||i3==i1||i3==i2) continue;
                        p[3] = vv[i3];
                        int i4 = 10-i0-i1-i2-i3;
                        p[4] = vv[i4];

                        /* Check if directed cycle */
                        int ok = 1;
                        for (int q = 0; q < 5; q++) {
                            if (!arc(T, p[q], p[(q+1)%5])) { ok=0; break; }
                        }
                        if (!ok) continue;

                        /* Canonical: among this cycle's rotations, the one starting at min.
                           Since p[0]=a=min, this IS the canonical form.
                           But wait: rotation (a,p1,p2,p3,p4) and (a,p4,p3,p2,p1) are
                           DIFFERENT directed cycles (reverse). Both start at a.
                           So we may find multiple canonical forms for same vertex set. */

                        /* Pack as unique ID: encode p[1..4] positions */
                        int id = i1*125 + i2*25 + i3*5 + i4;

                        /* Check if already seen (for this vertex set) */
                        int dup = 0;
                        for (int s = 0; s < nseen; s++) {
                            if (seen[s] == id) { dup = 1; break; }
                        }
                        if (!dup) {
                            seen[nseen++] = id;
                        }
                    }
                }
            }
        }
        for (int k = 0; k < nseen; k++)
            masks[count++] = mask;
    }
    return count;
}

int find_7cycles(unsigned int T, int *masks) {
    int count = 0;
    int vv[7], p[7];

    for (int a = 0; a < N; a++)
    for (int b = a+1; b < N; b++)
    for (int c = b+1; c < N; c++)
    for (int d = c+1; d < N; d++)
    for (int e = d+1; e < N; e++)
    for (int f = e+1; f < N; f++)
    for (int g = f+1; g < N; g++) {
        vv[0]=a; vv[1]=b; vv[2]=c; vv[3]=d; vv[4]=e; vv[5]=f; vv[6]=g;
        int mask = (1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e)|(1<<f)|(1<<g);

        int nseen = 0;
        int rest[6];
        p[0] = vv[0]; /* = a = min vertex */
        for (int x = 0; x < 6; x++) rest[x] = vv[x+1];

        for (int i1=0;i1<6;i1++) {
            p[1]=rest[i1];
            for (int i2=0;i2<6;i2++) { if(i2==i1) continue;
                p[2]=rest[i2];
            for (int i3=0;i3<6;i3++) { if(i3==i1||i3==i2) continue;
                p[3]=rest[i3];
            for (int i4=0;i4<6;i4++) { if(i4==i1||i4==i2||i4==i3) continue;
                p[4]=rest[i4];
            for (int i5=0;i5<6;i5++) { if(i5==i1||i5==i2||i5==i3||i5==i4) continue;
                p[5]=rest[i5];
                int i6 = 15-i1-i2-i3-i4-i5;
                p[6]=rest[i6];
                int ok = 1;
                for (int q = 0; q < 7; q++) {
                    if (!arc(T, p[q], p[(q+1)%7])) { ok=0; break; }
                }
                if (ok) nseen++;
            }}}}}
        for (int k = 0; k < nseen; k++)
            masks[count++] = mask;
    }
    return count;
}

unsigned int rand32() {
    return ((unsigned int)rand() << 17) ^ ((unsigned int)rand() << 2) ^ (unsigned int)rand();
}

/* Held-Karp for H(T) verification */
long long hamiltonian_paths(unsigned int T) {
    long long dp[1<<N][N];
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

int main(int argc, char **argv) {
    long long nsamples = 5000000;
    if (argc > 1) nsamples = atoll(argv[1]);

    srand(time(NULL) ^ getpid());

    int all_masks[2000];
    int c3_masks[200], c5_masks[500], c7_masks[500];

    long long count_8_1 = 0, count_10_0 = 0;
    long long total_a8 = 0, total_a10 = 0;
    int a8_i2_dist[100] = {0};
    int a10_i2_dist[100] = {0};
    long long h21_count = 0;

    for (long long s = 0; s < nsamples; s++) {
        unsigned int T = rand32() & ((1<<NBITS)-1);

        int t3 = find_3cycles(T, c3_masks);
        if (t3 > 10) continue;

        int t5 = find_5cycles(T, c5_masks);
        if (t3 + t5 > 12) continue;

        int t7 = find_7cycles(T, c7_masks);
        int alpha1 = t3 + t5 + t7;

        if (alpha1 != 8 && alpha1 != 10) continue;

        int ncycles = 0;
        for (int i = 0; i < t3; i++) all_masks[ncycles++] = c3_masks[i];
        for (int i = 0; i < t5; i++) all_masks[ncycles++] = c5_masks[i];
        for (int i = 0; i < t7; i++) all_masks[ncycles++] = c7_masks[i];

        int i2 = 0;
        for (int i = 0; i < ncycles; i++)
            for (int j = i+1; j < ncycles; j++)
                if ((all_masks[i] & all_masks[j]) == 0)
                    i2++;

        /* Compute H for verification when i2 matches target */
        int target = (alpha1 == 8 && i2 == 1) || (alpha1 == 10 && i2 == 0);

        if (alpha1 == 8) {
            total_a8++;
            if (i2 < 100) a8_i2_dist[i2]++;
            if (i2 == 1) {
                count_8_1++;
                long long H = hamiltonian_paths(T);
                printf("CHECK (8,1): T=%u t3=%d t5=%d t7=%d i2=%d H=%lld\n", T, t3, t5, t7, i2, H);
                if (H == 21) {
                    printf("*** H=21 FOUND! ***\n");
                    h21_count++;
                }
                fflush(stdout);
            }
        }
        if (alpha1 == 10) {
            total_a10++;
            if (i2 < 100) a10_i2_dist[i2]++;
            if (i2 == 0) {
                count_10_0++;
                long long H = hamiltonian_paths(T);
                printf("CHECK (10,0): T=%u t3=%d t5=%d t7=%d i2=%d H=%lld\n", T, t3, t5, t7, i2, H);
                if (H == 21) {
                    printf("*** H=21 FOUND! ***\n");
                    h21_count++;
                }
                fflush(stdout);
            }
        }

        if ((s+1) % 500000 == 0) {
            fprintf(stderr, "Progress: %lld/%lld. a8=%lld a10=%lld (8,1)=%lld (10,0)=%lld h21=%lld\n",
                    s+1, nsamples, total_a8, total_a10, count_8_1, count_10_0, h21_count);
        }
    }

    printf("\n=== RESULTS (n=8, %lld samples) ===\n", nsamples);
    printf("alpha_1=8: %lld tournaments\n", total_a8);
    printf("alpha_1=10: %lld tournaments\n", total_a10);
    printf("(8,1) found: %lld\n", count_8_1);
    printf("(10,0) found: %lld\n", count_10_0);
    printf("H=21 found: %lld\n", h21_count);

    printf("\nalpha_1=8 i_2 distribution:\n");
    for (int i = 0; i < 100; i++)
        if (a8_i2_dist[i] > 0)
            printf("  i_2=%d: %d\n", i, a8_i2_dist[i]);

    printf("\nalpha_1=10 i_2 distribution:\n");
    for (int i = 0; i < 100; i++)
        if (a10_i2_dist[i] > 0)
            printf("  i_2=%d: %d\n", i, a10_i2_dist[i]);

    return 0;
}
"""

def main():
    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, "h21_n8_v2.c")
    exe = os.path.join(tmpdir, "h21_n8_v2")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running 5M random n=8 tournaments with H verification...")
    sys.stdout.flush()

    proc = subprocess.run([exe, "5000000"], capture_output=True, text=True, timeout=600)
    print(proc.stdout[-3000:] if len(proc.stdout) > 3000 else proc.stdout)
    if proc.stderr:
        lines = proc.stderr.strip().split('\n')
        for l in lines[-3:]:
            print(l, file=sys.stderr)

if __name__ == "__main__":
    main()
