"""
h21_targeted_n8.py — Targeted search for H=21 at n=8

Specifically searches for:
1. alpha_1=10, i_2=0 (the (10,0) case)
2. alpha_1=8, i_2=1 (the (8,1) case)

Uses C for speed. Random sampling of tournaments.

Author: opus-2026-03-07-S41
"""

import ctypes, os, sys, tempfile, subprocess
from collections import Counter, defaultdict
import random

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

static inline int is_3cycle(unsigned int T, int a, int b, int c) {
    return arc(T,a,b) && arc(T,b,c) && arc(T,c,a);
}

// Find all directed 3-cycles, return count and fill masks
int find_3cycles(unsigned int T, int *masks) {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                if (is_3cycle(T, a, b, c) || is_3cycle(T, a, c, b)) {
                    masks[count++] = (1<<a)|(1<<b)|(1<<c);
                }
            }
    return count;
}

// Find directed 5-cycles (as vertex-set masks with multiplicity)
int find_5cycles(unsigned int T, int *masks) {
    int count = 0;
    int vv[5];
    for (int a = 0; a < N; a++)
    for (int b = a+1; b < N; b++)
    for (int c = b+1; c < N; c++)
    for (int d = c+1; d < N; d++)
    for (int e = d+1; e < N; e++) {
        vv[0]=a; vv[1]=b; vv[2]=c; vv[3]=d; vv[4]=e;
        int mask = (1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e);
        // Count distinct directed 5-cycles on this vertex set
        // Fix first vertex, try all orderings of the rest
        int ncyc = 0;
        int p[5];
        for (int i0 = 0; i0 < 5; i0++) {
            p[0] = vv[i0];
            for (int i1 = 0; i1 < 5; i1++) {
                if (i1==i0) continue;
                p[1] = vv[i1];
                if (p[0] > p[1]) continue; // canonical: first < second
                for (int i2 = 0; i2 < 5; i2++) {
                    if (i2==i0||i2==i1) continue;
                    p[2] = vv[i2];
                    for (int i3 = 0; i3 < 5; i3++) {
                        if (i3==i0||i3==i1||i3==i2) continue;
                        p[3] = vv[i3];
                        int i4 = 10-i0-i1-i2-i3;
                        p[4] = vv[i4];
                        int ok = 1;
                        for (int q = 0; q < 5; q++) {
                            if (!arc(T, p[q], p[(q+1)%5])) { ok=0; break; }
                        }
                        if (ok) ncyc++;
                    }
                }
            }
        }
        for (int k = 0; k < ncyc; k++)
            masks[count++] = mask;
    }
    return count;
}

// Find directed 7-cycles
int find_7cycles(unsigned int T, int *masks) {
    int count = 0;
    int vv[7];
    for (int a = 0; a < N; a++)
    for (int b = a+1; b < N; b++)
    for (int c = b+1; c < N; c++)
    for (int d = c+1; d < N; d++)
    for (int e = d+1; e < N; e++)
    for (int f = e+1; f < N; f++)
    for (int g = f+1; g < N; g++) {
        vv[0]=a; vv[1]=b; vv[2]=c; vv[3]=d; vv[4]=e; vv[5]=f; vv[6]=g;
        int mask = (1<<a)|(1<<b)|(1<<c)|(1<<d)|(1<<e)|(1<<f)|(1<<g);
        // Count directed 7-cycles: fix vertex 0, try all 720 orderings of rest
        int ncyc = 0;
        int rest[6], p[7];
        p[0] = vv[0];
        for (int x = 0; x < 6; x++) rest[x] = vv[x+1];
        // Generate all 720 permutations of rest
        // Use simple nested loops
        for (int i1=0;i1<6;i1++) {
            p[1]=rest[i1];
            if (p[0] > p[1]) continue; // canonical
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
                if (ok) ncyc++;
            }}}}}
        for (int k = 0; k < ncyc; k++)
            masks[count++] = mask;
    }
    return count;
}

// Random 32-bit int
unsigned int rand32() {
    return ((unsigned int)rand() << 17) ^ ((unsigned int)rand() << 2) ^ (unsigned int)rand();
}

int main(int argc, char **argv) {
    long long nsamples = 2000000;
    if (argc > 1) nsamples = atoll(argv[1]);

    srand(time(NULL));

    int all_masks[2000];
    int c3_masks[200], c5_masks[500], c7_masks[500];

    long long count_8_1 = 0, count_10_0 = 0;
    long long total_a8 = 0, total_a10 = 0;
    int a8_i2_dist[100] = {0};
    int a10_i2_dist[100] = {0};

    for (long long s = 0; s < nsamples; s++) {
        unsigned int T = rand32() & ((1<<NBITS)-1);

        int t3 = find_3cycles(T, c3_masks);

        // Quick filter
        if (t3 > 12) continue; // alpha_1 >= t3, and we want alpha_1 <= 10

        int t5 = find_5cycles(T, c5_masks);
        int alpha1_so_far = t3 + t5;

        if (alpha1_so_far > 12) continue;

        int t7 = find_7cycles(T, c7_masks);
        int alpha1 = alpha1_so_far + t7;

        if (alpha1 != 8 && alpha1 != 10) continue;

        // Compute i_2
        int ncycles = 0;
        for (int i = 0; i < t3; i++) all_masks[ncycles++] = c3_masks[i];
        for (int i = 0; i < t5; i++) all_masks[ncycles++] = c5_masks[i];
        for (int i = 0; i < t7; i++) all_masks[ncycles++] = c7_masks[i];

        int i2 = 0;
        for (int i = 0; i < ncycles; i++)
            for (int j = i+1; j < ncycles; j++)
                if ((all_masks[i] & all_masks[j]) == 0)
                    i2++;

        if (alpha1 == 8) {
            total_a8++;
            if (i2 < 100) a8_i2_dist[i2]++;
            if (i2 == 1) {
                count_8_1++;
                printf("FOUND (8,1): T=%u t3=%d t5=%d t7=%d i2=%d\n", T, t3, t5, t7, i2);
                fflush(stdout);
            }
        }
        if (alpha1 == 10) {
            total_a10++;
            if (i2 < 100) a10_i2_dist[i2]++;
            if (i2 == 0) {
                count_10_0++;
                printf("FOUND (10,0): T=%u t3=%d t5=%d t7=%d i2=%d\n", T, t3, t5, t7, i2);
                fflush(stdout);
            }
        }

        if ((s+1) % 100000 == 0) {
            fprintf(stderr, "Progress: %lld/%lld samples. a8=%lld, a10=%lld, (8,1)=%lld, (10,0)=%lld\n",
                    s+1, nsamples, total_a8, total_a10, count_8_1, count_10_0);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Total samples: %lld\n", nsamples);
    printf("alpha_1=8: %lld tournaments\n", total_a8);
    printf("alpha_1=10: %lld tournaments\n", total_a10);
    printf("(8,1) found: %lld\n", count_8_1);
    printf("(10,0) found: %lld\n", count_10_0);

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
    src = os.path.join(tmpdir, "h21_n8.c")
    exe = os.path.join(tmpdir, "h21_n8")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running 2M random n=8 tournaments...")
    sys.stdout.flush()

    proc = subprocess.run([exe, "2000000"], capture_output=True, text=True, timeout=600)
    print(proc.stdout)
    if proc.stderr:
        # Print last few progress lines
        lines = proc.stderr.strip().split('\n')
        for l in lines[-3:]:
            print(l, file=sys.stderr)

if __name__ == "__main__":
    main()
