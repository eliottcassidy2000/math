"""
h21_alpha_structure_n7.py — Exhaustive analysis of alpha_1=8 and alpha_1=10 at n=7

For the H=21 gap proof, we need to understand:
1. When alpha_1=10: why i_2=2 always (not 0)?
2. When alpha_1=8: why i_2=0 always (not 1)?

This script extracts the full Omega(T) structure for these cases.

Author: opus-2026-03-07-S41
"""

import ctypes, os, sys, tempfile, subprocess
from itertools import combinations
from collections import Counter, defaultdict

# Build C extension for fast tournament enumeration
C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Tournament on n=7 vertices, encoded as 21 bits
// Bit (i,j) for i<j at position i*(2*7-i-1)/2 + (j-i-1)

static inline int arc(int T, int n, int i, int j) {
    // Returns 1 if i->j, 0 if j->i
    if (i < j) {
        int pos = i*(2*n-i-1)/2 + (j-i-1);
        return (T >> pos) & 1;
    } else {
        int pos = j*(2*n-j-1)/2 + (i-j-1);
        return 1 - ((T >> pos) & 1);
    }
}

// Check if (a,b,c) is a directed 3-cycle a->b->c->a
static inline int is_3cycle(int T, int n, int a, int b, int c) {
    return arc(T,n,a,b) && arc(T,n,b,c) && arc(T,n,c,a);
}

// Check if vertices in perm[] form a directed cycle
static inline int is_directed_cycle(int T, int n, int *perm, int len) {
    for (int i = 0; i < len; i++) {
        if (!arc(T, n, perm[i], perm[(i+1)%len]))
            return 0;
    }
    return 1;
}

// Find all directed 3-cycles
int find_3cycles(int T, int n, int cycles[][3]) {
    int count = 0;
    for (int a = 0; a < n; a++)
        for (int b = a+1; b < n; b++)
            for (int c = b+1; c < n; c++) {
                if (is_3cycle(T, n, a, b, c)) {
                    cycles[count][0] = a; cycles[count][1] = b; cycles[count][2] = c;
                    count++;
                } else if (is_3cycle(T, n, a, c, b)) {
                    cycles[count][0] = a; cycles[count][1] = c; cycles[count][2] = b;
                    count++;
                }
            }
    return count;
}

// Find all directed 5-cycles (as vertex sets, not oriented)
// Returns count and fills vsets[][5] with sorted vertex sets
int find_5cycle_vsets(int T, int n, int vsets[][5]) {
    int count = 0;
    int perm[5];
    for (int a = 0; a < n; a++)
    for (int b = a+1; b < n; b++)
    for (int c = b+1; c < n; c++)
    for (int d = c+1; d < n; d++)
    for (int e = d+1; e < n; e++) {
        int vv[5] = {a,b,c,d,e};
        // Check all 12 distinct directed 5-cycles on these vertices
        // (5-1)!/2 = 12 distinct cyclic orderings
        int found = 0;
        // Generate all permutations of vv and check
        int p[5];
        for (int i0 = 0; i0 < 5 && !found; i0++)
        for (int i1 = 0; i1 < 5 && !found; i1++) {
            if (i1==i0) continue;
            for (int i2 = 0; i2 < 5 && !found; i2++) {
                if (i2==i0||i2==i1) continue;
                for (int i3 = 0; i3 < 5 && !found; i3++) {
                    if (i3==i0||i3==i1||i3==i2) continue;
                    int i4 = 10-i0-i1-i2-i3;
                    p[0]=vv[i0]; p[1]=vv[i1]; p[2]=vv[i2]; p[3]=vv[i3]; p[4]=vv[i4];
                    // Only check if p[0] < p[1] to avoid double counting
                    if (p[0] > p[1]) continue;
                    if (is_directed_cycle(T, n, p, 5)) {
                        found = 1;
                    }
                }
            }
        }
        if (found) {
            vsets[count][0]=a; vsets[count][1]=b; vsets[count][2]=c;
            vsets[count][3]=d; vsets[count][4]=e;
            count++;
        }
    }
    return count;
}

// Count directed 5-cycles (individual oriented cycles, not vertex sets)
int count_5cycles(int T, int n) {
    int count = 0;
    int perm[5];
    for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++) {
        if (b==a) continue;
        for (int c = 0; c < n; c++) {
            if (c==a||c==b) continue;
            for (int d = 0; d < n; d++) {
                if (d==a||d==b||d==c) continue;
                for (int e = 0; e < n; e++) {
                    if (e==a||e==b||e==c||e==d) continue;
                    perm[0]=a; perm[1]=b; perm[2]=c; perm[3]=d; perm[4]=e;
                    if (is_directed_cycle(T, n, perm, 5))
                        count++;
                }
            }
        }
    }
    return count / 5; // Each cycle counted 5 times (cyclic shifts)
}

// Main analysis: for each tournament, compute alpha_1 and i_2
// Output lines: "T alpha1 i2 t3 t5_vsets t7"
void analyze_all(int n) {
    int nbits = n*(n-1)/2;
    long long total = 1LL << nbits;

    int cycles3[100][3];
    int vsets5[100][5];

    for (long long T = 0; T < total; T++) {
        int t3 = find_3cycles(T, n, cycles3);

        // Quick filter: only look at alpha_1=8 or alpha_1=10 candidates
        // At n=7, t5+t7 add to alpha_1, so alpha_1 >= t3
        // For alpha_1=10, need t3 <= 10
        // For alpha_1=8, need t3 <= 8
        if (t3 > 10) continue;

        int t5v = find_5cycle_vsets(T, n, vsets5);

        // Count directed 5-cycles per vertex set
        // For now, count vertex sets that support at least one directed 5-cycle
        // Each such vertex set is ONE vertex in Omega
        // Actually, each DIRECTED cycle is a separate vertex in Omega
        // For simplicity, let's count directed cycles
        int t5 = count_5cycles(T, n);

        // Check 7-cycle (Hamiltonian cycle)
        int t7 = 0;
        {
            // Check all (7-1)!/2 = 360 distinct directed 7-cycles
            int p[7] = {0,1,2,3,4,5,6};
            // Fix vertex 0 first, permute rest
            int rest[6] = {1,2,3,4,5,6};
            // Generate all 720 permutations of rest
            // Check if 0->rest[0]->...->rest[5]->0 is a cycle
            // Only count if rest[0] < rest[5] (to avoid double counting)
            // Actually, just count all and divide by 7*2=14... no.
            // Directed 7-cycle counted once per starting vertex = 7 times per cycle.
            // So count all starting from 0 and that IS the count.
            // 0->p1->p2->p3->p4->p5->p6->0
            int pp[7];
            pp[0] = 0;
            for (int i1=1;i1<7;i1++)
            for (int i2=1;i2<7;i2++) { if(i2==i1) continue;
            for (int i3=1;i3<7;i3++) { if(i3==i1||i3==i2) continue;
            for (int i4=1;i4<7;i4++) { if(i4==i1||i4==i2||i4==i3) continue;
            for (int i5=1;i5<7;i5++) { if(i5==i1||i5==i2||i5==i3||i5==i4) continue;
            {
                int i6 = 21-i1-i2-i3-i4-i5;
                pp[1]=i1; pp[2]=i2; pp[3]=i3; pp[4]=i4; pp[5]=i5; pp[6]=i6;
                if (is_directed_cycle(T, n, pp, 7))
                    t7++;
            }}}}}
        }

        int alpha1 = t3 + t5 + t7;

        if (alpha1 != 8 && alpha1 != 10) continue;

        // Now compute i_2: number of pairs of vertex-disjoint odd cycles
        // Build list of all odd cycles (as vertex sets + cycle info)
        // For 3-cycles: vertex set from cycles3
        // For 5-cycles: vertex set from vsets5 (but need directed cycles)
        // For 7-cycle: vertex set = {0,...,6}

        // For i_2 computation, two cycles are "disjoint" if they share no vertex
        // In Omega, they are non-adjacent (independent set of size 2)

        // Build vertex bitmasks for each cycle
        int all_masks[200];
        int ncycles = 0;

        for (int i = 0; i < t3; i++) {
            all_masks[ncycles++] = (1<<cycles3[i][0]) | (1<<cycles3[i][1]) | (1<<cycles3[i][2]);
        }

        // For 5-cycles: need individual directed cycles, not just vertex sets
        // Two directed 5-cycles on the same vertex set have the same mask
        // They are always adjacent (sharing all vertices)
        // So for i_2, we can use vertex sets
        // But alpha_1 counts directed cycles, not vertex sets
        // However, for i_2 (independent pairs), cycles on same vertex set are never independent
        // So i_2 = number of pairs with disjoint VERTEX SETS

        // Actually, let me count directed 5-cycles individually but use masks
        // Multiple directed 5-cycles on same vertex set all have same mask

        // For 5-cycle vertex sets with multiplicity
        for (int i = 0; i < t5v; i++) {
            int mask = 0;
            for (int j = 0; j < 5; j++) mask |= (1 << vsets5[i][j]);
            // How many directed 5-cycles on this vertex set?
            int v[5] = {vsets5[i][0], vsets5[i][1], vsets5[i][2], vsets5[i][3], vsets5[i][4]};
            int cnt = 0;
            int pp[5];
            for (int i0=0;i0<5;i0++)
            for (int i1=0;i1<5;i1++) { if(i1==i0) continue;
            for (int i2=0;i2<5;i2++) { if(i2==i0||i2==i1) continue;
            for (int i3=0;i3<5;i3++) { if(i3==i0||i3==i1||i3==i2) continue;
            {
                int i4 = 10-i0-i1-i2-i3;
                pp[0]=v[i0]; pp[1]=v[i1]; pp[2]=v[i2]; pp[3]=v[i3]; pp[4]=v[i4];
                if (pp[0] > pp[1]) continue; // avoid double-counting
                if (is_directed_cycle(T, n, pp, 5)) cnt++;
            }}}}
            // cnt should equal the number of directed 5-cycles (each counted once)
            // Actually with the p[0]>p[1] filter, we get each cycle once
            for (int k = 0; k < cnt; k++)
                all_masks[ncycles++] = mask;
        }

        // 7-cycles
        for (int k = 0; k < t7; k++)
            all_masks[ncycles++] = 0x7F; // all 7 bits

        // Compute i_2
        int i2 = 0;
        for (int i = 0; i < ncycles; i++)
            for (int j = i+1; j < ncycles; j++)
                if ((all_masks[i] & all_masks[j]) == 0)
                    i2++;

        printf("%lld %d %d %d %d %d\n", T, alpha1, i2, t3, t5, t7);
    }
}

int main() {
    analyze_all(7);
    return 0;
}
"""

def main():
    # Compile and run C code
    tmpdir = tempfile.mkdtemp()
    src = os.path.join(tmpdir, "alpha_struct.c")
    exe = os.path.join(tmpdir, "alpha_struct")

    with open(src, 'w') as f:
        f.write(C_CODE)

    subprocess.run(["gcc", "-O3", "-o", exe, src], check=True)
    print("Compiled. Running exhaustive n=7 analysis (filtering alpha_1=8,10)...")
    print("This may take a few minutes...")
    sys.stdout.flush()

    proc = subprocess.Popen([exe], stdout=subprocess.PIPE, text=True)

    results = defaultdict(list)
    count = 0

    for line in proc.stdout:
        parts = line.strip().split()
        T, alpha1, i2, t3, t5, t7 = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
        key = (alpha1, i2, t3, t5, t7)
        results[key].append(T)
        count += 1
        if count % 10000 == 0:
            print(f"  Found {count} tournaments so far...", file=sys.stderr)

    proc.wait()

    print(f"\nTotal tournaments with alpha_1 in {{8,10}}: {count}")
    print(f"\nDistribution of (alpha_1, i_2, t3, t5, t7) -> count:")
    for key in sorted(results.keys()):
        alpha1, i2, t3, t5, t7 = key
        print(f"  alpha_1={alpha1}, i_2={i2}, t3={t3}, t5={t5}, t7={t7}: {len(results[key])} tournaments")

    # For alpha_1=10, show detailed structure
    print("\n=== alpha_1=10 analysis ===")
    a10_keys = [k for k in results if k[0]==10]
    for key in sorted(a10_keys):
        alpha1, i2, t3, t5, t7 = key
        print(f"\nalpha_1=10, i_2={i2}, t3={t3}, t5={t5}, t7={t7}: {len(results[key])} tournaments")
        # Show first example
        T = results[key][0]
        print(f"  Example tournament T={T} (binary: {bin(T)})")
        # Show adjacency matrix
        n = 7
        print("  Adjacency matrix:")
        for i in range(n):
            row = []
            for j in range(n):
                if i == j:
                    row.append('.')
                elif i < j:
                    pos = i*(2*n-i-1)//2 + (j-i-1)
                    row.append(str((T >> pos) & 1))
                else:
                    pos = j*(2*n-j-1)//2 + (i-j-1)
                    row.append(str(1 - ((T >> pos) & 1)))
            print(f"    {''.join(row)}")

    # For alpha_1=8, show structure
    print("\n=== alpha_1=8 analysis ===")
    a8_keys = [k for k in results if k[0]==8]
    for key in sorted(a8_keys):
        alpha1, i2, t3, t5, t7 = key
        print(f"\nalpha_1=8, i_2={i2}, t3={t3}, t5={t5}, t7={t7}: {len(results[key])} tournaments")
        T = results[key][0]
        print(f"  Example tournament T={T}")

if __name__ == "__main__":
    main()
