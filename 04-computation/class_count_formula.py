#!/usr/bin/env python3
"""
class_count_formula.py -- kind-pasteur-2026-03-13-S60

How many H-classes exist at each prime p?

Observed:
  p=7:  2 classes (sizes 2, 6)
  p=11: 4 classes (sizes 2, 10, 10, 10)
  p=13: 5 classes (sizes 4, 24, 12, 12, 12)

The H-class structure is determined by the eigenvalue type
(the multiset of |lambda_j|^2, or equivalently sum|lambda|^{2k}
for all k).

For circulant tournaments on Z_p with S = {one from each {j,p-j}},
the eigenvalue lambda_k = sum_{s in S} omega^{ks}.

Two orientations have the same eigenvalue type iff the
POWER SUMS S_{2k} = sum |lambda_j|^{2k} agree for all k.

This is equivalent to the D^2 multiset being the same (from HYP-682).
D^2_j = |lambda_j|^2 for j=1,...,(p-1)/2 (the non-trivial eigenvalues,
paired by complex conjugation).

The number of distinct D^2 multisets = number of orbits under the
symmetry group acting on the orientation choices.

At p=11: all 32 orientations are non-isomorphic (32 D^2 multisets),
but only 4 distinct sum_4 values. So it's NOT about isomorphism —
it's about a COARSER invariant.

Let me compute the class structure at p=5 and p=7 to see the pattern.
"""

import cmath
import math
from collections import defaultdict

def compute_classes(p):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)

    def bits_to_S(bits):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        return S

    # Group by D^2 multiset (rounded)
    d2_types = defaultdict(list)
    sum4_types = defaultdict(list)

    for bits in range(1 << m):
        S = bits_to_S(bits)
        eigs = [sum(omega**(j*s) for s in S) for j in range(p)]

        # D^2 multiset (non-trivial eigenvalues, sorted)
        d2 = tuple(sorted([round(abs(e)**2, 6) for e in eigs[1:]], reverse=True))
        d2_types[d2].append(bits)

        # sum_4 for coarser grouping
        s4 = round(sum(abs(e)**4 for e in eigs).real, 2)
        sum4_types[s4].append(bits)

    return d2_types, sum4_types

for p in [3, 5, 7, 11, 13, 17, 19]:
    m = (p - 1) // 2
    d2_types, sum4_types = compute_classes(p)

    print(f"\np={p}: m={m}, 2^m={1 << m} orientations")
    print(f"  D^2 types: {len(d2_types)}")
    print(f"  sum_4 types: {len(sum4_types)}")

    d2_sizes = sorted([len(v) for v in d2_types.values()], reverse=True)
    s4_sizes = sorted([len(v) for v in sum4_types.values()], reverse=True)

    print(f"  D^2 class sizes: {d2_sizes}")
    print(f"  sum_4 class sizes: {s4_sizes}")

    # Check if D^2 type = sum_4 type (are they the same partition?)
    same = len(d2_types) == len(sum4_types)
    print(f"  D^2 = sum_4 partition? {same}")

    # Count QR elements per class for sum_4
    S_qr = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

    if p <= 13:
        print(f"  QR = {sorted(S_qr)}")
        for s4, bits_list in sorted(sum4_types.items()):
            qr_counts = []
            for bits in bits_list:
                S = []
                for j in range(m):
                    if bits & (1 << j):
                        S.append(j + 1)
                    else:
                        S.append(p - (j + 1))
                n_qr = sum(1 for s in S if s in S_qr)
                qr_counts.append(n_qr)
            print(f"    sum_4={s4:.1f}: size={len(bits_list)}, "
                  f"#QR range={min(qr_counts)}-{max(qr_counts)}, "
                  f"dist={sorted(set(qr_counts))}")

# OEIS check: number of D^2 types
print(f"\n\nD^2 type counts by prime:")
for p in [3, 5, 7, 11, 13, 17, 19]:
    d2_types, _ = compute_classes(p)
    print(f"  p={p}: {len(d2_types)} types")

# Also: number of sum_4 types
print(f"\nsum_4 type counts by prime:")
for p in [3, 5, 7, 11, 13, 17, 19]:
    _, s4_types = compute_classes(p)
    print(f"  p={p}: {len(s4_types)} types")

print("\nDONE.")
