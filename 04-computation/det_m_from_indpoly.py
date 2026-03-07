#!/usr/bin/env python3
"""
Can det(M) be expressed in terms of independence polynomial I(Omega,x)?

We know: tr(M) = H = I(Omega, 2)
We know: at n=5, det(M) takes 8 distinct values across 10 isomorphism classes.
We know: I(Omega, x) at n=5 is always linear (I = 1 + alpha_1*x) since
         no two disjoint 3-cycles can fit in 5 vertices.

So at n=5: I(Omega, x) = 1 + c3*x where c3 = alpha_1 = number of 3-cycles.
H = I(2) = 1 + 2*c3, so c3 = (H-1)/2.

The I polynomial is determined by H alone at n=5.
But det(M) is NOT determined by H (e.g., H=5 gives det=9 or 17).
So det(M) carries ADDITIONAL information beyond the independence polynomial.

At n=7: I(Omega, x) can be quadratic (two disjoint 3-cycles possible).
I = 1 + alpha_1*x + alpha_2*x^2 where alpha_1 = c3, alpha_2 = # disjoint 3-cycle pairs.
H = 1 + 2*c3 + 4*alpha_2.

QUESTION: Does det(M) depend on the "cycle structure" beyond what I captures?

APPROACH: At n=7, sample tournaments with SAME I polynomial but different det(M).
If such exist, det(M) is strictly finer than I.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np
from collections import defaultdict
import random

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def count_3cycles(T, n):
    c3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check both cycle directions
        if T.get((i,j),0) and T.get((j,k),0) and T.get((k,i),0):
            c3 += 1
        if T.get((i,k),0) and T.get((k,j),0) and T.get((j,i),0):
            c3 += 1
    return c3

def count_disjoint_3cycle_pairs(T, n):
    # Find all 3-cycles
    cycles = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if T.get((i,j),0) and T.get((j,k),0) and T.get((k,i),0):
            cycles.append(frozenset(triple))
        elif T.get((i,k),0) and T.get((k,j),0) and T.get((j,i),0):
            cycles.append(frozenset(triple))
    # Count disjoint pairs
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_2

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        M[a, a] = val
        for b in range(a+1, n):
            U = [v for v in range(n) if v != a and v != b]
            val = 0
            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                sign = (-1)**len(S_list)
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})
                ea = 0
                if len(S_set) == 1:
                    ea = 1
                else:
                    for p in permutations(S_set):
                        if p[-1] != a: continue
                        prod = 1
                        for pk in range(len(p)-1):
                            prod *= T.get((p[pk], p[pk+1]), 0)
                        ea += prod
                bb2 = 0
                if len(R_set) == 1:
                    bb2 = 1
                else:
                    for p in permutations(R_set):
                        if p[0] != b: continue
                        prod = 1
                        for pk in range(len(p)-1):
                            prod *= T.get((p[pk], p[pk+1]), 0)
                        bb2 += prod
                val += sign * ea * bb2
            M[a, b] = val
            M[b, a] = val
    return M


# ============================================================
# n=5: det(M) vs (H, c3) — exhaustive
# ============================================================
print("=" * 70)
print("n=5: det(M) vs independence polynomial")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

seen = {}  # key=(H, c3, scores) -> det_M
for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    c3 = count_3cycles(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    key = (H, c3, scores)
    if key in seen:
        continue

    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    seen[key] = det_M

# Group by (H, c3)
hc3_to_dets = defaultdict(set)
for (H, c3, scores), det_M in seen.items():
    hc3_to_dets[(H, c3)].add(det_M)

print(f"\n  (H, c3) -> det(M):")
for (H, c3), dets in sorted(hc3_to_dets.items()):
    if len(dets) > 1:
        print(f"    H={H:>3}, c3={c3:>2}: det(M) = {sorted(dets)} -- NOT determined by (H,c3)!")
    else:
        print(f"    H={H:>3}, c3={c3:>2}: det(M) = {sorted(dets)[0]}")


# ============================================================
# n=5: det(M) vs (H, c3, scores) — is it determined by scores?
# ============================================================
print(f"\n  (H, c3, scores) -> det(M):")
scores_to_dets = defaultdict(set)
for (H, c3, scores), det_M in seen.items():
    scores_to_dets[(H, scores)].add(det_M)

for (H, scores), dets in sorted(scores_to_dets.items()):
    if len(dets) > 1:
        print(f"    H={H:>3}, scores={scores}: det(M) = {sorted(dets)} -- NOT determined!")
    else:
        print(f"    H={H:>3}, scores={scores}: det(M) = {sorted(dets)[0]}")


# ============================================================
# What DOES determine det(M)?
# ============================================================
print("\n" + "=" * 70)
print("det(M) classification at n=5")
print("=" * 70)

det_to_info = defaultdict(list)
for (H, c3, scores), det_M in seen.items():
    det_to_info[det_M].append((H, c3, scores))

for det_M, infos in sorted(det_to_info.items()):
    print(f"\n  det(M) = {det_M}:")
    for H, c3, scores in sorted(infos):
        print(f"    H={H}, c3={c3}, scores={scores}")


# ============================================================
# n=5: tr(M^2) as additional invariant
# ============================================================
print("\n" + "=" * 70)
print("n=5: tr(M^2) as invariant")
print("=" * 70)

seen2 = {}
for bits in range(1 << len(pairs)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    c3 = count_3cycles(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    key = (H, c3, scores)
    if key in seen2:
        continue

    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    tr_M2 = int(round(np.trace(M.astype(float) @ M.astype(float))))
    frobenius = int(round(np.sum(M.astype(float)**2)))

    seen2[key] = (det_M, tr_M2, frobenius)

print(f"\n  {'H':>3} {'c3':>3} {'scores':<20} {'det':>8} {'tr(M^2)':>8} {'||M||_F^2':>10}")
for (H, c3, scores), (det_M, tr_M2, frob) in sorted(seen2.items()):
    print(f"  {H:>3} {c3:>3} {str(scores):<20} {det_M:>8} {tr_M2:>8} {frob:>10}")


# ============================================================
# Does (H, tr(M^2)) determine det(M)?
# ============================================================
print("\n" + "=" * 70)
print("Does (H, tr(M^2)) determine det(M)?")
print("=" * 70)

h_trm2_to_det = defaultdict(set)
for (H, c3, scores), (det_M, tr_M2, frob) in seen2.items():
    h_trm2_to_det[(H, tr_M2)].add(det_M)

for (H, tr_M2), dets in sorted(h_trm2_to_det.items()):
    if len(dets) > 1:
        print(f"  H={H:>3}, tr(M^2)={tr_M2:>5}: det = {sorted(dets)} -- NOT determined!")
    else:
        print(f"  H={H:>3}, tr(M^2)={tr_M2:>5}: det = {sorted(dets)[0]}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
