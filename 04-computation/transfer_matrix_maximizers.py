#!/usr/bin/env python3
"""
QUESTION (opus-S11b): Does M = (H/n)*I for ALL H-maximizers, not just VT?

At n=5: All 24 regular tournaments have H=15, M=3I.
At n=7: VT tournaments have H=189, M=27I. But are there non-VT maximizers?
         And do ALL H-maximizers have M scalar?

Also explore: eigenvalue structure of M for general tournaments.

kind-pasteur-2026-03-06-S25b (continuation)
"""

from itertools import permutations
import numpy as np
import random
import sys

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_diagonal(T, n, a):
    """Compute M[a,a] = sum_P (-1)^pos(a,P) over Hamiltonian paths P."""
    val = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            pos = list(perm).index(a)
            val += (-1)**pos
    return val

def compute_M_offdiag(T, n, a, b):
    """Compute M[a,b] for a != b using subset convolution."""
    U = [v for v in range(n) if v != a and v != b]
    val = 0
    for mask in range(1 << len(U)):
        S = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S)
        S_set = sorted(set(S) | {a})
        R_set = sorted(set(R) | {b})
        ea = 0
        if len(S_set) == 1:
            ea = 1
        else:
            for p in permutations(S_set):
                if p[-1] != a: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod
        bb2 = 0
        if len(R_set) == 1:
            bb2 = 1
        else:
            for p in permutations(R_set):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod
        val += sign * ea * bb2
    return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_diagonal(T, n, a)
        for b in range(a+1, n):
            M[a, b] = compute_M_offdiag(T, n, a, b)
            M[b, a] = M[a, b]  # M is symmetric (THM-030)
    return M

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def tournament_from_bits(n, bits):
    T = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[idx]:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1
            idx += 1
    return T

def score_sequence(T, n):
    return tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))


# ============================================================
# n=5: Exhaustive check of all tournaments
# ============================================================
print("=" * 70)
print("n=5: ALL tournaments — H vs det(M) vs eigenvalues")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
num_pairs = len(pairs)

H_to_tournaments = {}
for bits in range(1 << num_pairs):
    b = [(bits >> k) & 1 for k in range(num_pairs)]
    T = tournament_from_bits(n, b)
    H = count_H(T, n)
    if H not in H_to_tournaments:
        H_to_tournaments[H] = []
    if len(H_to_tournaments[H]) < 5:  # Keep at most 5 per H value
        H_to_tournaments[H].append(T)

print(f"\n  H values: {sorted(H_to_tournaments.keys())}")

for H_val in sorted(H_to_tournaments.keys()):
    tournaments = H_to_tournaments[H_val]
    print(f"\n  H={H_val}: testing {len(tournaments)} tournaments...")

    all_scalar = True
    for T in tournaments:
        M = compute_M(T, n)
        evals = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]
        is_scalar = np.allclose(M, M[0,0] * np.eye(n))
        scores = score_sequence(T, n)

        if not is_scalar:
            all_scalar = False
            det_M = int(round(np.linalg.det(M.astype(float))))
            print(f"    scores={scores}: M diag={[M[i,i] for i in range(n)]}, det={det_M}")
            print(f"      evals = {[round(e, 4) for e in evals]}")

    if all_scalar:
        M = compute_M(tournaments[0], n)
        print(f"    ALL scalar: M = {M[0,0]}*I (tested {len(tournaments)})")


# ============================================================
# n=7: Sample maximizers and near-maximizers
# ============================================================
print("\n" + "=" * 70)
print("n=7: H-maximizers — is M always scalar?")
print("=" * 70)

n = 7
pairs_7 = [(i,j) for i in range(n) for j in range(i+1, n)]
num_pairs_7 = len(pairs_7)

# Known VT: Paley QR7 = {1,2,4}
paley7 = make_circulant(7, {1,2,4})
H_paley = count_H(paley7, 7)
M_paley = compute_M(paley7, 7)
print(f"\n  Paley T_7: H={H_paley}")
print(f"    M = {M_paley[0,0]}*I: {np.array_equal(M_paley, M_paley[0,0]*np.eye(7, dtype=int))}")

# Other circulant tournaments on 7 vertices
print("\n  Other circulants on 7 vertices:")
from itertools import combinations as combs
for S_set in combs(range(1, 7), 3):
    # Check if S is a valid tournament generating set (S and 7-S partition {1..6})
    S_set = set(S_set)
    complement = {(7-s) % 7 for s in S_set}
    if S_set & complement:
        continue
    if S_set | complement != set(range(1,7)):
        continue

    T = make_circulant(7, S_set)
    H = count_H(T, 7)
    M = compute_M(T, 7)
    is_scalar = np.array_equal(M, M[0,0]*np.eye(7, dtype=int))
    det_M = int(round(np.linalg.det(M.astype(float))))
    print(f"    S={sorted(S_set)}: H={H}, M scalar={is_scalar}" + (f" ({M[0,0]}*I)" if is_scalar else f" det={det_M}"))

# Random search for H-maximizers at n=7
print("\n  Random search for n=7 maximizers (H=189):")
random.seed(2026)
max_H = 0
maximizer_count = 0
non_vt_maximizer = False

for trial in range(50000):
    T = {}
    for (i,j) in pairs_7:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    H = count_H(T, 7)
    if H > max_H:
        max_H = H
        print(f"    New max: H={H}")

    if H == 189:
        maximizer_count += 1
        scores = score_sequence(T, 7)
        if scores != (3,3,3,3,3,3,3):
            non_vt_maximizer = True
            M = compute_M(T, 7)
            is_scalar = np.array_equal(M, M[0,0]*np.eye(7, dtype=int))
            print(f"    NON-REGULAR maximizer! scores={scores}, M scalar={is_scalar}")
            if not is_scalar:
                print(f"      M diag = {[M[i,i] for i in range(7)]}")

        if maximizer_count >= 10:
            break

if maximizer_count > 0:
    print(f"    Found {maximizer_count} maximizers (H=189)")
    if not non_vt_maximizer:
        print(f"    All maximizers are regular (3,3,3,3,3,3,3)")

print(f"    Max H seen: {max_H}")


# ============================================================
# Eigenvalue patterns of M for various tournament families
# ============================================================
print("\n" + "=" * 70)
print("Eigenvalue structure of M")
print("=" * 70)

# At n=5, classify eigenvalue patterns
n = 5
eig_patterns = {}
for bits in range(1 << num_pairs):
    b = [(bits >> k) & 1 for k in range(num_pairs)]
    T = tournament_from_bits(n, b)
    M = compute_M(T, n)
    evals = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]
    evals_rounded = tuple(round(e, 3) for e in evals)

    H = count_H(T, n)
    scores = score_sequence(T, n)

    key = evals_rounded
    if key not in eig_patterns:
        eig_patterns[key] = []
    eig_patterns[key].append((H, scores))

print(f"\n  n=5: {len(eig_patterns)} distinct eigenvalue patterns")
for evals, entries in sorted(eig_patterns.items(), key=lambda x: -len(x[1])):
    H_vals = sorted(set(h for h, s in entries))
    score_vals = sorted(set(s for h, s in entries))
    det_val = 1
    for e in evals:
        det_val *= e
    print(f"    evals={evals}: count={len(entries)}, H={H_vals}, det={round(det_val)}")


# ============================================================
# Is det(M) always odd? (Comprehensive n=5 test)
# ============================================================
print("\n" + "=" * 70)
print("det(M) parity at n=5 (exhaustive)")
print("=" * 70)

n = 5
odd_count = 0
even_count = 0
det_values = set()

for bits in range(1 << num_pairs):
    b = [(bits >> k) & 1 for k in range(num_pairs)]
    T = tournament_from_bits(n, b)
    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    det_values.add(det_M)
    if det_M % 2 == 0:
        even_count += 1
    else:
        odd_count += 1

print(f"  Total: {odd_count + even_count} tournaments")
print(f"  Odd det: {odd_count}, Even det: {even_count}")
print(f"  det values: {sorted(det_values)}")

# Check: is det(M) always = H (mod 2)?
print("\n  Checking det(M) = H (mod 2):")
mismatch = 0
for bits in range(1 << num_pairs):
    b = [(bits >> k) & 1 for k in range(num_pairs)]
    T = tournament_from_bits(n, b)
    H = count_H(T, n)
    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    if det_M % 2 != H % 2:
        mismatch += 1
        print(f"    MISMATCH: H={H}, det={det_M}")
        if mismatch >= 5:
            break

if mismatch == 0:
    print("  det(M) = H (mod 2) for ALL n=5 tournaments!")
else:
    print(f"  {mismatch} mismatches found")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
