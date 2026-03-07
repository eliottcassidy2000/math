#!/usr/bin/env python3
"""
WHY does W(r) = tr(M(r)) at odd n but NOT at even n?

W(r) = sum over all n! permutations P of prod_{edges e in P} (r + s_e)
tr(M(r)) = sum_a M[a,a](r) via IE decomposition

At r=1/2:
  W(1/2) = H (number of directed Ham paths)
  tr(M(1/2)) = sum_a inshat(a) or similar via IE

But W(1/2) should always equal H, and H = tr(M(1/2))... right?

Actually at n=4, the test showed tr(M(r)) = 0 for ALL r. This means
sum_a M[a,a](r) = 0 for all r, including r=1/2. But H > 0 at n=4!

So EITHER:
  (a) tr(M(1/2)) ≠ H at even n, or
  (b) the IE formula is computing something different from what I think

Let me investigate by computing M[a,a](r) individually for each vertex a.

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
import numpy as np

def _count_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_entry_at_r(A, a, r_val):
    """Compute M[a,a](r) via IE formula."""
    n = len(A)
    U = [v for v in range(n) if v != a]
    total = 0.0
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            S_verts = sorted(list(S) + [a])
            R_verts = sorted(R + [a])
            ea = _count_weighted(A, S_verts, r_val, end=a)
            ba = _count_weighted(A, R_verts, r_val, start=a)
            total += ((-1)**k) * ea * ba
    return total

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

# =====================================================================
import random
print("=" * 70)
print("W(r) vs tr(M(r)): DETAILED ANALYSIS")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n--- n={n} ---")
    random.seed(n * 100)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count(A)
    print(f"  H = {H}")
    print(f"  Scores: {[sum(A[i]) for i in range(n)]}")

    # Compute M[a,a](r) for each vertex a at several r values
    r_vals = [0.0, 0.25, 0.5, 1.0]
    for r in r_vals:
        entries = [transfer_entry_at_r(A, a, r) for a in range(n)]
        tr_M = sum(entries)

        # W(r)
        W = 0.0
        for p in permutations(range(n)):
            w = 1.0
            for i in range(n-1):
                w *= r + (A[p[i]][p[i+1]] - 0.5)
            W += w

        print(f"  r={r:.2f}: M[a,a]={[f'{e:.2f}' for e in entries]}, "
              f"tr(M)={tr_M:.4f}, W={W:.4f}, diff={W-tr_M:.4f}")

# =====================================================================
# KEY: Check what M[a,a] equals at r=1/2
# =====================================================================
print("\n" + "=" * 70)
print("WHAT IS M[a,a](1/2)?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n  n={n}:")
    random.seed(n * 100)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count(A)

    for a in range(n):
        M_aa = transfer_entry_at_r(A, a, 0.5)

        # Count: how many Ham paths end at a? How many start at a?
        end_count = sum(1 for p in permutations(range(n))
                       if p[-1] == a and all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))
        start_count = sum(1 for p in permutations(range(n))
                         if p[0] == a and all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

        # inshat(a) from the paper
        # At r=1/2, IE formula: M[a,a] = sum_S (-1)^|S| * E_a(S∪{a}) * B_a(R∪{a})
        # E_a(S∪{a}) = number of directed paths ending at a through S∪{a}
        # B_a(R∪{a}) = number of directed paths starting at a through R∪{a}

        print(f"    a={a}: M[{a},{a}](1/2)={M_aa:.4f}, "
              f"#paths ending at {a}={end_count}, starting at {a}={start_count}")

    print(f"    Sum of M[a,a](1/2) = {sum(transfer_entry_at_r(A, a, 0.5) for a in range(n)):.4f}")
    print(f"    H = {H}")

# =====================================================================
# Check: does the full transfer matrix M[a,b](1/2) recover H?
# =====================================================================
print("\n" + "=" * 70)
print("FULL TRANSFER MATRIX M[a,b](1/2)")
print("=" * 70)

for n in [3, 4]:
    print(f"\n  n={n}:")
    random.seed(n * 100)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count(A)

    # Compute full M[a,b] at r=1/2
    M = np.zeros((n,n))
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            if a == b:
                U = [v for v in range(n) if v != a]
            total = 0.0
            # IE formula for M[a,b]
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    if a == b:
                        S_verts = sorted(list(S) + [a])
                        R_verts = sorted(R + [a])
                    else:
                        S_verts = sorted(list(S) + [a])
                        R_verts = sorted(R + [b])
                    ea = _count_weighted(A, S_verts, 0.5, end=a)
                    ba_val = _count_weighted(A, R_verts, 0.5, start=b if a != b else a)
                    total += ((-1)**k) * ea * ba_val
            M[a][b] = total

    print(f"    M[a,b] at r=1/2:")
    for row in M:
        print(f"      {[f'{x:.2f}' for x in row]}")
    print(f"    tr(M) = {np.trace(M):.4f}")
    print(f"    sum(M) = {M.sum():.4f}")
    print(f"    H = {H}")

    # What about sum of entries?
    # H = sum_{a,b} M[a,b]?
    # Or H = sum_a M[a,a]?
    # Or something else?

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
