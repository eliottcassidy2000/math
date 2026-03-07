#!/usr/bin/env python3
"""
Connection between even-r polynomial and independence polynomial.

We know:
  H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...  (OCF)
  H(T) = tr(M) = tr(c_0) + tr(c_2)/4 + tr(c_4)/16 + ...    (even-r at r=1/2)

At n=5: tr(c_4) = 120 (universal), so tr(c_4)/16 = 7.5
  H = tr(c_0) + tr(c_2)/4 + 7.5
Also: H = 1 + 2*alpha_1 + 4*alpha_2

Therefore: tr(c_0) + tr(c_2)/4 = H - 7.5 = 2*alpha_1 + 4*alpha_2 - 6.5

At n=5: tr(c_2) = 12*t_3 - 30, and alpha_1 = t_3 (number of 3-cycles).
So tr(c_2)/4 = 3*t_3 - 7.5 = 3*alpha_1 - 7.5.
And tr(c_0) = H - 7.5 - 3*alpha_1 + 7.5 = H - 3*alpha_1
            = 1 + 2*alpha_1 + 4*alpha_2 - 3*alpha_1
            = 1 - alpha_1 + 4*alpha_2

FORMULA: tr(c_0) = 1 - alpha_1 + 4*alpha_2 at n=5

This connects the constant term of the even-r polynomial to the independence
polynomial of the odd-cycle graph!

Let's verify and look for patterns at n=7.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_trace_r(A, r_val):
    n = len(A)
    total = 0.0
    for a in range(n):
        U = [v for v in range(n) if v != a]
        val = 0.0
        for k in range(len(U)+1):
            for S in combinations(U, k):
                S_set = set(S)
                R = [v for v in U if v not in S_set]
                S_verts = sorted(list(S) + [a])
                R_verts = sorted(R + [a])
                ea = count_paths_weighted(A, S_verts, r_val, end=a)
                bb = count_paths_weighted(A, R_verts, r_val, start=a)
                val += ((-1)**k) * ea * bb
        total += val
    return total

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if (A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and
            A[perm[2]][perm[3]] and A[perm[3]][perm[4]] and A[perm[4]][perm[0]]):
            count += 1
    return count // 5

def independence_poly(A):
    """Compute independence polynomial of Omega(T).
    Odd cycles = independent sets in cycle conflict graph.
    I(x) = 1 + alpha_1*x + alpha_2*x^2 + ..."""
    n = len(A)

    # Find all odd directed cycles
    cycles = []

    # 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]:
                    cycles.append(frozenset([i,j,k]))
                if A[i][k]*A[k][j]*A[j][i]:
                    cycles.append(frozenset([i,j,k]))

    # 5-cycles (if n >= 5)
    if n >= 5:
        seen_5 = set()
        for perm in permutations(range(n), 5):
            if (A[perm[0]][perm[1]] and A[perm[1]][perm[2]] and
                A[perm[2]][perm[3]] and A[perm[3]][perm[4]] and A[perm[4]][perm[0]]):
                vs = frozenset(perm)
                if vs not in seen_5:
                    seen_5.add(vs)
                    cycles.append(vs)

    # 7-cycles (if n >= 7)
    if n >= 7:
        seen_7 = set()
        for perm in permutations(range(n), 7):
            if all(A[perm[i]][perm[(i+1)%7]] for i in range(7)):
                vs = frozenset(perm)
                if vs not in seen_7:
                    seen_7.add(vs)
                    cycles.append(vs)

    # Independence polynomial: count collections of vertex-disjoint cycles
    # alpha_k = number of collections of k pairwise vertex-disjoint odd cycles
    m = len(cycles)
    max_k = n // 3

    alphas = [0] * (max_k + 1)
    alphas[0] = 1

    # For small m, enumerate subsets
    for mask in range(1 << m):
        subset = [cycles[i] for i in range(m) if mask & (1 << i)]
        # Check pairwise vertex-disjoint
        all_verts = set()
        disjoint = True
        for c in subset:
            if c & all_verts:
                disjoint = False
                break
            all_verts |= c
        if disjoint and len(subset) > 0:
            alphas[len(subset)] += 1

    return alphas

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

# =====================================================================
print("=" * 70)
print("n=5: EVEN-r COEFFICIENTS vs INDEPENDENCE POLYNOMIAL")
print("=" * 70)

n = 5

# Get all iso classes
iso_data = {}
for bits in range(1024):
    A = [[0]*n for _ in range(n)]
    edge_idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> edge_idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            edge_idx += 1

    canon = tournament_canonical(A)
    if canon in iso_data:
        continue

    H = ham_path_count(A)
    t3 = count_3cycles(A)
    alphas = independence_poly(A)

    # Even-r polynomial: tr(M(r)) = tr(c_0) + tr(c_2)*u + tr(c_4)*u^2
    # where u = r^2. Sample at u = 0, 0.04, 0.16
    u_samples = [0.0, 0.04, 0.16]
    trace_vals = [transfer_trace_r(A, np.sqrt(u)) for u in u_samples]
    poly = np.polyfit(u_samples, trace_vals, 2)
    tr_c4 = poly[0]
    tr_c2 = poly[1]
    tr_c0 = poly[2]

    # OCF: H = 1 + 2*a1 + 4*a2
    # Expected: tr(c_0) = 1 - a1 + 4*a2
    a1 = alphas[1] if len(alphas) > 1 else 0
    a2 = alphas[2] if len(alphas) > 2 else 0

    expected_c0 = 1 - a1 + 4*a2
    scores = tuple(sorted(sum(row) for row in A))

    iso_data[canon] = {
        'H': H, 't3': t3, 'a1': a1, 'a2': a2,
        'tr_c0': tr_c0, 'tr_c2': tr_c2, 'tr_c4': tr_c4,
        'expected_c0': expected_c0, 'scores': scores
    }

print(f"\n  {'Scores':<20s} | H | a1 | a2 | tr(c0) | expected | tr(c2) | tr(c4)")
print("  " + "-" * 85)
for canon, d in sorted(iso_data.items(), key=lambda x: x[1]['H']):
    match = abs(d['tr_c0'] - d['expected_c0']) < 0.01
    mark = "OK" if match else "FAIL"
    print(f"  {str(d['scores']):<20s} | {d['H']:2d} | {d['a1']:2d} | {d['a2']:1d} | {d['tr_c0']:6.2f} | {d['expected_c0']:8.2f} | {d['tr_c2']:7.2f} | {d['tr_c4']:6.1f} | {mark}")

# =====================================================================
# Derive the formula algebraically
# =====================================================================
print("\n" + "=" * 70)
print("ALGEBRAIC DERIVATION AT n=5")
print("=" * 70)

print("""
  H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16       [even-r at r=1/2]
  H = 1 + 2*alpha_1 + 4*alpha_2                [OCF]
  tr(c_4) = 120 = 5!                           [universal]
  tr(c_2) = 12*alpha_1 - 30                    [verified at n=5]

  Substituting:
  1 + 2*a1 + 4*a2 = tr(c_0) + (12*a1 - 30)/4 + 120/16
                   = tr(c_0) + 3*a1 - 7.5 + 7.5
                   = tr(c_0) + 3*a1

  Therefore: tr(c_0) = 1 + 2*a1 + 4*a2 - 3*a1
                      = 1 - a1 + 4*a2

  FORMULA: tr(c_0) = 1 - alpha_1 + 4*alpha_2  at n=5

  This means:
    c_0 encodes the "excess" from disjoint pairs of 3-cycles
    c_2 encodes the linear 3-cycle count
    c_4 = (n-1)! * I is universal (encodes full permutation count)

  The decomposition H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16 is a
  REFINEMENT of the OCF decomposition H = 1 + 2*a1 + 4*a2!

  Specifically:
    tr(c_4)/16 = 120/16 = 7.5         [constant part]
    tr(c_2)/4  = 3*a1 - 7.5           [linear in a1]
    tr(c_0)    = 1 - a1 + 4*a2        [quadratic in a1, a2]

  Adding: 7.5 + 3*a1 - 7.5 + 1 - a1 + 4*a2 = 1 + 2*a1 + 4*a2 = H  CHECK!
""")

# =====================================================================
# What happens at n=7?
# =====================================================================
print("=" * 70)
print("n=7: EVEN-r vs INDEPENDENCE POLYNOMIAL (CIRCULANT)")
print("=" * 70)

n = 7

# Circulant tournaments
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

u_samples_7 = np.array([0.0, 0.04, 0.16, 0.36])
r_samples_7 = np.sqrt(u_samples_7)

for gs in sorted(gen_sets):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in gs:
                A[i][j] = 1

    H = ham_path_count(A)
    t3 = count_3cycles(A)
    t5 = count_5cycles(A)
    alphas = independence_poly(A)

    # Even-r: tr(M(r)) = tr(c_0) + tr(c_2)*u + tr(c_4)*u^2 + tr(c_6)*u^3
    trace_vals = [transfer_trace_r(A, rv) for rv in r_samples_7]
    poly = np.polyfit(u_samples_7, trace_vals, 3)
    tr_c6 = poly[0]
    tr_c4 = poly[1]
    tr_c2 = poly[2]
    tr_c0 = poly[3]

    a1 = alphas[1] if len(alphas) > 1 else 0
    a2 = alphas[2] if len(alphas) > 2 else 0

    # OCF: H = 1 + 2*a1 + 4*a2
    H_ocf = 1 + 2*a1 + 4*a2

    # Reconstruction: H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16 + tr(c_6)/64
    H_recon = tr_c0 + tr_c2/4 + tr_c4/16 + tr_c6/64

    print(f"\n  gen={sorted(gs)}: H={H}, a1={a1}, a2={a2}, H_ocf={H_ocf}")
    print(f"    tr(c_0)={tr_c0:.2f}, tr(c_2)={tr_c2:.2f}, tr(c_4)={tr_c4:.2f}, tr(c_6)={tr_c6:.2f}")
    print(f"    H_recon={H_recon:.2f}")
    print(f"    tr(c_6)/64 = {tr_c6/64:.4f}")
    print(f"    tr(c_4)/16 = {tr_c4/16:.4f}")
    print(f"    tr(c_2)/4  = {tr_c2/4:.4f}")
    print(f"    tr(c_0)    = {tr_c0:.4f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
