#!/usr/bin/env python3
"""
Connection between even-r polynomial and OCF independence polynomial.

DERIVED FORMULA at n=5:
  tr(c_0) = 1 - alpha_1 + 4*alpha_2
  tr(c_2) = 12*alpha_1 - 30
  tr(c_4) = 120 = 5!

where alpha_k = #collections of k vertex-disjoint odd cycles.
At n=5, alpha_1 = #(3-cycles) = t_3 since 5-cycles use all vertices
and can't be in a collection with anything else.

Verification at n=5 using known iso class data, then extension to n=7.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
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
    """Compute tr(M(r)) efficiently — only diagonal entries."""
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

def count_disjoint_3cycle_pairs(A):
    """Count pairs of vertex-disjoint 3-cycles."""
    n = len(A)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]:
                    cycles.append(frozenset([i,j,k]))
                if A[i][k]*A[k][j]*A[j][i]:
                    cycles.append(frozenset([i,j,k]))
    # Count disjoint pairs
    count = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                count += 1
    return count

def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            count += 1
    return count // 5

# =====================================================================
print("=" * 70)
print("n=5: VERIFY tr(c_0) = 1 - alpha_1 + 4*alpha_2")
print("=" * 70)

# Representative tournaments from each iso class at n=5
# (using known score sequences and specific constructions)
n = 5

# Build a set of non-isomorphic representatives
reps = []

# Score (0,1,2,3,4) - transitive tournament
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(i+1, 5):
        A[i][j] = 1
reps.append(("transitive", A))

# Score (2,2,2,2,2) - Paley: gen={1,2} on Z/5
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i)%5 in {1,2}:
            A[i][j] = 1
reps.append(("Paley", A))

# Score (2,2,2,2,2) - C5: gen={1,3}
A = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i)%5 in {1,3}:
            A[i][j] = 1
reps.append(("C5_circulant", A))

# Score (1,2,2,2,3) - Non-VT PU (cone over C3)
reps.append(("cone_C3", [
    [0, 0, 0, 0, 1],
    [1, 0, 0, 1, 0],
    [1, 1, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [0, 1, 1, 1, 0],
]))

# Other specific tournaments for coverage
# Score (0,2,2,2,4)
reps.append(("score_02224", [
    [0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0],
    [1, 0, 0, 1, 0],
    [1, 1, 0, 0, 0],
    [1, 1, 1, 1, 0],
]))

# Score (1,1,2,3,3) - two variants
reps.append(("score_11233a", [
    [0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0],
    [1, 1, 0, 0, 0],
    [0, 1, 1, 0, 1],
    [1, 1, 1, 0, 0],
]))

u_samples = np.array([0.0, 0.04, 0.16])
r_samples = np.sqrt(u_samples)

print(f"\n  {'Name':<16s} | H | a1 | a2 | tr(c0) | 1-a1+4a2 | tr(c2) | 12a1-30 | tr(c4)")
print("  " + "-" * 95)

for name, A in reps:
    H = ham_path_count(A)
    t3 = count_3cycles(A)
    a2 = count_disjoint_3cycle_pairs(A)

    # 5-cycles contribute to alpha_1 but not alpha_2 at n=5
    t5 = count_5cycles(A)
    a1 = t3 + t5  # alpha_1 = total odd cycles

    # Even-r trace polynomial
    trace_vals = [transfer_trace_r(A, rv) for rv in r_samples]
    poly = np.polyfit(u_samples, trace_vals, 2)
    tr_c4 = poly[0]
    tr_c2 = poly[1]
    tr_c0 = poly[2]

    expected_c0 = 1 - a1 + 4*a2
    expected_c2 = 12*a1 - 30

    c0_match = abs(tr_c0 - expected_c0) < 0.5
    c2_match = abs(tr_c2 - expected_c2) < 0.5

    print(f"  {name:<16s} | {H:2d} | {a1:2d} | {a2:1d} | {tr_c0:6.1f} | {expected_c0:8.1f} | {tr_c2:6.1f} | {expected_c2:7.1f} | {tr_c4:5.0f} | {'OK' if c0_match and c2_match else 'FAIL'}")

# =====================================================================
# n=7 circulant: what are the formulas?
# =====================================================================
print("\n" + "=" * 70)
print("n=7: EVEN-r vs OCF FOR CIRCULANT TOURNAMENTS")
print("=" * 70)

n = 7
half = list(range(1, (n+1)//2))
gen_sets_7 = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets_7.add(frozenset(gs))

u_samples_7 = np.array([0.0, 0.04, 0.16, 0.36])
r_samples_7 = np.sqrt(u_samples_7)

print(f"\n  {'gen':<12s} | H | a1 | a2 | tr(c0) | tr(c2) | tr(c4) | tr(c6)")
print("  " + "-" * 80)

n7_data = []
for gs in sorted(gen_sets_7):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in gs:
                A[i][j] = 1

    H = ham_path_count(A)
    t3 = count_3cycles(A)
    t5 = count_5cycles(A)

    # 7-cycles
    t7 = 0
    for perm in permutations(range(7)):
        if all(A[perm[i]][perm[(i+1)%7]] for i in range(7)):
            t7 += 1
    t7 //= 7

    a1 = t3 + t5 + t7

    # Disjoint pairs (only 3+3 possible at n=7, or 3+5... wait, 3+5=8>7, impossible)
    # So alpha_2 = #disjoint pairs of 3-cycles (needing 6 vertices out of 7)
    a2 = count_disjoint_3cycle_pairs(A)

    trace_vals = [transfer_trace_r(A, rv) for rv in r_samples_7]
    poly = np.polyfit(u_samples_7, trace_vals, 3)
    tr_c6 = poly[0]
    tr_c4 = poly[1]
    tr_c2 = poly[2]
    tr_c0 = poly[3]

    H_recon = tr_c0 + tr_c2/4 + tr_c4/16 + tr_c6/64

    print(f"  {str(sorted(gs)):<12s} | {H:3d} | {a1:2d} | {a2:1d} | {tr_c0:6.2f} | {tr_c2:7.2f} | {tr_c4:7.1f} | {tr_c6:7.1f}")

    n7_data.append({
        'gs': sorted(gs), 'H': H, 't3': t3, 't5': t5, 't7': t7,
        'a1': a1, 'a2': a2,
        'tr_c0': tr_c0, 'tr_c2': tr_c2, 'tr_c4': tr_c4, 'tr_c6': tr_c6
    })

# Try to find formulas
print("\n  Looking for pattern: tr(c_k) = f(a1, a2, ...)")
for d in n7_data:
    # Expected: tr(c_6) = 5040 = 7!
    # tr(c_4) = ?
    # tr(c_2) = ?
    # tr(c_0) = ?
    expected_c0 = d['H'] - d['tr_c2']/4 - d['tr_c4']/16 - d['tr_c6']/64
    print(f"  gen={d['gs']}: H={d['H']}, a1={d['a1']}, a2={d['a2']}, "
          f"t3={d['t3']}, t5={d['t5']}, t7={d['t7']}, tr_c0={d['tr_c0']:.2f}")

# Check: tr(c_2) = f(a1)?
print("\n  tr(c_2) vs a1:")
for d in n7_data:
    ratio = d['tr_c2'] / d['a1'] if d['a1'] > 0 else float('inf')
    print(f"    a1={d['a1']}, tr_c2={d['tr_c2']:.2f}, ratio={ratio:.4f}")

# Check: tr(c_2) = c*t3 + d*t5 + e?
# Build linear system
print("\n  Fitting tr(c_2) = a*t3 + b*t5 + c*t7 + d:")
X = np.array([[d['t3'], d['t5'], d['t7'], 1] for d in n7_data])
y = np.array([d['tr_c2'] for d in n7_data])
coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
print(f"    tr(c_2) = {coeffs[0]:.4f}*t3 + {coeffs[1]:.4f}*t5 + {coeffs[2]:.4f}*t7 + {coeffs[3]:.4f}")
print(f"    Residuals: {residuals}")

# Check: tr(c_4) = ?
print("\n  Fitting tr(c_4) = a*t3 + b*t5 + c*t7 + d:")
y4 = np.array([d['tr_c4'] for d in n7_data])
coeffs4, res4, _, _ = np.linalg.lstsq(X, y4, rcond=None)
print(f"    tr(c_4) = {coeffs4[0]:.4f}*t3 + {coeffs4[1]:.4f}*t5 + {coeffs4[2]:.4f}*t7 + {coeffs4[3]:.4f}")
print(f"    Residuals: {res4}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
