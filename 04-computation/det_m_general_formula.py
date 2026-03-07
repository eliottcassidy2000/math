#!/usr/bin/env python3
"""
det(M) for general tournaments: look for patterns.

Known:
  VT odd n: det(M) = (H/n)^n
  Transitive: det(M) = (-1)^{n(n-1)/2} F(n+1)

Questions:
  1. Is det(M) determined by H? (No — H=5 gives det=9 and det=17)
  2. Is det(M) determined by the isomorphism class? (Yes at n=5)
  3. Can we relate det(M) to cycle structure?

Also: look at det(M) mod small primes.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np
from collections import defaultdict

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        # Diagonal
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        M[a, a] = val
        # Off-diagonal
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
            M[a, b] = val
            M[b, a] = val
    return M

def count_3cycles(T, n):
    count = 0
    from itertools import combinations
    for a, b, c in combinations(range(n), 3):
        if T.get((a,b),0) and T.get((b,c),0) and T.get((c,a),0):
            count += 1
        if T.get((a,c),0) and T.get((c,b),0) and T.get((b,a),0):
            count += 1
    return count

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


# ============================================================
# n=5: EXHAUSTIVE det(M) analysis with cycle counts
# ============================================================
print("=" * 70)
print("n=5: det(M) vs H vs c3")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

# Classify by isomorphism class (approximated by (scores, H, c3))
classes = {}

for bits in range(1 << len(pairs)):
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = tournament_from_bits(n, b_list)

    H = count_H(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    c3 = count_3cycles(T, n)

    key = (scores, H, c3)
    if key in classes:
        classes[key]['count'] += 1
        continue

    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    tr_M = int(np.trace(M))

    classes[key] = {
        'count': 1, 'H': H, 'scores': scores, 'c3': c3,
        'det': det_M, 'tr': tr_M
    }

print(f"\n  {len(classes)} isomorphism classes")
print(f"\n  {'scores':<20} {'H':>3} {'c3':>3} {'det(M)':>8} {'tr(M)':>6} {'count':>6}")
print(f"  {'-----':<20} {'---':>3} {'---':>3} {'------':>8} {'----':>6} {'-----':>6}")

for key in sorted(classes.keys(), key=lambda x: (x[1], x[0])):
    c = classes[key]
    print(f"  {str(c['scores']):<20} {c['H']:>3} {c['c3']:>3} {c['det']:>8} {c['tr']:>6} {c['count']:>6}")


# ============================================================
# Patterns in det(M)
# ============================================================
print("\n" + "=" * 70)
print("det(M) patterns")
print("=" * 70)

det_values = sorted(set(c['det'] for c in classes.values()))
print(f"\n  det values: {det_values}")

# Factor each det
print("\n  Factorizations:")
for d in det_values:
    if d == 0:
        print(f"    {d} = 0")
        continue
    abs_d = abs(d)
    factors = []
    temp = abs_d
    for p in [2, 3, 5, 7, 11, 13]:
        while temp % p == 0:
            factors.append(p)
            temp //= p
    if temp > 1:
        factors.append(temp)
    sign = "-" if d < 0 else ""
    print(f"    {d:>8} = {sign}{' * '.join(str(f) for f in factors) if factors else '1'}")

# Is det(M) determined by (H, c3)?
print("\n  det(M) determined by (H, c3)?")
hc3_to_det = defaultdict(set)
for key, c in classes.items():
    hc3_to_det[(c['H'], c['c3'])].add(c['det'])

for (h, c3), dets in sorted(hc3_to_det.items()):
    if len(dets) > 1:
        print(f"    (H={h}, c3={c3}): det = {sorted(dets)} -- NOT determined!")
    else:
        print(f"    (H={h}, c3={c3}): det = {sorted(dets)[0]}")


# ============================================================
# Relation: det(M) and independence polynomial I(Omega, x)
# ============================================================
print("\n" + "=" * 70)
print("det(M) vs I(Omega, x) coefficients")
print("=" * 70)

from itertools import combinations

def compute_odd_cycles(T, n):
    """Return list of vertex sets of directed odd cycles."""
    cycles = set()
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if perm[0] != min(verts): continue
                is_cycle = True
                for k in range(length):
                    if T.get((perm[k], perm[(k+1) % length]), 0) != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(verts))
                    break
    return cycles

def independence_polynomial(cycles, n):
    """Compute I(Omega, x) where Omega has cycles as vertices."""
    cycle_list = list(cycles)
    m = len(cycle_list)
    # Find independent sets
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        subset = [cycle_list[k] for k in range(m) if mask & (1 << k)]
        # Check pairwise vertex-disjoint
        independent = True
        vertices = set()
        for c in subset:
            if vertices & c:
                independent = False
                break
            vertices |= c
        if independent:
            coeffs[len(subset)] += 1
    return coeffs

for key in sorted(classes.keys(), key=lambda x: (x[1], x[0])):
    c = classes[key]
    T = tournament_from_bits(n, None)  # Need to reconstruct...

# Actually, let me just compute for representatives
print("\n  Computing for representative tournaments...")

for bits in range(0, 1 << len(pairs), 37):  # Sample
    b_list = [(bits >> k) & 1 for k in range(len(pairs))]
    T = tournament_from_bits(n, b_list)
    H = count_H(T, n)
    c3 = count_3cycles(T, n)
    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))

    cycles = compute_odd_cycles(T, n)
    coeffs = independence_polynomial(cycles, n)
    I_at_2 = sum(c * 2**k for k, c in enumerate(coeffs))

    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    print(f"  scores={scores}: H={H}, c3={c3}, I(Omega,2)={I_at_2}, det(M)={det_M}")
    print(f"    I coeffs = {coeffs}")

    if I_at_2 != H:
        print(f"    WARNING: I(Omega,2) != H!")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
