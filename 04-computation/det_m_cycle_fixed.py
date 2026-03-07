#!/usr/bin/env python3
"""
Fixed cycle counting: count ALL directed odd cycles, not just one per vertex set.

A vertex set of size k can support MULTIPLE directed k-cycles.
Each directed cycle is a separate vertex of Omega(T).

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
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

def find_all_directed_odd_cycles(T, n):
    """Find ALL directed odd cycles. Each cycle = frozenset of (vertex, position) pairs.
    We represent cycles as (vertex_set, canonical_ordering) pairs.
    Two cycles are the same if they traverse the same edges.
    """
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            # Find ALL directed cycles on this vertex set
            min_v = min(verts)
            for perm in permutations(verts):
                if perm[0] != min_v:
                    continue  # canonical: start at smallest vertex
                # Also require perm[1] < perm[-1] to avoid reverse counting
                if perm[1] > perm[-1]:
                    continue
                is_cycle = True
                for k in range(length):
                    if T.get((perm[k], perm[(k+1) % length]), 0) != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(frozenset(verts))
    return cycles

def independence_polynomial(cycles):
    """Compute I(Omega, x) coefficients."""
    m = len(cycles)
    if m == 0:
        return [1]
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        subset = [cycles[k] for k in range(m) if mask & (1 << k)]
        independent = True
        vertices = set()
        for c in subset:
            if vertices & c:
                independent = False
                break
            vertices |= c
        if independent:
            coeffs[len(subset)] += 1
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def tournament_from_bits(n, bits_int):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits_int >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# n=5: Fixed cycle counting
# ============================================================
print("=" * 70)
print("n=5: det(M) vs I(Omega, x) with CORRECT cycle counting")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
num_pairs = len(pairs)

seen = set()
results = []

for bits in range(1 << num_pairs):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    cycles = find_all_directed_odd_cycles(T, n)
    alpha = independence_polynomial(cycles)
    I_at_2 = sum(a * 2**k for k, a in enumerate(alpha))

    c3 = sum(1 for c in cycles if len(c) == 3)
    c5 = sum(1 for c in cycles if len(c) == 5)

    key = (scores, H, c3, c5)
    if key in seen:
        continue
    seen.add(key)

    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))

    results.append({
        'scores': scores, 'H': H, 'c3': c3, 'c5': c5,
        'alpha': alpha, 'I2': I_at_2, 'det': det_M,
        'ncycles': len(cycles)
    })

print(f"\n  {'scores':<20} {'H':>3} {'c3':>3} {'c5':>3} {'#cyc':>4} {'alpha':>15} {'I(2)':>5} {'det':>8} {'OCF':>4}")
print(f"  {'-----':<20} {'---':>3} {'---':>3} {'---':>3} {'----':>4} {'-----':>15} {'----':>5} {'---':>8} {'---':>4}")

for r in sorted(results, key=lambda x: x['H']):
    ocf_ok = "OK" if r['I2'] == r['H'] else "FAIL"
    print(f"  {str(r['scores']):<20} {r['H']:>3} {r['c3']:>3} {r['c5']:>3} {r['ncycles']:>4} {str(r['alpha']):>15} {r['I2']:>5} {r['det']:>8} {ocf_ok:>4}")


# ============================================================
# det(M) vs I coefficients
# ============================================================
print("\n" + "=" * 70)
print("Is det(M) determined by alpha?")
print("=" * 70)

alpha_to_det = defaultdict(set)
alpha_to_H = defaultdict(set)
for r in results:
    alpha_to_det[tuple(r['alpha'])].add(r['det'])
    alpha_to_H[tuple(r['alpha'])].add(r['H'])

for alpha, dets in sorted(alpha_to_det.items()):
    H_vals = sorted(alpha_to_H[alpha])
    if len(dets) > 1:
        print(f"  alpha={list(alpha)}: H={H_vals}, det = {sorted(dets)} — NOT determined!")
    else:
        print(f"  alpha={list(alpha)}: H={H_vals}, det = {sorted(dets)[0]}")


# ============================================================
# Compute I(Omega, x) at various x values
# ============================================================
print("\n" + "=" * 70)
print("det(M) vs I(Omega, x) at special x")
print("=" * 70)

for r in sorted(results, key=lambda x: x['H']):
    alpha = r['alpha']
    det_M = r['det']
    H = r['H']

    I_vals = {}
    for x_val in range(-5, 6):
        I_vals[x_val] = sum(a * x_val**k for k, a in enumerate(alpha))

    # Check: is det(M) = product of I at roots?
    # Or: is det(M) = I(x) for some x?
    matches = [x for x, v in I_vals.items() if v == det_M]

    # I(Omega, sqrt(2))^n? No, need integer...
    # Try: det(M) = I(Omega, -1)^? or det(M) * something = ...

    print(f"  H={H:>3}, alpha={r['alpha']}, det={det_M:>8}")
    if matches:
        print(f"    det = I({matches[0]})")

    # Try: is det(M) = prod over cycles C of (-1)^{|C|/2} ?
    # Or: det(M) = (-1)^{c5} * ...?


# ============================================================
# Try: det(M) = resultant or discriminant of I(Omega, x)?
# ============================================================
print("\n" + "=" * 70)
print("Algebraic invariants of I(Omega, x)")
print("=" * 70)

for r in sorted(results, key=lambda x: x['H']):
    alpha = r['alpha']
    det_M = r['det']

    # I(Omega, x) as polynomial
    degree = len(alpha) - 1
    print(f"  H={r['H']:>3}: I(x) = {' + '.join(f'{a}*x^{k}' for k, a in enumerate(alpha) if a)}")
    print(f"    deg={degree}, det(M)={det_M}")

    if degree >= 1:
        # Discriminant of I(Omega, x) (for quadratics)
        if degree == 1:
            # I = a0 + a1*x, discriminant = a1^2 (trivially)
            disc = alpha[1]**2
        elif degree == 2:
            # I = a0 + a1*x + a2*x^2
            disc = alpha[1]**2 - 4*alpha[0]*alpha[2]
        else:
            disc = None

        if disc is not None:
            print(f"    disc(I) = {disc}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
