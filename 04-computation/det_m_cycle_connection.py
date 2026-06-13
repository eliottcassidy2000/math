#!/usr/bin/env python3
"""
Can det(M) be expressed in terms of the independence polynomial I(Omega, x)?

At n=5:
  H=15, c3=5: det=243=3^5, I=[1,5,5], I(2)=15
  H=15, c3=4: det=243=3^5, I=[1,4,2], I(2)=15
  H=13, c3=4: det=16=2^4, I=[1,4,2], I(2)=13  ← SAME I coeffs, DIFFERENT det!

Wait, can H=13 and H=15 have the same I coefficients? That would violate OCF!
H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2. If alpha=[4,2] then H=1+8+8=17. Hmm.

Let me actually COMPUTE I(Omega, x) for each isomorphism class.

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

def find_odd_cycles(T, n):
    """Find all directed odd cycles (as frozensets of vertices)."""
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

def independence_polynomial_coeffs(cycles):
    """Compute alpha_k = number of independent sets of size k."""
    cycle_list = list(cycles)
    m = len(cycle_list)
    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        subset = [cycle_list[k] for k in range(m) if mask & (1 << k)]
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
# n=5: Complete analysis with I(Omega, x) and det(M)
# ============================================================
print("=" * 70)
print("n=5: det(M) vs I(Omega, x) — exhaustive")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
num_pairs = len(pairs)

# Process one representative per isomorphism class
seen = set()
results = []

for bits in range(1 << num_pairs):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))

    cycles = find_odd_cycles(T, n)
    alpha = independence_polynomial_coeffs(cycles)
    alpha_key = tuple(alpha)
    I_at_2 = sum(a * 2**k for k, a in enumerate(alpha))

    key = (scores, H, alpha_key)
    if key in seen:
        continue
    seen.add(key)

    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))

    # Cycle details
    c3 = sum(1 for c in cycles if len(c) == 3)
    c5 = sum(1 for c in cycles if len(c) == 5)

    results.append({
        'scores': scores, 'H': H, 'c3': c3, 'c5': c5,
        'alpha': alpha, 'I2': I_at_2, 'det': det_M
    })

print(f"\n  {'scores':<20} {'H':>3} {'c3':>3} {'c5':>3} {'alpha':>15} {'I(2)':>5} {'det(M)':>8}")
print(f"  {'-----':<20} {'---':>3} {'---':>3} {'---':>3} {'-----':>15} {'----':>5} {'------':>8}")

for r in sorted(results, key=lambda x: x['H']):
    print(f"  {str(r['scores']):<20} {r['H']:>3} {r['c3']:>3} {r['c5']:>3} {str(r['alpha']):>15} {r['I2']:>5} {r['det']:>8}")
    assert r['I2'] == r['H'], f"OCF VIOLATION: I(2)={r['I2']} != H={r['H']}"

print("\n  OCF verified for all isomorphism classes!")

# ============================================================
# Is det(M) determined by alpha (= I coefficients)?
# ============================================================
print("\n" + "=" * 70)
print("Is det(M) determined by I(Omega, x)?")
print("=" * 70)

alpha_to_det = defaultdict(set)
for r in results:
    alpha_to_det[tuple(r['alpha'])].add(r['det'])

for alpha, dets in sorted(alpha_to_det.items()):
    H = sum(a * 2**k for k, a in enumerate(alpha))
    if len(dets) > 1:
        print(f"  alpha={list(alpha)}: H={H}, det = {sorted(dets)} — NOT determined!")
    else:
        print(f"  alpha={list(alpha)}: H={H}, det = {sorted(dets)[0]}")


# ============================================================
# Try: is det(M) related to I(Omega, x) evaluated at special points?
# ============================================================
print("\n" + "=" * 70)
print("det(M) vs I(Omega, x) at various x")
print("=" * 70)

for r in sorted(results, key=lambda x: x['H']):
    alpha = r['alpha']
    det_M = r['det']

    I_vals = {}
    for x in [-1, 0, 1, 2, 3, -2, -3]:
        I_vals[x] = sum(a * x**k for k, a in enumerate(alpha))

    # Check: is det(M) = I(Omega, x) for some x?
    matches = [x for x, v in I_vals.items() if v == det_M]

    # Check: is det(M) = I(Omega, x)^k for some x, k?
    print(f"  H={r['H']}, alpha={r['alpha']}, det={det_M}")
    print(f"    I(-1)={I_vals[-1]}, I(0)={I_vals[0]}, I(1)={I_vals[1]}, I(2)={I_vals[2]}, I(3)={I_vals[3]}")

    if matches:
        print(f"    det = I(x) at x = {matches}")


# ============================================================
# Direct formula attempt: det(M) from (H, c3, c5, alpha)
# ============================================================
print("\n" + "=" * 70)
print("Pattern search: det(M) formula")
print("=" * 70)

for r in sorted(results, key=lambda x: x['H']):
    H = r['H']
    c3 = r['c3']
    c5 = r['c5']
    det_M = r['det']
    alpha = r['alpha']

    # Try various formulas
    formulas = {
        'H^2 - n*H + 1': H**2 - 5*H + 1,
        '(-1)^c3 * H': (-1)**c3 * H,
        '(-1)^c3 * 3^(c3-1) * H/5': (-1)**c3 * 3**(c3-1) * H // 5 if c3 > 0 else None,
        'H*(H-8)/n + 8': H*(H-8)//5 + 8 if (H*(H-8)) % 5 == 0 else None,
        'prod eigenvalues': None  # Already known
    }

    print(f"  H={H:>3}, c3={c3}, det={det_M:>8}", end="")

    # Check H^n / n^n
    if H % n == 0:
        ratio = (H // n) ** n
        if ratio == det_M:
            print(f"  = (H/n)^n", end="")

    print()


# ============================================================
# DEEPER: characteristic polynomial of M
# ============================================================
print("\n" + "=" * 70)
print("Characteristic polynomial of M")
print("=" * 70)

for r in sorted(results, key=lambda x: x['H']):
    H = r['H']
    T = tournament_from_bits(n, None)  # need to find representative...

# Let me just compute for specific representatives
for bits in [0, 1, 100, 200, 500, 700, 900, 1000, 1023]:
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    M = compute_M(T, n)
    det_M = int(round(np.linalg.det(M.astype(float))))
    evals = sorted(np.linalg.eigvalsh(M.astype(float)))[::-1]

    # Characteristic polynomial coefficients
    # p(x) = det(xI - M) = x^5 - tr*x^4 + ... - det
    # Using Newton's identities
    s1 = sum(evals)
    s2 = sum(e**2 for e in evals)
    s3 = sum(e**3 for e in evals)
    s4 = sum(e**4 for e in evals)
    s5 = sum(e**5 for e in evals)

    # p1 = s1 = tr(M) = H
    p1 = round(s1)
    # p2 = (s1*p1 - s2)/2
    p2 = round((s1*p1 - s2) / 2)
    # p3 = (s1*p2 - s2*p1 + s3) / 3  (not quite right... use numpy)

    char_coeffs = np.round(np.poly(evals)).astype(int)
    scores = tuple(sorted(sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)))
    print(f"  bits={bits:>4}, H={H:>3}, scores={scores}, char poly = {list(char_coeffs)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
