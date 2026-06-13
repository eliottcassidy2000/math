#!/usr/bin/env python3
"""
CLOSED-FORM moment formulas at general n — derivation and verification.

PROVED so far:
  m2 = n!*(3n^2-5n+4)/12 + 4*(n-2)!*t3

From v1 output, m3_t3/(n-2)! = 12*(n-1)/2 per step... actually 12,24,36,48,60,72
at n=3,5,7,9,11,13 → these are 12*(n-1)/2 = 6*(n-1). So m3_t3 = 6*(n-1)*(n-2)! = 6*(n-1)!.

m3_const/n!: 2, 11, 33, 74, 140, 237. Third differences are 6 → cubic in n.
Fitting: (n-1)(n^2-n+2)/8.

So m3 = n!*(n-1)(n^2-n+2)/8 + 6*(n-1)!*t3.

Now derive m4 for general n including t5 and bc contributions.

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
from collections import Counter
import random

def count_position_patterns(n, k):
    if k == 0: return {(): 1}
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps = []
        comp = [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        pat = tuple(sorted(comps, reverse=True))
        patterns[pat] += 1
    return dict(patterns)

def sigma_extended(pattern, n):
    """Return (const, t3_coeff, t5_coeff, bc_coeff) for sigma(pattern) at general n."""
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0:
        return (0, 0, 0, 0)

    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)

    def pair_product(start_remaining):
        p = 1
        r = start_remaining
        for i in range(num_size1):
            p *= comb(r, 2)
            r -= 2
        return p

    F = factorial(free)

    if len(big_sizes) == 0:
        k = len(sizes)
        return (factorial(n) // (2**k), 0, 0, 0)

    if big_sizes == [2]:
        pp = pair_product(n - 3)
        return (F * comb(n, 3) * pp, F * 2 * pp, 0, 0)

    if big_sizes == [3]:
        pp = pair_product(n - 4)
        return (F * comb(n, 4) * pp, F * 2 * (n-3) * pp, 0, 0)

    if big_sizes == [4]:
        pp = pair_product(n - 5)
        return (F * comb(n, 5) * pp, F * 2 * comb(n-3, 2) * pp, F * 2 * pp, 0)

    if big_sizes == [2, 2]:
        pp = pair_product(n - 6)
        # sum_{ordered (G1,G2) disjoint 3-sets} (1+2*cyc1)(1+2*cyc2)
        # = C(n,3)*C(n-3,3) + 4*C(n-3,3)*t3 + 8*bc
        base_const = comb(n, 3) * comb(n-3, 3)
        base_t3 = 4 * comb(n-3, 3)  # 2*C(n-3,3) from cyc1 + 2*C(n-3,3) from cyc2
        base_bc = 8  # 4 * 2*bc (ordered pairs of disjoint cyclic triples)
        return (F * base_const * pp, F * base_t3 * pp, 0, F * base_bc * pp)

    if big_sizes == [5]:
        pp = pair_product(n - 6)
        # H(6-set) via OCF at n=6: 1 + 2*c3 + 2*c5 + 4*bc
        # sum_{6-subsets} = C(n,6) + 2*C(n-3,3)*t3 + 2*(n-5)*t5 + 4*bc
        return (F * comb(n, 6) * pp,
                F * 2 * comb(n-3, 3) * pp,
                F * 2 * (n-5) * pp,
                F * 4 * pp)

    return None

# =====================================================================
# VERIFY m3 CLOSED FORM
# =====================================================================
print("=" * 70)
print("m3 = n!*(n-1)(n^2-n+2)/8 + 6*(n-1)!*t3")
print("=" * 70)

for n in [3, 5, 7, 9, 11, 13, 15]:
    SIGMA_1 = (n-1) * factorial(n) // 2

    SIGMA_2_c = SIGMA_2_t3 = 0
    for pat, cnt in count_position_patterns(n, 2).items():
        r = sigma_extended(pat, n)
        SIGMA_2_c += cnt * r[0]
        SIGMA_2_t3 += cnt * r[1]

    SIGMA_3_c = SIGMA_3_t3 = 0
    for pat, cnt in count_position_patterns(n, 3).items():
        r = sigma_extended(pat, n)
        if r is None: continue
        SIGMA_3_c += cnt * r[0]
        SIGMA_3_t3 += cnt * r[1]

    m3_c = SIGMA_1 + 6*SIGMA_2_c + 6*SIGMA_3_c
    m3_t3 = 6*SIGMA_2_t3 + 6*SIGMA_3_t3

    pred_c = factorial(n) * (n-1) * (n*n - n + 2) // 8
    pred_t3 = 6 * factorial(n-1)

    ok = (m3_c == pred_c and m3_t3 == pred_t3)
    print(f"  n={n:2d}: m3 = {m3_c} + {m3_t3}*t3  closed: {pred_c} + {pred_t3}*t3  {'OK' if ok else 'FAIL'}")

# =====================================================================
# m4 COMPLETE FORMULA at general n
# =====================================================================
print(f"\n{'='*70}")
print("m4 = SIGMA_1 + 14*SIGMA_2 + 36*SIGMA_3 + 24*SIGMA_4")
print("  (S(4,1)=1, S(4,2)=7, S(4,3)=6, S(4,4)=1)")
print(f"{'='*70}")

for n in [5, 7, 9, 11, 13]:
    coefficients = [0, 0, 0, 0]  # const, t3, t5, bc

    for k in range(5):
        # Stirling weights: S(4,k)*k!
        # S(4,k)*k!: S(4,0)=0, S(4,1)=1, S(4,2)=7, S(4,3)=6, S(4,4)=1
        stirling_wts = {0: 0, 1: 1, 2: 14, 3: 36, 4: 24}
        wt = stirling_wts[k]
        if wt == 0: continue

        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_extended(pat, n)
            if r is None:
                print(f"  n={n}, k={k}: UNHANDLED pattern {pat} x {cnt}")
                continue
            for i in range(4):
                coefficients[i] += wt * cnt * r[i]

    print(f"  n={n:2d}: m4 = {coefficients[0]} + {coefficients[1]}*t3 + {coefficients[2]}*t5 + {coefficients[3]}*bc")

    # Analyze ratios
    if coefficients[2] != 0:
        print(f"         t5_coeff = {coefficients[2]}, t5/(n-2)! = {coefficients[2]/factorial(n-2):.4f}")
    if coefficients[3] != 0:
        print(f"         bc_coeff = {coefficients[3]}, bc/(n-2)! = {coefficients[3]/factorial(n-2):.4f}")

# =====================================================================
# Fit m4 coefficients
# =====================================================================
print(f"\n{'='*70}")
print("m4 closed form fitting")
print(f"{'='*70}")

ns = [5, 7, 9, 11, 13]
m4_data = {}

for n in ns:
    coefficients = [0, 0, 0, 0]
    for k in range(5):
        stirling_wts = {0: 0, 1: 1, 2: 14, 3: 36, 4: 24}
        wt = stirling_wts[k]
        if wt == 0: continue
        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_extended(pat, n)
            if r is None: continue
            for i in range(4):
                coefficients[i] += wt * cnt * r[i]
    m4_data[n] = coefficients

# t3 coefficient
print("\n  t3 coefficient analysis:")
for n in ns:
    r = m4_data[n][1] / factorial(n-2)
    print(f"    n={n}: t3_coeff/(n-2)! = {r:.4f}")

# Try ratios for t3:
print("\n  t3_coeff/(n-2)! differences:")
vals = [m4_data[n][1] / factorial(n-2) for n in ns]
for i in range(1, len(vals)):
    print(f"    {ns[i]}-{ns[i-1]}: {vals[i]-vals[i-1]:.4f}")

# t5 coefficient
print("\n  t5 coefficient analysis:")
for n in ns:
    if m4_data[n][2] != 0:
        r = m4_data[n][2] / factorial(n-2)
        print(f"    n={n}: t5_coeff/(n-2)! = {r:.4f}")

# bc coefficient
print("\n  bc coefficient analysis:")
for n in ns:
    if m4_data[n][3] != 0:
        r = m4_data[n][3] / factorial(n-2)
        print(f"    n={n}: bc_coeff/(n-2)! = {r:.4f}")

# const
print("\n  constant/(n!) analysis:")
for n in ns:
    r = m4_data[n][0] / factorial(n)
    print(f"    n={n}: const/n! = {r:.6f}")

# =====================================================================
# Brute-force verify at n=5 and n=7
# =====================================================================
print(f"\n{'='*70}")
print("Brute-force m4 verification")
print(f"{'='*70}")

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_5_cycles(A, n):
    count = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%5]] for i in range(5)):
                count += 1
    return count // 5  # each 5-cycle counted 5 times (rotations)

def count_bc(A, n):
    """Count pairs of vertex-disjoint 3-cycles summed over 6-subsets."""
    total = 0
    for S in combinations(range(n), 6):
        cyc_triples = []
        for triple in combinations(S, 3):
            a, b, c = triple
            if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
                cyc_triples.append(set(triple))
        # Count disjoint pairs
        for i in range(len(cyc_triples)):
            for j in range(i+1, len(cyc_triples)):
                if cyc_triples[i].isdisjoint(cyc_triples[j]):
                    total += 1
    return total

for n in [5, 7]:
    print(f"\n  n={n}:")
    for trial in range(5):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        t3 = count_3_cycles(A, n)
        t5 = count_5_cycles(A, n)
        bc = count_bc(A, n)

        m4_actual = sum(
            sum(1 for i in range(n-1) if A[p[i]][p[i+1]])**4
            for p in permutations(range(n))
        )

        c = m4_data[n]
        m4_pred = c[0] + c[1]*t3 + c[2]*t5 + c[3]*bc

        print(f"    Trial {trial}: t3={t3}, t5={t5}, bc={bc}, "
              f"m4={m4_actual}, pred={m4_pred}, {'OK' if m4_actual == m4_pred else 'FAIL'}")

# =====================================================================
# Summary of results
# =====================================================================
print(f"\n{'='*70}")
print("SUMMARY OF CLOSED FORMS")
print(f"{'='*70}")
print("""
  m0 = n!  (universal)
  m1 = n! * (n-1)/2  (universal)
  m2 = n!*(3n^2-5n+4)/12 + 4*(n-2)!*t3
  m3 = n!*(n-1)(n^2-n+2)/8 + 6*(n-1)!*t3

  For m4: depends on (t3, t5, bc). Closed forms TBD from fitting above.
""")

print("DONE")
