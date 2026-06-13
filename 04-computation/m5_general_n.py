#!/usr/bin/env python3
"""
Compute m5 at general n using extended sigma formulas.

m5 = SIGMA_1 + 30*SIGMA_2 + 150*SIGMA_3 + 240*SIGMA_4 + 120*SIGMA_5
where S(5,k)*k! = {1,30,150,240,120} for k=1,...,5.

Patterns at k=5:
  (1,1,1,1,1): 5 isolated pairs - universal
  (2,1,1,1): one 3-set + three 2-sets - depends on t3
  (3,1,1): one 4-set + two 2-sets - depends on t3
  (2,2,1): two 3-sets + one 2-set - depends on t3, bc
  (4,1): one 5-set + one 2-set - depends on t3, t5
  (3,2): one 4-set + one 3-set - depends on t3, bc
  (5,): one 6-set - depends on t3, t5, bc

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
from collections import Counter
import random
import numpy as np

def count_position_patterns(n, k):
    if k == 0: return {(): 1}
    patterns = Counter()
    for S in combinations(range(n-1), k):
        pos = sorted(S)
        comps, comp = [], [pos[0]]
        for i in range(1, len(pos)):
            if pos[i] == comp[-1] + 1:
                comp.append(pos[i])
            else:
                comps.append(len(comp))
                comp = [pos[i]]
        comps.append(len(comp))
        patterns[tuple(sorted(comps, reverse=True))] += 1
    return dict(patterns)

def sigma_full(pattern, n):
    """
    Return (const, t3_coeff, t5_coeff, bc_coeff) for sigma(pattern) at general n.
    Now handles ALL patterns up to total size 11 (needed for k=5).
    """
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0: return (0, 0, 0, 0)

    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)
    F = factorial(free)

    def pair_prod(start):
        """Product of C(start-2i, 2) for i=0,...,num_size1-1."""
        p, r = 1, start
        for i in range(num_size1):
            p *= comb(r, 2)
            r -= 2
        return p

    # === ALL-ISOLATED: sigma = n!/2^k ===
    if len(big_sizes) == 0:
        return (factorial(n) // (2**len(sizes)), 0, 0, 0)

    # === SINGLE BIG COMPONENT ===
    if len(big_sizes) == 1:
        s = big_sizes[0]
        verts_big = s + 1

        if s == 2:  # 3-set: H = 1 + 2*cyc
            pp = pair_prod(n - 3)
            return (F*comb(n,3)*pp, F*2*pp, 0, 0)

        if s == 3:  # 4-set: H = 1 + 2*c3
            pp = pair_prod(n - 4)
            # sum_{4-sets} H = C(n,4) + 2*(n-3)*t3
            return (F*comb(n,4)*pp, F*2*(n-3)*pp, 0, 0)

        if s == 4:  # 5-set: H = 1 + 2*(c3+c5) via OCF at n=5
            pp = pair_prod(n - 5)
            # sum_{5-sets} = C(n,5) + 2*C(n-3,2)*t3 + 2*t5
            return (F*comb(n,5)*pp, F*2*comb(n-3,2)*pp, F*2*pp, 0)

        if s == 5:  # 6-set: H = 1 + 2*c3 + 2*c5 + 4*bc via OCF at n=6
            pp = pair_prod(n - 6)
            # sum_{6-sets} = C(n,6) + 2*C(n-3,3)*t3 + 2*(n-5)*t5 + 4*bc
            return (F*comb(n,6)*pp, F*2*comb(n-3,3)*pp, F*2*(n-5)*pp, F*4*pp)

    # === TWO BIG COMPONENTS ===
    if big_sizes == [2, 2]:  # Two 3-sets
        pp = pair_prod(n - 6)
        # sum_{ordered (G1,G2) disjoint 3-sets} (1+2*cyc1)(1+2*cyc2)
        # = C(n,3)*C(n-3,3) + 4*C(n-3,3)*t3 + 8*bc
        return (F*comb(n,3)*comb(n-3,3)*pp, F*4*comb(n-3,3)*pp, 0, F*8*pp)

    if big_sizes == [3, 2]:  # 4-set + 3-set
        pp = pair_prod(n - 7)
        # sum_{(G1 of 4, G2 of 3) disjoint} (1+2*c3(G1))(1+2*cyc(G2))
        # const: C(n,4)*C(n-4,3)
        # t3: 2*(n-3)*C(n-4,3) + 2*C(n-3,4)
        # bc: 8*(n-6)  [from cross-term 4*2*(n-6)*bc]
        const = comb(n,4)*comb(n-4,3)
        t3_c = 2*(n-3)*comb(n-4,3) + 2*comb(n-3,4)
        bc_c = 8*(n-6) if n >= 6 else 0
        return (F*const*pp, F*t3_c*pp, 0, F*bc_c*pp)

    if big_sizes == [4, 2]:  # 5-set + 3-set
        # NOT needed for k<=5 unless pattern is (4,2) which has total_verts=5+3=8
        # Actually pattern (4,2) has sizes [4,2], total_verts = 5+3 = 8, but this
        # is a k=6 pattern. Not needed for m5.
        pp = pair_prod(n - 8)
        # sum_{(G1 of 5, G2 of 3) disjoint} H5(G1)*H3(G2)
        # H5 = 1+2*c3+2*c5, H3 = 1+2*cyc
        # This involves a cross-term: c3(G1)*cyc(G2) and c5(G1)*cyc(G2)
        # For now, not implementing this
        return None

    if big_sizes == [3, 3]:  # Two 4-sets
        # NOT needed for k<=5
        return None

    if big_sizes == [2, 2, 2]:  # Three 3-sets
        # NOT needed for k<=5 (would be k=6)
        return None

    return None

# =====================================================================
# Compute m5 at various n
# =====================================================================
print("=" * 70)
print("m5 = SIGMA_1 + 30*SIGMA_2 + 150*SIGMA_3 + 240*SIGMA_4 + 120*SIGMA_5")
print("S(5,1)=1, S(5,2)=15, S(5,3)=25, S(5,4)=10, S(5,5)=1")
print("=" * 70)

# Stirling wts: S(j,k)*k!
stirling5 = {0: 0, 1: 1, 2: 30, 3: 150, 4: 240, 5: 120}

ns = [7, 9, 11, 13, 15, 17]
m5_data = {}

for n in ns:
    coeffs = [0, 0, 0, 0]  # const, t3, t5, bc
    complete = True
    for k in range(6):
        wt = stirling5.get(k, 0)
        if wt == 0: continue
        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_full(pat, n)
            if r is None:
                print(f"  n={n}, k={k}: UNHANDLED {pat} x {cnt}")
                complete = False
                continue
            for i in range(4):
                coeffs[i] += wt * cnt * r[i]

    m5_data[n] = coeffs
    print(f"\n  n={n:2d}: m5 = {coeffs[0]} + {coeffs[1]}*t3 + {coeffs[2]}*t5 + {coeffs[3]}*bc")
    if not complete:
        print(f"         (INCOMPLETE — some patterns unhandled)")

# =====================================================================
# Analyze m5 coefficients
# =====================================================================
print(f"\n{'='*70}")
print("m5 coefficient analysis")
print(f"{'='*70}")

print("\n  t3 coefficient:")
for n in ns:
    c = m5_data[n]
    r = c[1] / factorial(n-2) if c[1] != 0 else 0
    print(f"    n={n}: {c[1]}, /(n-2)! = {r:.4f}")

# Try: does t3_coeff relate to some polynomial?
print("\n  t5 coefficient:")
for n in ns:
    c = m5_data[n]
    if c[2] != 0:
        r = c[2] / factorial(n-4)
        print(f"    n={n}: {c[2]}, /(n-4)! = {r:.4f}")

print("\n  bc coefficient:")
for n in ns:
    c = m5_data[n]
    if c[3] != 0:
        r = c[3] / factorial(n-4)
        print(f"    n={n}: {c[3]}, /(n-4)! = {r:.4f}")

print("\n  bc/t5 ratio:")
for n in ns:
    c = m5_data[n]
    if c[2] != 0:
        print(f"    n={n}: bc/t5 = {c[3]/c[2]:.4f}")

# =====================================================================
# Brute-force verify at n=7
# =====================================================================
print(f"\n{'='*70}")
print("Brute-force verification at n=7")
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
    return count // 5

def count_bc(A, n):
    cyc_triples = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            cyc_triples.append(set(triple))
    total = 0
    for i in range(len(cyc_triples)):
        for j in range(i+1, len(cyc_triples)):
            if cyc_triples[i].isdisjoint(cyc_triples[j]):
                total += 1
    return total

n = 7
c = m5_data[7]
for trial in range(10):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles(A, n)
    bc = count_bc(A, n)

    m5_actual = sum(
        sum(1 for i in range(n-1) if A[p[i]][p[i+1]])**5
        for p in permutations(range(n))
    )
    m5_pred = c[0] + c[1]*t3 + c[2]*t5 + c[3]*bc
    ok = m5_actual == m5_pred
    print(f"  T{trial}: t3={t3:2d}, t5={t5:2d}, bc={bc:2d}, m5={m5_actual}, pred={m5_pred}, {'OK' if ok else 'FAIL'}")

# =====================================================================
# Brute-force verify at n=9 via DP
# =====================================================================
print(f"\n{'='*70}")
print("DP verification at n=9")
print(f"{'='*70}")

n = 9
c = m5_data[9]

def compute_moments_dp(A, n, max_j):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for f in range(n):
                cnt = dp.get((mask, v, f), 0)
                if cnt == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    key = (mask | (1 << u), u, f + A[v][u])
                    dp[key] = dp.get(key, 0) + cnt
    full = (1 << n) - 1
    moments = [0] * (max_j + 1)
    for v in range(n):
        for f in range(n):
            cnt = dp.get((full, v, f), 0)
            if cnt == 0: continue
            for j in range(max_j + 1):
                moments[j] += cnt * f**j
    return moments

for trial in range(10):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles(A, n)
    bc = count_bc(A, n)
    moments = compute_moments_dp(A, n, 5)

    m5_pred = c[0] + c[1]*t3 + c[2]*t5 + c[3]*bc
    ok = moments[5] == m5_pred
    print(f"  T{trial}: t3={t3:3d}, t5={t5:4d}, bc={bc:3d}, m5={moments[5]}, pred={m5_pred}, {'OK' if ok else 'FAIL'}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
