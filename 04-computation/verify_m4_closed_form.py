#!/usr/bin/env python3
"""
Verify the complete m4 closed form at general n.

m4 = n!(15n^4-30n^3+65n^2-82n+48)/240
     + 2(n-2)!(3n^2-5n+4) * t3
     + 48(n-4)! * t5
     + 96(n-4)! * bc

Also brute-force verify at n=5,7,9.

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

def sigma_extended(pattern, n):
    sizes = list(pattern)
    total_verts = sum(s + 1 for s in sizes)
    free = n - total_verts
    if free < 0: return (0, 0, 0, 0)
    num_size1 = sizes.count(1)
    big_sizes = sorted([s for s in sizes if s > 1], reverse=True)
    F = factorial(free)
    def pair_product(start_remaining):
        p, r = 1, start_remaining
        for i in range(num_size1):
            p *= comb(r, 2); r -= 2
        return p
    if len(big_sizes) == 0: return (factorial(n) // (2**len(sizes)), 0, 0, 0)
    if big_sizes == [2]:
        pp = pair_product(n-3); return (F*comb(n,3)*pp, F*2*pp, 0, 0)
    if big_sizes == [3]:
        pp = pair_product(n-4); return (F*comb(n,4)*pp, F*2*(n-3)*pp, 0, 0)
    if big_sizes == [4]:
        pp = pair_product(n-5); return (F*comb(n,5)*pp, F*2*comb(n-3,2)*pp, F*2*pp, 0)
    if big_sizes == [2, 2]:
        pp = pair_product(n-6)
        return (F*comb(n,3)*comb(n-3,3)*pp, F*4*comb(n-3,3)*pp, 0, F*8*pp)
    return None

# Compute m4 from sigma decomposition
def compute_m4_sigma(n):
    wts = {0: 0, 1: 1, 2: 14, 3: 36, 4: 24}
    coeffs = [0, 0, 0, 0]
    for k in range(5):
        if wts[k] == 0: continue
        for pat, cnt in count_position_patterns(n, k).items():
            r = sigma_extended(pat, n)
            if r is None: return None
            for i in range(4): coeffs[i] += wts[k] * cnt * r[i]
    return coeffs

# Closed form
def m4_closed(n):
    const = factorial(n) * (15*n**4 - 30*n**3 + 65*n**2 - 82*n + 48) // 240
    t3 = 2 * factorial(n-2) * (3*n**2 - 5*n + 4)
    t5 = 48 * factorial(n-4) if n >= 4 else 0
    bc = 96 * factorial(n-4) if n >= 4 else 0
    return [const, t3, t5, bc]

# =====================================================================
# Part 1: Verify closed form matches sigma computation
# =====================================================================
print("=" * 70)
print("Part 1: Closed form vs sigma computation")
print("=" * 70)

for n in range(5, 22, 2):
    sigma_c = compute_m4_sigma(n)
    closed_c = m4_closed(n)
    if sigma_c is None:
        print(f"  n={n}: sigma incomplete")
        continue
    ok = all(sigma_c[i] == closed_c[i] for i in range(4))
    if not ok:
        for i, name in enumerate(["const", "t3", "t5", "bc"]):
            if sigma_c[i] != closed_c[i]:
                print(f"  n={n}: {name} FAIL: sigma={sigma_c[i]}, closed={closed_c[i]}")
    else:
        print(f"  n={n:2d}: ALL OK")

# =====================================================================
# Part 2: Brute-force verification at n=5,7
# =====================================================================
print(f"\n{'='*70}")
print("Part 2: Brute-force verification")
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
    total = 0
    for S in combinations(range(n), 6):
        cyc_triples = []
        for triple in combinations(S, 3):
            a, b, c = triple
            if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
                cyc_triples.append(set(triple))
        for i in range(len(cyc_triples)):
            for j in range(i+1, len(cyc_triples)):
                if cyc_triples[i].isdisjoint(cyc_triples[j]):
                    total += 1
    return total

for n in [5, 7]:
    print(f"\n  n={n}:")
    c = m4_closed(n)
    for trial in range(8):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5: A[i][j] = 1
                else: A[j][i] = 1

        t3 = count_3_cycles(A, n)
        t5 = count_5_cycles(A, n)
        bc = count_bc(A, n)

        m4_actual = sum(
            sum(1 for i in range(n-1) if A[p[i]][p[i+1]])**4
            for p in permutations(range(n))
        )
        m4_pred = c[0] + c[1]*t3 + c[2]*t5 + c[3]*bc
        ok = m4_actual == m4_pred
        print(f"    T{trial}: t3={t3:2d}, t5={t5:2d}, bc={bc:2d}, m4={m4_actual}, pred={m4_pred}, {'OK' if ok else 'FAIL'}")

# =====================================================================
# Part 3: Verify at n=9 with DP-based moment computation
# =====================================================================
print(f"\n{'='*70}")
print("Part 3: n=9 verification via DP")
print(f"{'='*70}")

n = 9

def compute_moments_dp(A, n, max_j):
    """Compute sum_P f^j for j=0..max_j using DP on bitmask."""
    # dp[mask][v] = sum over Ham paths ending at v using vertices in mask of f^power...
    # Actually easier: dp[mask][v] = number of paths * weight tracking
    # Just compute sum_P f^j directly via bitmask DP tracking f distribution.

    # dp[mask][v] = list of (count, f) pairs? No, too complex.
    # Better: dp[mask][v][f] = #{paths from some start to v using mask with f forward arcs}
    # But f can be up to n-1=8, and mask up to 2^9=512, v up to 9. So dp is 512*9*9 ~ 40K. Fine.

    dp = {}  # dp[(mask, v, f)] = count
    for v in range(n):
        dp[(1 << v, v, 0)] = 1  # single vertex, 0 arcs

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for f in range(n):
                key = (mask, v, f)
                if key not in dp or dp[key] == 0: continue
                cnt = dp[key]
                for u in range(n):
                    if mask & (1 << u): continue
                    new_mask = mask | (1 << u)
                    new_f = f + A[v][u]
                    new_key = (new_mask, u, new_f)
                    dp[new_key] = dp.get(new_key, 0) + cnt

    full = (1 << n) - 1
    moments = [0] * (max_j + 1)
    for v in range(n):
        for f in range(n):
            key = (full, v, f)
            cnt = dp.get(key, 0)
            if cnt == 0: continue
            for j in range(max_j + 1):
                moments[j] += cnt * f**j
    return moments

c = m4_closed(n)
print(f"  Formula: m4 = {c[0]} + {c[1]}*t3 + {c[2]}*t5 + {c[3]}*bc")

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

    moments = compute_moments_dp(A, n, 4)
    m4_actual = moments[4]
    m4_pred = c[0] + c[1]*t3 + c[2]*t5 + c[3]*bc
    ok = m4_actual == m4_pred
    print(f"  T{trial:2d}: t3={t3:3d}, t5={t5:4d}, bc={bc:3d}, m4={m4_actual}, pred={m4_pred}, {'OK' if ok else 'FAIL'}")

# Also verify m2 and m3
print(f"\n  Also checking m2 and m3 at n=9:")
for trial in range(5):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    t3 = count_3_cycles(A, n)
    moments = compute_moments_dp(A, n, 3)

    m2_pred = factorial(n)*(3*n**2-5*n+4)//12 + 4*factorial(n-2)*t3
    m3_pred = factorial(n)*(n-1)*(n**2-n+2)//8 + 6*factorial(n-1)*t3

    ok2 = moments[2] == m2_pred
    ok3 = moments[3] == m3_pred
    print(f"  T{trial}: t3={t3}, m2 {'OK' if ok2 else 'FAIL'}, m3 {'OK' if ok3 else 'FAIL'}")

print(f"\n{'='*70}")
print("ALL DONE")
print(f"{'='*70}")
