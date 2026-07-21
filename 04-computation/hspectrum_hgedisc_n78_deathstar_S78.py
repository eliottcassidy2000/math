#!/usr/bin/env python3
"""
hspectrum_hgedisc_n78_deathstar_S78.py  (HYP-8636)

Two n=7,8 tests from the spectral-Graffiti miner:
 (1) H-SPECTRUM: at n<=6 the missing odd H-values are {7,21,35,39}, not just {7,21}
     (S70). Do 35 and 39 appear at n=7,8? (tests S70's 'odds \ {7,21}' claim; a
     marginal-threshold check on my own result.)
 (2) H >= disc(T) := |det(I+K)|/2^{n-1}  (survived n<=6, tight at transitive):
     does it survive n=7,8, and is equality UNIQUELY the transitive tournament?
Random sampling (n=7: 2^21, n=8: 2^28 too large to enumerate); samples find
achievable values / counterexamples, not completeness.
"""
import random
from itertools import combinations
random.seed(20260721)
try:
    import numpy as np
except Exception:
    np = None

def rand_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if random.random() < 0.5: A[i][j] = 1
        else: A[j][i] = 1
    return A

def ham_paths(A, n):
    full = (1 << n) - 1
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                ov = out[v]
                for w in range(n):
                    if not (mask >> w & 1) and (ov >> w & 1):
                        dp[mask | 1 << w][w] += c
    return sum(dp[full][v] for v in range(n))

def disc(A, n):
    K = np.array(A, float) - np.array(A, float).T
    return round(abs(np.linalg.det(np.eye(n) + K))) / 2**(n-1)

def is_transitive(A, n):
    s = sorted(sum(A[i]) for i in range(n))
    return s == list(range(n))

def run(n, R):
    Hset = set(); viol = 0; eq_nontrans = 0; eqcount = 0
    for _ in range(R):
        A = rand_tournament(n)
        H = ham_paths(A, n)
        Hset.add(H)
        if np is not None:
            d = disc(A, n)
            if H < d - 1e-9: viol += 1
            if abs(H - d) < 1e-9:
                eqcount += 1
                if not is_transitive(A, n): eq_nontrans += 1
    return Hset, viol, eqcount, eq_nontrans

for n, R in [(7, 400000), (8, 120000)]:
    Hset, viol, eqc, eqnt = run(n, R)
    hs = sorted(Hset)
    checks = {v: (v in Hset) for v in (7, 21, 35, 39, 49, 63)}
    print(f"\nn={n} ({R} random samples):")
    print(f"  #distinct H seen = {len(hs)}; min={hs[0]}, max={hs[-1]}")
    print(f"  small odd values present? {checks}")
    if np is not None:
        print(f"  H >= disc violations: {viol}   equality cases: {eqc} (of which NON-transitive: {eqnt})")
    # which odds up to 45 still missing (compare to n<=6 gap {7,21,35,39})
    missing = [k for k in range(1, 46, 2) if k not in Hset]
    print(f"  odd values in [1,45] still MISSING at n={n}: {missing}")

print("\n(S70 claim: H-spectrum = odds minus {7,21}. If 35,39 appear here, consistent;")
print(" if they persist as missing, S70's forbidden set is larger than {7,21}.)")
