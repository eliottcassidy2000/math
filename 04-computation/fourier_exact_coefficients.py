"""
fourier_exact_coefficients.py -- kind-pasteur-2026-03-14-S75
Compute EXACT Fourier coefficients and find closed-form formulas.

From S73:
- All odd-level coefficients = 0 (proved: H(T) = H(T^op))
- Level-2 coefficients at n=3: +-0.5, n=4: +-0.5, n=5: +-0.75
- Level-4 coefficients at n=5: +-0.125

GOAL: Find closed-form for the level-2 Fourier coefficient in terms of n
and the arc pair structure.

KEY OBSERVATION: The level-2 coefficient H_hat({e1, e2}) depends on
the GRAPH STRUCTURE of the arc pair:
- Case 1: e1 and e2 share a vertex (adjacent arcs)
  Sub-case 1a: they form a "V" (both from or both to the shared vertex)
  Sub-case 1b: they form a "path" (one in, one out of shared vertex)
- Case 2: e1 and e2 are disjoint (non-adjacent arcs)

From S73 data:
n=3: all level-2 are +-0.5 (all arcs share a vertex at n=3)
n=4: all level-2 are +-0.5 (12 nonzero, all +-0.5)
n=5: all level-2 are +-0.75 (30 nonzero, all +-0.75)

QUESTION: Does the magnitude depend on n, or on n and the adjacency type?
"""

import numpy as np
from collections import Counter, defaultdict
import sys, math

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("EXACT FOURIER COEFFICIENTS OF H(T)")
    print("kind-pasteur-2026-03-14-S75")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        m = n * (n - 1) // 2
        N = 2**m
        arcs = [(i, j) for i in range(n) for j in range(i+1, n)]

        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m}")
        print(f"{'='*70}")

        if n > 5:
            print(f"  (n=6 too large for exact Fourier, skipping)")
            continue

        # Compute all H values
        H_values = np.zeros(N, dtype=float)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H_values[bits] = compute_H_dp(A, n)

        # Fast Walsh-Hadamard
        H_hat = H_values.copy()
        for i in range(m):
            step = 1 << (i + 1)
            half = 1 << i
            for j in range(0, N, step):
                for k in range(half):
                    u, v = H_hat[j+k], H_hat[j+k+half]
                    H_hat[j+k], H_hat[j+k+half] = u+v, u-v
        H_hat /= N

        # Level-0
        print(f"\n  Level 0: H_hat(empty) = {H_hat[0]}")
        print(f"    = n!/2^(n-1) = {math.factorial(n)/2**(n-1)}")

        # Level-1 (all should be 0)
        level1 = [H_hat[1 << i] for i in range(m)]
        print(f"\n  Level 1: all zero? {all(abs(c) < 1e-10 for c in level1)}")

        # Level-2: classify by arc pair type
        print(f"\n  Level 2 coefficients by arc pair type:")

        # Types: adjacent (share vertex) or disjoint
        # For adjacent: which vertex is shared, and what's the orientation?
        type_coeffs = defaultdict(list)

        for i in range(m):
            for j in range(i+1, m):
                S = (1 << i) | (1 << j)
                coeff = H_hat[S]
                if abs(coeff) < 1e-10:
                    continue

                arc1 = arcs[i]
                arc2 = arcs[j]
                shared = set(arc1) & set(arc2)

                if shared:
                    v_shared = shared.pop()
                    # Check if shared vertex is "low" or "high" in each arc
                    pos1 = 0 if arc1[0] == v_shared else 1
                    pos2 = 0 if arc2[0] == v_shared else 1
                    # pos=0 means shared is the lower vertex (arc goes FROM shared)
                    # pos=1 means shared is the higher vertex (arc goes TO shared)
                    type_key = f"adjacent({pos1},{pos2})"
                else:
                    type_key = "disjoint"

                type_coeffs[type_key].append(coeff)

        for type_key, coeffs in sorted(type_coeffs.items()):
            vals = sorted(set(round(c, 10) for c in coeffs))
            print(f"    {type_key}: {len(coeffs)} pairs, coefficients = {vals}")

        # Exact rational values
        # From the data: all level-2 coefficients have the same magnitude
        nonzero_level2 = [H_hat[S] for S in range(N)
                          if bin(S).count('1') == 2 and abs(H_hat[S]) > 1e-10]
        if nonzero_level2:
            magnitude = abs(nonzero_level2[0])
            all_same_mag = all(abs(abs(c) - magnitude) < 1e-10 for c in nonzero_level2)
            print(f"\n  All level-2 same magnitude? {all_same_mag}")
            print(f"  Magnitude = {magnitude}")

            # Try to identify as rational
            # magnitude * 2^n = ?
            for denom in range(1, 100):
                numer = round(magnitude * denom)
                if abs(magnitude - numer/denom) < 1e-10:
                    print(f"  Magnitude = {numer}/{denom}")
                    break

            # Check: is magnitude = (n-2)!/2^(n-1)?
            candidate = math.factorial(n-2) / 2**(n-1)
            print(f"  (n-2)!/2^(n-1) = {candidate}")
            print(f"  Match? {abs(magnitude - candidate) < 1e-10}")

        # Level-4
        level4_coeffs = [H_hat[S] for S in range(N)
                         if bin(S).count('1') == 4 and abs(H_hat[S]) > 1e-10]
        if level4_coeffs:
            magnitude4 = abs(level4_coeffs[0])
            all_same_mag4 = all(abs(abs(c) - magnitude4) < 1e-10 for c in level4_coeffs)
            print(f"\n  Level-4 coefficients: {len(level4_coeffs)} nonzero")
            print(f"  All same magnitude? {all_same_mag4}")
            if all_same_mag4:
                print(f"  Magnitude = {magnitude4}")
                for denom in range(1, 200):
                    numer = round(magnitude4 * denom)
                    if abs(magnitude4 - numer/denom) < 1e-10:
                        print(f"  Magnitude = {numer}/{denom}")
                        break

        # SIGN PATTERN of level-2 coefficients
        print(f"\n  Level-2 SIGN PATTERN:")
        for i in range(m):
            for j in range(i+1, m):
                S = (1 << i) | (1 << j)
                coeff = H_hat[S]
                if abs(coeff) < 1e-10:
                    continue
                sign = '+' if coeff > 0 else '-'
                arc1 = arcs[i]
                arc2 = arcs[j]
                shared = set(arc1) & set(arc2)

                if shared:
                    v = shared.pop()
                    # Is v the "middle" of a 2-path?
                    # arc1 = (a,v) or (v,b), arc2 = (c,v) or (v,d)
                    # 2-path: one arc ends at v, other starts at v
                    is_path = (arc1[1] == v and arc2[0] == v) or (arc2[1] == v and arc1[0] == v)
                    # V-shape: both arcs start at v or both end at v
                    is_v = (arc1[0] == v and arc2[0] == v) or (arc1[1] == v and arc2[1] == v)

                    if is_path:
                        type_str = "2-path"
                    elif is_v:
                        type_str = "V-shape"
                    else:
                        type_str = "other"
                else:
                    type_str = "disjoint"

                if n <= 4 or (i < 3 and j < 5):  # limit output
                    print(f"    {arcs[i]}-{arcs[j]} [{type_str}]: {sign}{abs(coeff):.4f}")

        # Summary of sign pattern
        path_signs = []
        v_signs = []
        disj_signs = []
        for i in range(m):
            for j in range(i+1, m):
                S = (1 << i) | (1 << j)
                coeff = H_hat[S]
                if abs(coeff) < 1e-10:
                    continue
                arc1 = arcs[i]
                arc2 = arcs[j]
                shared = set(arc1) & set(arc2)
                if shared:
                    v = shared.pop()
                    is_path = (arc1[1] == v and arc2[0] == v) or (arc2[1] == v and arc1[0] == v)
                    if is_path:
                        path_signs.append('+' if coeff > 0 else '-')
                    else:
                        v_signs.append('+' if coeff > 0 else '-')
                else:
                    disj_signs.append('+' if coeff > 0 else '-')

        print(f"\n  Sign pattern summary:")
        print(f"    2-path pairs: {Counter(path_signs)} (expect all +?)")
        print(f"    V-shape pairs: {Counter(v_signs)} (expect all -?)")
        print(f"    Disjoint pairs: {Counter(disj_signs)}")

    # FORMULA CONJECTURE
    print(f"\n{'='*70}")
    print("FORMULA CONJECTURE")
    print(f"{'='*70}")
    print(f"""
  Level-2 Fourier coefficient H_hat({{e1, e2}}):

  Magnitude = (n-2)! / 2^(n-1) for all n?
    n=3: (1)!/4 = 0.25? But actual = 0.5. NO.

  Let me check another formula:
    n=3: 0.5 = 1/2
    n=4: 0.5 = 1/2
    n=5: 0.75 = 3/4

  Pattern: 1/2, 1/2, 3/4, ???

  Trying: f(n) = (n-1)/2^(floor(n/2))
    n=3: 2/2 = 1. NO.

  Trying: f(n) = n!/2^(2n-3)
    n=3: 6/8 = 3/4. NO (actual 1/2).

  Trying: f(n) = (n-2)!/2^(n-2)
    n=3: 1/2. YES!
    n=4: 2/4 = 1/2. YES!
    n=5: 6/8 = 3/4. YES!

  CONJECTURE: |H_hat(S)| = (n-2)!/2^(n-2) for all level-2 S.

  Sign: H_hat(S) > 0 iff arcs form a 2-path (shared vertex is "middle")
        H_hat(S) < 0 iff arcs form a V-shape (shared vertex is "source" or "sink")
        H_hat(S) = 0 for disjoint arcs at n=3,4 (not enough vertices)
        H_hat(S) = ??? for disjoint arcs at n>=5
""")

    # Verify the formula conjecture
    for n in [3, 4, 5]:
        predicted = math.factorial(n-2) / 2**(n-2)
        m = n * (n - 1) // 2
        N = 2**m
        H_values = np.zeros(N, dtype=float)
        for bits in range(N):
            A = bits_to_adj(bits, n)
            H_values[bits] = compute_H_dp(A, n)
        H_hat = H_values.copy()
        for i in range(m):
            step = 1 << (i + 1)
            half = 1 << i
            for j in range(0, N, step):
                for k in range(half):
                    u, v = H_hat[j+k], H_hat[j+k+half]
                    H_hat[j+k], H_hat[j+k+half] = u+v, u-v
        H_hat /= N

        nonzero_level2 = [abs(H_hat[S]) for S in range(N) if bin(S).count('1') == 2 and abs(H_hat[S]) > 1e-10]
        actual = nonzero_level2[0] if nonzero_level2 else 0
        print(f"  n={n}: predicted = {predicted}, actual = {actual}, match = {abs(predicted - actual) < 1e-10}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
