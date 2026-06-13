#!/usr/bin/env python3
"""
disj_fraction_analysis.py -- kind-pasteur-2026-03-13-S60

KEY QUESTION: What fraction of (k,k')-cycle pairs are disjoint?

disj_frac(k,k') = disj(k,k') / total_pairs(k,k')
                 = 1 - conflict_ratio

At p=7 Paley: disj_frac(3,3) = 7/91 = 1/13
              For uniform random: C(4,3)/C(7,3) = 4/35

So the disjoint fraction for CYCLES is different from uniform random subsets.
This makes sense: cycles are not uniformly distributed over vertex sets.

DEEPER QUESTION: Is disj(k,k') = (c_k * c_k' / something) * [combinatorial factor]?

Actually, since H = linear(c_k), and H = 1 + 2*alpha_1 + 4*alpha_2 + ...,
and alpha_2 = sum_{k<=k'} disj(k,k'), we know:
  alpha_2 = quadratic in c_k (product of cycle counts)
  But H = linear in c_k
  So the quadratic terms must cancel or reduce to linear!

HOW? If disj(k,k') = A(k,k') * c_k * c_k' + B(k) * c_k + B(k') * c_k' + C(k,k')
then for alpha_2 to be linear in c_k, the A(k,k') terms must cancel in the sum.

Actually, from the data:
  At p=7: disj(3,3) = 7 or 14 (for Paley vs Interval)
  c3 = 14 for both (constant)
  So disj(3,3) varies even though c3 is constant!

This means disj(3,3) is NOT just a function of c_3.
But it IS a function of S4 (from THM-156).

Let me compute: disj(k,k') / c_k for each (k,k') pair.
If disj(k,k') = L(k,k') * c_k + R(k,k'), then disj/c_k = L + R/c_k.
"""

from itertools import combinations
from fractions import Fraction


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def analyze_disj_structure(p):
    """Compute disjoint pair counts for all orientations and test linearity in c_k."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"DISJOINT PAIR STRUCTURE at p={p}, m={m}, N={N}")
    print(f"{'='*70}")

    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)

        # Enumerate cycles
        cycles_by_k = {}
        for k in range(3, p + 1, 2):
            cyc_list = []
            for subset in combinations(range(p), k):
                nc = count_ham_cycles(A, list(subset))
                for _ in range(nc):
                    cyc_list.append(frozenset(subset))
            if cyc_list:
                cycles_by_k[k] = cyc_list

        c_k = {k: len(v) for k, v in cycles_by_k.items()}

        # Count disjoint pairs for each (k1, k2) with k1 <= k2
        disj = {}
        for k1 in sorted(cycles_by_k):
            for k2 in sorted(cycles_by_k):
                if k1 > k2:
                    continue
                if k1 + k2 > p:
                    disj[(k1, k2)] = 0  # pigeonhole
                    continue

                cyc1 = cycles_by_k[k1]
                cyc2 = cycles_by_k[k2]
                count = 0
                if k1 == k2:
                    for i in range(len(cyc1)):
                        for j in range(i + 1, len(cyc1)):
                            if not (cyc1[i] & cyc1[j]):
                                count += 1
                else:
                    for c1 in cyc1:
                        for c2 in cyc2:
                            if not (c1 & c2):
                                count += 1
                disj[(k1, k2)] = count

        data.append({
            'bits': bits, 'c_k': c_k, 'disj': disj,
            'alpha_2': sum(disj.values())
        })

    # Display
    cycle_lengths = sorted(data[0]['c_k'].keys())
    pair_types = sorted(set(k for d in data for k in d['disj'].keys()))

    # Which c_k vary?
    varying = [k for k in cycle_lengths
               if len(set(d['c_k'][k] for d in data)) > 1]
    constant = [(k, data[0]['c_k'][k]) for k in cycle_lengths
                if len(set(d['c_k'][k] for d in data)) == 1]

    print(f"\n  Constant c_k: {constant}")
    print(f"  Varying c_k: {varying}")

    # For each pair type with non-zero disjoint count, test linearity in c_k
    print(f"\n  Disjoint pair counts by type:")
    print(f"  {'bits':>4} " + " ".join(f"disj{k1}{k2}" for k1, k2 in pair_types if any(d['disj'][(k1,k2)] > 0 for d in data)) + f"  alpha_2")

    active_types = [(k1, k2) for k1, k2 in pair_types
                    if any(d['disj'][(k1, k2)] > 0 for d in data)]

    for d in data:
        vals = " ".join(f"{d['disj'][t]:>8}" for t in active_types)
        ck_str = " ".join(f"c{k}={d['c_k'][k]}" for k in varying)
        print(f"  {d['bits']:>4} {vals}  a2={d['alpha_2']:>6}  {ck_str}")

    # Test: is each disj(k1,k2) a LINEAR function of the varying c_k?
    print(f"\n  LINEARITY TEST for each disj(k1,k2):")
    import numpy as np

    for k1, k2 in active_types:
        y = np.array([d['disj'][(k1, k2)] for d in data], dtype=float)
        if len(set(y)) <= 1:
            print(f"    disj({k1},{k2}) = CONSTANT = {int(y[0])}")
            continue

        # Fit: disj = sum a_k * c_k + const
        X_cols = [np.array([d['c_k'][k] for d in data], dtype=float) for k in varying]
        X = np.column_stack(X_cols + [np.ones(len(data))])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        pred = X @ coeffs
        max_err = np.max(np.abs(y - pred))

        coeff_str = " + ".join(f"{c:.6f}*c{k}" for c, k in zip(coeffs[:-1], varying))
        coeff_str += f" + {coeffs[-1]:.2f}"
        exact = max_err < 0.5
        print(f"    disj({k1},{k2}) = {coeff_str}  max_err={max_err:.4f} {'EXACT' if exact else ''}")

    # Test: alpha_2 = f(c_k)?
    print(f"\n  ALPHA_2 LINEARITY:")
    y = np.array([d['alpha_2'] for d in data], dtype=float)
    X_cols = [np.array([d['c_k'][k] for d in data], dtype=float) for k in varying]
    X = np.column_stack(X_cols + [np.ones(len(data))])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    max_err = np.max(np.abs(y - pred))
    coeff_str = " + ".join(f"{c:.6f}*c{k}" for c, k in zip(coeffs[:-1], varying))
    coeff_str += f" + {coeffs[-1]:.2f}"
    print(f"    alpha_2 = {coeff_str}  max_err={max_err:.4f}")

    # KEY: Test disj(k1,k2) as function of c_k1 * c_k2 AND c_k1 and c_k2
    print(f"\n  PRODUCT TEST: disj(k1,k2) = a*c_k1*c_k2 + b*c_k1 + c*c_k2 + d?")
    for k1, k2 in active_types:
        y = np.array([d['disj'][(k1, k2)] for d in data], dtype=float)
        if len(set(y)) <= 1:
            continue

        if k1 == k2:
            ck_vals = np.array([d['c_k'][k1] for d in data], dtype=float)
            X = np.column_stack([
                ck_vals * (ck_vals - 1) / 2,  # C(c_k, 2)
                ck_vals,
                np.ones(len(data))
            ])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(y - pred))
            print(f"    disj({k1},{k1}) = {coeffs[0]:.6f}*C(c{k1},2) + {coeffs[1]:.6f}*c{k1} + {coeffs[2]:.2f}  err={max_err:.4f}")
        else:
            ck1_vals = np.array([d['c_k'][k1] for d in data], dtype=float)
            ck2_vals = np.array([d['c_k'][k2] for d in data], dtype=float)
            X = np.column_stack([
                ck1_vals * ck2_vals,
                ck1_vals,
                ck2_vals,
                np.ones(len(data))
            ])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(y - pred))
            print(f"    disj({k1},{k2}) = {coeffs[0]:.6f}*c{k1}*c{k2} + {coeffs[1]:.6f}*c{k1} + {coeffs[2]:.6f}*c{k2} + {coeffs[3]:.2f}  err={max_err:.4f}")

    return data


for p_val in [7, 11]:
    analyze_disj_structure(p_val)

print("\nDONE.")
