#!/usr/bin/env python3
"""
H_maximization_mechanism.py -- kind-pasteur-2026-03-13-S60

THE KEY QUESTION: Why does Paley maximize H = 1 + 2*N + 4*alpha_2 + 8*alpha_3?

At p=11, the term dominance analysis showed:
- Paley has MORE total cycles (N=21169 vs 18397-19629)
- Paley has FEWER disjoint pairs (alpha_2=10879 vs 10912-11220)
- Paley has FEWER disjoint triples (alpha_3=1155 vs 1155-1474)

The linear term 2*N dominates over the reductions in 4*alpha_2 and 8*alpha_3.

But WHY are more cycles and fewer disjoint pairs linked?
Hypothesis: More cycles => more overlapping => fewer disjoint pairs.
This is the OVERLAP DENSITY mechanism.

Let's quantify: compute the overlap density D(T) = alpha_2 / C(N,2)
and the packing fraction alpha_3 / C(N,3).

Also: compare the "virtual linear model" H_lin = 1 + c*N (best fit)
with actual H values.
"""

from itertools import combinations
from collections import defaultdict
import math

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

def compute_H_heldkarp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


for p in [7, 11]:
    m = (p - 1) // 2
    print(f"\n{'='*70}")
    print(f"  H-MAXIMIZATION MECHANISM AT p={p}")
    print(f"{'='*70}")

    results = []
    for bits in range(1 << m):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        H = compute_H_heldkarp(A, p)

        active = []
        for k in range(3, p+1, 2):
            for subset in combinations(range(p), k):
                fs = frozenset(subset)
                nc = count_ham_cycles(A, list(subset))
                if nc > 0:
                    active.append((fs, k, nc))

        N = sum(nc for _, _, nc in active)

        alpha_2 = 0
        disj = []
        for i in range(len(active)):
            for j in range(i+1, len(active)):
                if not (active[i][0] & active[j][0]):
                    alpha_2 += active[i][2] * active[j][2]
                    disj.append((i, j))

        alpha_3 = 0
        if p <= 11:
            for i, j in disj:
                V1, k1, n1 = active[i]
                V2, k2, n2 = active[j]
                used = V1 | V2
                for l in range(j+1, len(active)):
                    V3, k3, n3 = active[l]
                    if not (used & V3):
                        alpha_3 += n1 * n2 * n3

        # Count active vertex sets and total disjoint SET pairs (unweighted)
        n_sets = len(active)
        n_disj_sets = len(disj)

        # Overlap density
        n_choose_2 = N * (N - 1) // 2
        overlap_frac = 1 - alpha_2 / n_choose_2 if n_choose_2 > 0 else 0

        results.append({
            'bits': bits, 'H': H, 'N': N,
            'alpha_2': alpha_2, 'alpha_3': alpha_3,
            'n_sets': n_sets, 'n_disj_sets': n_disj_sets,
            'overlap_frac': overlap_frac
        })

    results.sort(key=lambda r: r['H'], reverse=True)

    # Compute mean overlap weight
    print(f"\n  {'bits':>5} {'H':>8} {'N':>6} {'a2':>8} {'a3':>6} {'ovlp%':>7} "
          f"{'H/N':>6} {'#sets':>5} {'#disj':>5}")
    print(f"  {'-'*60}")

    seen = set()
    for r in results:
        if r['H'] in seen:
            continue
        seen.add(r['H'])
        hn = r['H'] / r['N'] if r['N'] > 0 else 0
        ovlp_pct = 100 * r['overlap_frac']
        print(f"  {r['bits']:>5} {r['H']:>8} {r['N']:>6} {r['alpha_2']:>8} {r['alpha_3']:>6} "
              f"{ovlp_pct:>6.3f}% {hn:>6.3f} {r['n_sets']:>5} {r['n_disj_sets']:>5}")

    # Linear regression: H = a + b*N
    Ns = [r['N'] for r in results]
    Hs = [r['H'] for r in results]
    n_data = len(Ns)
    mx = sum(Ns) / n_data
    my = sum(Hs) / n_data
    cov_xy = sum((x-mx)*(y-my) for x,y in zip(Ns, Hs)) / n_data
    var_x = sum((x-mx)**2 for x in Ns) / n_data
    if var_x > 0:
        b = cov_xy / var_x
        a = my - b * mx
        residuals = [h - (a + b*n) for h, n in zip(Hs, Ns)]
        max_res = max(abs(r) for r in residuals)
        r_sq = 1 - sum(r**2 for r in residuals) / sum((h-my)**2 for h in Hs)
    else:
        b, a, max_res, r_sq = 0, my, 0, 1

    print(f"\n  Linear model: H = {a:.2f} + {b:.4f} * N")
    print(f"  R^2 = {r_sq:.6f}")
    print(f"  Max residual = {max_res:.2f}")
    print(f"  Note: exact formula coefficient is 2.0 for the N term")
    print(f"  The slope b={b:.4f} > 2 absorbs part of the alpha_2, alpha_3 contributions")

    # Is there a SINGLE function of N that gives H exactly?
    # H = 1 + 2N + 4*a2 + 8*a3
    # If a2 = f(N) and a3 = g(N), then H = h(N)
    # Let's check: is a2 a function of N?
    n_to_a2 = defaultdict(set)
    n_to_a3 = defaultdict(set)
    for r in results:
        n_to_a2[r['N']].add(r['alpha_2'])
        n_to_a3[r['N']].add(r['alpha_3'])

    a2_functional = all(len(v) == 1 for v in n_to_a2.values())
    a3_functional = all(len(v) == 1 for v in n_to_a3.values())

    print(f"\n  Is alpha_2 a function of N alone? {a2_functional}")
    if a2_functional:
        print(f"  N -> alpha_2 map: {dict((k, list(v)[0]) for k, v in sorted(n_to_a2.items()))}")
    else:
        for n_val, a2_vals in sorted(n_to_a2.items()):
            if len(a2_vals) > 1:
                print(f"  N={n_val}: alpha_2 takes {len(a2_vals)} values: {sorted(a2_vals)}")

    print(f"  Is alpha_3 a function of N alone? {a3_functional}")

    # The DEEPER question: what determines H? Let's check if H is a function of c_k tuple.
    ck_to_H = defaultdict(set)
    for r in results:
        # Need to recompute c_k for this
        pass

    # N/alpha_2 TRADE-OFF FORMULA
    # Define: f = alpha_2 / alpha_2_max where alpha_2_max = C(N,2)/something
    # Or better: look at the RATIO of the marginal gain from more cycles
    # vs the marginal loss from fewer disjoint pairs.
    #
    # dH/dN (holding alpha_2 constant) = 2
    # dH/d(alpha_2) (holding N constant) = 4
    #
    # So each disjoint pair is worth 4/2 = 2 "effective cycles".
    # But increasing N by 1 typically decreases the number of disjoint pairs.
    # The net effect depends on the OVERLAP COEFFICIENT:
    #   dN = +delta_N
    #   d(alpha_2) = -overlap_coeff * delta_N
    #   dH = 2*delta_N - 4*overlap_coeff*delta_N = (2 - 4*overlap_coeff)*delta_N
    #
    # So increasing N helps H iff overlap_coeff < 1/2.

    # Estimate overlap coefficient from data
    a2s = [r['alpha_2'] for r in results]

    if var_x > 0:
        cov_na2 = sum((n-mx)*(a-sum(a2s)/n_data) for n, a in zip(Ns, a2s)) / n_data
        overlap_coeff = -cov_na2 / var_x  # negative because anti-correlated
        print(f"\n  Overlap coefficient: d(alpha_2)/dN ~ {-overlap_coeff:.4f}")
        print(f"  Critical threshold: 1/2 = 0.5")
        print(f"  Since |overlap_coeff| = {abs(overlap_coeff):.4f} {'< 0.5: MORE cycles always helps' if abs(overlap_coeff) < 0.5 else '>= 0.5: trade-off ambiguous'}")

    # For p=7, H = 1 + 2N + 4*alpha_2 (no alpha_3 at p=7)
    # This is simpler: just 2 terms.
    if p == 7:
        print(f"\n  At p=7: H = 1 + 2N + 4*a2 (alpha_3 = 0)")
        print(f"  Net dH/dN = 2 - 4*(da2/dN)")
        # Direct computation of da2/dN from consecutive pairs
        sorted_results = sorted(results, key=lambda r: r['N'])
        for i in range(1, len(sorted_results)):
            prev = sorted_results[i-1]
            curr = sorted_results[i]
            if curr['N'] != prev['N']:
                dN = curr['N'] - prev['N']
                da2 = curr['alpha_2'] - prev['alpha_2']
                dH = curr['H'] - prev['H']
                print(f"    N: {prev['N']}->{curr['N']}: dN={dN}, da2={da2}, "
                      f"da2/dN={da2/dN:.3f}, dH={dH}, predicted_dH={2*dN+4*da2}")

print("\nDONE.")
