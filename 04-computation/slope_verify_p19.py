#!/usr/bin/env python3
"""
slope_verify_p19.py -- Verify binomial slope formula b_k = C(m-(k-1)/2, (k-1)/2) at p=19

Conjecture (HYP-524): For Interval S={1,...,m}, co_occ_k(d) = a_k + b_k * min(d, p-d)
where b_k = C(m - (k-1)/2, (k-1)/2) for k >= 5, and b_3 = 1 (anomalous).

Predictions for p=19, m=9:
  k=3:  b = 1 (always)
  k=5:  C(7, 2) = 21
  k=7:  C(6, 3) = 20
  k=9:  C(5, 4) = 5
  k=11: C(4, 5) = 0
  k=13: C(3, 6) = 0

Also tests p=17, m=8:
  k=3:  b = 1
  k=5:  C(6, 2) = 15
  k=7:  C(5, 3) = 10
  k=9:  C(4, 4) = 1

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations
from math import comb


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def has_ham_cycle(A, verts):
    """Check if vertex set supports at least one directed Ham cycle (DP)."""
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a] + A[a][c] * A[c][b] * A[b][a]) > 0

    dp = {}
    dp[(1 << 0, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    if nkey not in dp:
                        dp[nkey] = 1  # only need existence

    full = (1 << k) - 1
    for v in range(1, k):
        key = (full, v)
        if key in dp and A[verts[v]][verts[0]]:
            return True
    return False


def compute_co_occ_profile(A, p, k):
    """Compute co_occ_k(d) = #{k-cycle vertex sets containing both 0 and d}."""
    profile = [0] * p
    n_sets = 0
    n_0 = 0

    for subset in combinations(range(p), k):
        if has_ham_cycle(A, list(subset)):
            n_sets += 1
            if 0 in subset:
                n_0 += 1
                for v in subset:
                    if v != 0:
                        profile[v] += 1

    return profile, n_sets, n_0


def fit_linear(profile, p):
    """Fit co_occ(d) = a + b*d for d=1,...,m where m=(p-1)/2."""
    m = (p - 1) // 2
    vals = [profile[d] for d in range(1, m + 1)]
    n = len(vals)
    sum_x = sum(range(1, n + 1))
    sum_y = sum(vals)
    sum_xx = sum(d*d for d in range(1, n + 1))
    sum_xy = sum(d * vals[d-1] for d in range(1, n + 1))
    b = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x**2)
    a = (sum_y - b * sum_x) / n
    residuals = [abs(vals[d-1] - (a + b*d)) for d in range(1, n + 1)]
    max_res = max(residuals)
    return a, b, max_res


def main():
    print("=" * 70)
    print("BINOMIAL SLOPE VERIFICATION: b_k = C(m-(k-1)/2, (k-1)/2)")
    print("=" * 70)

    for p in [17, 19]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, Interval S={S_int}")
        print(f"{'='*70}")

        # Predictions
        predictions = {}
        for k in range(3, min(p, 14), 2):
            j = (k - 1) // 2
            if k == 3:
                predictions[k] = 1  # anomalous
            else:
                predictions[k] = comb(m - j, j)

        print(f"\nPredictions:")
        for k, pred in predictions.items():
            j = (k-1)//2
            if k == 3:
                print(f"  k={k}: b = 1 (THM-141, anomalous)")
            else:
                print(f"  k={k}: b = C({m-j},{j}) = {pred}")

        for k in range(3, min(p, 14), 2):
            t0 = time.time()
            n_subsets = comb(p, k)
            print(f"\n  k={k}: checking {n_subsets} subsets...", end=" ", flush=True)

            # Skip if too many subsets
            if n_subsets > 500000:
                print(f"SKIPPED (>{500000} subsets)")
                continue

            profile, n_sets, n_0 = compute_co_occ_profile(A, p, k)
            t1 = time.time()

            a, b, max_res = fit_linear(profile, p)
            is_linear = max_res < 0.01

            predicted = predictions[k]
            match = abs(b - predicted) < 0.01

            sym = all(profile[d] == profile[p-d] for d in range(1, (p+1)//2))

            print(f"{t1-t0:.1f}s")
            print(f"    {n_sets} vertex sets, {n_0} through 0")
            print(f"    Profile: {[profile[d] for d in range(1, m+1)]}")
            print(f"    Linear: {is_linear} (res={max_res:.6f})")
            print(f"    Slope = {b:.4f}, Intercept = {a:.4f}")
            print(f"    Predicted slope = {predicted}")
            print(f"    MATCH: {match}")

    # Also verify the intercept pattern
    print(f"\n{'='*70}")
    print("INTERCEPT ANALYSIS")
    print("=" * 70)

    all_data = {}
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        A = build_adj(p, S_int)
        all_data[p] = {}
        for k in range(3, min(p, 10), 2):
            profile, n_sets, n_0 = compute_co_occ_profile(A, p, k)
            a, b, max_res = fit_linear(profile, p)
            all_data[p][k] = {'a': round(a), 'b': round(b), 'n_sets': n_sets, 'n_0': n_0}

    print(f"\nIntercept table (a_k):")
    print(f"  {'k':>3} | {'p=7':>8} | {'p=11':>8} | {'p=13':>8}")
    print(f"  {'-'*3}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}")
    for k in [3, 5, 7, 9]:
        row = f"  {k:>3} |"
        for p in [7, 11, 13]:
            if k in all_data.get(p, {}):
                row += f" {all_data[p][k]['a']:>8} |"
            else:
                row += f" {'---':>8} |"
        print(row)

    # Check intercept formula candidates
    print(f"\nIntercept pattern search:")
    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")
        for k in sorted(all_data.get(p, {}).keys()):
            j = (k-1)//2
            a_val = all_data[p][k]['a']
            n_sets = all_data[p][k]['n_sets']
            n_0 = all_data[p][k]['n_0']
            # n_0 = C(p-1, k-1) * c_k(S, 0) / C(p, k) ... no, n_0 = n_sets * k / p
            ratio = n_0 * p / (n_sets * k) if n_sets > 0 else 0
            # Various formulas to check
            check_comb = comb(p-2, k-2)  # choosing k-2 from remaining p-2
            print(f"    k={k}: a={a_val}, n_sets={n_sets}, n_0={n_0}, "
                  f"C(p-2,k-2)={check_comb}, n_0/n_sets={n_0/n_sets:.4f}, k/p={k/p:.4f}")

    # k=3 special analysis: WHY is slope 1 instead of C(m-1,1)=m-1?
    print(f"\n{'='*70}")
    print("k=3 ANOMALY ANALYSIS")
    print("=" * 70)
    print("\nFor k=3, the formula C(m-1,1) = m-1 would predict slopes:")
    for p in [7, 11, 13, 17, 19]:
        m = (p-1)//2
        print(f"  p={p}: C({m-1},1) = {m-1}, but actual slope = 1")
    print("\nk=3 vertex sets = directed triangles.")
    print("co_occ_3(d) counts triangle vertex sets containing {0, d}.")
    print("The third vertex v must have 0->v->d->0 or 0->d->v->0 (one of these).")
    print("For Interval S={1,...,m}, 0->v requires v in S = {1,...,m}.")
    print("v->d requires d-v in S. d->0 requires p-d in S, i.e., d <= m.")
    print("If 0->d (i.e., d in S), then third vertex v needs:")
    print("  v in S AND (d-v) in S AND (p-d+v)...wait, not quite.")
    print("Actually for a 3-cycle {0,d,v}, exactly ONE cyclic ordering works.")
    print("THM-141 proof: co_occ_3(d) = |S intersect (d-S)| = min(d, p-d) for Interval.")
    print("\nFor k>=5, the 'third vertex' argument doesn't apply;")
    print("instead we're choosing k-2 more vertices, and the constraint is global (Ham cycle).")
    print("The binomial coefficient C(m-j, j) with j=(k-1)/2 counts something about")
    print("the number of ways to increase the feasible region as d grows.")


if __name__ == '__main__':
    main()
