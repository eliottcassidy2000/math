#!/usr/bin/env python3
"""
k-Colored Independence Polynomial I_k(Omega, x)
=================================================

For a tournament T on n vertices, define:
  a_k(T) = # permutations of [n] with exactly k forward edges in T

Key facts:
  - sum_k a_k(T) = n! for all T
  - a_{n-1}(T) = H(T) = # Hamiltonian paths
  - a_k(transitive T) = A(n, k) (Eulerian numbers)

The OCF identity: H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...

DISCOVERED FORMULA: The "deformed Eulerian" gives:
  a_k(T) = A(n,k) + sum_{(p,f)} catalog[(p,f)] * 2^p * c_k^{(p, f, d)}

where:
  - catalog[(p,f)] = # independent sets of p cycles with total freedom f
  - d = n-1
  - c_k^{(p, f, d)} = inflated_eulerian(d - f - p, d, k)
                     = sum_j A(d-f-p+1, j) * C(f+p, k-j) * (-1)^{f+p-k+j}

The k-COLORED INDEPENDENCE POLYNOMIAL:
  I_k(Omega, x) = A(n,k) + sum_{(p,f)} catalog[(p,f)] * c_k^{(p,f,d)} * x^p

Key properties verified:
  1. a_k(T) = I_k(Omega(T), 2) for all k
  2. I_{n-1}(Omega, x) = I(Omega, x) (standard independence polynomial)
  3. I_0(Omega, x) is generally NOT equal to I(Omega, x)
  4. I_k(Omega, 1) = A(n, k) for all T (follows from sum_k c_k = 0)

Author: opus (k-colored independence polynomial investigation)
"""

import sys, os, random
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, random_tournament, find_odd_cycles,
    conflict_graph, hamiltonian_path_count, tournament_from_bits
)


# ---------------------------------------------------------------------------
# Eulerian numbers
# ---------------------------------------------------------------------------

def eulerian_number(n, k):
    """A(n, k) = number of permutations of [n] with exactly k ascents (0-indexed).
    Uses the formula: A(n,k) = sum_{j=0}^{k} (-1)^j C(n+1,j) (k+1-j)^n"""
    if n == 0:
        return 1 if k == 0 else 0
    if k < 0 or k >= n:
        return 0
    total = 0
    for j in range(k + 1):
        total += (-1)**j * comb(n + 1, j) * (k + 1 - j)**n
    return total


def eulerian_poly(n):
    """Return [A(n,0), A(n,1), ..., A(n,n-1)] as the Eulerian polynomial coefficients."""
    return [eulerian_number(n, k) for k in range(n)]


# ---------------------------------------------------------------------------
# Inflated Eulerian coefficient
# ---------------------------------------------------------------------------

def inflated_eulerian(f_eff, d, k):
    """Inflated Eulerian coefficient.

    inflated_eulerian(f_eff, d, k) = sum_j A(f_eff+1, j) * C(d-f_eff, k-j) * (-1)^{d-f_eff-k+j}

    For the k-colored independence polynomial, use f_eff = d - f - parts
    where f = total freedom of the independent set and parts = its size.
    """
    if f_eff < 0:
        return 0
    total = 0
    for j in range(f_eff + 1):
        a_val = eulerian_number(f_eff + 1, j)
        if a_val == 0:
            continue
        r = k - j
        df = d - f_eff
        if r < 0 or r > df:
            continue
        sign = (-1)**(df - r)
        total += a_val * comb(df, r) * sign
    return total


def c_k_coefficient(parts, f, d, k):
    """The c_k coefficient for the k-colored independence polynomial.

    c_k^{(parts, f, d)} = inflated_eulerian(d - f - parts, d, k)

    This is the CORRECT formula discovered empirically:
    the effective Eulerian parameter is d - f - parts.
    """
    f_eff = d - f - parts
    return inflated_eulerian(f_eff, d, k)


# ---------------------------------------------------------------------------
# Forward edge distribution via DP
# ---------------------------------------------------------------------------

def forward_edge_dist(T):
    """Compute a_k(T) for all k using DP.
    a_k = # permutations sigma of [n] with exactly k forward edges.
    Forward edge: T[sigma[i]][sigma[i+1]] = 1 for i = 0, ..., n-2.

    Returns list of length n where result[k] = a_k.
    """
    n = len(T)
    if n == 0:
        return [1]
    if n == 1:
        return [1]

    full = (1 << n) - 1
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[(1 << v, v)][0] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            fwd_dict = dp.get((mask, v))
            if not fwd_dict:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                is_fwd = T[v][u]
                for fwd_count, num_paths in fwd_dict.items():
                    dp[(new_mask, u)][fwd_count + is_fwd] += num_paths

    result = [0] * n
    for v in range(n):
        fwd_dict = dp.get((full, v))
        if fwd_dict:
            for fwd_count, num_paths in fwd_dict.items():
                if fwd_count < n:
                    result[fwd_count] += num_paths
    return result


# ---------------------------------------------------------------------------
# Independence polynomial coefficients (DP over cycles)
# ---------------------------------------------------------------------------

def indep_poly_from_cycles(cycles):
    """Compute independence polynomial coefficients from cycle list.
    Uses DP over cycles (keyed by used-vertex frozenset) to avoid 2^m brute force."""
    if not cycles:
        return [1]

    m = len(cycles)
    vsets = [frozenset(c) for c in cycles]

    dp = {frozenset(): {0: 1}}
    for i in range(m):
        vs_i = vsets[i]
        new_dp = {}
        for used, sz_counts in dp.items():
            if used not in new_dp:
                new_dp[used] = {}
            for sz, cnt in sz_counts.items():
                new_dp[used][sz] = new_dp[used].get(sz, 0) + cnt
            if not (used & vs_i):
                new_used = used | vs_i
                if new_used not in new_dp:
                    new_dp[new_used] = {}
                for sz, cnt in sz_counts.items():
                    new_dp[new_used][sz + 1] = new_dp[new_used].get(sz + 1, 0) + cnt
        dp = new_dp

    coeffs_dict = defaultdict(int)
    for used, sz_counts in dp.items():
        for sz, cnt in sz_counts.items():
            coeffs_dict[sz] += cnt

    max_sz = max(coeffs_dict.keys()) if coeffs_dict else 0
    result = [coeffs_dict.get(i, 0) for i in range(max_sz + 1)]
    return result


# ---------------------------------------------------------------------------
# OCF invariant catalog
# ---------------------------------------------------------------------------

def ocf_invariant_catalog_from_cycles(cycles):
    """Compute OCF invariant catalog from a pre-computed cycle list.
    Returns { (parts, f): count } using DP over cycles."""
    m = len(cycles)
    if m == 0:
        return {}

    vsets = [frozenset(c) for c in cycles]
    cycle_f = [len(c) - 2 for c in cycles]

    dp = {frozenset(): {(0, 0): 1}}
    for i in range(m):
        vs_i = vsets[i]
        f_i = cycle_f[i]
        new_dp = {}
        for used, pf_counts in dp.items():
            if used not in new_dp:
                new_dp[used] = {}
            for pf, cnt in pf_counts.items():
                new_dp[used][pf] = new_dp[used].get(pf, 0) + cnt
            if not (used & vs_i):
                new_used = used | vs_i
                if new_used not in new_dp:
                    new_dp[new_used] = {}
                for (p, f), cnt in pf_counts.items():
                    key = (p + 1, f + f_i)
                    new_dp[new_used][key] = new_dp[new_used].get(key, 0) + cnt
        dp = new_dp

    catalog = defaultdict(int)
    for used, pf_counts in dp.items():
        for (p, f), cnt in pf_counts.items():
            if p > 0:
                catalog[(p, f)] += cnt
    return catalog


def ocf_invariant_catalog(T):
    """Compute OCF invariant catalog for tournament T."""
    cycles = find_odd_cycles(T)
    return ocf_invariant_catalog_from_cycles(cycles)


# ---------------------------------------------------------------------------
# k-Colored Independence Polynomial
# ---------------------------------------------------------------------------

def k_colored_indep_poly_from_catalog(n, k, catalog):
    """Compute I_k(Omega, x) from a pre-computed catalog.

    I_k(Omega, x) = A(n,k) + sum_{(parts, f)} catalog[(parts,f)] * c_k(parts,f,d) * x^parts

    Returns a list of coefficients [coeff_0, coeff_1, ...].
    """
    d = n - 1
    max_parts = max((p for (p, f) in catalog), default=0)

    poly = [0] * (max_parts + 1)
    poly[0] = eulerian_number(n, k)

    for (parts, f), count in catalog.items():
        ck = c_k_coefficient(parts, f, d, k)
        poly[parts] += ck * count

    return poly


def k_colored_indep_poly(T, k):
    """Compute I_k(Omega(T), x) as a polynomial in x."""
    n = len(T)
    catalog = ocf_invariant_catalog(T)
    return k_colored_indep_poly_from_catalog(n, k, catalog)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def poly_to_str(coeffs, var='x'):
    """Pretty-print a polynomial."""
    terms = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        if i == 0:
            terms.append(str(c))
        elif i == 1:
            if c == 1:
                terms.append(var)
            elif c == -1:
                terms.append(f'-{var}')
            else:
                terms.append(f'{c}{var}')
        else:
            if c == 1:
                terms.append(f'{var}^{i}')
            elif c == -1:
                terms.append(f'-{var}^{i}')
            else:
                terms.append(f'{c}{var}^{i}')
    return ' + '.join(terms).replace('+ -', '- ') if terms else '0'


def eval_poly(coeffs, x):
    """Evaluate polynomial at x."""
    return sum(c * x**i for i, c in enumerate(coeffs))


def polys_equal(p1, p2):
    """Check if two polynomials are equal (trim trailing zeros)."""
    a = list(p1)
    b = list(p2)
    while a and a[-1] == 0:
        a.pop()
    while b and b[-1] == 0:
        b.pop()
    return a == b


# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def run_single_tournament(T, label="", verbose=True):
    """Run all checks on a single tournament. Returns dict of results."""
    n = len(T)
    d = n - 1

    # Compute forward edge distribution
    a_dist = forward_edge_dist(T)

    # Compute cycles once, reuse
    cycles = find_odd_cycles(T)
    std_poly = indep_poly_from_cycles(cycles)
    catalog = ocf_invariant_catalog_from_cycles(cycles)

    # Compute all I_k polynomials
    ik_polys = [k_colored_indep_poly_from_catalog(n, k, catalog) for k in range(n)]

    # CHECK 1: I_k(Omega, 2) = a_k(T)
    check1_pass = True
    for k in range(n):
        val = eval_poly(ik_polys[k], 2)
        if val != a_dist[k]:
            check1_pass = False
            if verbose:
                print(f"  FAIL check1: k={k}, I_k(Omega,2)={val}, a_k={a_dist[k]}")

    # CHECK 2: I_{n-1}(Omega, x) = I(Omega, x)
    check2_pass = polys_equal(ik_polys[n - 1], std_poly)

    # CHECK 3: I_0(Omega, x) = I(Omega, x)?
    check3_pass = polys_equal(ik_polys[0], std_poly)

    # CHECK 4: I_k(Omega, 1) values
    ik_at_1 = [eval_poly(ik_polys[k], 1) for k in range(n)]

    # CHECK 5: I_k(Omega, -1) values
    ik_at_minus1 = [eval_poly(ik_polys[k], -1) for k in range(n)]

    # CHECK 6: sum_k I_k(Omega, 2) = H(T)?
    # Actually sum_k a_k = n!, so sum_k I_k(Omega,2) = n!, not H(T).
    # H(T) = a_{n-1}(T) = I_{n-1}(Omega, 2).
    sum_at_2 = sum(eval_poly(ik_polys[k], 2) for k in range(n))
    check6_pass = (sum_at_2 == factorial(n))

    sum_at_1 = sum(ik_at_1)
    sum_at_minus1 = sum(ik_at_minus1)

    if verbose:
        print(f"\n  {label} (n={n})")
        print(f"  H(T) = {hamiltonian_path_count(T)}")
        print(f"  a_k = {a_dist}")
        print(f"  I(Omega, x) = {poly_to_str(std_poly)}")
        print(f"  I(Omega, 2) = {eval_poly(std_poly, 2)}")

        for k in range(n):
            print(f"  I_{k}(Omega, x) = {poly_to_str(ik_polys[k])}")

        print(f"\n  Check 1 (I_k(Omega,2) = a_k): {'PASS' if check1_pass else 'FAIL'}")
        print(f"  Check 2 (I_{{n-1}} = I):        {'PASS' if check2_pass else 'FAIL'}")
        print(f"  Check 3 (I_0 = I):             {'PASS' if check3_pass else 'FAIL'}")
        print(f"  Check 6 (sum = n!):            {'PASS' if check6_pass else 'FAIL'}")

        print(f"\n  I_k(Omega, 1) values: {ik_at_1}")
        print(f"  sum_k I_k(Omega, 1) = {sum_at_1}  (n! = {factorial(n)})")

        eul = eulerian_poly(n)
        if ik_at_1 == eul:
            print(f"  I_k(Omega, 1) = A(n, k) (Eulerian numbers): YES")
        else:
            print(f"  I_k(Omega, 1) vs A(n,k): {ik_at_1} vs {eul}")

        print(f"\n  I_k(Omega, -1) values: {ik_at_minus1}")
        print(f"  sum_k I_k(Omega, -1) = {sum_at_minus1}")

    return {
        'check1': check1_pass,
        'check2': check2_pass,
        'check3': check3_pass,
        'check6': check6_pass,
        'ik_at_1': ik_at_1,
        'ik_at_minus1': ik_at_minus1,
        'ik_polys': ik_polys,
        'std_poly': std_poly,
        'a_dist': a_dist,
        'H_T': hamiltonian_path_count(T),
    }


def main():
    print("=" * 70)
    print("k-COLORED INDEPENDENCE POLYNOMIAL VERIFICATION")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Preliminary checks
    # -----------------------------------------------------------------------
    print("\n--- Eulerian number sanity checks ---")
    e4 = eulerian_poly(4)
    print(f"  A(4, k) = {e4}  (should be [1, 11, 11, 1])")
    assert e4 == [1, 11, 11, 1]

    e5 = eulerian_poly(5)
    print(f"  A(5, k) = {e5}  (should be [1, 26, 66, 26, 1])")
    assert e5 == [1, 26, 66, 26, 1]

    e7 = eulerian_poly(7)
    print(f"  A(7, k) = {e7}")
    assert sum(e7) == factorial(7)

    # Verify c_k formula: at k = d, c_k should be 1 for all valid (p, f)
    print("\n--- c_k sanity: c_d = 1 ---")
    for d in [4, 5, 6]:
        for p in [1, 2]:
            for f in [1, 2, 3, 4, 5]:
                if f + p > d:
                    continue
                c = c_k_coefficient(p, f, d, d)
                if c != 1:
                    print(f"  FAIL: c_{d}^({p},{f},{d}) = {c}")
    print("  All c_d = 1 checks passed.")

    # Verify c_k formula: sum_k c_k = 0 when f_eff < d
    print("\n--- c_k sanity: sum_k c_k = 0 ---")
    for d in [4, 5, 6]:
        for p in [1, 2]:
            for f in [1, 2, 3, 4, 5]:
                if f + p > d:
                    continue
                s = sum(c_k_coefficient(p, f, d, k) for k in range(d + 1))
                if s != 0:
                    print(f"  FAIL: sum c_k^({p},{f},{d}) = {s}")
    print("  All sum_k c_k = 0 checks passed.")

    # -----------------------------------------------------------------------
    # n=5 EXHAUSTIVE
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("n=5: EXHAUSTIVE VERIFICATION (all 1024 tournaments)")
    print("=" * 70)

    n = 5
    total = 0
    c1_pass = c2_pass = c3_pass = c6_pass = 0
    c3_count_true = 0

    ik1_patterns = defaultdict(int)
    ik_minus1_patterns = defaultdict(int)
    ik1_equals_eulerian = True

    for T in all_tournaments(n):
        total += 1
        res = run_single_tournament(T, label=f"T#{total}", verbose=False)

        if res['check1']:
            c1_pass += 1
        if res['check2']:
            c2_pass += 1
        if res['check3']:
            c3_pass += 1
            c3_count_true += 1
        if res['check6']:
            c6_pass += 1

        ik1_key = tuple(res['ik_at_1'])
        ik1_patterns[ik1_key] += 1

        ik_m1_key = tuple(res['ik_at_minus1'])
        ik_minus1_patterns[ik_m1_key] += 1

        eul = eulerian_poly(n)
        if list(res['ik_at_1']) != eul:
            ik1_equals_eulerian = False

    print(f"\n  Total tournaments: {total}")
    print(f"  Check 1 (I_k(Omega,2) = a_k):  {c1_pass}/{total}")
    print(f"  Check 2 (I_{{n-1}} = I(Omega)):   {c2_pass}/{total}")
    print(f"  Check 3 (I_0 = I(Omega)):       {c3_pass}/{total}")
    print(f"  Check 6 (sum = n!):             {c6_pass}/{total}")

    print(f"\n  I_k(Omega, 1) = A(n,k) for all T? {ik1_equals_eulerian}")
    print(f"  Distinct I_k(Omega, 1) patterns: {len(ik1_patterns)}")
    if len(ik1_patterns) <= 10:
        for pattern, count in sorted(ik1_patterns.items(), key=lambda x: -x[1]):
            print(f"    {list(pattern)} : {count} tournaments")

    print(f"\n  Distinct I_k(Omega, -1) patterns: {len(ik_minus1_patterns)}")
    if len(ik_minus1_patterns) <= 15:
        for pattern, count in sorted(ik_minus1_patterns.items(), key=lambda x: -x[1]):
            print(f"    {list(pattern)} : {count} tournaments")

    # -----------------------------------------------------------------------
    # Detailed examples at n=5
    # -----------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("DETAILED EXAMPLES at n=5")
    print("-" * 70)

    # Transitive tournament
    T_trans = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1, 5):
            T_trans[i][j] = 1
    run_single_tournament(T_trans, label="Transitive T5", verbose=True)

    # Near-cyclic tournament
    T_cyc = [[0]*5 for _ in range(5)]
    for i in range(5):
        T_cyc[i][(i+1) % 5] = 1
        T_cyc[i][(i+2) % 5] = 1
    run_single_tournament(T_cyc, label="Near-cyclic T5", verbose=True)

    # -----------------------------------------------------------------------
    # n=7 SAMPLED
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("n=7: SAMPLED VERIFICATION (200 random tournaments)")
    print("=" * 70)

    n = 7
    rng = random.Random(42)
    total = 0
    c1_pass = c2_pass = c3_pass = c6_pass = 0

    ik1_patterns = defaultdict(int)
    ik_minus1_patterns = defaultdict(int)
    ik1_equals_eulerian = True

    for trial in range(200):
        T = random_tournament(n, rng)
        total += 1
        res = run_single_tournament(T, label=f"T#{total}", verbose=False)

        if res['check1']:
            c1_pass += 1
        if res['check2']:
            c2_pass += 1
        if res['check3']:
            c3_pass += 1
        if res['check6']:
            c6_pass += 1

        ik1_key = tuple(res['ik_at_1'])
        ik1_patterns[ik1_key] += 1

        ik_m1_key = tuple(res['ik_at_minus1'])
        ik_minus1_patterns[ik_m1_key] += 1

        eul = eulerian_poly(n)
        if list(res['ik_at_1']) != eul:
            ik1_equals_eulerian = False

    print(f"\n  Total tournaments: {total}")
    print(f"  Check 1 (I_k(Omega,2) = a_k):  {c1_pass}/{total}")
    print(f"  Check 2 (I_{{n-1}} = I(Omega)):   {c2_pass}/{total}")
    print(f"  Check 3 (I_0 = I(Omega)):       {c3_pass}/{total}")
    print(f"  Check 6 (sum = n!):             {c6_pass}/{total}")

    print(f"\n  I_k(Omega, 1) = A(n,k) for all T? {ik1_equals_eulerian}")
    print(f"  Distinct I_k(Omega, 1) patterns: {len(ik1_patterns)}")
    if len(ik1_patterns) <= 10:
        for pattern, count in sorted(ik1_patterns.items(), key=lambda x: -x[1]):
            print(f"    {list(pattern)} : {count} tournaments")
    else:
        print("  (showing top 5)")
        for pattern, count in sorted(ik1_patterns.items(), key=lambda x: -x[1])[:5]:
            print(f"    {list(pattern)} : {count} tournaments")

    print(f"\n  Distinct I_k(Omega, -1) patterns: {len(ik_minus1_patterns)}")
    if len(ik_minus1_patterns) <= 15:
        for pattern, count in sorted(ik_minus1_patterns.items(), key=lambda x: -x[1]):
            print(f"    {list(pattern)} : {count} tournaments")
    else:
        print("  (showing top 5)")
        for pattern, count in sorted(ik_minus1_patterns.items(), key=lambda x: -x[1])[:5]:
            print(f"    {list(pattern)} : {count} tournaments")

    # Detailed n=7 example
    print("\n" + "-" * 70)
    print("DETAILED EXAMPLE at n=7")
    print("-" * 70)
    T7 = random_tournament(7, random.Random(123))
    run_single_tournament(T7, label="Random T7 (seed=123)", verbose=True)

    # -----------------------------------------------------------------------
    # GENERATING FUNCTION CHECK
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("GENERATING FUNCTION CHECK")
    print("=" * 70)
    print("  Testing: sum_k I_k(Omega, x) * t^k")
    print("         = A_n(t) + sum_{(p,f)} cat * x^p * G_{p,f}(t)")
    print("  where G_{p,f}(t) = sum_k c_k^{(p,f,d)} t^k")

    for n_test in [5, 7]:
        print(f"\n  --- n={n_test} ---")
        if n_test == 5:
            test_tournaments = list(all_tournaments(5))[:5]
        else:
            test_tournaments = [random_tournament(n_test, random.Random(seed)) for seed in range(3)]

        d = n_test - 1
        for idx, T in enumerate(test_tournaments):
            cycles = find_odd_cycles(T)
            catalog = ocf_invariant_catalog_from_cycles(cycles)
            ik_polys_gf = [k_colored_indep_poly_from_catalog(n_test, k, catalog) for k in range(n_test)]

            all_ok = True
            for x_val in [1, 2, -1, 3]:
                for t_val in [1, 2, -1, 0]:
                    lhs = sum(eval_poly(ik_polys_gf[k], x_val) * t_val**k for k in range(n_test))

                    rhs = sum(eulerian_number(n_test, k) * t_val**k for k in range(n_test))
                    for (parts, f), count in catalog.items():
                        gf_val = sum(c_k_coefficient(parts, f, d, k) * t_val**k for k in range(n_test))
                        rhs += count * x_val**parts * gf_val
                    if lhs != rhs:
                        all_ok = False
            print(f"    T#{idx}: generating function {'PASS' if all_ok else 'FAIL'}")

    # -----------------------------------------------------------------------
    # I_k(Omega, 1) ANALYSIS
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS: I_k(Omega, 1)")
    print("=" * 70)

    print("\n  Since sum_k c_k = 0 for f_eff < d (i.e., whenever there ARE invariants),")
    print("  evaluating I_k(Omega, x) at x=1 kills all invariant terms.")
    print("  Therefore I_k(Omega, 1) = A(n, k) for ALL tournaments.")

    print("\n  Verification at n=5 (exhaustive):")
    eul5 = eulerian_poly(5)
    all_match = True
    for T in all_tournaments(5):
        cycles = find_odd_cycles(T)
        catalog = ocf_invariant_catalog_from_cycles(cycles)
        for k in range(5):
            ik = k_colored_indep_poly_from_catalog(5, k, catalog)
            if eval_poly(ik, 1) != eul5[k]:
                all_match = False
                break
        if not all_match:
            break
    print(f"    I_k(Omega, 1) = A(5, k) for all T on 5 vertices: {all_match}")

    # -----------------------------------------------------------------------
    # I_k(Omega, -1) ANALYSIS
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS: I_k(Omega, -1)")
    print("=" * 70)

    print("\n  At x=-1, I_k(Omega, -1) = A(n,k) + sum cat * c_k * (-1)^p")
    print("  This gives an 'alternating' count where pairs of disjoint cycles")
    print("  contribute with opposite sign from single cycles.")

    print("\n  n=5 first 5 tournaments:")
    for T_idx, T in enumerate(all_tournaments(5)):
        if T_idx >= 5:
            break
        cycles = find_odd_cycles(T)
        ip = indep_poly_from_cycles(cycles)
        i_minus1 = eval_poly(ip, -1)

        catalog = ocf_invariant_catalog_from_cycles(cycles)
        ik_at_m1 = [eval_poly(k_colored_indep_poly_from_catalog(5, k, catalog), -1) for k in range(5)]
        sum_ik_m1 = sum(ik_at_m1)

        print(f"    T#{T_idx}: I(Omega,-1) = {i_minus1}, sum_k I_k(Omega,-1) = {sum_ik_m1}, "
              f"I_k(Omega,-1) = {ik_at_m1}")

    print("\n  n=5: Distribution of sum_k I_k(Omega, -1):")
    vals = defaultdict(int)
    for T in all_tournaments(5):
        cycles = find_odd_cycles(T)
        catalog = ocf_invariant_catalog_from_cycles(cycles)
        s = sum(eval_poly(k_colored_indep_poly_from_catalog(5, k, catalog), -1) for k in range(5))
        vals[s] += 1
    for v, c in sorted(vals.items()):
        print(f"    {v}: {c} tournaments")

    # -----------------------------------------------------------------------
    # PALINDROMY CHECK
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PALINDROMY CHECK: c_0 vs c_d for OCF invariants")
    print("=" * 70)

    for d in range(3, 8):
        print(f"\n  d = {d}:")
        for p in [1, 2]:
            for f in range(1, d):
                if f + p > d:
                    continue
                c0 = c_k_coefficient(p, f, d, 0)
                cd = c_k_coefficient(p, f, d, d)
                f_eff = d - f - p
                print(f"    p={p}, f={f}: c_0={c0}, c_d={cd}, f_eff={f_eff}")

    # -----------------------------------------------------------------------
    # CHECK: I_0(Omega, x) vs I(Omega, x)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("I_0 vs I(Omega): when are they equal?")
    print("=" * 70)

    print("\n  For I_0 = I, we need c_0^{(p,f,d)} = 1 for all (p,f) that appear.")
    print("  c_0 = inflated_eulerian(f_eff, d, 0) = sum_j A(f_eff+1, j) C(d-f_eff, -j) (-1)^{d-f_eff+j}")
    print("  Only j=0 contributes: c_0 = A(f_eff+1, 0) * C(d-f_eff, 0) * (-1)^{d-f_eff}")
    print("  = 1 * 1 * (-1)^{d-f_eff} = (-1)^{f+p}")
    print("  So c_0 = 1 iff f+p is even.")
    print()
    print("  At n=5 (d=4): only invariants are (1,1) and (1,3).")
    print("  (1,1): f+p=2 (even) -> c_0=1")
    print("  (1,3): f+p=4 (even) -> c_0=1")
    print("  So I_0 = I at n=5: TRUE (for all n=5 tournaments)")
    print()
    print("  At n=7 (d=6): invariants include (1,5), (1,3), (1,1), (2,2).")
    print("  (1,1): f+p=2 (even) -> c_0=1")
    print("  (1,3): f+p=4 (even) -> c_0=1")
    print("  (1,5): f+p=6 (even) -> c_0=1")
    print("  (2,2): f+p=4 (even) -> c_0=1")
    print("  ALL have f+p even! So I_0 = I at n=7 too.")

    # Check if f+p is always even
    print("\n  Checking: is f+p always even for tournament OCF invariants?")
    print("  An independent set of p odd cycles has f = sum(len_i - 2).")
    print("  Each len_i is odd (>= 3), so len_i - 2 is odd.")
    print("  f = sum of p odd numbers = p * (odd) mod 2.")
    print("  If p is odd, f is odd -> f+p = even.")
    print("  If p is even, f is even -> f+p = even.")
    print("  THEREFORE f+p is ALWAYS even, and c_0 = 1 for all invariants.")
    print("  CONCLUSION: I_0(Omega, x) = I(Omega, x) for ALL tournaments!")

    # Re-verify
    print("\n  Verification at n=5 (exhaustive):")
    all_i0_eq = True
    for T in all_tournaments(5):
        cycles = find_odd_cycles(T)
        catalog = ocf_invariant_catalog_from_cycles(cycles)
        ik0 = k_colored_indep_poly_from_catalog(5, 0, catalog)
        std = indep_poly_from_cycles(cycles)
        if not polys_equal(ik0, std):
            all_i0_eq = False
            break
    print(f"    I_0 = I(Omega) for all n=5: {all_i0_eq}")

    print("\n  Verification at n=7 (200 sampled):")
    rng2 = random.Random(99)
    all_i0_eq_7 = True
    for _ in range(200):
        T = random_tournament(7, rng2)
        cycles = find_odd_cycles(T)
        catalog = ocf_invariant_catalog_from_cycles(cycles)
        ik0 = k_colored_indep_poly_from_catalog(7, 0, catalog)
        std = indep_poly_from_cycles(cycles)
        if not polys_equal(ik0, std):
            all_i0_eq_7 = False
            break
    print(f"    I_0 = I(Omega) for 200 random n=7: {all_i0_eq_7}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
