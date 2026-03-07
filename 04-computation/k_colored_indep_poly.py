#!/usr/bin/env python3
"""
k-Colored Independence Polynomial I_k(Omega, x)
=================================================

For a tournament T on n vertices, define:
  a_k(T) = # permutations of [n] with exactly k forward edges in T

The OCF identity says H(T) = I(Omega(T), 2) = sum_k a_k(T).

The "deformed Eulerian" formula (proved by Grinberg-Stanley) gives:
  a_k(T) = A(n,k) + sum_I  2^{parts(I)} * c_k^{(f_I, n-1)} * inv_I(T)

where the sum is over OCF invariants I, each with:
  - f_I = degrees of freedom (related to cycle lengths)
  - parts(I) = number of independent cycles in I
  - inv_I(T) = the invariant value for tournament T
  - c_k^{(f, d)} = sum_j A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}

NEW CONCEPT: Define the k-colored independence polynomial:
  I_k(Omega, x) = A(n,k) + sum_I c_k^{(f_I, n-1)} * x^{parts(I)} * inv_I(T)

KEY PREDICTIONS:
  1. a_k(T) = I_k(Omega(T), 2) for all k
  2. I_{n-1}(Omega, x) = I(Omega, x) = standard independence polynomial
  3. I_0(Omega, x) = I(Omega, x) (palindromy)
  4. Generating function: sum_k I_k(Omega, x) t^k = A_n(t) + sum_I x^parts * inv_I(T) * A_{f+1}(t) * (t-1)^{d-f}

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
# Inflated Eulerian coefficient c_k^{(f, d)}
# ---------------------------------------------------------------------------

def inflated_eulerian(f, d, k):
    """c_k^{(f, d)} = sum_{j} A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}

    This is the coefficient that appears in the deformed Eulerian formula.
    f = degrees of freedom, d = n-1 (total dimension).
    k = target number of forward edges (0-indexed: k forward edges among n-1 pairs).
    """
    total = 0
    for j in range(f + 1):  # A(f+1, j) is nonzero for 0 <= j <= f
        a_val = eulerian_number(f + 1, j)
        if a_val == 0:
            continue
        r = k - j
        df = d - f
        if r < 0 or r > df:
            continue
        sign = (-1)**(df - r)
        total += a_val * comb(df, r) * sign
    return total


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

    # dp[mask][v][fwd] = # paths through vertices in mask, ending at v, with fwd forward edges
    # This is expensive but correct
    full = (1 << n) - 1
    # Use dict for sparse storage: dp[(mask, v)] = {fwd_count: number_of_paths}
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
                is_fwd = T[v][u]  # 1 if forward edge v->u
                for fwd_count, num_paths in fwd_dict.items():
                    dp[(new_mask, u)][fwd_count + is_fwd] += num_paths

    result = [0] * n  # k ranges from 0 to n-1
    for v in range(n):
        fwd_dict = dp.get((full, v))
        if fwd_dict:
            for fwd_count, num_paths in fwd_dict.items():
                if fwd_count < n:
                    result[fwd_count] += num_paths
    return result


# ---------------------------------------------------------------------------
# Independence polynomial coefficients (from conflict graph)
# ---------------------------------------------------------------------------

def indep_poly_coefficients(adj):
    """Compute [alpha_0, alpha_1, ..., alpha_d] of I(G, x)."""
    m = len(adj)
    if m == 0:
        return [1]

    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j

    counts = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            counts[bin(mask).count('1')] += 1

    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


def indep_poly_from_cycles(cycles):
    """Compute independence polynomial coefficients from cycle list.
    Uses DP over cycles (keyed by used-vertex frozenset) to avoid 2^m brute force."""
    if not cycles:
        return [1]

    m = len(cycles)
    vsets = [frozenset(c) for c in cycles]

    # DP: dp[used_verts] = {size: count}
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

    # Aggregate over all used-vertex sets
    coeffs_dict = defaultdict(int)
    for used, sz_counts in dp.items():
        for sz, cnt in sz_counts.items():
            coeffs_dict[sz] += cnt

    max_sz = max(coeffs_dict.keys()) if coeffs_dict else 0
    result = [coeffs_dict.get(i, 0) for i in range(max_sz + 1)]
    return result


# ---------------------------------------------------------------------------
# OCF invariant catalog for a tournament
# ---------------------------------------------------------------------------

def ocf_invariant_catalog(T):
    """Compute all OCF invariants for tournament T.
    Wrapper that finds cycles and delegates to ocf_invariant_catalog_from_cycles."""
    cycles = find_odd_cycles(T)
    return ocf_invariant_catalog_from_cycles(cycles)


# ---------------------------------------------------------------------------
# k-Colored Independence Polynomial
# ---------------------------------------------------------------------------

def ocf_invariant_catalog_from_cycles(cycles):
    """Compute OCF invariant catalog from a pre-computed cycle list.
    Returns { (parts, f): count } using DP."""
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


def k_colored_indep_poly_from_catalog(n, k, catalog):
    """Compute I_k(Omega, x) from a pre-computed catalog.

    I_k(Omega, x) = A(n,k) + sum_{(parts, f)} c_k^{(f, n-1)} * catalog[(parts,f)] * x^parts

    Returns a list of coefficients [coeff_0, coeff_1, ...].
    """
    d = n - 1
    max_parts = max((p for (p, f) in catalog), default=0)

    poly = [0] * (max_parts + 1)
    poly[0] = eulerian_number(n, k)

    for (parts, f), count in catalog.items():
        ck = inflated_eulerian(f, d, k)
        poly[parts] += ck * count

    return poly


def k_colored_indep_poly(T, k):
    """Compute I_k(Omega(T), x) as a polynomial in x.

    I_k(Omega, x) = A(n,k) + sum_{(parts, f)} c_k^{(f, n-1)} * catalog[(parts,f)] * x^parts

    Returns a list of coefficients [coeff_0, coeff_1, ...] where
    result[p] is the coefficient of x^p.
    """
    n = len(T)
    catalog = ocf_invariant_catalog(T)
    return k_colored_indep_poly_from_catalog(n, k, catalog)


def standard_indep_poly(T):
    """Compute I(Omega(T), x) as polynomial coefficients."""
    cycles = find_odd_cycles(T)
    return indep_poly_from_cycles(cycles)


# ---------------------------------------------------------------------------
# Main verification
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
    while len(a) > 0 and a[-1] == 0:
        a.pop()
    while len(b) > 0 and b[-1] == 0:
        b.pop()
    return a == b


def run_single_tournament(T, label="", verbose=True):
    """Run all checks on a single tournament. Returns dict of results."""
    n = len(T)
    d = n - 1

    # Compute forward edge distribution
    a_dist = forward_edge_dist(T)

    # Compute cycles once, reuse for both standard and k-colored
    cycles = find_odd_cycles(T)
    std_poly = indep_poly_from_cycles(cycles)

    # Compute catalog once, reuse for all k
    catalog = ocf_invariant_catalog_from_cycles(cycles)

    # Compute all I_k polynomials using shared catalog
    ik_polys = []
    for k in range(n):
        ik = k_colored_indep_poly_from_catalog(n, k, catalog)
        ik_polys.append(ik)

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

    # CHECK 3: I_0(Omega, x) = I(Omega, x)
    check3_pass = polys_equal(ik_polys[0], std_poly)

    # CHECK 4: I_k(Omega, 1) values
    ik_at_1 = [eval_poly(ik_polys[k], 1) for k in range(n)]

    # CHECK 5: I_k(Omega, -1) values
    ik_at_minus1 = [eval_poly(ik_polys[k], -1) for k in range(n)]

    # CHECK 6: sum_k I_k(Omega, x) = sum_k a_k = H(T) when evaluated at x=2
    sum_at_2 = sum(eval_poly(ik_polys[k], 2) for k in range(n))
    H_T = hamiltonian_path_count(T)
    check6_pass = (sum_at_2 == H_T)

    # Check: sum_k I_k(Omega, 1) = ?
    sum_at_1 = sum(ik_at_1)

    # Check: sum_k I_k(Omega, -1) = ?
    sum_at_minus1 = sum(ik_at_minus1)

    if verbose:
        print(f"\n  {label} (n={n})")
        print(f"  H(T) = {H_T}")
        print(f"  a_k = {a_dist}")
        print(f"  I(Omega, x) = {poly_to_str(std_poly)}")
        print(f"  I(Omega, 2) = {eval_poly(std_poly, 2)}")

        for k in range(n):
            print(f"  I_{k}(Omega, x) = {poly_to_str(ik_polys[k])}")

        print(f"\n  Check 1 (I_k(Omega,2) = a_k): {'PASS' if check1_pass else 'FAIL'}")
        print(f"  Check 2 (I_{{n-1}} = I):        {'PASS' if check2_pass else 'FAIL'}")
        print(f"  Check 3 (I_0 = I):             {'PASS' if check3_pass else 'FAIL'}")
        print(f"  Check 6 (sum consistency):     {'PASS' if check6_pass else 'FAIL'}")

        print(f"\n  I_k(Omega, 1) values: {ik_at_1}")
        print(f"  sum_k I_k(Omega, 1) = {sum_at_1}  (n! = {factorial(n)})")

        print(f"\n  I_k(Omega, -1) values: {ik_at_minus1}")
        print(f"  sum_k I_k(Omega, -1) = {sum_at_minus1}")

        # Check if I_k(Omega, 1) = A(n, k) (Eulerian numbers)
        eul = eulerian_poly(n)
        if ik_at_1 == eul:
            print(f"  I_k(Omega, 1) = A(n, k) (Eulerian numbers): YES")
        else:
            print(f"  I_k(Omega, 1) vs A(n,k): {ik_at_1} vs {eul}")

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
        'H_T': H_T,
    }


def main():
    print("=" * 70)
    print("k-COLORED INDEPENDENCE POLYNOMIAL VERIFICATION")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Preliminary: verify Eulerian / inflated Eulerian basics
    # -----------------------------------------------------------------------
    print("\n--- Eulerian number sanity checks ---")
    # A(4, k) = [1, 11, 11, 1]
    e4 = eulerian_poly(4)
    print(f"  A(4, k) = {e4}  (should be [1, 11, 11, 1])")
    assert e4 == [1, 11, 11, 1]

    e5 = eulerian_poly(5)
    print(f"  A(5, k) = {e5}  (should be [1, 26, 66, 26, 1])")
    assert e5 == [1, 26, 66, 26, 1]

    e7 = eulerian_poly(7)
    print(f"  A(7, k) = {e7}")
    assert sum(e7) == factorial(7)

    # Check inflated Eulerian: c_k^{(f, d)} for some known values
    # When f = d (all free), c_k^{(d, d)} = A(d+1, k) (Eulerian of n=d+1)
    # Actually c_k^{(d, d)} = sum_j A(d+1, j) C(0, k-j) (-1)^{-k+j}
    #                        = A(d+1, k) * C(0, 0) * (-1)^0 = A(d+1, k)
    for d in [3, 4, 5, 6]:
        for k in range(d + 1):
            assert inflated_eulerian(d, d, k) == eulerian_number(d + 1, k), \
                f"inflated_eulerian({d},{d},{k}) failed"
    print("  Inflated Eulerian c_k^(d,d) = A(d+1,k): PASS")

    # When f = 1 (single 3-cycle), d = n-1:
    # c_k^{(1, d)} = sum_j A(2, j) C(d-1, k-j) (-1)^{d-1-k+j}
    # A(2, 0) = 1, A(2, 1) = 1
    # = C(d-1, k) (-1)^{d-1-k} + C(d-1, k-1) (-1)^{d-k}
    # = (-1)^{d-1-k} [C(d-1,k) - C(d-1,k-1)]
    print("  c_k^{(1, d)} check:")
    for d in [4, 5, 6]:
        vals = [inflated_eulerian(1, d, k) for k in range(d + 1)]
        print(f"    d={d}: {vals}")
        # Verify sum = 0 (since d-f > 0, the signed binomial sum vanishes)
        # Actually sum_k c_k^{(f,d)} = sum_k [what we get evaluating A_{f+1}(t)(t-1)^{d-f} at t=1]
        # = A_{f+1}(1) * 0^{d-f} = 0 when d > f
        s = sum(vals)
        print(f"    sum = {s}  (should be 0 when d > f={1})")

    # Key property: c_0^{(f,d)} and c_{d}^{(f,d)} when d-f is even
    print("\n  c_0 and c_d values (palindromy check):")
    for f in [1, 3, 5]:
        for d in [4, 5, 6]:
            if f > d:
                continue
            c0 = inflated_eulerian(f, d, 0)
            cd = inflated_eulerian(f, d, d)
            parity_df = (d - f) % 2
            print(f"    f={f}, d={d}: c_0={c0}, c_d={cd}, d-f={d-f} (parity={parity_df})")

    # -----------------------------------------------------------------------
    # n=5 EXHAUSTIVE
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("n=5: EXHAUSTIVE VERIFICATION (all 1024 tournaments)")
    print("=" * 70)

    n = 5
    total = 0
    c1_pass = c2_pass = c3_pass = c6_pass = 0
    c2_fail_examples = []
    c3_fail_examples = []

    # Track I_k(Omega, 1) patterns
    ik1_patterns = defaultdict(int)
    ik_minus1_patterns = defaultdict(int)

    # Track whether I_k(Omega, 1) = A(n, k) always
    ik1_equals_eulerian = True

    for T in all_tournaments(n):
        total += 1
        res = run_single_tournament(T, label=f"T#{total}", verbose=False)

        if res['check1']:
            c1_pass += 1
        if res['check2']:
            c2_pass += 1
        else:
            if len(c2_fail_examples) < 3:
                c2_fail_examples.append((total, res))
        if res['check3']:
            c3_pass += 1
        else:
            if len(c3_fail_examples) < 3:
                c3_fail_examples.append((total, res))
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
    print(f"  Check 6 (sum consistency):      {c6_pass}/{total}")

    if c2_fail_examples:
        print(f"\n  I_{{n-1}} != I(Omega) examples:")
        for idx, res in c2_fail_examples:
            print(f"    T#{idx}: I_{{n-1}} = {poly_to_str(res['ik_polys'][n-1])}")
            print(f"           I(Omega) = {poly_to_str(res['std_poly'])}")

    if c3_fail_examples:
        print(f"\n  I_0 != I(Omega) examples:")
        for idx, res in c3_fail_examples:
            print(f"    T#{idx}: I_0 = {poly_to_str(res['ik_polys'][0])}")
            print(f"           I(Omega) = {poly_to_str(res['std_poly'])}")

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
    # Show detailed output for a few interesting n=5 tournaments
    # -----------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("DETAILED EXAMPLES at n=5")
    print("-" * 70)

    # Transitive tournament
    T_trans = tournament_from_bits(5, (1 << 10) - 1)  # all arcs i->j for i<j
    run_single_tournament(T_trans, label="Transitive T5", verbose=True)

    # A tournament with many cycles
    # Try the cyclic tournament C5: 0->1->2->3->4->0
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
    c2_fail_examples = []
    c3_fail_examples = []

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
        else:
            if len(c2_fail_examples) < 3:
                c2_fail_examples.append((total, res))
        if res['check3']:
            c3_pass += 1
        else:
            if len(c3_fail_examples) < 3:
                c3_fail_examples.append((total, res))
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
    print(f"  Check 6 (sum consistency):      {c6_pass}/{total}")

    if c2_fail_examples:
        print(f"\n  I_{{n-1}} != I(Omega) examples:")
        for idx, res in c2_fail_examples:
            print(f"    T#{idx}: I_{{n-1}} = {poly_to_str(res['ik_polys'][n-1])}")
            print(f"           I(Omega) = {poly_to_str(res['std_poly'])}")

    if c3_fail_examples:
        print(f"\n  I_0 != I(Omega) examples:")
        for idx, res in c3_fail_examples:
            print(f"    T#{idx}: I_0 = {poly_to_str(res['ik_polys'][0])}")
            print(f"           I(Omega) = {poly_to_str(res['std_poly'])}")

    print(f"\n  I_k(Omega, 1) = A(n,k) for all T? {ik1_equals_eulerian}")
    print(f"  Distinct I_k(Omega, 1) patterns: {len(ik1_patterns)}")
    if len(ik1_patterns) <= 10:
        for pattern, count in sorted(ik1_patterns.items(), key=lambda x: -x[1]):
            print(f"    {list(pattern)} : {count} tournaments")
    else:
        # Show top 5
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

    # Show one detailed n=7 example
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
    print("  Prediction: sum_k I_k(Omega, x) * t^k = A_n(t) + sum_I x^parts * inv(T) * A_{f+1}(t) * (t-1)^{d-f}")

    # For n=5, verify the generating function at specific x, t values
    for n_test in [5, 7]:
        print(f"\n  --- n={n_test} ---")
        if n_test == 5:
            test_tournaments = list(all_tournaments(5))[:5]
        else:
            test_tournaments = [random_tournament(n_test, random.Random(seed)) for seed in range(3)]

        for idx, T in enumerate(test_tournaments):
            d = n_test - 1
            cycles = find_odd_cycles(T)
            catalog = ocf_invariant_catalog_from_cycles(cycles)

            # Pre-compute all I_k polys
            ik_polys_gf = [k_colored_indep_poly_from_catalog(n_test, k, catalog) for k in range(n_test)]

            # Just print pass/fail summary
            all_ok = True
            for x_val in [1, 2, -1, 3]:
                for t_val in [1, 2, -1, 0]:
                    lhs = sum(eval_poly(ik_polys_gf[k], x_val) * t_val**k for k in range(n_test))
                    rhs = sum(eulerian_number(n_test, k) * t_val**k for k in range(n_test))
                    for (parts, f), count in catalog.items():
                        a_f1_t = sum(eulerian_number(f + 1, j) * t_val**j for j in range(f + 1))
                        rhs += count * x_val**parts * a_f1_t * (t_val - 1)**(d - f)
                    if lhs != rhs:
                        all_ok = False
            print(f"    T#{idx}: generating function {'PASS' if all_ok else 'FAIL'} (tested x in {{1,2,-1,3}}, t in {{1,2,-1,0}})")

    # -----------------------------------------------------------------------
    # ADDITIONAL ANALYSIS: What does I_k(Omega, 1) mean?
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS: I_k(Omega, 1) INTERPRETATION")
    print("=" * 70)

    print("\n  From the generating function at x=1:")
    print("  sum_k I_k(Omega, 1) t^k = A_n(t) + sum_{(parts,f)} count * A_{f+1}(t) * (t-1)^{d-f}")
    print("  At t=1: sum_k I_k(Omega, 1) = A_n(1) = n!")
    print("  This is because (t-1)^{d-f} = 0 at t=1 for d > f.")
    print("  So I_k(Omega, 1) = A(n, k) for ALL tournaments!")

    # Verify
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
    # ADDITIONAL ANALYSIS: I_k(Omega, -1)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS: I_k(Omega, -1)")
    print("=" * 70)

    print("\n  From the generating function at x=-1:")
    print("  sum_k I_k(Omega, -1) t^k = A_n(t) + sum_{(parts,f)} count * (-1)^parts * A_{f+1}(t) * (t-1)^{d-f}")
    print("  This is I(Omega, -1) evaluated 'per k'.")
    print()

    # For n=5, compute I(Omega, -1) and compare
    print("  n=5 first 5 tournaments:")
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
    # PALINDROMY deep dive: check c_0 = c_{n-1} for OCF invariants
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PALINDROMY CHECK: c_0^{(f,d)} vs c_d^{(f,d)}")
    print("=" * 70)

    for d in range(2, 8):
        print(f"\n  d = {d}:")
        for f in range(1, d + 1, 2):  # OCF invariants have odd f (from odd cycles)
            # Actually f = sum (len_cycle - 2), and odd cycles have len >= 3
            # For a single 3-cycle: f = 1. For single 5-cycle: f = 3. For two disjoint 3-cycles: f = 2.
            # So f can be even or odd.
            c0 = inflated_eulerian(f, d, 0)
            cd = inflated_eulerian(f, d, d)
            df = d - f
            # The theory: c_0 = 1 when d-f is even (and = -1 when odd?)
            # Let's just check
            print(f"    f={f}: c_0={c0}, c_d={cd}, d-f={df}, d-f even={df%2==0}")

    # Also check with even f (disjoint cycle pairs)
    print("\n  Even f values (e.g., two disjoint 3-cycles, f=2):")
    for d in range(2, 8):
        for f in [2, 4]:
            if f > d:
                continue
            c0 = inflated_eulerian(f, d, 0)
            cd = inflated_eulerian(f, d, d)
            df = d - f
            print(f"    d={d}, f={f}: c_0={c0}, c_d={cd}, d-f={df}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
