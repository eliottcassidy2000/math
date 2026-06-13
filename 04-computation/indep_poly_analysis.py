#!/usr/bin/env python3
"""
Independence Polynomial Analysis of Omega(T)
=============================================
For tournaments T, build the conflict graph Omega(T) on odd cycles,
compute the full independence polynomial, and test:
  1. Real-rootedness
  2. Ultra-log-concavity
  3. Interlacing under vertex deletion

Also: Paley tournament analysis for p=3,7,11.
"""

import sys, os, random, time
from math import comb
from itertools import combinations

# Add parent so we can import tournament_lib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, random_tournament, find_odd_cycles,
    conflict_graph, tournament_from_bits
)

import numpy as np


# ---------------------------------------------------------------------------
# Independence polynomial: full coefficient computation
# ---------------------------------------------------------------------------

def indep_poly_coefficients(adj):
    """Compute independence polynomial coefficients [alpha_0, alpha_1, ..., alpha_d].
    adj is adjacency matrix (list of lists). alpha_k = number of independent sets of size k.
    Uses bitmask enumeration - practical for m <= ~22."""
    m = len(adj)
    if m == 0:
        return [1]

    # Precompute neighbor bitmasks
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j

    # Count independent sets by size
    max_possible = m
    counts = [0] * (max_possible + 1)

    for mask in range(1 << m):
        # Check independence
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

    # Trim trailing zeros
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


def indep_poly_coefficients_fast(cycles, max_indep_size=None):
    """Faster coefficient computation exploiting vertex-disjointness structure.
    Independent sets in Omega correspond to vertex-disjoint cycle collections."""
    m = len(cycles)
    if m == 0:
        return [1]

    if max_indep_size is None:
        total_verts = len(set(v for c in cycles for v in c))
        max_indep_size = total_verts // 3

    vsets = [frozenset(c) for c in cycles]
    adj_bits = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    # For small m, just do bitmask enumeration
    if m <= 22:
        nbr = adj_bits
        counts = [0] * (max_indep_size + 1)
        for mask in range(1 << m):
            ok = True
            seen = 0
            temp = mask
            sz = 0
            while temp:
                v = (temp & -temp).bit_length() - 1
                if nbr[v] & seen:
                    ok = False
                    break
                seen |= 1 << v
                temp &= temp - 1
                sz += 1
            if ok:
                if sz < len(counts):
                    counts[sz] += 1
                else:
                    counts.extend([0] * (sz - len(counts) + 1))
                    counts[sz] += 1
        while len(counts) > 1 and counts[-1] == 0:
            counts.pop()
        return counts

    # For larger m, enumerate by size up to max_indep_size
    counts = [1, m]  # size 0 and 1
    if max_indep_size >= 2:
        pairs = []
        ct2 = 0
        for a in range(m):
            for b in range(a + 1, m):
                if not (adj_bits[a] & (1 << b)):
                    pairs.append((a, b))
                    ct2 += 1
        counts.append(ct2)

        if max_indep_size >= 3:
            ct3 = 0
            triples = []
            for a, b in pairs:
                for c in range(b + 1, m):
                    if not (adj_bits[a] & (1 << c)) and not (adj_bits[b] & (1 << c)):
                        triples.append((a, b, c))
                        ct3 += 1
            counts.append(ct3)

            if max_indep_size >= 4:
                ct4 = 0
                for a, b, c in triples:
                    for d in range(c + 1, m):
                        if (not (adj_bits[a] & (1 << d)) and
                            not (adj_bits[b] & (1 << d)) and
                            not (adj_bits[c] & (1 << d))):
                            ct4 += 1
                counts.append(ct4)

    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def check_real_roots(coeffs):
    """Check if all roots of the polynomial are real.
    Returns (all_real, roots)."""
    if len(coeffs) <= 1:
        return True, np.array([])
    # numpy.roots wants highest degree first
    poly_rev = list(reversed(coeffs))
    roots = np.roots(poly_rev)
    all_real = all(abs(r.imag) < 1e-8 for r in roots)
    return all_real, roots


def check_ultra_log_concave(coeffs, m):
    """Check ultra-log-concavity: alpha_k / C(m,k) is log-concave.
    That means (a_k/C(m,k))^2 >= (a_{k-1}/C(m,k-1)) * (a_{k+1}/C(m,k+1)) for all valid k."""
    if len(coeffs) <= 2:
        return True, []
    violations = []
    for k in range(1, len(coeffs) - 1):
        if k > m or k - 1 > m or k + 1 > m:
            break
        c_k = comb(m, k)
        c_km1 = comb(m, k - 1)
        c_kp1 = comb(m, k + 1)
        if c_k == 0 or c_km1 == 0 or c_kp1 == 0:
            continue
        lhs = (coeffs[k] / c_k) ** 2
        rhs = (coeffs[k - 1] / c_km1) * (coeffs[k + 1] / c_kp1)
        if lhs < rhs - 1e-12:
            violations.append((k, lhs, rhs))
    return len(violations) == 0, violations


def check_interlacing(coeffs_g, coeffs_f):
    """Check if f interlaces g.
    Between consecutive roots of g there must be a root of f.
    f should have degree = deg(g) - 1 or deg(g)."""
    deg_g = len(coeffs_g) - 1
    deg_f = len(coeffs_f) - 1

    if deg_g <= 0 or deg_f <= 0:
        return True  # trivial

    if deg_f != deg_g and deg_f != deg_g - 1:
        return False  # wrong degree relationship

    roots_g = np.sort(np.roots(list(reversed(coeffs_g))).real)
    roots_f = np.sort(np.roots(list(reversed(coeffs_f))).real)

    # Check: between consecutive roots of g, there's a root of f
    for i in range(len(roots_g) - 1):
        lo, hi = roots_g[i], roots_g[i + 1]
        if hi - lo < 1e-10:
            continue  # repeated root
        found = any(lo - 1e-8 <= r <= hi + 1e-8 for r in roots_f)
        if not found:
            return False
    return True


def indep_poly_delete_vertex(adj, v):
    """Compute independence polynomial of G - v (delete vertex v from graph with adj matrix)."""
    m = len(adj)
    new_adj = []
    for i in range(m):
        if i == v:
            continue
        row = []
        for j in range(m):
            if j == v:
                continue
            row.append(adj[i][j])
        new_adj.append(row)
    return indep_poly_coefficients(new_adj)


# ---------------------------------------------------------------------------
# Paley tournaments
# ---------------------------------------------------------------------------

def paley_tournament(p, qr):
    """Build Paley tournament T_p. Arc i->j iff (j-i) mod p in qr."""
    T = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                T[i][j] = 1
    return T


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_tournament(T, label=""):
    """Full analysis of one tournament. Returns dict of results."""
    n = len(T)
    cycles = find_odd_cycles(T)
    m = len(cycles)

    if m == 0:
        return {
            'n': n, 'num_cycles': 0, 'coeffs': [1], 'max_indep': 0,
            'real_roots': True, 'ulc': True, 'interlacing': True,
            'label': label
        }

    # Build conflict graph
    cg = conflict_graph(cycles)

    # Compute coefficients - use fast method for large m, exact for small
    if m <= 20:
        coeffs = indep_poly_coefficients(cg)
    else:
        coeffs = indep_poly_coefficients_fast(cycles)

    max_indep = len(coeffs) - 1

    # Check real roots
    real_roots, roots = check_real_roots(coeffs)

    # Check ULC
    ulc, ulc_violations = check_ultra_log_concave(coeffs, m)

    # Check interlacing (pick a few vertices)
    interlacing_ok = True
    vertices_to_check = list(range(min(m, 5)))  # check first 5 vertices of Omega
    if m <= 20:
        for v in vertices_to_check:
            coeffs_minus_v = indep_poly_delete_vertex(cg, v)
            if not check_interlacing(coeffs, coeffs_minus_v):
                interlacing_ok = False
                break

    return {
        'n': n, 'num_cycles': m, 'coeffs': coeffs, 'max_indep': max_indep,
        'real_roots': real_roots, 'ulc': ulc, 'interlacing': interlacing_ok,
        'label': label, 'roots': roots if m > 0 else []
    }


def main():
    rng = random.Random(42)

    stats = {
        'total': 0, 'real_roots_pass': 0, 'real_roots_fail': 0,
        'ulc_pass': 0, 'ulc_fail': 0,
        'interlacing_pass': 0, 'interlacing_fail': 0,
        'trivial': 0,  # no cycles
    }
    fail_examples = {'real_roots': [], 'ulc': [], 'interlacing': []}

    # -----------------------------------------------------------------------
    # Part 1: All n=5 tournaments
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("PART 1: ALL n=5 TOURNAMENTS (1024 total)")
    print("=" * 70)
    t0 = time.time()
    n5_count = 0
    n5_max_cycles = 0
    n5_max_indep = 0

    for T in all_tournaments(5):
        n5_count += 1
        res = analyze_tournament(T, f"n5-{n5_count}")
        stats['total'] += 1

        if res['num_cycles'] == 0:
            stats['trivial'] += 1
            stats['real_roots_pass'] += 1
            stats['ulc_pass'] += 1
            stats['interlacing_pass'] += 1
            continue

        n5_max_cycles = max(n5_max_cycles, res['num_cycles'])
        n5_max_indep = max(n5_max_indep, res['max_indep'])

        if res['real_roots']:
            stats['real_roots_pass'] += 1
        else:
            stats['real_roots_fail'] += 1
            if len(fail_examples['real_roots']) < 3:
                fail_examples['real_roots'].append(res)

        if res['ulc']:
            stats['ulc_pass'] += 1
        else:
            stats['ulc_fail'] += 1
            if len(fail_examples['ulc']) < 3:
                fail_examples['ulc'].append(res)

        if res['interlacing']:
            stats['interlacing_pass'] += 1
        else:
            stats['interlacing_fail'] += 1
            if len(fail_examples['interlacing']) < 3:
                fail_examples['interlacing'].append(res)

    t1 = time.time()
    print(f"  Processed {n5_count} tournaments in {t1-t0:.1f}s")
    print(f"  Max cycles in Omega: {n5_max_cycles}")
    print(f"  Max independent set size: {n5_max_indep}")
    print(f"  Trivial (no cycles): {stats['trivial']}")
    print()

    # -----------------------------------------------------------------------
    # Part 1b: 200 random n=6 tournaments
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("PART 1b: 200 RANDOM n=6 TOURNAMENTS")
    print("=" * 70)
    t0 = time.time()
    n6_max_cycles = 0
    n6_max_indep = 0

    for i in range(200):
        T = random_tournament(6, rng)
        res = analyze_tournament(T, f"n6-{i}")
        stats['total'] += 1

        if res['num_cycles'] == 0:
            stats['trivial'] += 1
            stats['real_roots_pass'] += 1
            stats['ulc_pass'] += 1
            stats['interlacing_pass'] += 1
            continue

        n6_max_cycles = max(n6_max_cycles, res['num_cycles'])
        n6_max_indep = max(n6_max_indep, res['max_indep'])

        if res['real_roots']:
            stats['real_roots_pass'] += 1
        else:
            stats['real_roots_fail'] += 1
            if len(fail_examples['real_roots']) < 3:
                fail_examples['real_roots'].append(res)

        if res['ulc']:
            stats['ulc_pass'] += 1
        else:
            stats['ulc_fail'] += 1
            if len(fail_examples['ulc']) < 3:
                fail_examples['ulc'].append(res)

        if res['interlacing']:
            stats['interlacing_pass'] += 1
        else:
            stats['interlacing_fail'] += 1
            if len(fail_examples['interlacing']) < 3:
                fail_examples['interlacing'].append(res)

    t1 = time.time()
    print(f"  Processed 200 n=6 tournaments in {t1-t0:.1f}s")
    print(f"  Max cycles in Omega: {n6_max_cycles}")
    print(f"  Max independent set size: {n6_max_indep}")
    print()

    # -----------------------------------------------------------------------
    # Part 1c: 100 random n=7 tournaments
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("PART 1c: 100 RANDOM n=7 TOURNAMENTS")
    print("=" * 70)
    t0 = time.time()
    n7_max_cycles = 0
    n7_max_indep = 0

    for i in range(100):
        T = random_tournament(7, rng)
        res = analyze_tournament(T, f"n7-{i}")
        stats['total'] += 1

        if res['num_cycles'] == 0:
            stats['trivial'] += 1
            stats['real_roots_pass'] += 1
            stats['ulc_pass'] += 1
            stats['interlacing_pass'] += 1
            continue

        n7_max_cycles = max(n7_max_cycles, res['num_cycles'])
        n7_max_indep = max(n7_max_indep, res['max_indep'])

        if res['real_roots']:
            stats['real_roots_pass'] += 1
        else:
            stats['real_roots_fail'] += 1
            if len(fail_examples['real_roots']) < 3:
                fail_examples['real_roots'].append(res)

        if res['ulc']:
            stats['ulc_pass'] += 1
        else:
            stats['ulc_fail'] += 1
            if len(fail_examples['ulc']) < 3:
                fail_examples['ulc'].append(res)

        if res['interlacing']:
            stats['interlacing_pass'] += 1
        else:
            stats['interlacing_fail'] += 1
            if len(fail_examples['interlacing']) < 3:
                fail_examples['interlacing'].append(res)

    t1 = time.time()
    print(f"  Processed 100 n=7 tournaments in {t1-t0:.1f}s")
    print(f"  Max cycles in Omega: {n7_max_cycles}")
    print(f"  Max independent set size: {n7_max_indep}")
    print()

    # -----------------------------------------------------------------------
    # Part 2: Paley tournaments (3-cycles only for Omega_3)
    # -----------------------------------------------------------------------
    print("=" * 70)
    print("PART 2: PALEY TOURNAMENTS - Omega_3 (3-cycles only)")
    print("=" * 70)

    paley_data = [
        (3, {1}),
        (7, {1, 2, 4}),
        (11, {1, 3, 4, 5, 9}),
    ]

    for p, qr in paley_data:
        print(f"\n--- T_{p} (Paley tournament on {p} vertices) ---")
        T = paley_tournament(p, qr)

        # Find 3-cycles only
        n = len(T)
        three_cycles = []
        for verts in combinations(range(n), 3):
            a, b, c = verts
            # Check all cyclic orderings (2 per triple)
            if T[a][b] and T[b][c] and T[c][a]:
                three_cycles.append((a, b, c))
            elif T[a][c] and T[c][b] and T[b][a]:
                three_cycles.append((a, c, b))

        m3 = len(three_cycles)
        print(f"  Number of 3-cycles: {m3}")

        if m3 == 0:
            print("  No 3-cycles, I(Omega_3, x) = 1")
            continue

        # Build conflict graph on 3-cycles
        cg3 = conflict_graph(three_cycles)

        # Compute full coefficients
        if m3 <= 22:
            coeffs3 = indep_poly_coefficients(cg3)
        else:
            coeffs3 = indep_poly_coefficients_fast(three_cycles)

        print(f"  Independence polynomial coefficients: {coeffs3}")
        poly_str = " + ".join(
            f"{c}*x^{k}" if k > 0 else str(c)
            for k, c in enumerate(coeffs3) if c != 0
        )
        print(f"  I(Omega_3(T_{p}), x) = {poly_str}")

        # Compute and print roots
        if len(coeffs3) > 1:
            real_ok, roots = check_real_roots(coeffs3)
            print(f"  Roots: {np.sort(roots.real) if real_ok else roots}")
            print(f"  All roots real? {real_ok}")
        else:
            print("  Polynomial is constant (degree 0), no roots")

        # Also do full Omega (all odd cycles)
        print(f"\n  --- Full Omega(T_{p}) (all odd cycles) ---")
        all_cyc = find_odd_cycles(T)
        m_all = len(all_cyc)
        print(f"  Total odd cycles: {m_all}")
        if m_all > 0 and m_all <= 22:
            cg_all = conflict_graph(all_cyc)
            coeffs_all = indep_poly_coefficients(cg_all)
            print(f"  Full independence polynomial coefficients: {coeffs_all}")
            real_ok_all, roots_all = check_real_roots(coeffs_all)
            if len(coeffs_all) > 1:
                print(f"  All roots real? {real_ok_all}")
                print(f"  Roots: {np.sort(roots_all.real) if real_ok_all else roots_all}")
        elif m_all > 22:
            coeffs_all = indep_poly_coefficients_fast(all_cyc)
            print(f"  Full independence polynomial coefficients: {coeffs_all}")
            real_ok_all, roots_all = check_real_roots(coeffs_all)
            if len(coeffs_all) > 1:
                print(f"  All roots real? {real_ok_all}")
                if real_ok_all:
                    print(f"  Roots (real parts): {np.sort(roots_all.real)}")
                else:
                    print(f"  Roots: {roots_all}")

    # -----------------------------------------------------------------------
    # Part 3: Summary
    # -----------------------------------------------------------------------
    print()
    print("=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    print(f"  Total tournaments tested: {stats['total']}")
    print(f"  Trivial (no odd cycles):  {stats['trivial']}")
    print()
    print(f"  REAL ROOTS:")
    print(f"    Pass: {stats['real_roots_pass']}")
    print(f"    Fail: {stats['real_roots_fail']}")
    print()
    print(f"  ULTRA-LOG-CONCAVITY (alpha_k/C(m,k)):")
    print(f"    Pass: {stats['ulc_pass']}")
    print(f"    Fail: {stats['ulc_fail']}")
    print()
    print(f"  INTERLACING (vertex deletion):")
    print(f"    Pass: {stats['interlacing_pass']}")
    print(f"    Fail: {stats['interlacing_fail']}")

    if fail_examples['real_roots']:
        print()
        print("  REAL ROOTS FAILURE EXAMPLES:")
        for ex in fail_examples['real_roots'][:3]:
            print(f"    {ex['label']}: coeffs={ex['coeffs']}, m={ex['num_cycles']}")

    if fail_examples['ulc']:
        print()
        print("  ULC FAILURE EXAMPLES:")
        for ex in fail_examples['ulc'][:3]:
            print(f"    {ex['label']}: coeffs={ex['coeffs']}, m={ex['num_cycles']}")

    if fail_examples['interlacing']:
        print()
        print("  INTERLACING FAILURE EXAMPLES:")
        for ex in fail_examples['interlacing'][:3]:
            print(f"    {ex['label']}: coeffs={ex['coeffs']}, m={ex['num_cycles']}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
