#!/usr/bin/env python3
"""
Independence Polynomial Special Values for Omega(T)
====================================================
Instance: opus-2026-03-07-S1

We know I(Omega(T), 2) = H(T) (OCF identity).
This script explores I(Omega(T), x) at other special values of x:

1. I(Omega, -1) — the alternating independence count
2. I(Omega, x) at roots of unity (i, omega = e^{2pi i/3})
3. I'(Omega, 2) — the derivative at the evaluation point
4. Chromatic connection: I(Omega, -1) and acyclic orientations
5. Roots of I(Omega, x) — real-rootedness, distribution, bounds

Uses tournament_lib for core computations.
"""

import sys
import os
import random
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import cmath

# Add tournament_lib to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, random_tournament, find_odd_cycles,
    conflict_graph, hamiltonian_path_count
)


# ---------------------------------------------------------------------------
# Independence polynomial: full coefficient computation
# ---------------------------------------------------------------------------

def indep_poly_coefficients(cycles):
    """Compute independence polynomial coefficients [alpha_0, alpha_1, ..., alpha_d].

    Builds the conflict graph from the cycle list and enumerates all independent
    sets by bitmask. Returns list where coefficients[k] = number of independent
    sets of size k.

    Practical for |cycles| <= ~22.
    """
    m = len(cycles)
    if m == 0:
        return [1]

    # Build adjacency via vertex sharing
    vsets = [frozenset(c) for c in cycles]
    nbr = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if vsets[i] & vsets[j]:
                nbr[i] |= 1 << j
                nbr[j] |= 1 << i

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

    # Trim trailing zeros
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


def eval_poly(coeffs, x):
    """Evaluate polynomial sum(coeffs[k] * x^k) at point x."""
    return sum(c * x**k for k, c in enumerate(coeffs))


def eval_derivative(coeffs, x):
    """Evaluate I'(G, x) = sum(k * coeffs[k] * x^{k-1})."""
    return sum(k * c * x**(k-1) for k, c in enumerate(coeffs) if k > 0)


def eval_second_derivative(coeffs, x):
    """Evaluate I''(G, x) = sum(k*(k-1) * coeffs[k] * x^{k-2})."""
    return sum(k * (k-1) * c * x**(k-2) for k, c in enumerate(coeffs) if k > 1)


def find_roots(coeffs):
    """Find roots of the polynomial. Returns numpy array of complex roots."""
    if len(coeffs) <= 1:
        return np.array([])
    # numpy.roots wants highest degree first
    poly_rev = list(reversed(coeffs))
    return np.roots(poly_rev)


# ---------------------------------------------------------------------------
# Cycle finding using vertex sets (for conflict graph — no duplicate vertex sets)
# ---------------------------------------------------------------------------

def find_odd_cycle_vertex_sets(T):
    """Find all directed odd cycles as frozensets of vertices.

    Two cycles on the same vertex set are considered the SAME node in Omega(T).
    This is the correct definition: Omega(T) has one vertex per odd-cycle vertex set.
    """
    n = len(T)
    cycle_vsets = set()
    for length in range(3, n + 1, 2):
        for combo in combinations(range(n), length):
            vset = frozenset(combo)
            if vset in cycle_vsets:
                continue
            # Check if there's a directed Hamiltonian cycle on this vertex subset
            first = combo[0]
            found = False
            for perm in permutations(combo[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not T[path[i]][path[(i + 1) % length]]:
                        valid = False
                        break
                if valid:
                    found = True
                    break
            if found:
                cycle_vsets.add(vset)
    return list(cycle_vsets)


def build_conflict_adj_from_vsets(vsets):
    """Build adjacency (as neighbor bitmasks) from list of frozensets."""
    m = len(vsets)
    nbr = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if vsets[i] & vsets[j]:
                nbr[i] |= 1 << j
                nbr[j] |= 1 << i
    return nbr


# ---------------------------------------------------------------------------
# Section 1: I(Omega, -1) — alternating independence count
# ---------------------------------------------------------------------------

def section1_alternating_count():
    print("=" * 70)
    print("SECTION 1: I(Omega(T), -1) — Alternating Independence Count")
    print("=" * 70)
    print()
    print("I(G, -1) = sum_k (-1)^k alpha_k")
    print("  = (# even-sized independent sets) - (# odd-sized independent sets)")
    print()

    for n in [3, 4, 5, 6]:
        print(f"--- n = {n} ---")
        results = Counter()
        h_to_iminus1 = defaultdict(list)
        count = 0

        for T in all_tournaments(n):
            cycles = find_odd_cycle_vertex_sets(T)
            coeffs = indep_poly_coefficients(cycles)
            val = eval_poly(coeffs, -1)
            H = hamiltonian_path_count(T)
            results[val] += 1
            h_to_iminus1[H].append(val)
            count += 1

            if n <= 5 and count <= 4:
                print(f"  T#{count}: |Omega| = {len(cycles)}, coeffs = {coeffs}, "
                      f"I(Omega,-1) = {val}, H(T) = {H}")

        print(f"  Total tournaments: {count}")
        print(f"  I(Omega, -1) distribution: {dict(sorted(results.items()))}")

        # Check if I(Omega, -1) is always 0 or 1
        always_01 = all(v in [0, 1] for v in results.keys())
        always_pm1 = all(v in [-1, 0, 1] for v in results.keys())
        print(f"  Always in {{0,1}}: {always_01}")
        print(f"  Always in {{-1,0,1}}: {always_pm1}")

        # Check relationship to H mod 2
        h_mod2_match = True
        for H, vals in h_to_iminus1.items():
            for v in vals:
                if v % 2 != H % 2:
                    h_mod2_match = False
        print(f"  I(Omega,-1) ≡ H(T) mod 2: {h_mod2_match}")

        # Check if I(Omega,-1) = 1 always (for connected graphs, this means
        # the simplicial complex has Euler char = 1)
        all_one = all(v == 1 for v in results.keys())
        print(f"  I(Omega,-1) = 1 always: {all_one}")

        # Redei: H(T) is always odd. So I(Omega,-1) should be odd too if
        # it matches mod 2.
        all_odd = all(v % 2 == 1 for v in results.keys())
        print(f"  I(Omega,-1) always odd: {all_odd}")

        print()

    # Sample for n=6
    print("--- n = 6 (sample of 2000 random tournaments) ---")
    results6 = Counter()
    for _ in range(2000):
        T = random_tournament(6)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val = eval_poly(coeffs, -1)
        results6[val] += 1
    print(f"  I(Omega, -1) distribution: {dict(sorted(results6.items()))}")
    print()


# ---------------------------------------------------------------------------
# Section 2: I(Omega, x) at roots of unity
# ---------------------------------------------------------------------------

def section2_roots_of_unity():
    print("=" * 70)
    print("SECTION 2: I(Omega(T), x) at Roots of Unity")
    print("=" * 70)
    print()

    omega3 = cmath.exp(2j * cmath.pi / 3)  # primitive cube root of unity
    imag_i = 1j  # fourth root of unity

    for n in [3, 4, 5]:
        print(f"--- n = {n} ---")

        vals_i = []
        vals_omega = []
        vals_minus_i = []
        count = 0

        for T in all_tournaments(n):
            count += 1
            cycles = find_odd_cycle_vertex_sets(T)
            coeffs = indep_poly_coefficients(cycles)

            vi = eval_poly(coeffs, imag_i)
            vo = eval_poly(coeffs, omega3)
            vmi = eval_poly(coeffs, -imag_i)

            vals_i.append(vi)
            vals_omega.append(vo)
            vals_minus_i.append(vmi)

        # Analyze |I(Omega, i)|
        mods_i = [abs(v) for v in vals_i]
        mods_omega = [abs(v) for v in vals_omega]

        print(f"  {count} tournaments")
        print(f"  I(Omega, i):")
        print(f"    |I| range: [{min(mods_i):.4f}, {max(mods_i):.4f}]")
        print(f"    Mean |I|: {np.mean(mods_i):.4f}")

        # Check if I(Omega, i) is always real or purely imaginary
        real_parts = [v.real for v in vals_i]
        imag_parts = [v.imag for v in vals_i]
        print(f"    Real part range: [{min(real_parts):.4f}, {max(real_parts):.4f}]")
        print(f"    Imag part range: [{min(imag_parts):.4f}, {max(imag_parts):.4f}]")

        # Check if I(Omega, i) is always an integer (Gaussian integer)
        all_gaussian = all(abs(v.real - round(v.real)) < 1e-8 and
                          abs(v.imag - round(v.imag)) < 1e-8 for v in vals_i)
        print(f"    Always Gaussian integer: {all_gaussian}")
        if all_gaussian:
            gi_dist = Counter((int(round(v.real)), int(round(v.imag))) for v in vals_i)
            print(f"    Gaussian integer distribution: {dict(sorted(gi_dist.items()))}")

        print(f"  I(Omega, omega_3):")
        print(f"    |I| range: [{min(mods_omega):.4f}, {max(mods_omega):.4f}]")
        print(f"    Mean |I|: {np.mean(mods_omega):.4f}")

        # Check if always Eisenstein integer
        # omega_3 = (-1 + sqrt(3)i)/2
        all_eisenstein = True
        for v in vals_omega:
            # Eisenstein integers: a + b*omega_3 where a,b are integers
            # v = a + b*omega_3 = a + b*(-1/2 + sqrt(3)/2 * i)
            # real part = a - b/2, imag part = b*sqrt(3)/2
            b_cand = v.imag / (np.sqrt(3) / 2)
            a_cand = v.real + b_cand / 2
            if abs(a_cand - round(a_cand)) > 1e-6 or abs(b_cand - round(b_cand)) > 1e-6:
                all_eisenstein = False
                break
        print(f"    Always Eisenstein integer: {all_eisenstein}")

        # Check I(Omega, i) * I(Omega, -i) = |I(Omega, i)|^2
        # This is always true; more interesting: is the product always an integer?
        products = [v * vc for v, vc in zip(vals_i, vals_minus_i)]
        all_real_int = all(abs(p.imag) < 1e-8 and abs(p.real - round(p.real)) < 1e-8
                          for p in products)
        print(f"  I(Omega, i) * I(Omega, -i) always integer: {all_real_int}")
        if all_real_int:
            norm_dist = Counter(int(round(p.real)) for p in products)
            print(f"    Norm distribution: {dict(sorted(norm_dist.items()))}")

        print()

    # Sample n=6
    print("--- n = 6 (sample of 1000 random) ---")
    vals_i_6 = []
    vals_omega_6 = []
    for _ in range(1000):
        T = random_tournament(6)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        vals_i_6.append(eval_poly(coeffs, imag_i))
        vals_omega_6.append(eval_poly(coeffs, omega3))
    mods_i_6 = [abs(v) for v in vals_i_6]
    mods_o_6 = [abs(v) for v in vals_omega_6]
    print(f"  I(Omega, i): |I| range [{min(mods_i_6):.4f}, {max(mods_i_6):.4f}], "
          f"mean {np.mean(mods_i_6):.4f}")
    all_gi_6 = all(abs(v.real - round(v.real)) < 1e-6 and
                   abs(v.imag - round(v.imag)) < 1e-6 for v in vals_i_6)
    print(f"  Always Gaussian integer: {all_gi_6}")
    if all_gi_6:
        gi6 = Counter((int(round(v.real)), int(round(v.imag))) for v in vals_i_6)
        top5 = sorted(gi6.items(), key=lambda x: -x[1])[:10]
        print(f"  Top Gaussian int values: {top5}")
    print(f"  I(Omega, omega_3): |I| range [{min(mods_o_6):.4f}, {max(mods_o_6):.4f}], "
          f"mean {np.mean(mods_o_6):.4f}")
    print()


# ---------------------------------------------------------------------------
# Section 3: I'(Omega, 2) — derivative at the evaluation point
# ---------------------------------------------------------------------------

def section3_derivative():
    print("=" * 70)
    print("SECTION 3: I'(Omega(T), 2) — Derivative at the OCF Point")
    print("=" * 70)
    print()
    print("I'(G, x) = sum_k k * alpha_k * x^{k-1}")
    print("At x=2: I'(Omega, 2) = alpha_1 + 4*alpha_2 + 12*alpha_3 + 32*alpha_4 + ...")
    print()

    for n in [3, 4, 5, 6]:
        print(f"--- n = {n} ---")
        deriv_vals = []
        H_vals = []
        ratio_vals = []
        second_deriv_vals = []
        count = 0

        gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(2000))
        label = "exhaustive" if n <= 5 else "sample of 2000"

        for T in gen:
            count += 1
            cycles = find_odd_cycle_vertex_sets(T)
            coeffs = indep_poly_coefficients(cycles)
            H = eval_poly(coeffs, 2)  # = H(T) by OCF
            dH = eval_derivative(coeffs, 2)
            d2H = eval_second_derivative(coeffs, 2)

            deriv_vals.append(dH)
            H_vals.append(H)
            second_deriv_vals.append(d2H)
            if H > 0:
                ratio_vals.append(dH / H)

            if n <= 4 and count <= 8:
                print(f"  T#{count}: coeffs = {coeffs}, H = {H}, "
                      f"I'(Omega,2) = {dH}, ratio = {dH/H:.4f}" if H > 0 else "")

        print(f"  {count} tournaments ({label})")
        print(f"  I'(Omega, 2) range: [{min(deriv_vals)}, {max(deriv_vals)}]")
        print(f"  Mean I'(Omega, 2): {np.mean(deriv_vals):.2f}")
        print(f"  H(T) range: [{min(H_vals)}, {max(H_vals)}]")

        if ratio_vals:
            print(f"  Ratio I'/H range: [{min(ratio_vals):.4f}, {max(ratio_vals):.4f}]")
            print(f"  Mean I'/H: {np.mean(ratio_vals):.4f}")

        # Check if I'(Omega, 2) is always odd
        all_int = all(isinstance(d, int) or abs(d - round(d)) < 1e-8 for d in deriv_vals)
        if all_int:
            int_vals = [int(round(d)) for d in deriv_vals]
            all_odd = all(v % 2 == 1 for v in int_vals)
            print(f"  I'(Omega, 2) always integer: True")
            print(f"  I'(Omega, 2) always odd: {all_odd}")
            parity_dist = Counter(v % 2 for v in int_vals)
            print(f"  Parity distribution: {dict(parity_dist)}")

        # Check if I'(Omega, 2) relates to sum of H over single-vertex deletions
        # For n<=5, we can verify
        if n <= 5:
            print("  Checking combinatorial interpretations...")
            mismatch_count = 0
            for T in all_tournaments(n):
                cycles = find_odd_cycle_vertex_sets(T)
                coeffs = indep_poly_coefficients(cycles)
                dH = eval_derivative(coeffs, 2)

                # Compare with sum of H(T-v)
                sum_H_del = 0
                for v in range(n):
                    Tv = [[T[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
                    sum_H_del += hamiltonian_path_count(Tv)

                if count <= 4 and n <= 4:
                    print(f"    I'(Omega,2) = {dH}, sum_v H(T-v) = {sum_H_del}")

                if dH != sum_H_del:
                    mismatch_count += 1

            if mismatch_count == 0:
                print(f"  *** I'(Omega,2) = sum_v H(T-v) for ALL tournaments! ***")
            else:
                print(f"  I'(Omega,2) != sum_v H(T-v) for {mismatch_count} tournaments")

        # Logarithmic derivative: I'/I = d/dx log I at x=2
        if ratio_vals:
            print(f"  Log derivative d/dx log I(Omega,x)|_{{x=2}}:")
            print(f"    Range: [{min(ratio_vals):.6f}, {max(ratio_vals):.6f}]")

        print()


# ---------------------------------------------------------------------------
# Section 4: Chromatic / Acyclic Orientation Connection
# ---------------------------------------------------------------------------

def section4_chromatic():
    print("=" * 70)
    print("SECTION 4: Chromatic Connection — I(Omega, -1) and Acyclic Orientations")
    print("=" * 70)
    print()
    print("For a graph G with m vertices:")
    print("  I(G, -1) = (-1)^m * (# acyclic orientations of complement G_bar)")
    print("  (Stanley's theorem, when G is the incomparability graph)")
    print()
    print("For Omega(T), the complement is the 'compatibility graph':")
    print("  cycles are compatible iff they DON'T share vertices.")
    print()

    for n in [3, 4, 5]:
        print(f"--- n = {n} ---")
        count = 0

        for T in all_tournaments(n):
            count += 1
            cycles = find_odd_cycle_vertex_sets(T)
            m = len(cycles)
            coeffs = indep_poly_coefficients(cycles)
            val = eval_poly(coeffs, -1)

            # Count acyclic orientations of Omega(T) itself
            # An acyclic orientation of an undirected graph: orient each edge
            # such that no directed cycle exists.
            # For small graphs we can enumerate.
            if m <= 12:
                # Build edge list for Omega(T)
                edges = []
                for i in range(m):
                    for j in range(i + 1, m):
                        if cycles[i] & cycles[j]:
                            edges.append((i, j))

                # Count acyclic orientations by brute force
                ne = len(edges)
                acyclic_count = 0
                if ne <= 18:
                    for bits in range(1 << ne):
                        # Orient edges
                        out = [set() for _ in range(m)]
                        for k, (u, v) in enumerate(edges):
                            if bits & (1 << k):
                                out[u].add(v)
                            else:
                                out[v].add(u)

                        # Check acyclicity via topological sort
                        in_deg = [0] * m
                        for u in range(m):
                            for v in out[u]:
                                in_deg[v] += 1
                        queue = [u for u in range(m) if in_deg[u] == 0]
                        visited = 0
                        while queue:
                            u = queue.pop()
                            visited += 1
                            for v in out[u]:
                                in_deg[v] -= 1
                                if in_deg[v] == 0:
                                    queue.append(v)
                        if visited == m:
                            acyclic_count += 1

                    sign = (-1) ** m
                    prediction = sign * acyclic_count

                    if count <= 8 and n <= 4:
                        print(f"  T#{count}: m={m}, |E|={ne}, I(Omega,-1)={val}, "
                              f"acyclic_orient={acyclic_count}, (-1)^m * ao = {prediction}")

        # Summary statistics
        print(f"  Total tournaments checked: {count}")

        # Collect for aggregate
        vals_list = []
        ao_list = []
        for T in all_tournaments(n):
            cycles = find_odd_cycle_vertex_sets(T)
            m = len(cycles)
            coeffs = indep_poly_coefficients(cycles)
            val = eval_poly(coeffs, -1)
            vals_list.append(val)

            # Build complement graph edges
            edges = []
            for i in range(m):
                for j in range(i + 1, m):
                    if not (cycles[i] & cycles[j]):  # complement: NON-adjacent
                        edges.append((i, j))
            ne = len(edges)

            # Count acyclic orientations of COMPLEMENT
            ao_comp = 0
            if ne <= 18 and m > 0:
                for bits in range(1 << ne):
                    out = [set() for _ in range(m)]
                    for k, (u, v) in enumerate(edges):
                        if bits & (1 << k):
                            out[u].add(v)
                        else:
                            out[v].add(u)
                    in_deg = [0] * m
                    for u in range(m):
                        for v in out[u]:
                            in_deg[v] += 1
                    queue = [u for u in range(m) if in_deg[u] == 0]
                    visited = 0
                    while queue:
                        u = queue.pop()
                        visited += 1
                        for v in out[u]:
                            in_deg[v] -= 1
                            if in_deg[v] == 0:
                                queue.append(v)
                    if visited == m:
                        ao_comp += 1
            elif m == 0:
                ao_comp = 1
            ao_list.append(ao_comp)

        match_omega = sum(1 for v, a, m2 in zip(vals_list, ao_list,
                         [len(find_odd_cycle_vertex_sets(T)) for T in all_tournaments(n)])
                         if v == (-1)**m2 * a)
        match_comp = sum(1 for v, a in zip(vals_list, ao_list)
                         if abs(v) == a or v == a)
        total = len(vals_list)
        print(f"  I(Omega,-1) = (-1)^m * AO(Omega) match rate: {match_omega}/{total}")
        print(f"  |I(Omega,-1)| = AO(complement) match rate: {match_comp}/{total}")
        print()


# ---------------------------------------------------------------------------
# Section 5: Roots of I(Omega(T), x)
# ---------------------------------------------------------------------------

def section5_roots():
    print("=" * 70)
    print("SECTION 5: Roots of I(Omega(T), x)")
    print("=" * 70)
    print()
    print("Known: I(G, x) has only real roots when G is claw-free (Chudnovsky-Seymour).")
    print("Omega(T) is claw-free for n <= 8 (verified elsewhere).")
    print()

    all_roots_real = True
    largest_real_roots = []
    smallest_real_roots = []
    root_data_by_n = {}

    for n in [3, 4, 5]:
        print(f"--- n = {n} (exhaustive) ---")
        n_real_rooted = 0
        n_total = 0
        max_roots = []
        min_roots = []
        degree_dist = Counter()
        non_real_examples = []

        for T in all_tournaments(n):
            n_total += 1
            cycles = find_odd_cycle_vertex_sets(T)
            coeffs = indep_poly_coefficients(cycles)
            degree_dist[len(coeffs) - 1] += 1

            if len(coeffs) <= 1:
                n_real_rooted += 1
                continue

            roots = find_roots(coeffs)
            real_check = all(abs(r.imag) < 1e-8 for r in roots)
            if real_check:
                n_real_rooted += 1
                real_roots = sorted([r.real for r in roots])
                max_roots.append(real_roots[-1])
                min_roots.append(real_roots[0])
            else:
                all_roots_real = False
                non_real_examples.append((n_total, coeffs, roots))

        print(f"  {n_total} tournaments")
        print(f"  Degree distribution: {dict(sorted(degree_dist.items()))}")
        print(f"  Real-rooted: {n_real_rooted}/{n_total}")

        if non_real_examples:
            print(f"  Non-real-rooted examples:")
            for idx, c, r in non_real_examples[:3]:
                print(f"    T#{idx}: coeffs={c}")
                for root in r:
                    if abs(root.imag) > 1e-8:
                        print(f"      root: {root}")

        if max_roots:
            print(f"  Largest real root range: [{min(max_roots):.6f}, {max(max_roots):.6f}]")
            print(f"  Mean largest real root: {np.mean(max_roots):.6f}")
            largest_real_roots.extend(max_roots)
        if min_roots:
            print(f"  Smallest real root range: [{min(min_roots):.6f}, {max(min_roots):.6f}]")
            smallest_real_roots.extend(min_roots)

        root_data_by_n[n] = {
            'total': n_total, 'real_rooted': n_real_rooted,
            'max_roots': max_roots, 'min_roots': min_roots
        }
        print()

    # n=6: sample
    print("--- n = 6 (sample of 2000 random tournaments) ---")
    n_real_6 = 0
    n_total_6 = 0
    max_roots_6 = []
    min_roots_6 = []
    degree_dist_6 = Counter()
    non_real_6 = []

    for _ in range(2000):
        T = random_tournament(6)
        n_total_6 += 1
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        degree_dist_6[len(coeffs) - 1] += 1

        if len(coeffs) <= 1:
            n_real_6 += 1
            continue

        roots = find_roots(coeffs)
        real_check = all(abs(r.imag) < 1e-8 for r in roots)
        if real_check:
            n_real_6 += 1
            real_roots = sorted([r.real for r in roots])
            max_roots_6.append(real_roots[-1])
            min_roots_6.append(real_roots[0])
        else:
            if len(non_real_6) < 5:
                non_real_6.append((n_total_6, coeffs, roots))

    print(f"  {n_total_6} tournaments sampled")
    print(f"  Degree distribution: {dict(sorted(degree_dist_6.items()))}")
    print(f"  Real-rooted: {n_real_6}/{n_total_6}")

    if non_real_6:
        print(f"  Non-real-rooted examples found!")
        for idx, c, r in non_real_6[:3]:
            print(f"    T#{idx}: coeffs={c}")
            for root in r:
                if abs(root.imag) > 1e-8:
                    print(f"      root: {root}")
    else:
        print(f"  ALL sampled tournaments are real-rooted!")

    if max_roots_6:
        print(f"  Largest real root range: [{min(max_roots_6):.6f}, {max(max_roots_6):.6f}]")
        print(f"  Mean largest real root: {np.mean(max_roots_6):.6f}")
    if min_roots_6:
        print(f"  Smallest real root range: [{min(min_roots_6):.6f}, {max(min_roots_6):.6f}]")

    print()

    # Summary
    print("=" * 70)
    print("ROOT ANALYSIS SUMMARY")
    print("=" * 70)
    all_real = all(root_data_by_n[n]['real_rooted'] == root_data_by_n[n]['total']
                   for n in root_data_by_n) and n_real_6 == n_total_6
    print(f"All roots real for n<=5 (exhaustive): "
          f"{all(root_data_by_n[n]['real_rooted'] == root_data_by_n[n]['total'] for n in root_data_by_n)}")
    print(f"All roots real for n=6 (sample): {n_real_6 == n_total_6}")

    if largest_real_roots:
        print(f"Overall largest real root seen (n<=5): {max(largest_real_roots):.6f}")
    if max_roots_6:
        print(f"Largest real root seen (n=6 sample): {max(max_roots_6):.6f}")

    # Universal bound: are all roots > -1/2?
    all_above = True
    bound = None
    for n in root_data_by_n:
        if root_data_by_n[n]['min_roots']:
            mn = min(root_data_by_n[n]['min_roots'])
            if bound is None or mn < bound:
                bound = mn
    if min_roots_6:
        mn6 = min(min_roots_6)
        if bound is None or mn6 < bound:
            bound = mn6
    print(f"Most negative root seen: {bound:.6f}" if bound else "No roots found")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Independence Polynomial Special Values for Omega(T)")
    print("Instance: opus-2026-03-07-S1")
    print()

    section1_alternating_count()
    section2_roots_of_unity()
    section3_derivative()
    section4_chromatic()
    section5_roots()

    print("=" * 70)
    print("COMPUTATION COMPLETE")
    print("=" * 70)
