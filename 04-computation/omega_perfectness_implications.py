#!/usr/bin/env python3
"""
Omega(T) Perfectness Implications for OCF

Investigates structural consequences of Omega(T) being perfect for the
independence polynomial I(Omega(T), x) and its relationship to H(T).

Uses tournament_lib for correct directed cycle enumeration.
Omega(T) vertices = directed odd cycles (tuples), edges = shared vertex.

Key findings explored:
  1. OCF verification: I(Omega(T), 2) = H(T)
  2. Roots of I(Omega(T), x): real? negative?
  3. Perfectness and chordality of Omega(T)
  4. Log-concavity and unimodality of coefficients
  5. Modular arithmetic: H(T) mod 2^k from OCF
  6. Component factorization
  7. alpha/omega of Omega vs tournament invariants
"""

import sys
import os
import itertools
import random
from collections import Counter, defaultdict
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count,
    find_odd_cycles, conflict_graph, independence_poly_at,
    independence_poly_at_fast
)

# ─────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────

def random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def score_sequence(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def num_3_cycles(T):
    n = len(T)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (T[i][j] and T[j][k] and T[k][i]) or \
                   (T[i][k] and T[k][j] and T[j][i]):
                    count += 1
    return count

# ─────────────────────────────────────────────────
# Independence polynomial coefficients (for small m)
# ─────────────────────────────────────────────────

def independence_polynomial_coeffs(cg_adj, m):
    """Full [alpha_0, ..., alpha_k] using bitmask. Only for m <= 22."""
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if cg_adj[i][j]:
                nbr[i] |= (1 << j)

    coeffs = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        bits = mask
        while bits:
            v = (bits & -bits).bit_length() - 1
            if nbr[v] & mask & ((1 << v) - 1):
                ok = False
                break
            bits &= bits - 1
        if ok:
            coeffs[bin(mask).count('1')] += 1

    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

# ─────────────────────────────────────────────────
# Graph properties
# ─────────────────────────────────────────────────

def adj_to_sets(cg, m):
    return [set(j for j in range(m) if cg[i][j]) for i in range(m)]

def max_clique(adj_sets, m):
    if m == 0:
        return 0
    comp_nbr = [0]*m
    for i in range(m):
        for j in range(m):
            if j != i and j not in adj_sets[i]:
                comp_nbr[i] |= (1<<j)
    best = 1
    for mask in range(1, 1<<m):
        sz = bin(mask).count('1')
        if sz <= best:
            continue
        ok = True
        bits = mask
        while bits:
            v = (bits & -bits).bit_length() - 1
            if comp_nbr[v] & mask:
                ok = False
                break
            bits &= bits - 1
        if ok:
            best = sz
    return best

def max_indep(adj_sets, m):
    if m == 0:
        return 0
    nbr = [0]*m
    for i in range(m):
        for j in adj_sets[i]:
            nbr[i] |= (1<<j)
    best = 0
    for mask in range(1<<m):
        sz = bin(mask).count('1')
        if sz <= best:
            continue
        ok = True
        bits = mask
        while bits:
            v = (bits & -bits).bit_length() - 1
            if nbr[v] & mask & ((1<<v)-1):
                ok = False
                break
            bits &= bits - 1
        if ok:
            best = sz
    return best

def connected_comps(adj_sets, m):
    visited = [False]*m
    comps = []
    for s in range(m):
        if visited[s]:
            continue
        c = []
        stack = [s]
        visited[s] = True
        while stack:
            u = stack.pop()
            c.append(u)
            for v in adj_sets[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
        comps.append(c)
    return comps

def has_odd_hole(adj_sets, m):
    if m < 5:
        return False
    for length in range(5, min(m+1, 10), 2):
        for verts in itertools.combinations(range(m), length):
            ok = True
            for v in verts:
                if sum(1 for u in verts if u != v and u in adj_sets[v]) != 2:
                    ok = False
                    break
            if not ok:
                continue
            vis = {verts[0]}
            stk = [verts[0]]
            while stk:
                u = stk.pop()
                for w in verts:
                    if w not in vis and w in adj_sets[u]:
                        vis.add(w)
                        stk.append(w)
            if len(vis) == length:
                return True
    return False

def has_induced_cycle_geq4(adj_sets, m):
    if m < 4:
        return False
    for length in range(4, min(m+1, 9)):
        for verts in itertools.combinations(range(m), length):
            ok = True
            for v in verts:
                if sum(1 for u in verts if u != v and u in adj_sets[v]) != 2:
                    ok = False
                    break
            if not ok:
                continue
            vis = {verts[0]}
            stk = [verts[0]]
            while stk:
                u = stk.pop()
                for w in verts:
                    if w not in vis and w in adj_sets[u]:
                        vis.add(w)
                        stk.append(w)
            if len(vis) == length:
                return True
    return False

def is_perfect(adj_sets, m):
    if m <= 4:
        return True
    if has_odd_hole(adj_sets, m):
        return False
    comp = [set() for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if j not in adj_sets[i]:
                comp[i].add(j)
                comp[j].add(i)
    return not has_odd_hole(comp, m)

def find_roots(coeffs):
    try:
        import numpy as np
        if len(coeffs) <= 1:
            return []
        return np.roots(list(reversed(coeffs)))
    except ImportError:
        return None

# ─────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────

def sep(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")

def main():
    print("Omega(T) Perfectness Implications for OCF")
    print("="*70)
    print("Using directed odd cycles from tournament_lib")

    numpy_ok = True
    try:
        import numpy as np
    except ImportError:
        numpy_ok = False
        print("WARNING: numpy not available")

    # ═══════════════════════════════════════════════
    # n=5 exhaustive
    # ═══════════════════════════════════════════════
    sep("SECTION 1: Exhaustive n=5")

    n5_polys = Counter()
    n5_ocf_fail = 0
    n5_perf_fail = 0
    n5_chord_fail = 0
    n5_alpha_dist = Counter()
    n5_omega_dist = Counter()
    n5_comp_dist = Counter()
    n5_total = 0
    n5_detailed = []  # store (H, poly, alpha, omega, score_seq) for each

    for T in all_tournaments(5):
        n5_total += 1
        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        ss = score_sequence(T)

        if m == 0:
            poly = [1]
            I2 = 1
            alpha_v, omega_v = 0, 0
            is_perf, is_chord = True, True
            ncomp = 0
            comp_sizes = ()
        else:
            cg = conflict_graph(cycles)
            m = len(cg)
            poly = independence_polynomial_coeffs(cg, m)
            I2 = eval_poly(poly, 2)
            adj_s = adj_to_sets(cg, m)
            alpha_v = max_indep(adj_s, m)
            omega_v = max_clique(adj_s, m)
            is_perf = is_perfect(adj_s, m)
            is_chord = not has_induced_cycle_geq4(adj_s, m)
            comps = connected_comps(adj_s, m)
            ncomp = len(comps)
            comp_sizes = tuple(sorted([len(c) for c in comps], reverse=True))

        if H != I2:
            n5_ocf_fail += 1
        if not is_perf:
            n5_perf_fail += 1
        if not is_chord:
            n5_chord_fail += 1
        n5_polys[tuple(poly)] += 1
        n5_alpha_dist[alpha_v] += 1
        n5_omega_dist[omega_v] += 1
        n5_comp_dist[comp_sizes] += 1
        n5_detailed.append((H, tuple(poly), alpha_v, omega_v, ss, m))

    print(f"Total: {n5_total}")
    print(f"OCF failures: {n5_ocf_fail}")
    print(f"Non-perfect: {n5_perf_fail}")
    print(f"Non-chordal: {n5_chord_fail}")
    print(f"\nDistinct I(Omega,x): {len(n5_polys)}")
    for p, cnt in n5_polys.most_common(15):
        s = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(p) if c > 0)
        print(f"  [{cnt:4d}x] I(x) = {s}  => I(2) = {eval_poly(list(p), 2)}")
    print(f"\nalpha: {dict(sorted(n5_alpha_dist.items()))}")
    print(f"omega: {dict(sorted(n5_omega_dist.items()))}")
    print(f"\nComponent sizes: {dict(n5_comp_dist.most_common(10))}")

    # ═══════════════════════════════════════════════
    # n=6 exhaustive (OCF only uses fast evaluation)
    # ═══════════════════════════════════════════════
    sep("SECTION 2: Exhaustive n=6 (OCF verification + stats)")

    n6_ocf_fail = 0
    n6_total = 0
    n6_cycle_count_dist = Counter()
    n6_H_by_score = defaultdict(list)
    n6_alpha_by_score = defaultdict(set)
    n6_cycle_by_score = defaultdict(set)
    n6_poly_samples = {}  # poly_key -> (H, cycles, cg) for distinct polys
    n6_m_distribution = Counter()

    for T in all_tournaments(6):
        n6_total += 1
        if n6_total % 5000 == 0:
            print(f"  ... {n6_total}", flush=True)

        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        m_cyc = len(cycles)
        ss = score_sequence(T)

        if m_cyc == 0:
            I2 = 1
        else:
            I2 = independence_poly_at_fast(cycles, 2)

        if H != I2:
            n6_ocf_fail += 1
            if n6_ocf_fail <= 3:
                print(f"  FAIL #{n6_ocf_fail}: H={H}, I(2)={I2}, m={m_cyc}")

        n6_cycle_count_dist[m_cyc] += 1
        n6_H_by_score[ss].append(H)
        n6_m_distribution[m_cyc] += 1

    print(f"Total: {n6_total}")
    print(f"OCF failures: {n6_ocf_fail}")
    print(f"\n#directed_cycles distribution:")
    for m_val in sorted(n6_cycle_count_dist.keys()):
        print(f"  m={m_val:3d}: {n6_cycle_count_dist[m_val]:5d} tournaments")

    # ═══════════════════════════════════════════════
    # n=6 detailed analysis on a SAMPLE (for expensive properties)
    # ═══════════════════════════════════════════════
    sep("SECTION 3: n=6 detailed sample (first 2000 + all interesting)")

    # Sample: first 2000 tournaments + stratified by score sequence
    n6_sample_polys = Counter()
    n6_sample_perf = 0
    n6_sample_perf_fail = 0
    n6_sample_chord_fail = 0
    n6_sample_alpha = Counter()
    n6_sample_omega = Counter()
    n6_sample_total = 0
    n6_sample_lc_fail = 0
    n6_sample_uni_fail = 0
    n6_sample_comp_dist = Counter()
    n6_sample_detailed = []

    for T in all_tournaments(6):
        n6_sample_total += 1
        if n6_sample_total > 2000:
            break

        H = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        m = len(cycles)
        ss = score_sequence(T)

        if m == 0:
            poly = [1]
            alpha_v, omega_v = 0, 0
            n6_sample_perf += 1
            ncomp = 0
            comp_sizes = ()
            is_chord = True
        else:
            cg = conflict_graph(cycles)
            m = len(cg)
            if m <= 22:
                poly = independence_polynomial_coeffs(cg, m)
            else:
                poly = None
            adj_s = adj_to_sets(cg, m)
            if m <= 20:
                alpha_v = max_indep(adj_s, m)
                omega_v = max_clique(adj_s, m)
            else:
                alpha_v = None
                omega_v = None
            if m <= 12:
                perf = is_perfect(adj_s, m)
                is_chord = not has_induced_cycle_geq4(adj_s, m)
                if perf:
                    n6_sample_perf += 1
                else:
                    n6_sample_perf_fail += 1
                if not is_chord:
                    n6_sample_chord_fail += 1
            else:
                perf = None
                is_chord = None
            comps = connected_comps(adj_s, m)
            ncomp = len(comps)
            comp_sizes = tuple(sorted([len(c) for c in comps], reverse=True))

        if poly is not None:
            pk = tuple(poly)
            n6_sample_polys[pk] += 1
            # Log-concavity
            if len(poly) > 2:
                for i in range(1, len(poly)-1):
                    if poly[i]**2 < poly[i-1]*poly[i+1]:
                        n6_sample_lc_fail += 1
                        break
                peaked = False
                for i in range(1, len(poly)):
                    if poly[i] < poly[i-1]:
                        peaked = True
                    elif poly[i] > poly[i-1] and peaked:
                        n6_sample_uni_fail += 1
                        break

        if alpha_v is not None:
            n6_sample_alpha[alpha_v] += 1
        if omega_v is not None:
            n6_sample_omega[omega_v] += 1
        n6_sample_comp_dist[comp_sizes] += 1

        n6_sample_detailed.append((H, poly, alpha_v, omega_v, ss, m))

    print(f"Sample size: {n6_sample_total}")
    print(f"Non-perfect: {n6_sample_perf_fail}")
    print(f"Non-chordal: {n6_sample_chord_fail}")
    print(f"Log-concavity failures: {n6_sample_lc_fail}")
    print(f"Unimodality failures: {n6_sample_uni_fail}")

    print(f"\nDistinct I(Omega,x): {len(n6_sample_polys)}")
    for p, cnt in n6_sample_polys.most_common(20):
        s = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(p) if c > 0)
        print(f"  [{cnt:4d}x] I(x) = {s}  => I(2) = {eval_poly(list(p), 2)}")

    print(f"\nalpha: {dict(sorted(n6_sample_alpha.items()))}")
    print(f"omega: {dict(sorted(n6_sample_omega.items()))}")
    print(f"\nComponent sizes (top 10):")
    for cs, cnt in n6_sample_comp_dist.most_common(10):
        print(f"  {cs}: {cnt}")

    # alpha*omega >= m
    violations = sum(1 for H, poly, a, o, ss, m in n6_sample_detailed
                     if a is not None and o is not None and m > 0 and a*o < m)
    checked = sum(1 for H, poly, a, o, ss, m in n6_sample_detailed
                  if a is not None and o is not None and m > 0)
    print(f"\nalpha*omega >= m violations: {violations}/{checked}")

    # ═══════════════════════════════════════════════
    # Root analysis
    # ═══════════════════════════════════════════════
    sep("SECTION 4: Roots of I(Omega(T), x)")

    if numpy_ok:
        for label, polys_counter in [("n=5", n5_polys), ("n=6", n6_sample_polys)]:
            print(f"\n--- {label} ---")
            all_real = True
            all_neg = True
            mags = []
            neg_reals = []

            for pk in polys_counter:
                if len(pk) <= 1:
                    continue
                roots = find_roots(list(pk))
                if roots is None:
                    continue
                for r in roots:
                    mags.append(abs(r))
                    if abs(r.imag) > 1e-8:
                        all_real = False
                    if r.real > 1e-8:
                        all_neg = False
                    if abs(r.imag) < 1e-8 and r.real < 0:
                        neg_reals.append(r.real)

            print(f"  Distinct polynomials: {len(polys_counter)}")
            print(f"  All roots real: {all_real}")
            print(f"  All roots real & negative: {all_neg}")
            if mags:
                print(f"  |root| range: [{min(mags):.6f}, {max(mags):.6f}]")
            if neg_reals:
                print(f"  Neg real roots: closest to 0 = {max(neg_reals):.6f}, "
                      f"farthest = {min(neg_reals):.6f}")
                print(f"  All |neg root| > 2: {all(abs(r)>2 for r in neg_reals)}")

        # Detail for n=5
        print(f"\n--- n=5 root details ---")
        for pk in sorted(n5_polys.keys(), key=len):
            if len(pk) <= 1:
                continue
            roots = find_roots(list(pk))
            if roots is None or len(roots) == 0:
                continue
            s = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(pk) if c > 0)
            rs = []
            for r in sorted(roots, key=lambda r: r.real):
                if abs(r.imag) < 1e-8:
                    rs.append(f"{r.real:.6f}")
                else:
                    rs.append(f"({r.real:.4f}{r.imag:+.4f}i)")
            print(f"  I(x) = {s}")
            print(f"    Roots: {rs}")

        # n=6 interesting cases
        print(f"\n--- n=6 root details (degree >= 2 polynomials) ---")
        shown = 0
        for pk in sorted(n6_sample_polys.keys(), key=len, reverse=True):
            if len(pk) <= 2:
                continue
            roots = find_roots(list(pk))
            if roots is None or len(roots) == 0:
                continue
            s = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(pk) if c > 0)
            rs = []
            for r in sorted(roots, key=lambda r: r.real):
                if abs(r.imag) < 1e-8:
                    rs.append(f"{r.real:.6f}")
                else:
                    rs.append(f"({r.real:.4f}{r.imag:+.4f}i)")
            print(f"  I(x) = {s}  [{n6_sample_polys[pk]}x]")
            print(f"    Roots: {rs}")
            shown += 1
            if shown >= 15:
                break
    else:
        print("  (skipped)")

    # ═══════════════════════════════════════════════
    # Modular arithmetic
    # ═══════════════════════════════════════════════
    sep("SECTION 5: H(T) modular structure via OCF")

    for label, detailed in [("n=5", n5_detailed), ("n=6 sample", n6_sample_detailed)]:
        print(f"\n{label}:")
        Hs = [d[0] for d in detailed]
        polys = [d[1] for d in detailed if d[1] is not None]

        mod2 = Counter(H % 2 for H in Hs)
        print(f"  H mod 2: {dict(mod2)} (Redei)")

        mod4 = Counter(H % 4 for H in Hs)
        print(f"  H mod 4: {dict(sorted(mod4.items()))}")

        # H mod 4 = (1 + 2*alpha_1) mod 4?
        check4 = True
        for H, poly, *_ in detailed:
            if poly is None:
                continue
            a1 = poly[1] if len(poly) > 1 else 0
            if H % 4 != (1 + 2*a1) % 4:
                check4 = False
                break
        print(f"  H mod 4 = 1 + 2*alpha_1 mod 4: {check4}")

        # H mod 8
        mod8 = Counter(H % 8 for H in Hs)
        print(f"  H mod 8: {dict(sorted(mod8.items()))}")
        check8 = True
        for H, poly, *_ in detailed:
            if poly is None:
                continue
            a1 = poly[1] if len(poly) > 1 else 0
            a2 = poly[2] if len(poly) > 2 else 0
            if H % 8 != (1 + 2*a1 + 4*a2) % 8:
                check8 = False
                break
        print(f"  H mod 8 = 1 + 2*a1 + 4*a2 mod 8: {check8}")

    # ═══════════════════════════════════════════════
    # Score sequence correlations (n=5)
    # ═══════════════════════════════════════════════
    sep("SECTION 6: Score sequence vs Omega invariants")

    print("\nn=5:")
    by_ss5 = defaultdict(list)
    for H, poly, alpha, omega, ss, m in n5_detailed:
        by_ss5[ss].append((H, poly, alpha, omega, m))

    for ss in sorted(by_ss5.keys()):
        items = by_ss5[ss]
        Hs = [i[0] for i in items]
        ms = [i[4] for i in items]
        alphas = [i[2] for i in items]
        omegas = [i[3] for i in items]
        a1s = [i[1][1] if len(i[1]) > 1 else 0 for i in items]
        print(f"  {ss} ({len(items)}x): H={sorted(set(Hs))}, m={sorted(set(ms))}, "
              f"alpha={sorted(set(alphas))}, omega={sorted(set(omegas))}, "
              f"alpha_1={sorted(set(a1s))}")

    print("\nn=6 (sample, by score seq):")
    by_ss6 = defaultdict(list)
    for H, poly, alpha, omega, ss, m in n6_sample_detailed:
        by_ss6[ss].append((H, alpha, omega, m))

    for ss in sorted(by_ss6.keys()):
        items = by_ss6[ss]
        Hs = [i[0] for i in items]
        alphas = [i[1] for i in items if i[1] is not None]
        omegas = [i[2] for i in items if i[2] is not None]
        ms = [i[3] for i in items]
        print(f"  {ss} ({len(items)}x): "
              f"H=[{min(Hs)},{max(Hs)}] "
              f"m=[{min(ms)},{max(ms)}] "
              f"alpha={sorted(set(alphas)) if alphas else '?'} "
              f"omega=[{min(omegas) if omegas else '?'},{max(omegas) if omegas else '?'}]")

    # ═══════════════════════════════════════════════
    # I at special values
    # ═══════════════════════════════════════════════
    sep("SECTION 7: I(Omega, x) at special values")

    for label, detailed in [("n=5", n5_detailed)]:
        print(f"\n{label}:")
        for x_val in [-1, 1, 3]:
            vals = [eval_poly(list(poly), x_val) for H, poly, *_ in detailed if poly]
            vc = Counter(vals)
            print(f"  I(Omega, {x_val:2d}): {dict(sorted(vc.items()))}")

    # For n=6 sample
    print(f"\nn=6 sample:")
    for x_val in [-1, 1, 3]:
        vals = [eval_poly(list(poly), x_val) for H, poly, *_ in n6_sample_detailed
                if poly is not None]
        if vals:
            vc = Counter(vals)
            print(f"  I(Omega, {x_val:2d}): range=[{min(vals)},{max(vals)}], "
                  f"mean={sum(vals)/len(vals):.1f}, top: {vc.most_common(5)}")

    # ═══════════════════════════════════════════════
    # Dominant term
    # ═══════════════════════════════════════════════
    sep("SECTION 8: Dominant term alpha_k * 2^k")

    for label, detailed in [("n=5", n5_detailed), ("n=6 sample", n6_sample_detailed)]:
        dom = Counter()
        for H, poly, *_ in detailed:
            if poly is None:
                continue
            terms = [(k, poly[k] * 2**k) for k in range(len(poly)) if poly[k] > 0]
            if terms:
                dom[max(terms, key=lambda t: t[1])[0]] += 1
        print(f"\n{label}: dominant k: {dict(sorted(dom.items()))}")

    # ═══════════════════════════════════════════════
    # NEW: Coefficient ratios for n=6
    # ═══════════════════════════════════════════════
    sep("SECTION 9: Coefficient density alpha_k / C(m,k)")

    for target_m in sorted(set(m for *_, m in n6_sample_detailed if m > 0)):
        ratios = defaultdict(list)
        for H, poly, *_, m in n6_sample_detailed:
            if poly is None or m != target_m:
                continue
            for k in range(1, len(poly)):
                c_mk = comb(m, k)
                if c_mk > 0:
                    ratios[k].append(poly[k] / c_mk)
        if ratios and target_m <= 20:
            print(f"\n  m={target_m}:")
            for k in sorted(ratios.keys()):
                if ratios[k]:
                    v = ratios[k]
                    print(f"    k={k}: {sum(v)/len(v):.4f} [{min(v):.4f},{max(v):.4f}] "
                          f"(n={len(v)})")

    # ═══════════════════════════════════════════════
    # Max independent set patterns
    # ═══════════════════════════════════════════════
    sep("SECTION 10: Max independent set cycle-size patterns (n=5)")

    mis_patterns = Counter()
    for T in all_tournaments(5):
        cycles = find_odd_cycles(T)
        m = len(cycles)
        if m <= 1:
            continue
        cg = conflict_graph(cycles)
        adj_s = adj_to_sets(cg, m)
        nbr = [0]*m
        for i in range(m):
            for j in range(m):
                if cg[i][j]:
                    nbr[i] |= (1<<j)

        alpha = 0
        max_sets = []
        for mask in range(1<<m):
            sz = bin(mask).count('1')
            if sz < alpha:
                continue
            ok = True
            bits = mask
            while bits:
                v = (bits & -bits).bit_length() - 1
                if nbr[v] & mask & ((1<<v)-1):
                    ok = False
                    break
                bits &= bits-1
            if ok:
                if sz > alpha:
                    alpha = sz
                    max_sets = []
                max_sets.append(mask)

        for mask in max_sets:
            sizes = tuple(sorted([len(set(cycles[i])) for i in range(m)
                                  if mask & (1<<i)], reverse=True))
            mis_patterns[sizes] += 1

    print("Cycle sizes in maximum independent sets:")
    for pat, cnt in mis_patterns.most_common(15):
        print(f"  {pat}: {cnt}")

    # ═══════════════════════════════════════════════
    # Summary
    # ═══════════════════════════════════════════════
    sep("SUMMARY OF FINDINGS")
    print(f"""
METHODOLOGY: Uses tournament_lib with DIRECTED odd cycles (not vertex sets).
Each directed cycle is a separate vertex of Omega(T).

1. OCF: H(T) = I(Omega(T), 2)
   n=5: {n5_ocf_fail}/{n5_total} failures
   n=6: {n6_ocf_fail}/{n6_total} failures

2. PERFECTNESS of Omega(T):
   n=5: {n5_perf_fail} non-perfect out of {n5_total}
   n=6 sample: {n6_sample_perf_fail} non-perfect

3. CHORDALITY:
   n=5: {n5_chord_fail} non-chordal
   n=6 sample: {n6_sample_chord_fail} non-chordal

4. ROOTS of I(Omega(T), x): See Section 4.
   If all roots are real & negative, then I(x) > 0 for x > 0,
   giving a "positivity proof" of H(T) >= 1 (Redei).

5. LOG-CONCAVITY: n=6 failures: {n6_sample_lc_fail}
   UNIMODALITY: n=6 failures: {n6_sample_uni_fail}

6. MODULAR ARITHMETIC:
   H(T) mod 2 = 1 (Redei) follows from alpha_0 = 1.
   H(T) mod 4 = 1 + 2*(#directed_odd_cycles) mod 4.
   These give a "modular cascade" proof strategy for OCF.

7. For PERFECT Omega(T):
   - clique_cover = alpha (min # cliques to cover all cycles)
   - chi = omega (chromatic number = max clique)
   - These equalities give tight bounds relating cycle structure to counting.
""")


if __name__ == "__main__":
    main()
