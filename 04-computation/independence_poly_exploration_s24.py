#!/usr/bin/env python3
"""
opus-2026-04-05-S24: Deep exploration of I(Ω(T), x) as a polynomial.

The independence polynomial of the odd-cycle conflict graph, evaluated at x=2,
gives H(T). But the POLYNOMIAL itself — its zeros, factorizations, coefficients —
contains much more information. This script explores that structure.
"""

import numpy as np
from itertools import combinations, permutations
from collections import defaultdict, Counter
import sys

# ─── Tournament utilities ───

def all_tournaments(n):
    """Generate all labeled tournaments on n vertices."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        T = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield [row[:] for row in T]

def tournament_to_canon(T):
    """Canonical form under vertex relabeling."""
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def find_all_directed_odd_cycles(T):
    """Find ALL directed odd cycles in tournament T.
    Returns list of cycles, each as a tuple of vertices in cycle order.
    Canonical form: start from minimum vertex in the cycle.
    """
    n = len(T)
    cycles = set()

    for length in range(3, n+1, 2):  # odd lengths only
        for verts in combinations(range(n), length):
            # For these vertices, find all directed Hamiltonian cycles
            # A cycle (v0, v1, ..., v_{k-1}) means v0→v1→...→v_{k-1}→v0
            # Fix v0 = verts[0], try all permutations of the rest
            v0 = verts[0]
            rest = verts[1:]
            for perm in permutations(rest):
                path = (v0,) + perm
                # Check if directed cycle
                is_cycle = True
                for i in range(length):
                    if T[path[i]][path[(i+1) % length]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    # Canonical form: start from minimum vertex
                    min_v = min(path)
                    min_idx = path.index(min_v)
                    canon = tuple(path[(min_idx + j) % length] for j in range(length))
                    cycles.add(canon)

    return list(cycles)

def build_conflict_graph(cycles):
    """Build conflict graph: cycles adjacent iff they share a vertex."""
    k = len(cycles)
    adj = [[False]*k for _ in range(k)]
    vert_sets = [set(c) for c in cycles]
    for i in range(k):
        for j in range(i+1, k):
            if vert_sets[i] & vert_sets[j]:
                adj[i][j] = adj[j][i] = True
    return adj

def independence_polynomial_coeffs(adj):
    """Compute independence polynomial I(G, x) = sum_k alpha_k x^k."""
    k = len(adj)
    if k == 0:
        return [1]

    coeffs = [0] * (k + 1)
    for mask in range(2**k):
        verts = [i for i in range(k) if (mask >> i) & 1]
        is_indep = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            coeffs[len(verts)] += 1

    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def hamiltonian_path_count(T):
    """Count Hamiltonian paths via DP (bitmask)."""
    n = len(T)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def score_sequence(T):
    n = len(T)
    return tuple(sorted(sum(T[i]) for i in range(n)))

# ─── Main exploration ───

def explore_n(n, verbose=True):
    print(f"\n{'='*70}")
    print(f"  EXPLORING n = {n}")
    print(f"{'='*70}")

    iso_classes = {}
    for T in all_tournaments(n):
        canon = tournament_to_canon(T)
        if canon not in iso_classes:
            iso_classes[canon] = T

    print(f"Isomorphism classes: {len(iso_classes)}")

    results = []
    for idx, (canon, T) in enumerate(iso_classes.items()):
        H = hamiltonian_path_count(T)
        cycles = find_all_directed_odd_cycles(T)
        adj = build_conflict_graph(cycles)
        coeffs = independence_polynomial_coeffs(adj)

        H_ocf = eval_poly(coeffs, 2)
        if H != H_ocf:
            print(f"  OCF MISMATCH at class {idx}: H={H}, I(Ω,2)={H_ocf}")
            print(f"    scores={score_sequence(T)}, #cycles={len(cycles)}")
            print(f"    cycle lengths: {sorted(len(c) for c in cycles)}")
            print(f"    coeffs={coeffs}")
            continue

        val_at = {}
        for x_val in [-3, -2, -1, 0, 1, 2, 3, 4, 5]:
            val_at[x_val] = eval_poly(coeffs, x_val)

        if len(coeffs) > 1:
            np_coeffs = list(reversed(coeffs))
            zeros = np.roots(np_coeffs)
        else:
            zeros = np.array([])

        all_real = all(abs(z.imag) < 1e-10 for z in zeros) if len(zeros) > 0 else True
        all_neg_real = all(abs(z.imag) < 1e-10 and z.real < 1e-10 for z in zeros) if len(zeros) > 0 else True

        scores = score_sequence(T)
        num_3cycles = sum(1 for c in cycles if len(c) == 3)
        num_5cycles = sum(1 for c in cycles if len(c) == 5)
        num_7cycles = sum(1 for c in cycles if len(c) == 7)

        results.append({
            'T': T, 'H': H, 'scores': scores,
            'num_cycles': len(cycles), 'c3': num_3cycles,
            'c5': num_5cycles, 'c7': num_7cycles,
            'coeffs': coeffs, 'degree': len(coeffs) - 1,
            'vals': val_at, 'zeros': zeros,
            'all_real': all_real, 'all_neg_real': all_neg_real,
        })

    if not results:
        return results

    # ── Reports ──
    degrees = Counter(r['degree'] for r in results)
    print(f"\nPolynomial degrees: {dict(sorted(degrees.items()))}")

    n_real = sum(1 for r in results if r['all_real'])
    n_neg = sum(1 for r in results if r['all_neg_real'])
    print(f"All zeros real: {n_real}/{len(results)}")
    print(f"All zeros negative real: {n_neg}/{len(results)}")

    all_zeros = []
    for r in results:
        all_zeros.extend(r['zeros'])
    if all_zeros:
        real_zeros = [z.real for z in all_zeros if abs(z.imag) < 1e-10]
        complex_zeros = [z for z in all_zeros if abs(z.imag) >= 1e-10]
        print(f"\nTotal zeros: {len(all_zeros)} ({len(real_zeros)} real, {len(complex_zeros)} complex)")
        if real_zeros:
            print(f"  Real zero range: [{min(real_zeros):.6f}, {max(real_zeros):.6f}]")

    # Special values
    print(f"\n── Special values ──")
    for x_val in [-3, -2, -1, 1, 3, 4, 5]:
        vals = sorted(set(r['vals'][x_val] for r in results))
        print(f"  I(Ω, {x_val:+d}): {vals[:20]}{'...' if len(vals)>20 else ''}")

    # Achievable H
    achievable_H = sorted(set(r['H'] for r in results))
    print(f"\n── H values at n={n}: {achievable_H}")
    max_H = max(achievable_H)
    all_odd = set(range(1, max_H + 1, 2))
    missing = sorted(all_odd - set(achievable_H))
    print(f"  Missing odd: {missing}")

    # Polynomial as invariant
    poly_to_classes = defaultdict(list)
    for r in results:
        poly_to_classes[tuple(r['coeffs'])].append(r)
    print(f"\n── Polynomial as invariant: {len(poly_to_classes)} distinct out of {len(results)} classes")

    h_groups = defaultdict(list)
    for r in results:
        h_groups[r['H']].append(r)
    for H_val, group in sorted(h_groups.items()):
        if len(group) > 1:
            polys = set(tuple(r['coeffs']) for r in group)
            if len(polys) > 1:
                print(f"  H={H_val}: {len(polys)} distinct polynomials!")
                for r in group:
                    print(f"    scores={r['scores']}, I(x)={r['coeffs']}")

    # Detailed output for small n
    if verbose and n <= 6:
        print(f"\n── Full polynomial table ──")
        for r in sorted(results, key=lambda r: r['H']):
            poly_str = ' + '.join(
                f"{c}x^{k}" if k > 0 else str(c)
                for k, c in enumerate(r['coeffs']) if c
            )
            zeros_str = ""
            if len(r['zeros']) > 0:
                zeros_real = sorted([z.real for z in r['zeros'] if abs(z.imag) < 1e-10])
                zeros_str = f"  zeros(real)={[round(z,4) for z in zeros_real]}"
            print(f"  H={r['H']:3d} c3={r['c3']:2d} scores={r['scores']} "
                  f"I(x)={poly_str}{zeros_str}")

    # I(Ω, -1) — Euler char
    print(f"\n── I(Ω, -1) values: {sorted(Counter(r['vals'][-1] for r in results).items())}")

    # I'(Ω, 2) — derivative at fugacity 2
    print(f"\n── Derivative I'(Ω, 2) ──")
    for r in results:
        if len(r['coeffs']) > 1:
            r['deriv_at_2'] = sum(k * c * 2**(k-1) for k, c in enumerate(r['coeffs']) if k > 0)
        else:
            r['deriv_at_2'] = 0
    if verbose and n <= 6:
        for r in sorted(results, key=lambda r: r['H']):
            ratio = r['deriv_at_2'] / r['H'] if r['H'] > 0 else 0
            print(f"  H={r['H']:3d}, I'(2)={r['deriv_at_2']:5d}, I'(2)/H={ratio:.4f}")

    # NEW: Check I(Ω, 3) mod something
    print(f"\n── I(Ω, 3) values: {sorted(set(r['vals'][3] for r in results))}")
    print(f"  I(Ω,3) mod 2: {sorted(set(r['vals'][3] % 2 for r in results))}")
    print(f"  I(Ω,3) mod 3: {sorted(set(r['vals'][3] % 3 for r in results))}")

    return results

# ─── Run ───
all_results = {}
for n in [3, 4, 5, 6]:
    all_results[n] = explore_n(n, verbose=(n <= 6))
    sys.stdout.flush()

# ─── Cross-n analysis ───
print(f"\n{'='*70}")
print(f"  CROSS-n ANALYSIS")
print(f"{'='*70}")

print(f"\nAchievable H values:")
for n in [3, 4, 5, 6]:
    H_vals = sorted(set(r['H'] for r in all_results[n]))
    print(f"  n={n}: {H_vals}")

print(f"\nMissing odd H values:")
for n in [3, 4, 5, 6]:
    H_vals = set(r['H'] for r in all_results[n])
    max_H = max(H_vals)
    missing = sorted(set(range(1, max_H+1, 2)) - H_vals)
    print(f"  n={n}: {missing}")

all_achieved = set()
for n in [3, 4, 5, 6]:
    all_achieved |= set(r['H'] for r in all_results[n])
never_achieved = sorted(set(range(1, max(all_achieved)+1, 2)) - all_achieved)
print(f"\nNever achieved at any n=3..6: {never_achieved[:30]}")

# Zero structure
print(f"\n── Zero structure ──")
for n in [3, 4, 5, 6]:
    all_z = []
    for r in all_results[n]:
        all_z.extend(r['zeros'])
    real_z = [z.real for z in all_z if abs(z.imag) < 1e-10]
    n_real = sum(1 for r in all_results[n] if r['all_real'])
    n_neg = sum(1 for r in all_results[n] if r['all_neg_real'])
    total = len(all_results[n])
    rng = f"[{min(real_z):.4f}, {max(real_z):.4f}]" if real_z else "[]"
    print(f"  n={n}: {n_real}/{total} all-real, {n_neg}/{total} all-neg-real, range={rng}")

# Ratios I(3)/I(2)
print(f"\n── I(Ω,3)/I(Ω,2) ratios ──")
for n in [3, 4, 5, 6]:
    ratios = sorted(set(round(r['vals'][3]/r['H'], 4) for r in all_results[n]))
    print(f"  n={n}: {ratios[:15]}")

# NEW: Look at I(Ω, x) mod 2 for various x
print(f"\n── I(Ω, x) mod 2 universality ──")
for x_val in [1, 3, 4, 5]:
    for n in [3, 4, 5, 6]:
        vals_mod2 = set(r['vals'][x_val] % 2 for r in all_results[n])
        print(f"  n={n}, x={x_val}: I mod 2 ∈ {vals_mod2}", end="")
    print()

# NEW: Is I(Ω, 3) always odd? (generalized Rédei?)
print(f"\n── Is I(Ω(T), 3) always odd? ──")
for n in [3, 4, 5, 6]:
    all_odd = all(r['vals'][3] % 2 == 1 for r in all_results[n])
    print(f"  n={n}: {'YES' if all_odd else 'NO'}")

print(f"\n── Is I(Ω(T), x) always odd for odd x? ──")
for x_val in [1, 3, 5]:
    for n in [3, 4, 5, 6]:
        all_odd = all(r['vals'][x_val] % 2 == 1 for r in all_results[n])
        print(f"  n={n}, x={x_val}: {'YES' if all_odd else 'NO'}")

# NEW: What about I(Ω, -1)? Connection to matching polynomial?
print(f"\n── I(Ω, -1) distribution ──")
for n in [3, 4, 5, 6]:
    dist = Counter(r['vals'][-1] for r in all_results[n])
    print(f"  n={n}: {dict(sorted(dist.items()))}")

# NEW: Factorization over Z
print(f"\n── Polynomial factorization check ──")
for n in [3, 4, 5, 6]:
    for r in all_results[n]:
        coeffs = r['coeffs']
        if len(coeffs) <= 2:
            continue
        # Check if the polynomial factors as product of linear terms
        zeros = r['zeros']
        # Check if product of (x - z_i) for real rational zeros
        rational_zeros = []
        for z in zeros:
            if abs(z.imag) < 1e-10:
                zr = z.real
                # Check if rational with small denominator
                for d in range(1, 20):
                    if abs(zr * d - round(zr * d)) < 1e-8:
                        rational_zeros.append((int(round(zr * d)), d))
                        break
        if len(rational_zeros) == len(zeros) and len(zeros) > 1:
            print(f"  n={n} H={r['H']}: ALL RATIONAL zeros! {rational_zeros}")

print("\n\nDONE.")
