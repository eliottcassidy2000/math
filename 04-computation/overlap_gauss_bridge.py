#!/usr/bin/env python3
"""
overlap_gauss_bridge.py -- Connect overlap weights to Gauss layering theorem

KEY QUESTIONS:
1. Does disj3 = -c5_dir/2 + const hold at ALL primes? (found at p=11)
2. How do Walsh degree-2 of alpha_2 relate to Walsh degree-2 of c5?
3. What is the 4-fold Gauss sum structure encoding about overlaps?
4. Is there a simple formula for disjoint 3-cycle pair count?

Author: kind-pasteur-2026-03-12-S60
"""

import math
import cmath
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ck(A, p, k):
    """Count directed k-cycles."""
    total = 0
    for subset in combinations(range(p), k):
        verts = list(subset)
        n = len(verts)
        if n == 3:
            a, b, c = verts
            total += A[a][b] * A[b][c] * A[c][a]
            total += A[a][c] * A[c][b] * A[b][a]
        else:
            start = 0
            dp = {(1 << start, start): 1}
            for mask in range(1, 1 << n):
                if not (mask & (1 << start)):
                    continue
                for v in range(n):
                    if not (mask & (1 << v)):
                        continue
                    key = (mask, v)
                    if key not in dp or dp[key] == 0:
                        continue
                    cnt = dp[key]
                    for w in range(n):
                        if mask & (1 << w):
                            continue
                        if A[verts[v]][verts[w]]:
                            nkey = (mask | (1 << w), w)
                            dp[nkey] = dp.get(nkey, 0) + cnt
            full = (1 << n) - 1
            for v in range(n):
                if v == start:
                    continue
                key = (full, v)
                if key in dp and dp[key] > 0:
                    if A[verts[v]][verts[start]]:
                        total += dp[key]
    return total


def count_c3_vertex_sets(A, p):
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    return list(set(c3_sets))


def disjoint_pair_counts(cycle_sets):
    n = len(cycle_sets)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            if not (cycle_sets[i] & cycle_sets[j]):
                count += 1
    return count


def all_orientations_circ(p):
    m = (p - 1) // 2
    orientations = []
    for bits in range(1 << m):
        S = []
        for i in range(m):
            chord = i + 1
            if (bits >> i) & 1:
                S.append(chord)
            else:
                S.append(p - chord)
        orientations.append((bits, S))
    return orientations


def walsh_degree2(values, m):
    """Compute degree-2 Walsh coefficients h_hat[{a,b}]."""
    h2 = {}
    N = 1 << m
    for a in range(m):
        for b in range(a+1, m):
            total = 0
            for sigma_bits in range(N):
                sigma = tuple(1 if (sigma_bits >> i) & 1 else -1 for i in range(m))
                val = values.get(sigma, 0)
                total += val * sigma[a] * sigma[b]
            h2[(a, b)] = total / N
    return h2


def walsh_degree4(values, m):
    """Compute degree-4 Walsh coefficients h_hat[{a,b,c,d}]."""
    h4 = {}
    N = 1 << m
    for a in range(m):
        for b in range(a+1, m):
            for c in range(b+1, m):
                for d in range(c+1, m):
                    total = 0
                    for sigma_bits in range(N):
                        sigma = tuple(1 if (sigma_bits >> i) & 1 else -1 for i in range(m))
                        total += values[sigma] * sigma[a] * sigma[b] * sigma[c] * sigma[d]
                    h4[(a,b,c,d)] = total / N
    return h4


def compute_H_p7(A, p):
    """Compute H via full cycle enumeration + independence polynomial for p=7."""
    cycles = []
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n = len(verts)
            start = 0
            dp = {(1 << start, start): 1}
            for mask in range(1, 1 << n):
                if not (mask & (1 << start)):
                    continue
                for v in range(n):
                    if not (mask & (1 << v)):
                        continue
                    key = (mask, v)
                    if key not in dp or dp[key] == 0:
                        continue
                    cnt = dp[key]
                    for w in range(n):
                        if mask & (1 << w):
                            continue
                        if A[verts[v]][verts[w]]:
                            nkey = (mask | (1 << w), w)
                            dp[nkey] = dp.get(nkey, 0) + cnt
            full = (1 << n) - 1
            for v in range(n):
                if v == start:
                    continue
                key = (full, v)
                if key in dp and dp[key] > 0:
                    if A[verts[v]][verts[start]]:
                        cycles.append(frozenset(verts))

    n_cyc = len(cycles)
    # Build conflict graph adjacency
    nbr = [0] * n_cyc
    for i in range(n_cyc):
        for j in range(i+1, n_cyc):
            if cycles[i] & cycles[j]:
                nbr[i] |= (1 << j)
                nbr[j] |= (1 << i)

    alpha = [0] * (n_cyc + 1)
    def backtrack(v, mask, size):
        alpha[size] += 1
        for w in range(v + 1, n_cyc):
            if not (mask & (1 << w)):
                new_mask = mask | nbr[w]
                backtrack(w, new_mask, size + 1)
    backtrack(-1, 0, 0)

    H = sum(alpha[j] * (2**j) for j in range(len(alpha)))
    return H, alpha, cycles


def main():
    print("=" * 70)
    print("OVERLAP-GAUSS BRIDGE: Connecting overlap weights to Walsh layers")
    print("=" * 70)

    # ====== PART 1: Test disj3 = -c5_dir/2 + const at all primes ======
    print("\n" + "=" * 60)
    print("PART 1: disj3(3-cycle) vs c5_dir relationship")
    print("=" * 60)

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2
        n_orient = 1 << m
        print(f"\n  p={p}, m={m}, {n_orient} orientations")

        data = []
        for bits, S in all_orientations_circ(p):
            A = build_adj(p, S)
            c3_dir = count_ck(A, p, 3)
            c5_dir = count_ck(A, p, 5)
            c3_sets = count_c3_vertex_sets(A, p)
            n3 = len(c3_sets)
            disj3 = disjoint_pair_counts(c3_sets)

            sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
            data.append({
                'bits': bits, 'sigma': sigma, 'S': S,
                'c3_dir': c3_dir, 'c5_dir': c5_dir,
                'c3_vsets': n3, 'disj3': disj3
            })

        # Check c3 constancy
        c3_vals = sorted(set(d['c3_vsets'] for d in data))
        c5_vals = sorted(set(d['c5_dir'] for d in data))
        disj3_vals = sorted(set(d['disj3'] for d in data))
        print(f"  c3 vertex sets: {c3_vals[:5]}{'...' if len(c3_vals)>5 else ''}")
        print(f"  c5_dir range: [{c5_vals[0]}, {c5_vals[-1]}], {len(c5_vals)} distinct values")
        print(f"  disj3 range: [{disj3_vals[0]}, {disj3_vals[-1]}], {len(disj3_vals)} distinct values")

        # Is disj3 a function of c5_dir?
        disj_by_c5 = defaultdict(set)
        for d in data:
            disj_by_c5[d['c5_dir']].add(d['disj3'])

        is_function = all(len(v) == 1 for v in disj_by_c5.values())
        print(f"  disj3 is function of c5_dir: {is_function}")

        if is_function and len(c5_vals) >= 2:
            x_vals = c5_vals
            y_vals = [list(disj_by_c5[c5])[0] for c5 in c5_vals]
            slope = (y_vals[-1] - y_vals[0]) / (x_vals[-1] - x_vals[0])
            intercept = y_vals[0] - slope * x_vals[0]
            max_err = max(abs(y - (slope * x + intercept)) for x, y in zip(x_vals, y_vals))
            print(f"  LINEAR FIT: disj3 = {slope:.6f} * c5_dir + {intercept:.2f}")
            print(f"  Max error = {max_err:.10f}")

            if abs(slope + 0.5) < 0.001 and max_err < 0.01:
                print(f"  *** CONFIRMED: disj3 = -c5_dir/2 + {intercept:.0f} ***")
            elif max_err < 0.01:
                print(f"  *** EXACT LINEAR but slope = {slope:.6f} ***")
            else:
                print(f"  *** NOT exactly linear ***")

        elif not is_function:
            # Check how badly it fails
            multi = sum(1 for v in disj_by_c5.values() if len(v) > 1)
            print(f"  WARNING: {multi}/{len(disj_by_c5)} c5 values have multiple disj3 values")
            for c5 in sorted(disj_by_c5):
                vals = disj_by_c5[c5]
                if len(vals) > 1:
                    print(f"    c5={c5}: disj3 in {sorted(vals)}")
                    if len(list(sorted(vals))) <= 5:
                        break  # just show a few

            # Try fitting disj3 = a * c5 + b * c3 + const
            # Actually at p=13, c3 can vary (26 for Paley, 91 for Interval)
            # So maybe disj3 = f(c3, c5)?
            if len(c3_vals) > 1:
                print(f"  c3 varies: trying disj3 = a*c5 + b*c3 + const")
                # Simple 2D regression
                from numpy import array as _  # skip if unavailable
                # Do it manually: 3 unknowns a,b,c, minimize sum (disj3 - a*c5 - b*c3 - c)^2
                n = len(data)
                sx5 = sum(d['c5_dir'] for d in data)
                sx3 = sum(d['c3_vsets'] for d in data)
                sy = sum(d['disj3'] for d in data)
                sx5x5 = sum(d['c5_dir']**2 for d in data)
                sx3x3 = sum(d['c3_vsets']**2 for d in data)
                sx5x3 = sum(d['c5_dir']*d['c3_vsets'] for d in data)
                sx5y = sum(d['c5_dir']*d['disj3'] for d in data)
                sx3y = sum(d['c3_vsets']*d['disj3'] for d in data)

                # Solve 3x3 system
                # [sx5x5, sx5x3, sx5] [a]   [sx5y]
                # [sx5x3, sx3x3, sx3] [b] = [sx3y]
                # [sx5,   sx3,   n  ] [c]   [sy  ]
                det = sx5x5*(sx3x3*n - sx3**2) - sx5x3*(sx5x3*n - sx3*sx5) + sx5*(sx5x3*sx3 - sx3x3*sx5)
                if abs(det) > 0.001:
                    a = (sx5y*(sx3x3*n - sx3**2) - sx5x3*(sx3y*n - sx3*sy) + sx5*(sx3y*sx3 - sx3x3*sy)) / det
                    b = (sx5x5*(sx3y*n - sx3*sy) - sx5y*(sx5x3*n - sx3*sx5) + sx5*(sx5x3*sy - sx3y*sx5)) / det
                    c = (sx5x5*(sx3x3*sy - sx3*sx3y) - sx5x3*(sx5x3*sy - sx3y*sx5) + sx5y*(sx5x3*sx3 - sx3x3*sx5)) / det

                    max_err2 = max(abs(d['disj3'] - a*d['c5_dir'] - b*d['c3_vsets'] - c) for d in data)
                    print(f"  FIT: disj3 = {a:.6f}*c5 + {b:.6f}*c3 + {c:.2f}")
                    print(f"  Max residual = {max_err2:.6f}")

    # ====== PART 2: Walsh degree-2 of disj3, c5, c3 at all primes ======
    print("\n" + "=" * 60)
    print("PART 2: Walsh degree-2 of disj3 and c5")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}")

        disj3_vals = {}
        c5_vals_dict = {}
        c3_vals_dict = {}

        for bits, S in all_orientations_circ(p):
            A = build_adj(p, S)
            c3_dir = count_ck(A, p, 3)
            c5_dir = count_ck(A, p, 5)
            c3_sets = count_c3_vertex_sets(A, p)
            disj3 = disjoint_pair_counts(c3_sets)

            sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
            disj3_vals[sigma] = disj3
            c5_vals_dict[sigma] = c5_dir
            c3_vals_dict[sigma] = c3_dir

        h2_disj3 = walsh_degree2(disj3_vals, m)
        h2_c5 = walsh_degree2(c5_vals_dict, m)
        h2_c3 = walsh_degree2(c3_vals_dict, m)

        print(f"\n  Walsh degree-2 coefficients:")
        print(f"  {'pair':>8} {'h2[disj3]':>12} {'h2[c5]':>12} {'h2[c3]':>10} {'ratio d/c5':>12}")

        for (a, b) in sorted(h2_disj3.keys()):
            hd = h2_disj3[(a, b)]
            hc5 = h2_c5[(a, b)]
            hc3 = h2_c3[(a, b)]
            ratio = hd / hc5 if abs(hc5) > 0.001 else float('nan')
            if abs(hd) > 0.001 or abs(hc5) > 0.001:
                print(f"  ({a},{b})  {hd:>12.4f} {hc5:>12.4f} {hc3:>10.4f} {ratio:>12.6f}")

    # ====== PART 3: Full alpha decomposition + H at p=7 ======
    print("\n" + "=" * 60)
    print("PART 3: Full H decomposition at p=7 (all 8 orientations)")
    print("=" * 60)

    p = 7
    m = 3
    H_vals = {}
    alpha1_vals = {}
    alpha2_vals = {}
    disj3_only_vals = {}

    for bits, S in all_orientations_circ(p):
        A = build_adj(p, S)
        H, alpha, cycles = compute_H_p7(A, p)
        sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
        H_vals[sigma] = H
        alpha1_vals[sigma] = alpha[1]
        alpha2_vals[sigma] = alpha[2]

        c3_sets = count_c3_vertex_sets(A, p)
        disj3_only_vals[sigma] = disjoint_pair_counts(c3_sets)

        by_len = defaultdict(int)
        for c in cycles:
            by_len[len(c)] += 1

        print(f"  bits={bits:03b} S={S}: c3={by_len[3]} c5={by_len[5]} c7={by_len[7]} "
              f"alpha_1={alpha[1]} alpha_2={alpha[2]} disj3={disj3_only_vals[sigma]} H={H}")

    # Walsh of everything
    h2_H = walsh_degree2(H_vals, m)
    h2_a1 = walsh_degree2(alpha1_vals, m)
    h2_a2 = walsh_degree2(alpha2_vals, m)
    h2_d3 = walsh_degree2(disj3_only_vals, m)

    print(f"\n  Walsh degree-2:")
    print(f"  {'pair':>8} {'h2[H]':>10} {'h2[a1]':>10} {'h2[a2]':>10} {'h2[d3]':>10} {'2*h2[a1]+4*h2[a2]':>18}")
    for (a, b) in sorted(h2_H.keys()):
        hH = h2_H[(a,b)]
        ha1 = h2_a1[(a,b)]
        ha2 = h2_a2[(a,b)]
        hd3 = h2_d3[(a,b)]
        combo = 2*ha1 + 4*ha2
        if abs(hH) > 0.001:
            print(f"  ({a},{b})  {hH:>10.4f} {ha1:>10.4f} {ha2:>10.4f} {hd3:>10.4f} {combo:>18.4f}")

    # ====== PART 4: Degree-4 Walsh of H at p=7 and 4-fold Gauss sums ======
    print("\n" + "=" * 60)
    print("PART 4: Degree-4 Walsh of H at p=7")
    print("=" * 60)

    # Only 1 degree-4 term at m=3: (0,1,2)... wait, degree 4 needs 4 indices but m=3
    # So no degree-4 Walsh at p=7! The max is degree 3.
    # Let me compute degree-3 instead
    h3 = {}
    N = 1 << m
    for a in range(m):
        for b in range(a+1, m):
            for c in range(b+1, m):
                total = 0
                for sigma_bits in range(N):
                    sigma = tuple(1 if (sigma_bits >> i) & 1 else -1 for i in range(m))
                    total += H_vals[sigma] * sigma[a] * sigma[b] * sigma[c]
                h3[(a,b,c)] = total / N

    print(f"  Degree-3 Walsh of H at p=7:")
    for (a,b,c), val in sorted(h3.items()):
        if abs(val) > 0.001:
            print(f"    ({a},{b},{c}): h3 = {val:.4f}")

    # Full Walsh spectrum at p=7
    print(f"\n  Full Walsh spectrum of H at p=7:")
    h0 = sum(H_vals[sigma] for sigma in H_vals) / N
    print(f"    degree 0: {h0:.4f}")

    h1 = {}
    for a in range(m):
        total = 0
        for sigma_bits in range(N):
            sigma = tuple(1 if (sigma_bits >> i) & 1 else -1 for i in range(m))
            total += H_vals[sigma] * sigma[a]
        h1[a] = total / N
        if abs(h1[a]) > 0.001:
            print(f"    degree 1, chord {a+1}: {h1[a]:.4f}")

    for (a, b) in sorted(h2_H.keys()):
        if abs(h2_H[(a,b)]) > 0.001:
            print(f"    degree 2, ({a+1},{b+1}): {h2_H[(a,b)]:.4f}")

    for (a,b,c), val in sorted(h3.items()):
        if abs(val) > 0.001:
            print(f"    degree 3, ({a+1},{b+1},{c+1}): {val:.4f}")

    # Verify Parseval: sum of h_S^2 * 2^m should equal sum of H^2 / 2^m
    all_h = [h0] + list(h1.values()) + list(h2_H.values()) + list(h3.values())
    parseval_sum = sum(h**2 for h in all_h)
    H_sq_sum = sum(H**2 for H in H_vals.values()) / N
    print(f"\n  Parseval check: sum(h_S^2) = {parseval_sum:.4f}, sum(H^2)/N = {H_sq_sum:.4f}")

    # ====== PART 5: Degree-4 Walsh at p=11 ======
    print("\n" + "=" * 60)
    print("PART 5: Degree-4 Walsh of disj3 and c5 at p=11")
    print("=" * 60)

    p = 11
    m = 5
    print(f"  p={p}, m={m}")

    disj3_vals_11 = {}
    c5_vals_11 = {}
    c3_vals_11 = {}

    for bits, S in all_orientations_circ(p):
        A = build_adj(p, S)
        c5_dir = count_ck(A, p, 5)
        c3_sets = count_c3_vertex_sets(A, p)
        disj3 = disjoint_pair_counts(c3_sets)

        sigma = tuple(1 if (bits >> i) & 1 else -1 for i in range(m))
        disj3_vals_11[sigma] = disj3
        c5_vals_11[sigma] = c5_dir

    h4_disj3 = walsh_degree4(disj3_vals_11, m)
    h4_c5 = walsh_degree4(c5_vals_11, m)

    print(f"\n  Degree-4 Walsh coefficients:")
    print(f"  {'(a,b,c,d)':>12} {'h4[disj3]':>12} {'h4[c5]':>12} {'ratio':>10}")

    for key in sorted(h4_disj3.keys()):
        hd = h4_disj3[key]
        hc = h4_c5[key]
        ratio = hd / hc if abs(hc) > 0.001 else float('nan')
        if abs(hd) > 0.01 or abs(hc) > 0.01:
            print(f"  {key}  {hd:>12.4f} {hc:>12.4f} {ratio:>10.6f}")

    # 4-fold Gauss sums at p=11
    print(f"\n  4-fold Gauss sums at p=11:")
    def gauss4(ga, gb, gc, gd, p):
        total = 0
        for t in range(1, p):
            val = 1
            for g in [ga, gb, gc, gd]:
                val *= math.sin(2 * math.pi * g * t / p)
            total += val
        return total

    for key in sorted(h4_c5.keys()):
        hc = h4_c5[key]
        if abs(hc) > 0.01:
            a, b, c, d = key
            ga, gb, gc, gd = a+1, b+1, c+1, d+1
            g4 = gauss4(ga, gb, gc, gd, p)
            ratio = hc / g4 if abs(g4) > 0.001 else float('nan')

            # Zero-sum index
            zero_sums = 0
            for s1 in [1, -1]:
                for s2 in [1, -1]:
                    for s3 in [1, -1]:
                        if (ga + s1*gb + s2*gc + s3*gd) % p == 0:
                            zero_sums += 1
                        if (-ga + s1*gb + s2*gc + s3*gd) % p == 0:
                            zero_sums += 1

            print(f"  {key} g=({ga},{gb},{gc},{gd}): h4[c5]={hc:.4f}, G4={g4:.4f}, "
                  f"h4/G4={ratio:.4f}, zero_sums={zero_sums}")

    # ====== PART 6: Disjointness formula derivation ======
    print("\n" + "=" * 60)
    print("PART 6: Deriving disjointness formula")
    print("=" * 60)

    # For 3-cycles on Z_p: two 3-cycles {a,b,c} and {d,e,f} are disjoint
    # iff they use 6 distinct vertices.
    # c3 vertex sets = p * zero_sum_triples / (3 * 2) for unique vertex sets
    # No: c3_vsets = p * (# zero-sum unordered triples {g1,g2,g3} in S with g1+g2+g3=0) / 3
    # Actually for circulant: each 3-cycle vertex set {v, v+g1, v+g1+g2} is determined by
    # its gap triple (g1, g2, g3=-g1-g2) with all in S, and the starting vertex v.
    # There are p choices for v, giving p cycles per gap triple.
    # But {v, v+g1, v+g1+g2} = {v', v'+g1', v'+g1'+g2'} iff same vertex set.

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            S_set = set(S)

            # Count zero-sum triples (ordered)
            ordered_zst = 0
            for a in S:
                for b in S:
                    c = (p - a - b) % p
                    if c in S_set:
                        ordered_zst += 1

            # Each unordered triple {g1,g2,g3} with g1+g2+g3=0 contributes:
            # - 6 ordered if all distinct
            # - 3 ordered if exactly two equal (e.g., g1=g2, g3=-2g1)
            # - 1 ordered if all equal (g1=g2=g3, so 3g=0, g=p/3, only if 3|p)

            # Count the number of distinct gap triples
            gap_triples = set()
            for a in S:
                for b in S:
                    c = (p - a - b) % p
                    if c in S_set:
                        gap_triples.add(tuple(sorted([a, b, c])))

            # Each gap triple {g1,g2,g3} gives p vertex sets via translation
            # But some vertex sets may be the same if the gap triple has symmetry
            # For a 3-cycle (v, v+g1, v+g1+g2), the 3 rotations give gap triples:
            # (g1, g2, g3), (g2, g3, g1), (g3, g1, g2)
            # These are the SAME unordered triple!
            # So each unordered gap triple gives p vertex sets / 1 (no further identification)
            # Wait: different starting vertices give different vertex sets (translation by 1)
            # and the 3 rotations give the same vertex set (just different starting point)

            # So c3_vsets = p * |gap_triples| / 3... no
            # Each ordered gap triple (g1,g2,g3) gives p vertex sets.
            # ordered_zst triples, each starting from p vertices, but:
            # - Each vertex set is counted 3 times (3 starting vertices for same set)
            # - Each unordered triple has 6, 3, or 1 ordered forms
            # So c3_vsets = p * ordered_zst / (3 * p)... no, that's ordered_zst/3

            # Actually: each directed 3-cycle on vertices {a,b,c} corresponds to
            # ONE ordered gap triple (g1,g2,g3). There are p translations.
            # c3_directed = p * ordered_zst / 3 (divide by 3 for vertex rotations)
            # c3_vsets = c3_directed / 2 if each set has 2 directed cycles? No, for tournaments
            # each 3-vertex set has exactly 0 or 2 directed cycles (the cycle and its reverse)
            # For tournaments: if {a,b,c} forms a 3-cycle, there are 2 directed Ham cycles.
            # But in our counting, we get the directed count.

            # Let me just verify
            A = build_adj(p, S)
            c3_sets = count_c3_vertex_sets(A, p)
            n3 = len(c3_sets)
            c3_dir = count_ck(A, p, 3)
            disj3 = disjoint_pair_counts(c3_sets)

            print(f"\n  p={p}, {name}:")
            print(f"    |gap_triples| = {len(gap_triples)}")
            print(f"    ordered_zst = {ordered_zst}")
            print(f"    c3_vsets (actual) = {n3}")
            print(f"    c3_directed (actual) = {c3_dir}")
            print(f"    Prediction c3_vsets = p * |gap_triples| = {p * len(gap_triples)} "
                  f"(ratio: {n3 / (p * len(gap_triples)) if len(gap_triples) > 0 else 'N/A'})")
            print(f"    disj3 = {disj3}")

            # Disjointness: two 3-cycles {v,v+g1,v+g1+g2} and {w,w+h1,w+h1+h2} are disjoint
            # iff their vertex sets don't overlap.
            # For circulant: by symmetry disj3 = p * (# disjoint pairs with one cycle at v=0) / 2?
            # No, it's more subtle because of the orbit structure.

            # Direct formula: fix cycle C0 = {0, g1, g1+g2}.
            # Count cycles disjoint from C0.
            # By symmetry all c3 vertex sets have same number of disjoint partners.
            # So disj3 = n3 * D / 2 where D = # disjoint partners of any fixed cycle.
            if n3 > 0:
                # Count disjoint partners of first cycle
                c0 = list(c3_sets)[0]
                D = sum(1 for c in c3_sets if not (c & c0))
                print(f"    D (disjoint from fixed cycle) = {D}")
                print(f"    n3 * D / 2 = {n3 * D // 2} (should = {disj3})")

                # Count by gap from c0
                # For each vertex set disjoint from c0, what is min distance?
                c0_verts = sorted(c0)
                forbidden = set(c0_verts)

                # Count gap-triple compositions avoiding c0's vertices
                avail = set(range(p)) - forbidden
                avail_diffs = set()
                for u in avail:
                    for v in avail:
                        if u != v:
                            avail_diffs.add((v - u) % p)

                print(f"    Forbidden vertices: {sorted(forbidden)}")
                print(f"    Available vertices: {sorted(avail)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
