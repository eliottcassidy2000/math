#!/usr/bin/env python3
"""
ck_walsh_signs.py -- Walsh sign patterns for individual cycle counts c_k

At p=11, the degree-2 Walsh coefficients of c_k show different sign laws:
  c_3=55 CONSTANT (all Walsh = 0)
  c_5: anti-product law (sign * chi(ab) = -1 for some pairs)
  c_7, c_9: mixed?
  c_11 = c_p: product law (sign * chi(ab) = +1)

Investigate WHY different cycle lengths obey different sign rules.

The key: c_k counts directed k-cycles in T_sigma. When we flip orientation j,
we change the connection set S by swapping gap g_j <-> p-g_j. This:
  - Changes some k-cycles (those using arc with gap g_j) -> contribution to Walsh
  - The k-cycle count depends on which arcs use gap j vs p-j

The Fourier approach: c_k(sigma) = (1/p) sum_t S_hat(t)^k where S_hat(t) is DFT of 1_S.
Flipping gap j changes S_hat(t) by adding omega^{-g_j*t} - omega^{g_j*t} = -2i sin(2pi g_j t/p).

So the Walsh coefficient of c_k w.r.t. pair {a,b} involves the interaction of
sin terms in the k-th power of S_hat.

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def held_karp_cycles_by_length(A, p):
    """Count all directed k-cycles for each odd k, using Held-Karp DP."""
    c = defaultdict(int)
    for k in range(3, p + 1, 2):
        for subset in combinations(range(p), k):
            verts = list(subset)
            n_cyc = _count_ham_cycles(A, verts)
            c[k] += n_cyc
    return dict(c)


def _count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        cnt = 0
        if A[a][b] and A[b][c] and A[c][a]:
            cnt += 1
        if A[a][c] and A[c][b] and A[b][a]:
            cnt += 1
        return cnt
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
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


def fourier_ck(S, p, k):
    """Compute c_k via Fourier: c_k = (1/p) sum_t S_hat(t)^k."""
    omega = cmath.exp(2j * cmath.pi / p)
    total = 0
    for t in range(p):
        s_hat = sum(omega**(t*s) for s in S)
        total += s_hat**k
    return total / p


def main():
    p = 11
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    n_orient = 1 << m

    print("=" * 70)
    print(f"c_k WALSH SIGN PATTERNS at p={p}")
    print("=" * 70)

    # Compute c_k for all 2^m orientations
    # This is expensive, so use Fourier for c_k (exact for all k!)
    # Actually, Fourier gives EXACT c_k for circulant tournaments:
    # c_k = (1/p) sum_t |S_hat(t)|^k is WRONG for directed cycles
    # The correct formula for directed k-cycles is:
    # c_k = (1/p) sum_t S_hat(t)^k where S_hat(t) = sum_{s in S} omega^{st}
    # This works because the trace of the adjacency matrix A^k counts
    # closed walks, but for SIMPLE cycles we need inclusion-exclusion.
    # However, (1/p) sum S_hat(t)^k counts the number of directed CLOSED WALKS
    # of length k on the circulant, not simple cycles.
    # For k=3 (triangles): closed 3-walks = 2*c_3 (each triangle counted twice by direction)
    # Actually no, for the circulant adjacency matrix A with A[v][(v+s)%p]=1 for s in S:
    # tr(A^k) = (1/p) sum_t eigenvalue(t)^k * p = sum_t eigenvalue(t)^k
    # eigenvalue(t) = S_hat(t) = sum_{s in S} omega^{st}
    # tr(A^k) = sum_t S_hat(t)^k = # closed walks of length k starting at any vertex
    # For k=3: tr(A^3) = 6 * (# undirected triangles) = 2 * c_3 (directed 3-cycles)
    # Wait: each directed 3-cycle contributes p to the trace (p starting vertices).
    # No: tr(A^3) = sum_v [A^3]_{v,v} = sum of all closed walks v->...->v of length 3.
    # A directed 3-cycle {v1,v2,v3} with v1->v2->v3->v1 contributes 3 to the trace
    # (starting at v1, v2, or v3). The reverse cycle contributes 3 more.
    # So tr(A^3) = 3 * c_3 (directed 3-cycles counted with direction).
    # Wait, c_3 counts DIRECTED 3-cycles. Each unordered triangle gives 2 directed cycles.
    # Each directed cycle gives 3 closed walks (3 starting vertices).
    # So tr(A^3) = 3 * c_3.
    # For circulant: tr(A^k) = sum_t S_hat(t)^k.
    # So c_3 = (1/3) sum_t S_hat(t)^3.

    # For k=p (Hamiltonian cycles):
    # Each directed Hamiltonian cycle gives p closed walks (p starting vertices).
    # But tr(A^p) also includes NON-simple closed walks of length p.
    # So Fourier only gives SIMPLE c_k when k is small enough that all closed
    # k-walks are simple. At k=3 on complete circulant with p>=7, all closed 3-walks
    # are simple (no repeated vertices for 3-walks on > 4 vertices).
    # At k=5 with p>=11, not all 5-walks are simple.

    # So for EXACT c_k, we need Held-Karp DP (for simple cycles) or
    # inclusion-exclusion on the Fourier formula.

    # Let's use Held-Karp for exact values, but only for select orientations
    # to characterize the Walsh structure.

    # Strategy: compute c_k for all 32 orientations, then compute Walsh transforms.
    # At p=11 with 32 orientations and cycle enumeration up to k=11...
    # This is expensive but feasible (we already did it for the full alpha computation).

    # Actually, for Walsh sign analysis, we only need the degree-2 Walsh coefficients.
    # With 4 distinct H-orbits, the degree-2 Walsh only depends on the orbit values.
    # We need c_k per orbit. Let me compute for one representative per orbit.

    # 4 extended orbits at p=11: sizes 10, 10, 2, 10
    # Representatives: bits = 00000, 00010, 01011, 01001 (need to find correct reps)

    # Let me just compute all 32 orientations with c_k and then Walsh transform.
    # For p=11, the cycle enumeration takes ~90 seconds per orientation for full alpha.
    # But just c_k (without alpha) is much faster.

    print("\nComputing c_k for all 32 orientations...")
    all_ck = {}  # bits -> {k: count}

    for bits in range(n_orient):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        A = [[0]*p for _ in range(p)]
        for v in range(p):
            for s in S:
                A[v][(v + s) % p] = 1

        ck = held_karp_cycles_by_length(A, p)
        all_ck[bits] = ck
        if bits % 8 == 0:
            print(f"  bits={bits:05b} ({bits+1}/32) c_k = {ck}")

    # Compute degree-2 Walsh coefficients of each c_k
    print(f"\n{'='*70}")
    print("DEGREE-2 WALSH DECOMPOSITION of c_k")
    print("=" * 70)

    chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]
    sigma = {}  # bits -> sigma vector (+1/-1)^m
    for bits in range(n_orient):
        sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

    for k in range(3, p + 1, 2):
        print(f"\n  c_{k}:")

        # Values across orientations
        vals = [all_ck[bits].get(k, 0) for bits in range(n_orient)]
        val_set = sorted(set(vals))
        print(f"    Distinct values: {val_set}")
        print(f"    Range: {min(vals)} to {max(vals)}")

        # Mean (= degree-0 Walsh coefficient)
        mean = sum(vals) / n_orient
        print(f"    Mean (h_hat_0): {mean:.4f}")

        # Degree-2 Walsh coefficients
        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            chi_ab = legendre(ga * gb, p)

            h_hat = 0
            for bits in range(n_orient):
                sig = sigma[bits]
                h_hat += vals[bits] * sig[a] * sig[b]
            h_hat /= n_orient

            if abs(h_hat) < 0.001:
                continue

            sign = 1 if h_hat > 0 else -1
            sign_chi = sign * chi_ab
            print(f"    ({a},{b}) gaps=({ga},{gb}) q={q}: "
                  f"h_hat={h_hat:>10.4f}, |h_hat|={abs(h_hat):>10.4f}, "
                  f"sign*chi(ab)={sign_chi:+d}")

    # Also compute degree-1 Walsh coefficients
    print(f"\n{'='*70}")
    print("DEGREE-1 WALSH COEFFICIENTS of c_k")
    print("=" * 70)

    for k in range(3, p + 1, 2):
        vals = [all_ck[bits].get(k, 0) for bits in range(n_orient)]
        has_deg1 = False
        for a in range(m):
            h_hat = 0
            for bits in range(n_orient):
                sig = sigma[bits]
                h_hat += vals[bits] * sig[a]
            h_hat /= n_orient
            if abs(h_hat) > 0.001:
                if not has_deg1:
                    print(f"\n  c_{k}:")
                    has_deg1 = True
                print(f"    chord {a} (gap {a+1}): h_hat_1 = {h_hat:.4f}")
        if not has_deg1:
            print(f"\n  c_{k}: all degree-1 coefficients = 0")

    # Degree-4 Walsh coefficients (for comparison with H structure)
    print(f"\n{'='*70}")
    print("DEGREE-4 WALSH COEFFICIENTS of c_k")
    print("=" * 70)

    for k in range(3, p + 1, 2):
        vals = [all_ck[bits].get(k, 0) for bits in range(n_orient)]
        print(f"\n  c_{k}:")
        # Sample a few degree-4 sets
        for abcd in [(0,1,2,3), (0,1,2,4), (0,1,3,4), (0,2,3,4), (1,2,3,4)]:
            a, b, c, d = abcd
            ga, gb, gc, gd = a+1, b+1, c+1, d+1
            chi_prod = legendre(ga * gb * gc * gd, p)

            h_hat = 0
            for bits in range(n_orient):
                sig = sigma[bits]
                h_hat += vals[bits] * sig[a] * sig[b] * sig[c] * sig[d]
            h_hat /= n_orient

            sign = 1 if h_hat > 0 else (-1 if h_hat < 0 else 0)
            if abs(h_hat) > 0.001:
                sign_chi = sign * chi_prod
                print(f"    {abcd} gaps=({ga},{gb},{gc},{gd}): "
                      f"h_hat={h_hat:>12.4f}, chi(prod)={chi_prod:+d}, "
                      f"sign*chi={sign_chi:+d}")
            else:
                print(f"    {abcd}: h_hat = 0")

    # Summary: sign law classification
    print(f"\n{'='*70}")
    print("SIGN LAW CLASSIFICATION")
    print("=" * 70)

    for k in range(3, p + 1, 2):
        vals = [all_ck[bits].get(k, 0) for bits in range(n_orient)]
        product_count = 0
        anti_count = 0
        zero_count = 0
        total = 0
        magnitudes = set()
        by_q = defaultdict(lambda: {'product': 0, 'anti': 0, 'zero': 0, 'mags': set()})

        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            chi_ab = legendre(ga * gb, p)

            h_hat = 0
            for bits in range(n_orient):
                sig = sigma[bits]
                h_hat += vals[bits] * sig[a] * sig[b]
            h_hat /= n_orient

            total += 1
            mag = abs(h_hat)
            magnitudes.add(round(mag, 4))
            by_q[q]['mags'].add(round(mag, 4))

            if abs(h_hat) < 0.001:
                zero_count += 1
                by_q[q]['zero'] += 1
            else:
                sign = 1 if h_hat > 0 else -1
                if sign * chi_ab > 0:
                    product_count += 1
                    by_q[q]['product'] += 1
                else:
                    anti_count += 1
                    by_q[q]['anti'] += 1

        law = "ZERO" if zero_count == total else \
              "PRODUCT" if product_count == total else \
              "ANTI-PRODUCT" if anti_count == total else \
              "MIXED"

        print(f"\n  c_{k}: {law}")
        print(f"    product={product_count}, anti={anti_count}, zero={zero_count}/{total}")
        print(f"    magnitudes: {sorted(magnitudes)}")
        for q in sorted(by_q.keys()):
            d = by_q[q]
            print(f"    q={q}: product={d['product']}, anti={d['anti']}, "
                  f"zero={d['zero']}, mags={sorted(d['mags'])}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
