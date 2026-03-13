#!/usr/bin/env python3
"""
k6_walsh_decomposition.py -- Separate c_6 Walsh into simple-cycle vs N_6 correction

At k=6, tr(A^6) = 6*c_6 + N_6 where N_6 counts non-simple 6-walks.
Non-simple 6-walks decompose into two 3-cycles sharing a vertex.

This script:
1. Computes exact c_6 and N_6/6 for all orientations (p=7,11)
2. Walsh-decomposes c_6 and N_6/6 separately
3. Analytically evaluates the j=4 eigenvalue contribution to separate
   the "eigenvalue-predicted" vs "residual" parts of c_6 Walsh
4. Checks if the q=5 content comes from c_6, N_6, or both

Also: derives j=4 Walsh analytically from D^4 expansion.

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


def resonance_sign(ga, gb, p, q):
    """Sign eps for q resonance."""
    if (q * ga - gb) % p == 0:
        return -1
    elif (q * ga + gb) % p == 0:
        return +1
    elif (ga - q * gb) % p == 0:
        return -1
    elif (ga + q * gb) % p == 0:
        return +1
    return 0


def held_karp(A, verts):
    k = len(verts)
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


def main():
    print("=" * 70)
    print("k=6 WALSH DECOMPOSITION: c_6 vs N_6/6")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        n_orient = 1 << m
        pairs = [(j, p - j) for j in range(1, m + 1)]
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n{'='*60}")
        print(f"p={p}, m={m}")
        print("=" * 60)

        data = []
        for bits in range(n_orient):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            A = [[0]*p for _ in range(p)]
            for v in range(p):
                for s in S:
                    A[v][(v + s) % p] = 1

            # Exact cycle counts
            c3 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 3))
            c6 = sum(held_karp(A, list(sub)) for sub in combinations(range(p), 6))

            # Trace
            tr6 = sum(sum(omega**(s*t) for s in S)**6 for t in range(p)).real
            N6 = round(tr6 - 6 * c6)

            # N6 should = c3 * (c3 - 1) * something ... two 3-cycles sharing vertex
            # Actually: for each vertex v, count pairs of 3-cycles through v.
            # 3-cycles through v: for each pair (a,b) with a->b, b->v, v->a or v->a, a->b, b->v
            # In circulant: c3 is total directed 3-cycles.
            # Each 3-cycle has 3 vertices. Expected N6 = sum_v (c3_v choose 2) * 2! * ... hmm

            data.append({
                'bits': bits, 'S': S, 'c3': c3, 'c6': c6,
                'tr6': tr6, 'N6': N6, 'tr6_over_6': tr6 / 6
            })

        sigma = {}
        for bits in range(n_orient):
            sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

        chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

        # Summary
        print(f"  c3 values: {sorted(set(d['c3'] for d in data))}")
        print(f"  c6 values: {sorted(set(d['c6'] for d in data))}")
        print(f"  N6 values: {sorted(set(d['N6'] for d in data))}")

        # Walsh decomposition of c6, N6/6, tr6/6
        print(f"\n  Degree-2 Walsh coefficients:")
        print(f"  {'pair':>8} {'q':>3} {'h_c6':>12} {'h_N6/6':>12} {'h_tr/6':>12} {'h_c3':>12}")

        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)

            h_c6 = sum(data[bits]['c6'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_N6 = sum(data[bits]['N6'] / 6 * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_tr6 = sum(data[bits]['tr6'] / 6 * sigma[bits][a] * sigma[bits][b]
                        for bits in range(n_orient)) / n_orient
            h_c3 = sum(data[bits]['c3'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient

            if abs(h_c6) > 0.01 or abs(h_N6) > 0.01:
                print(f"  ({a},{b}) q={q:>3}: {h_c6:>12.4f} {h_N6:>12.4f} "
                      f"{h_tr6:>12.4f} {h_c3:>12.4f}")

        # Now: analytical j=2 and j=4 predictions for tr(A^6)/6
        # lambda_t = C + D, C = -1/2, D = sum sigma_i * i*sin(2*pi*g_i*t/p)
        # lambda^6 = sum_{j=0}^6 C(6,j) C^{6-j} D^j
        # Degree-2 Walsh of lambda^6:
        # j=2: C(6,2)*C^4 * [D^2]_{a,b} = 15*(1/16)*2*d_a*d_b
        # j=4: C(6,4)*C^2 * [D^4]_{a,b} = 15*(1/4)*[D^4]
        # j=6: C(6,6)*1 * [D^6]_{a,b} = 1*[D^6]

        print(f"\n  Analytical eigenvalue predictions for tr(A^6)/6:")
        print(f"  {'pair':>8} {'q':>3} {'j=2':>12} {'j=4':>12} {'j=6':>12} {'total':>12} {'actual':>12} {'residual':>12}")

        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)

            # j=2 contribution: sum_t [15/16 * 2 d_a d_b]
            # d_a = i*sin(alpha_a), d_b = i*sin(alpha_b)
            # d_a d_b = -sin(alpha_a)*sin(alpha_b)
            # sum_t sin(alpha_a)*sin(alpha_b) = sum_t sin(2pi ga t/p)*sin(2pi gb t/p)
            #   = eps(a,b)*p/2 if q=3, 0 otherwise (at degree-2 level)
            # Contribution: 15/16 * 2 * (-1) * [sum sin sin] = -15/8 * sum sin sin

            alpha_a = [2 * math.pi * ga * t / p for t in range(1, p)]
            alpha_b = [2 * math.pi * gb * t / p for t in range(1, p)]

            sin_a = [math.sin(x) for x in alpha_a]
            sin_b = [math.sin(x) for x in alpha_b]
            cos_a = [math.cos(x) for x in alpha_a]
            cos_b = [math.cos(x) for x in alpha_b]

            # Gauss sums
            G_ab = sum(sin_a[i] * sin_b[i] for i in range(p-1))
            G_a2b = sum(sin_a[i]**2 * sin_b[i] for i in range(p-1))
            G_ab2 = sum(sin_a[i] * sin_b[i]**2 for i in range(p-1))
            G_a3b = sum(sin_a[i]**3 * sin_b[i] for i in range(p-1))
            G_ab3 = sum(sin_a[i] * sin_b[i]**3 for i in range(p-1))
            G_a2b2 = sum(sin_a[i]**2 * sin_b[i]**2 for i in range(p-1))

            # Sum of sin^2 for all chords
            S2 = -(p - 1) / 4  # This is Sigma d_c^2 = -S_2(t) summed over t? No.
            # S2(t) = sum_c sin^2(2pi g_c t / p) = (p-1)/4 for each t.
            # So sum_c d_c^2 = sum_c (-sin^2) = -(p-1)/4 per t.

            # j=2: C(6,2)*(-1/2)^4 * 2*d_a*d_b summed over t, divided by 6 (for c=tr/6)
            # = 15 * (1/16) * (-2) * sum_t sin(a)*sin(b) / 6
            # = -15/(16*3) * G_ab = -5/16 * G_ab
            j2_per_t_coeff = 15 * (1/16) * (-2)  # = -15/8
            j2 = j2_per_t_coeff * G_ab / 6  # /6 for tr/6 -> c_k

            # j=4: C(6,4)*(-1/2)^2 * [D^4]_{a,b} summed over t / 6
            # [D^4]_{a,b} at each t = 4*d_a*d_b*(3*S_2_t - 2*d_a^2 - 2*d_b^2)
            # where S_2_t = sum_c d_c^2 = -sum_c sin^2(2pi g_c t/p) = -(p-1)/4
            # d_a^2 = -sin^2(a_t), d_b^2 = -sin^2(b_t)
            # So 3*S_2_t - 2*d_a^2 - 2*d_b^2 = -3(p-1)/4 + 2*sin^2(a) + 2*sin^2(b)

            D4_sum = 0
            for i in range(p-1):
                sa, sb = sin_a[i], sin_b[i]
                da_db = -sa * sb
                S2t = -(p-1)/4
                da2 = -sa**2
                db2 = -sb**2
                bracket = 3*S2t - 2*da2 - 2*db2
                D4_sum += 4 * da_db * bracket

            j4_coeff = 15 * (1/4)  # C(6,4) * (-1/2)^2 = 15 * 1/4
            j4 = j4_coeff * D4_sum / 6

            # j=6: [D^6]_{a,b} summed over t / 6
            # Complex to derive analytically. Compute numerically.
            D6_sum = 0
            for i in range(p-1):
                # Compute D^6 at this t, extract coefficient of sigma_a * sigma_b
                # D = sum_c sigma_c * d_c where d_c = i*sin(alpha_c_t)
                # For degree-2 at (a,b), we need:
                # (5,1): 6*d_a^5*d_b + 6*d_a*d_b^5
                # (3,3): 20*d_a^3*d_b^3
                # (3,1,2): 60*d_a^3*d_b*sum_{c!=a,b} d_c^2 + 60*d_a*d_b^3*sum d_c^2
                # (1,1,4): 30*d_a*d_b*sum_{c!=a,b} d_c^4
                # (1,1,2,2): 180*d_a*d_b*sum_{c<d,c!=a,b,d!=a,b} d_c^2*d_d^2

                da = 1j * sin_a[i]
                db = 1j * sin_b[i]
                # all chords
                d_all = []
                for c in range(m):
                    gc = c + 1
                    alpha_c = 2 * math.pi * gc * (i+1) / p  # t = i+1
                    d_all.append(1j * math.sin(alpha_c))

                da_v = d_all[a]
                db_v = d_all[b]
                d_other = [d_all[c] for c in range(m) if c != a and c != b]

                sum_dc2 = sum(dc**2 for dc in d_other)
                sum_dc4 = sum(dc**4 for dc in d_other)
                sum_dc2_sq = sum_dc2**2  # (sum d_c^2)^2
                sum_dc2dc2 = (sum_dc2_sq - sum_dc4) / 2  # sum_{c<d} dc^2*dd^2

                term_51 = 6 * (da_v**5 * db_v + da_v * db_v**5)
                term_33 = 20 * da_v**3 * db_v**3
                term_312 = 60 * (da_v**3 * db_v + da_v * db_v**3) * sum_dc2
                term_114 = 30 * da_v * db_v * sum_dc4
                term_1122 = 180 * da_v * db_v * sum_dc2dc2

                D6_t = term_51 + term_33 + term_312 + term_114 + term_1122
                D6_sum += D6_t.real  # should be real after sum over t

            j6 = D6_sum / 6  # C(6,6)*(-1/2)^0 = 1, /6 for tr/6

            total = j2 + j4 + j6

            # Actual Walsh of tr/6
            h_tr6 = sum(data[bits]['tr6'] / 6 * sigma[bits][a] * sigma[bits][b]
                        for bits in range(n_orient)) / n_orient

            if abs(h_tr6) > 0.01 or abs(total) > 0.01:
                resid = h_tr6 - total
                print(f"  ({a},{b}) q={q:>3}: {j2:>12.4f} {j4:>12.4f} {j6:>12.4f} "
                      f"{total:>12.4f} {h_tr6:>12.4f} {resid:>12.4f}")

        # Check: does q=5 content come from c_6 or N_6?
        print(f"\n  Q=5 content analysis:")
        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            if q != 5:
                continue

            h_c6 = sum(data[bits]['c6'] * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_N6 = sum(data[bits]['N6'] / 6 * sigma[bits][a] * sigma[bits][b]
                       for bits in range(n_orient)) / n_orient
            h_tr6 = sum(data[bits]['tr6'] / 6 * sigma[bits][a] * sigma[bits][b]
                        for bits in range(n_orient)) / n_orient

            print(f"    ({a},{b}) q=5: h_c6={h_c6:.4f}, h_N6/6={h_N6:.4f}, "
                  f"h_tr/6={h_tr6:.4f}")

        # N6 structure: N6 = 6 * sum_v C(c3_v, 2) * 2 where c3_v = directed 3-cycles through v
        # Actually N6 = sum_v (count of ordered pairs of directed 3-cycles through v)
        # For circulant: c3_v = c3/p (constant by homogeneity)... wait, c3 counts DIRECTED 3-cycles total.
        # Each has 3 vertices, so sum_v c3_v = 3*c3. So c3_v = 3*c3/p.
        # N6 = sum_v c3_v*(c3_v - 1) * ... hmm, ordered pairs including direction.

        # Let me compute N6 directly:
        c3_vals = [d['c3'] for d in data]
        N6_vals = [d['N6'] for d in data]

        # Test N6 = alpha * c3^2 + beta * c3 + gamma
        # Since c3 is constant for p=7, this can't be tested there.
        print(f"\n  N6 vs c3 relationship:")
        for d in data[:8]:
            print(f"    bits={d['bits']:>4d}: c3={d['c3']:>6d}, c6={d['c6']:>6d}, "
                  f"N6={d['N6']:>6d}, N6/c3={d['N6']/d['c3'] if d['c3'] > 0 else 'N/A':>10}")

        # For p=7, c3 is constant, so let's check if N6 varies
        if len(set(c3_vals)) == 1:
            print(f"    c3 is CONSTANT = {c3_vals[0]}")
            print(f"    N6 values: {sorted(set(N6_vals))}")
            print(f"    N6 is {'CONSTANT' if len(set(N6_vals)) == 1 else 'VARIABLE'}")

            if len(set(N6_vals)) > 1:
                # Walsh decompose N6
                print(f"\n    N6 Walsh coefficients:")
                for a2, b2 in chord_pairs:
                    ga2, gb2 = a2 + 1, b2 + 1
                    q2 = resonance_level(ga2, gb2, p)
                    h = sum(N6_vals[bits] * sigma[bits][a2] * sigma[bits][b2]
                            for bits in range(n_orient)) / n_orient
                    if abs(h) > 0.01:
                        eps3 = resonance_sign(ga2, gb2, p, 3) if q2 == 3 else 0
                        print(f"      ({a2},{b2}) q={q2}: h_N6 = {h:.4f}, "
                              f"h_N6/(p*eps) = {h/(p*eps3) if eps3 != 0 else 'N/A':.6f}")

    # PART 2: Analytical j=4 contribution derivation
    print(f"\n{'='*60}")
    print(f"ANALYTICAL j=4 CONTRIBUTION")
    print("=" * 60)

    print("""
    lambda_t^k = (C + D)^k where C=-1/2, D = sum sigma_i d_i

    Degree-2 Walsh of lambda^k at pair (a,b):
    = sum_{j even} C(k,j) * C^{k-j} * [D^j]_{a,b}

    j=2: [D^2]_{a,b} = 2*d_a*d_b = -2*sin(alpha_a)*sin(alpha_b)
    j=4: [D^4]_{a,b} = 4*d_a*d_b*(3*S2 - 2*d_a^2 - 2*d_b^2)
         where S2 = sum_c d_c^2 = -(p-1)/4

    Expanding j=4:
    [D^4]_{a,b} = 4*(-sin_a*sin_b)*(-3(p-1)/4 + 2*sin^2_a + 2*sin^2_b)
                = sin_a*sin_b * (3(p-1) - 8*sin^2_a - 8*sin^2_b)

    Summed over t with Gauss machinery:
    - sum sin_a*sin_b = eps*p/2 (q=3)
    - sum sin_a*sin_b*sin^2_a = sum sin^3_a*sin_b
      = (3/4)*sum sin_a*sin_b - (1/4)*sum sin(3ga*t)*sin(gb*t)
      The sin(3ga*t) is a Gauss sum at tripled frequency!
    """)

    # Compute the tripled-frequency Gauss sums
    for p in [7, 11, 13, 17, 23]:
        m = (p - 1) // 2
        print(f"\n  p={p}: Tripled-frequency Gauss sums")
        chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)
            if q != 3:
                continue

            # sum sin(3*ga*t)*sin(gb*t) for t=1..p-1
            G3ab = sum(math.sin(2*math.pi*3*ga*t/p) * math.sin(2*math.pi*gb*t/p)
                       for t in range(1, p))
            G_a3b = sum(math.sin(2*math.pi*ga*t/p)**3 * math.sin(2*math.pi*gb*t/p)
                        for t in range(1, p))
            Gab = sum(math.sin(2*math.pi*ga*t/p) * math.sin(2*math.pi*gb*t/p)
                      for t in range(1, p))

            # Check: G_a3b = (3/4)*Gab - (1/4)*G3ab
            pred = 3/4 * Gab - 1/4 * G3ab
            eps = resonance_sign(ga, gb, p, 3)

            # 3*ga mod p: what chord does it correspond to?
            ga3 = (3 * ga) % p
            if ga3 > p // 2:
                ga3 = p - ga3  # map to chord index
            q_3a_b = resonance_level(ga3, gb, p)

            print(f"    ({a},{b}): 3*ga={3*ga}=>{ga3} mod {p}, q(3ga,gb)={q_3a_b}, "
                  f"Gab={Gab:.4f}, G3ab={G3ab:.4f}, G_a3b={G_a3b:.4f}, "
                  f"check={(abs(G_a3b - pred) < 0.001)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
