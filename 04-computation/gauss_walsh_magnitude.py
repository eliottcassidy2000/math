#!/usr/bin/env python3
"""
gauss_walsh_magnitude.py -- Connect Gauss sum structure to Walsh magnitudes

KEY INSIGHT: For Paley T_p at p=3 mod 4, the Fourier transform of the
connection set S = QR is:
    S_hat(t) = (-1 + chi(t) * i*sqrt(p)) / 2    for t != 0
    S_hat(0) = m = (p-1)/2

This gives EXACT cycle count formulas:
    c_k = (1/k)(m^k + (p-1) * r^k * cos(k*theta))
where r = sqrt((p+1)/4) and theta = pi - arctan(sqrt(p)).

QUESTION: Can we derive |h_hat[{a,b}]| from this?

The Walsh coefficient h_hat[{a,b}] measures how H(sigma) varies when
chords a and b are simultaneously flipped. For circulant tournaments,
flipping chord i changes the connection set by swapping gap (i+1) <-> p-(i+1).

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def exact_cycle_count(p, k):
    """Exact k-cycle count for Paley T_p using Gauss sum formula."""
    m = (p - 1) // 2
    z = complex(-1, math.sqrt(p)) / 2  # (-1 + i*sqrt(p))/2
    # S_hat(0) = m
    # For t != 0 with chi(t)=+1: S_hat(t) = z
    # For t != 0 with chi(t)=-1: S_hat(t) = conj(z)
    # Number of each: (p-1)/2

    total = complex(m, 0)**k  # t=0 contribution
    total += m * z**k          # (p-1)/2 terms with S_hat = z
    total += m * z.conjugate()**k  # (p-1)/2 terms with S_hat = conj(z)

    return total.real / k


def verify_cycle_counts(p):
    """Verify the exact formula against direct enumeration."""
    print(f"\n  Cycle count verification at p={p}:")
    from itertools import combinations

    m = (p - 1) // 2
    QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

    # Build adjacency
    A = [[0]*p for _ in range(p)]
    for v in range(p):
        for s in QR:
            A[v][(v + s) % p] = 1

    # Direct count for small k
    for k in [3, 5, 7]:
        if k > p:
            break
        # Count via gap sequences summing to 0 mod p
        direct = 0
        def count_seqs(remaining, partial_sum, seq):
            nonlocal direct
            if remaining == 0:
                if partial_sum % p == 0:
                    direct += 1
                return
            for g in QR:
                count_seqs(remaining - 1, partial_sum + g, seq + [g])
        count_seqs(k, 0, [])
        direct_ck = direct // k  # directed cycles (not gap sequences)

        formula_ck = exact_cycle_count(p, k)
        print(f"    c_{k}: direct = {direct_ck}, formula = {formula_ck:.1f}, "
              f"match = {'YES' if abs(direct_ck - formula_ck) < 0.5 else 'NO'}")


def orientation_walsh_analysis(p):
    """Compute degree-2 Walsh of H over orientations using direct computation.
    Then compare with Gauss-sum-based predictions."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]

    # Compute H for all orientations
    H_vals = {}
    for bits in range(1 << m):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        adj_list = [[] for _ in range(p)]
        for v in range(p):
            for s in S:
                adj_list[v].append((v + s) % p)
        H_vals[bits] = count_ham_paths_fast(adj_list, p)

    # Compute degree-2 Walsh coefficients
    deg2 = {}
    for a in range(m):
        for b in range(a + 1, m):
            total = 0
            for bits in range(1 << m):
                sa = 1 if bits & (1 << a) else -1
                sb = 1 if bits & (1 << b) else -1
                total += H_vals[bits] * sa * sb
            deg2[(a, b)] = total / (1 << m)

    # Compute degree-4 Walsh coefficients
    deg4 = {}
    for a in range(m):
        for b in range(a+1, m):
            for c in range(b+1, m):
                for d in range(c+1, m):
                    total = 0
                    for bits in range(1 << m):
                        prod = 1
                        for idx in [a, b, c, d]:
                            prod *= (1 if bits & (1 << idx) else -1)
                        total += H_vals[bits] * prod
                    val = total / (1 << m)
                    if abs(val) > 0.001:
                        deg4[(a, b, c, d)] = val

    return H_vals, deg2, deg4


def count_ham_paths_fast(adj_list, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp[mask][v]
            if cnt == 0:
                continue
            for w in adj_list[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += cnt
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def gauss_sum_prediction(p, q):
    """Try to predict |h_hat(q)| from p and q using Gauss sum structure.

    The D^{2n} Walsh contribution at resonance level q has onset n=(q+1)/2.
    The Gauss sum z = (-1+i*sqrt(p))/2 has |z|^2 = (p+1)/4.

    D^{2n} for chord pair (a,b) involves:
    (1/p) sum_t S_hat(t)^{2n} * sin(2*pi*t*gap_a/p) * sin(2*pi*t*gap_b/p)

    At resonance qa = +/- b mod p, there's a special t value where the sine
    product resonates.
    """
    m = (p - 1) // 2
    z = complex(-1, math.sqrt(p)) / 2
    r_sq = abs(z)**2  # = (p+1)/4

    onset = (q + 1) // 2

    # The key quantity seems to be r^{2*onset} = ((p+1)/4)^onset
    r_onset = r_sq ** onset

    return r_onset


def main():
    print("=" * 70)
    print("GAUSS SUM <-> WALSH MAGNITUDE CONNECTION")
    print("=" * 70)

    # Part 1: Verify exact cycle count formula
    print("\n" + "=" * 60)
    print("PART 1: EXACT CYCLE COUNT FORMULA VERIFICATION")
    print("=" * 60)
    for p in [7, 11, 19]:
        verify_cycle_counts(p)

    # Part 2: Gauss sum Fourier analysis
    print("\n" + "=" * 60)
    print("PART 2: GAUSS SUM FOURIER STRUCTURE")
    print("=" * 60)

    for p in [7, 11, 19]:
        m = (p - 1) // 2
        z = complex(-1, math.sqrt(p)) / 2
        zbar = z.conjugate()

        print(f"\n  p={p}:")
        print(f"    z = (-1 + i*sqrt({p}))/2 = {z}")
        print(f"    |z|^2 = {abs(z)**2:.6f} = (p+1)/4 = {(p+1)/4}")
        print(f"    arg(z) = {cmath.phase(z):.6f} rad = {cmath.phase(z)*180/math.pi:.2f} deg")

        theta = cmath.phase(z)
        r = abs(z)
        print(f"    r = {r:.6f}, theta = {theta:.6f}")

        # z^k for small k
        for k in [2, 3, 4, 5, 6]:
            zk = z**k
            print(f"    z^{k} = {zk.real:>12.4f} + {zk.imag:>12.4f}i, "
                  f"|z^{k}| = {abs(zk):>10.4f}, Re(z^{k}) = {zk.real:>12.4f}")

    # Part 3: Walsh magnitude analysis
    print("\n" + "=" * 60)
    print("PART 3: WALSH MAGNITUDE FORMULA SEARCH")
    print("=" * 60)

    for p in [7, 11]:
        m = (p - 1) // 2
        z = complex(-1, math.sqrt(p)) / 2

        print(f"\n{'='*50}")
        print(f"p = {p}, m = {m}")
        print(f"{'='*50}")

        H_vals, deg2, deg4 = orientation_walsh_analysis(p)

        # Group deg-2 by q
        q_groups = defaultdict(list)
        for (a, b), val in deg2.items():
            gap_a, gap_b = a + 1, b + 1
            res = classify_resonance(gap_a, gap_b, p)
            min_q = min(qq for qq, t in res) if res else 'inf'
            q_groups[min_q].append(((a, b), val))

        for q in sorted(q_groups.keys()):
            items = q_groups[q]
            mag = abs(items[0][1])
            onset = (q + 1) // 2
            chi_q = legendre(q, p)

            print(f"\n  q={q}, onset=D^{2*onset}, chi(q)={chi_q:+d}")
            print(f"    |h_hat| = {mag:.4f}")
            print(f"    4*|h_hat| = {4*mag:.4f}")

            # Check various formulas
            r_sq = (p + 1) / 4
            print(f"    r^2 = (p+1)/4 = {r_sq:.4f}")
            print(f"    r^{2*onset} = {r_sq**onset:.4f}")
            print(f"    (p-1)/2 * r^{2*onset} = {(p-1)/2 * r_sq**onset:.4f}")

            # The onset D^{2*onset} has a specific structure
            # For resonance q, the t-sum has a dominant term from t where
            # q*gap_a = +/- gap_b mod p

            # Try: |h_hat| = C * r^{2*onset} * sin_factor
            # At p=7, q=3: |h_hat|=3.5, onset=2, r^4 = 4
            #   3.5/4 = 0.875
            # At p=11, q=3: |h_hat|=272.25, onset=2, r^4 = 9
            #   272.25/9 = 30.25
            # At p=11, q=5: |h_hat|=8.25, onset=3, r^6 = 27
            #   8.25/27 = 0.30556

            print(f"    |h_hat| / r^{2*onset} = {mag / r_sq**onset:.6f}")

            # Try p-dependent factors
            print(f"    |h_hat| / (m * r^{2*onset}) = {mag / (m * r_sq**onset):.6f}")
            print(f"    |h_hat| * 4 / r^{2*onset} = {4 * mag / r_sq**onset:.6f}")

            # Key: check if |h_hat| = (p-1)/4 * r^{2*(onset-1)}
            if onset >= 1:
                candidate = (p - 1) / 4 * r_sq**(onset - 1)
                print(f"    (p-1)/4 * r^{2*(onset-1)} = {candidate:.6f}")

            # Check: |h_hat| * 4 / p = ?
            print(f"    4*|h_hat|/p = {4*mag/p:.6f}")
            print(f"    4*|h_hat|/p^2 = {4*mag/p**2:.6f}")

            # Check ratio to Gauss sum powers
            g_sq = p  # |g|^2 = p for the full Gauss sum
            print(f"    |h_hat| / p^{onset-1} = {mag / p**(onset-1):.6f}")
            print(f"    |h_hat| * 4 / (p+1)^{onset} = {4*mag / (p+1)**onset:.6f}")

        # Degree-4 coefficients
        if deg4:
            print(f"\n  Degree-4 Walsh coefficients:")
            q4_groups = defaultdict(list)
            for (a, b, c, d), val in deg4.items():
                gaps = (a+1, b+1, c+1, d+1)
                chi_prod = legendre(math.prod(gaps), p)
                # What's the "resonance level" for a 4-element set?
                # For now, use min pairwise q
                pairwise_qs = []
                for i in range(4):
                    for j in range(i+1, 4):
                        res = classify_resonance(gaps[i], gaps[j], p)
                        if res:
                            pairwise_qs.append(min(qq for qq, t in res))
                min_pair_q = min(pairwise_qs) if pairwise_qs else 'inf'
                q4_groups[min_pair_q].append(((a,b,c,d), val, chi_prod, gaps))

            for q in sorted(q4_groups.keys()):
                items = q4_groups[q]
                mags = sorted(set(abs(v) for _, v, _, _ in items))
                print(f"    min_pair_q={q}: |h_hat_4| = {mags}, count={len(items)}")
                for S, v, chi_prod, gaps in items:
                    sign_ok = (1 if v > 0 else -1) == chi_prod
                    print(f"      {S} gaps={gaps}: h_hat={v:>10.4f}, "
                          f"chi(prod)={chi_prod:+d}, sign={'OK' if sign_ok else 'FAIL'}")

    # Part 4: What determines |h_hat(q)|?
    print("\n\n" + "=" * 60)
    print("PART 4: MAGNITUDE FORMULA CANDIDATES")
    print("=" * 60)

    # Collect all known |h_hat(q)| values
    data = []
    for p in [7, 11]:
        m = (p - 1) // 2
        H_vals, deg2, deg4 = orientation_walsh_analysis(p)

        q_groups = defaultdict(list)
        for (a, b), val in deg2.items():
            gap_a, gap_b = a + 1, b + 1
            res = classify_resonance(gap_a, gap_b, p)
            min_q = min(qq for qq, t in res) if res else 'inf'
            q_groups[min_q].append(abs(val))

        for q in sorted(q_groups.keys()):
            mag = q_groups[q][0]
            data.append((p, q, mag))

    print("\n  Known data points:")
    for p, q, mag in data:
        onset = (q + 1) // 2
        r_sq = (p + 1) / 4
        print(f"    p={p:>2}, q={q}, onset={onset}, |h_hat|={mag:>12.4f}, "
              f"r^{2*onset}={r_sq**onset:>10.4f}, "
              f"|h_hat|/r^{2*onset}={mag/r_sq**onset:.6f}")

    # Try to find pattern in |h_hat|/r^{2*onset}
    print("\n  Looking for pattern in |h_hat| / r^{2*onset}:")
    for p, q, mag in data:
        onset = (q + 1) // 2
        r_sq = (p + 1) / 4
        ratio = mag / r_sq**onset
        print(f"    p={p}, q={q}: ratio = {ratio:.6f}")
        # Factor out (p-1)/4
        factor = ratio / ((p - 1) / 4)
        print(f"      ratio / ((p-1)/4) = {factor:.6f}")
        # Factor out 1/(2^m)
        factor2 = ratio * (1 << ((p-1)//2))
        print(f"      ratio * 2^m = {factor2:.6f}")

    # Part 5: Gauss sum angle decomposition
    print("\n\n" + "=" * 60)
    print("PART 5: ANGLE DECOMPOSITION")
    print("=" * 60)

    for p in [7, 11]:
        m = (p - 1) // 2
        z = complex(-1, math.sqrt(p)) / 2
        theta = cmath.phase(z)  # angle of z

        print(f"\n  p={p}, theta = {theta:.6f} rad")

        # The angle theta determines how cycle counts oscillate
        # For the Walsh coefficient, flipping chords a and b changes
        # S_hat(t) at frequency t by a factor related to sin(2*pi*t*(a+1)/p)

        # Key: when sigma_a flips, gap (a+1) swaps to p-(a+1).
        # S_hat(t; sigma) changes by omega^{-t*(a+1)} - omega^{t*(a+1)}
        # = -2i sin(2*pi*t*(a+1)/p)
        # Times sigma_a (the sign)

        # For deg-2 Walsh, we need the correlation of H with sigma_a*sigma_b
        # This involves the product sin(2*pi*t*(a+1)/p) * sin(2*pi*t*(b+1)/p)

        # At resonance q: q*(a+1) = +/-(b+1) mod p
        # So the product becomes sin^2(2*pi*t*(a+1)/p) at the resonant t

        # Let's compute the "resonant sine product" for each q
        for q in [3, 5, 7, 9]:
            if q >= p:
                continue
            # Find representative (a, b) with qa = +/- b mod p
            found = False
            for a in range(1, m+1):
                for b in range(a+1, m+1):
                    res = classify_resonance(a, b, p)
                    if res and min(qq for qq, t in res) == q:
                        # Compute sine product sum
                        omega = cmath.exp(2j * cmath.pi / p)
                        sine_sum = 0
                        for t in range(1, p):
                            sa = math.sin(2 * math.pi * t * a / p)
                            sb = math.sin(2 * math.pi * t * b / p)
                            sine_sum += sa * sb
                        print(f"    q={q}, pair=({a},{b}): sine_prod_sum = {sine_sum:.4f}")
                        print(f"      sine_sum / (p/2) = {sine_sum / (p/2):.6f}")
                        print(f"      sine_sum / (p-1) = {sine_sum / (p-1):.6f}")

                        # Also compute the character-weighted sine product
                        chi_sine_sum = 0
                        for t in range(1, p):
                            chi_t = legendre(t, p)
                            sa = math.sin(2 * math.pi * t * a / p)
                            sb = math.sin(2 * math.pi * t * b / p)
                            chi_sine_sum += chi_t * sa * sb
                        print(f"      chi-weighted sine sum = {chi_sine_sum:.4f}")
                        print(f"      chi_sine / sqrt(p) = {chi_sine_sum / math.sqrt(p):.6f}")

                        found = True
                        break
                if found:
                    break

    # Part 6: Direct connection
    print("\n\n" + "=" * 60)
    print("PART 6: DIRECT |h_hat| FORMULA")
    print("=" * 60)

    # The key insight: H(sigma) involves ALL cycle lengths, not just c_3.
    # The Walsh coefficient h_hat[{a,b}] involves how flipping sigma_a and sigma_b
    # SIMULTANEOUSLY changes H.
    #
    # For a circulant tournament, flipping sigma_a changes gap (a+1) to p-(a+1).
    # This changes c_k by an amount that depends on the Fourier structure.
    #
    # But H = I(Omega, 2) is a NONLINEAR function of cycle counts.
    #
    # The simplest decomposition: H(sigma) = sum_{sigma'>=sigma} mu(sigma',sigma) * H(sigma')
    # which is just the Walsh-Hadamard transform itself.
    #
    # Let's try to express h_hat in terms of the Gauss sum directly.
    # h_hat[{a,b}] = (1/2^m) sum_bits H(bits) * sa * sb
    #              = (1/2^m) sum_bits I(Omega(bits), 2) * sa * sb
    #
    # Since I(Omega, 2) = sum over independent sets of 2^|S|,
    # and independent sets correspond to vertex-disjoint odd cycle collections,
    # h_hat involves a SUM over all disjoint-cycle-collections weighted by 2^j.
    #
    # This is hard to compute directly. But empirically:

    print("\n  EMPIRICAL MAGNITUDE TABLE:")
    print(f"  {'p':>3} {'q':>3} {'onset':>5} {'|h_hat|':>12} {'4|h|/p':>10} "
          f"{'4|h|/p^2':>10} {'r^{2*on}':>10} {'|h|/r^{2on}':>12}")

    for p, q, mag in data:
        onset = (q + 1) // 2
        r_sq = (p + 1) / 4
        print(f"  {p:>3} {q:>3} {onset:>5} {mag:>12.4f} {4*mag/p:>10.4f} "
              f"{4*mag/p**2:>10.6f} {r_sq**onset:>10.4f} {mag/r_sq**onset:>12.6f}")

    # KEY OBSERVATION:
    # p=7, q=3: 4|h|/p = 2.0000, |h|/r^4 = 0.875
    # p=11, q=3: 4|h|/p^2 = 9.0000, |h|/r^4 = 30.25
    # p=11, q=5: 4|h|/p = 3.0000, |h|/r^6 = 0.305556

    # Hmm... 4|h|/p is integer-like for q at the "boundary"
    # and 4|h|/p^2 is integer-like for the lowest q

    # Let me check if |h_hat| relates to the DERIVATIVE of cycle counts wrt orientation

    # For the Fourier expansion: when we flip sigma_a (gap_a <-> p-gap_a),
    # S_hat(t) changes by delta(t) = 2 * sigma_a * i * sin(2*pi*t*gap_a/p)
    # (from omega^{t*gap} - omega^{-t*gap} = 2i*sin)

    # The second-order Walsh (flipping both a and b) involves the "interaction":
    # delta_a(t) * delta_b(t) = -4 * sin(2*pi*t*a/p) * sin(2*pi*t*b/p)

    # For the k-th cycle count:
    # d^2 c_k / (d_sigma_a d_sigma_b) involves sums of
    # sin(ta) * sin(tb) * S_hat(t)^{k-2}

    # The Walsh coefficient of c_k at degree 2 is:
    # h_hat_ck[{a,b}] = (1/k) sum_{t=1}^{p-1} (-4) sin(ta) sin(tb) S_hat(t)^{k-2} * C(k,2)
    # No, this isn't right either since the Walsh is over discrete orientations.

    # Let me just compute the "c_k Walsh coefficients" directly for verification.
    print("\n\n  CYCLE-COUNT WALSH DECOMPOSITION:")
    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")
        pairs = [(s, p - s) for s in range(1, m + 1)]

        # For each orientation, compute cycle counts
        ck_vals = defaultdict(dict)  # ck_vals[k][bits] = c_k value
        QR_set = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        for bits in range(1 << m):
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
            S_set = set(S)

            # Compute c_3 by gap triples summing to 0 mod p
            c3 = 0
            for a in S:
                for b in S:
                    c_val = (p - a - b) % p
                    if c_val in S_set:
                        c3 += 1
            c3 = c3 * p // 3  # p translations, /3 for rotations
            # Actually c3 is the number of DIRECTED 3-cycles
            # A directed 3-cycle on vertices {v, v+g1, v+g1+g2} with g1,g2,g3=-(g1+g2) in S
            # For each ordered triple (g1,g2,g3) with sum=0 and all in S, there are p cycles
            # But each cycle has 3 rotations, so c3 = p * zero_sum_triples / 3

            ck_vals[3][bits] = c3

        # Walsh of c_3
        print(f"    Degree-2 Walsh of c_3:")
        for a in range(m):
            for b in range(a + 1, m):
                total = 0
                for bits in range(1 << m):
                    sa = 1 if bits & (1 << a) else -1
                    sb = 1 if bits & (1 << b) else -1
                    total += ck_vals[3][bits] * sa * sb
                h_c3 = total / (1 << m)
                if abs(h_c3) > 0.001:
                    gap_a, gap_b = a + 1, b + 1
                    res = classify_resonance(gap_a, gap_b, p)
                    min_q = min(qq for qq, t in res) if res else 'inf'
                    chi_ab = legendre(gap_a * gap_b, p)
                    print(f"      ({a},{b}): gaps=({gap_a},{gap_b}), q={min_q}, "
                          f"h_hat_c3={h_c3:>8.2f}, chi(ab)={chi_ab:+d}")


if __name__ == '__main__':
    main()
