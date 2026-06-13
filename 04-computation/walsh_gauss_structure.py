#!/usr/bin/env python3
"""
walsh_gauss_structure.py -- Analyze how Walsh coefficients decompose
via Gauss sum structure at p=11 (all 32 orientations tractable).

KEY INSIGHT: For Paley at p=3 mod 4:
  - sigma_Paley = (chi(1), chi(2), ..., chi(m)) where chord i has gap i+1
  - H(sigma_Paley) = max H (Paley maximizer for p=3 mod 4)
  - The sign law says h_hat[S] * chi(prod gaps in S) = |h_hat[S]|

QUESTION: Can we express |h_hat[{a,b}]| as a function of the resonance
level q alone? At p=11 we know |h_hat| takes only 2 values (8.25 and 272.25)
corresponding to q=5 and q=3. What formula gives these?

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


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


def count_ham_paths_fast(adj_list, n):
    """Bitmask DP Hamiltonian path count."""
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


def compute_all_H(p):
    """Compute H for all 2^m orientations at prime p."""
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    half = 1 << m
    H_vals = {}

    for bits in range(half):
        sigma = []
        for i in range(m):
            sigma.append(1 if bits & (1 << i) else -1)

        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        adj_list = [[] for _ in range(p)]
        for v in range(p):
            for s in S:
                adj_list[v].append((v + s) % p)

        H_vals[bits] = count_ham_paths_fast(adj_list, p)

    return H_vals, m


def walsh_coefficients(H_vals, m):
    """Compute ALL Walsh coefficients from H values."""
    half = 1 << m
    coeffs = {}

    for S_mask in range(half):
        S = [i for i in range(m) if S_mask & (1 << i)]
        h_hat = 0
        for bits in range(half):
            sigma_S = 1
            for i in S:
                sigma_S *= (1 if bits & (1 << i) else -1)
            # H(sigma) = H(-sigma), so each bits contributes 2x
            h_hat += H_vals[bits] * sigma_S * 2

        h_hat /= half  # normalize: 1/2^m * sum (with 2x factor already)
        if abs(h_hat) > 0.001:
            coeffs[tuple(S)] = h_hat

    return coeffs


def main():
    print("=" * 70)
    print("WALSH-GAUSS STRUCTURE ANALYSIS")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n{'='*60}")
        print(f"p = {p}, m = {m}, QR = {QR}")
        print(f"{'='*60}")

        # Compute all H values
        H_vals, m = compute_all_H(p)

        # Show H distribution
        H_dist = defaultdict(int)
        for bits, H in H_vals.items():
            H_dist[H] += 1
        print(f"\n  H distribution (counts symmetric orientations):")
        for H, count in sorted(H_dist.items()):
            print(f"    H = {H}: {count} orientations (x2 with symmetry = {2*count})")

        # Compute Walsh coefficients
        coeffs = walsh_coefficients(H_vals, m)

        # Show degree-2 coefficients
        print(f"\n  Degree-2 Walsh coefficients:")
        deg2 = {S: v for S, v in coeffs.items() if len(S) == 2}

        # Group by resonance level
        q_groups = defaultdict(list)
        for S, v in sorted(deg2.items()):
            a, b = S  # chord indices (0-based)
            gap_a, gap_b = a + 1, b + 1  # gap values
            res = classify_resonance(gap_a, gap_b, p)
            min_q = min(qq for qq, t in res) if res else 'inf'
            chi_ab = legendre(gap_a * gap_b, p)
            chi_prod = legendre(gap_a * gap_b, p)

            # Sign law prediction
            sign_pred = chi_prod
            sign_actual = 1 if v > 0 else -1

            q_groups[min_q].append((S, v, chi_prod, sign_pred == sign_actual))

            print(f"    {S}: gaps=({gap_a},{gap_b}), q={min_q}, "
                  f"h_hat={v:>12.4f}, |h_hat|={abs(v):>12.4f}, "
                  f"chi(prod)={chi_prod:+d}, sign={'OK' if sign_pred==sign_actual else 'FAIL'}")

        # Analyze by q-level
        print(f"\n  |h_hat| by resonance level:")
        for q in sorted(q_groups.keys()):
            items = q_groups[q]
            magnitudes = set(abs(v) for S, v, c, m in items)
            chi_q = legendre(q, p)
            print(f"    q={q} (chi(q)={chi_q:+d}): |h_hat| = {sorted(magnitudes)}")

        # Show degree-4 coefficients
        if m >= 4:
            print(f"\n  Degree-4 Walsh coefficients:")
            deg4 = {S: v for S, v in coeffs.items() if len(S) == 4}
            for S, v in sorted(deg4.items()):
                gaps = tuple(i+1 for i in S)
                chi_prod = legendre(math.prod(gaps), p)
                sign_actual = 1 if v > 0 else -1
                print(f"    {S}: gaps={gaps}, prod={math.prod(gaps)}, "
                      f"h_hat={v:>12.4f}, chi(prod)={chi_prod:+d}, "
                      f"sign={'OK' if chi_prod==sign_actual else 'FAIL'}")

        # PART 2: Ratio analysis
        print(f"\n  --- MAGNITUDE RATIOS ---")
        if len(set(abs(v) for v in deg2.values())) > 1:
            vals = sorted(set(abs(v) for v in deg2.values()))
            print(f"    Distinct |h_hat_2| values: {vals}")
            if len(vals) == 2:
                ratio = vals[1] / vals[0]
                print(f"    Ratio large/small: {ratio:.6f}")
                print(f"    = {ratio:.1f}")
                # Check if ratio = (p+1)/4 * something
                print(f"    Ratio / p = {ratio/p:.6f}")
                print(f"    Ratio / (p-1) = {ratio/(p-1):.6f}")
                print(f"    Ratio / (p+1) = {ratio/(p+1):.6f}")

        # PART 3: Express h_hat in terms of D^{2n} series
        print(f"\n  --- D^{{2n}} DECOMPOSITION ---")
        print(f"  h_hat[{{a,b}}] = sum_n D^{{2n}}_ab / 2^m")
        print(f"  For q-resonant pair: D^{{2n}}_ab nonzero only at n >= (q+1)/2")
        print(f"  Sign of D^{{2n}} at onset: chi(q) * chi(ab)")
        print()

        for q in sorted(q_groups.keys()):
            items = q_groups[q]
            chi_q = legendre(q, p)
            onset = (q + 1) // 2
            rep = items[0]  # representative
            S, v, chi_ab, _ = rep
            gap_a, gap_b = S[0]+1, S[1]+1

            print(f"    q={q}: onset at D^{2*onset}, sign(D^{2*onset}) = "
                  f"chi(q)*chi(ab) = ({chi_q:+d})*({chi_ab:+d}) = {chi_q*chi_ab:+d}")
            print(f"    h_hat = {v:>12.4f}, sign(h_hat) = {1 if v>0 else -1:+d}, "
                  f"chi(ab) = {chi_ab:+d}")

            if chi_q * chi_ab == chi_ab:
                print(f"    => D^{2*onset} has CORRECT sign for product law")
            else:
                print(f"    => D^{2*onset} has WRONG sign for product law!")
                print(f"       But h_hat still has correct sign => "
                      f"higher D^{{2n}} or OCF nonlinearity corrects it")
            print()

    # PART 4: Explicit formula for h_hat at p=7
    print("\n\n" + "=" * 60)
    print("PART 4: EXPLICIT h_hat FORMULA AT p=7")
    print("=" * 60)

    p = 7
    m = 3
    H_vals, _ = compute_all_H(p)

    print(f"\n  All H values (sigma_0=+1 half):")
    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]
        pairs = [(s, p-s) for s in range(1, m+1)]
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        print(f"    sigma={sigma}, S={S}, H={H_vals[bits]}")

    print(f"\n  h_hat[{{0,1}}] = (1/2^m) * sum H(sigma) * sigma_0 * sigma_1")
    print(f"  = (2/2^m) * sum_(sigma_0=+1) H(sigma) * sigma_0 * sigma_1")
    print(f"  = (2/8) * sum_(sigma_0=+1) H(sigma) * sigma_1")
    total = 0
    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]
        contrib = H_vals[bits] * sigma[0] * sigma[1]
        total += contrib
        print(f"    sigma={sigma}, H={H_vals[bits]}, s0*s1={sigma[0]*sigma[1]:+d}, "
              f"contrib={contrib:+d}")
    print(f"  Total = {total}, h_hat = 2*{total}/{1<<m} = {2*total/(1<<m):.4f}")

    # PART 5: H(Paley) - H(non-Paley) analysis
    print(f"\n  H(Paley) = 189, H(non-Paley) = 175")
    print(f"  Difference = 14 = 2*p = 2*7")
    print(f"  h_hat = 2*(189-175*3)/(2^3) ... no, need full Walsh")

    # Check: at p=7, only one q-level (q=3), so all deg-2 h_hat are +/-3.5 = +/-p/2
    print(f"  |h_hat_2| = p/2 = {p/2} for ALL degree-2 at p=7")


if __name__ == '__main__':
    main()
