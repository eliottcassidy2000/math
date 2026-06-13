#!/usr/bin/env python3
"""
d2n_sign_analysis.py -- Does D^{2n} degree-2 Walsh have sign = chi(ab) for EACH n?

From sign_mechanism.py:
- At p=7: D^4 alone gives all degree-2 Walsh (D^4 dominates)
- At p=11: D^6 dominates D^4 by factor ~19x, but BOTH have sign aligned with chi(ab)
- Magnitude ratio for resonant vs non-resonant pairs = 3p

KEY DISCOVERY SOUGHT: if EVERY D^{2n} has sign(deg-2 Walsh) = chi(ab),
then any positive linear combination also has sign = chi(ab), proving
the product law structurally.

This script also derives the ALGEBRAIC structure:
- D^4 deg-2 = T_pure + T_mixed = T_pure - 3*T_pure = -2*T_pure
- T_pure uses only 3-resonance (3a = +/-b mod p)
- The Gauss sum g_hat(a) = chi(a)*sqrt(p) factorization at p=3 mod 4

Author: kind-pasteur-2026-03-12-S60
"""

import math
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def D2n_deg2(a, b, p, n_power):
    """Degree-2 Walsh of sum_{k=1}^{p-1} D_k^{2n} at position (a,b).

    Uses the extraction formula:
    coeff_xy = (1/4)[f(1,1) - f(1,-1) - f(-1,1) + f(-1,-1)]
    where f(sx,sy) = (sx*Sa + sy*Sb + R)^{2n}.
    """
    m = (p - 1) // 2
    total = 0
    for k in range(1, p):
        Sa = math.sin(2 * math.pi * k * a / p)
        Sb = math.sin(2 * math.pi * k * b / p)
        R = sum(math.sin(2 * math.pi * k * c / p)
                for c in range(1, m + 1) if c != a and c != b)

        vals = {}
        for sx in [1, -1]:
            for sy in [1, -1]:
                D_val = sx * Sa + sy * Sb + R
                vals[(sx, sy)] = D_val ** (2 * n_power)

        coeff = (vals[(1,1)] - vals[(1,-1)] - vals[(-1,1)] + vals[(-1,-1)]) / 4.0
        total += coeff
    return total


def main():
    print("=" * 70)
    print("D^{2n} DEGREE-2 WALSH: SIGN vs CHI(ab) BY POWER")
    print("=" * 70)

    # PART 1: Each D^{2n} separately
    print("\n--- PART 1: SIGN ANALYSIS ---")

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        chi_neg1 = legendre(-1, p)
        print(f"\n  p={p} ({p%4} mod 4), m={m}, chi(-1)={chi_neg1:+d}:")
        max_n = 4 if p <= 13 else 3

        # Header
        hdr = f"  {'(a,b)':>8} {'chi':>5}"
        for n in range(1, max_n+1):
            hdr += f"  {'D^'+str(2*n):>12} {'sgn':>4}"
        print(hdr)

        match_count = {n: [0, 0] for n in range(1, max_n+1)}  # [matches, nonzero]
        total_pairs = 0

        for a in range(1, m+1):
            for b in range(a+1, m+1):
                chi_ab = legendre(a*b, p)
                total_pairs += 1
                line = f"  ({a},{b}){' ':>3} {chi_ab:>+5d}"

                for n in range(1, max_n+1):
                    val = D2n_deg2(a, b, p, n)
                    sign_val = 1 if val > 0.01 else (-1 if val < -0.01 else 0)
                    if sign_val != 0:
                        match_count[n][1] += 1
                        if sign_val == chi_ab:
                            match_count[n][0] += 1
                    line += f"  {val:>12.3f} {sign_val:>+4d}"
                print(line)

        print(f"\n  Summary for p={p}:")
        for n in range(1, max_n+1):
            matches, nonzero = match_count[n]
            print(f"    D^{2*n}: sign=chi in {matches}/{nonzero} nonzero "
                  f"(of {total_pairs} pairs)")

    # PART 2: Gauss sum connection
    print("\n" + "=" * 70)
    print("PART 2: GAUSS SUM STRUCTURE")
    print("=" * 70)
    print("g(a) = sum_{k=1}^{p-1} chi(k)*sin(2pi*ka/p)")
    print("At p=3 mod 4: g(a) = chi(a)*sqrt(p)")

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        print(f"\n  p={p}:")

        g = {}
        for a in range(1, m+1):
            g[a] = sum(
                (1 if k in QR else -1) * math.sin(2 * math.pi * k * a / p)
                for k in range(1, p)
            )
            expected = legendre(a, p) * math.sqrt(p)
            print(f"    g({a}) = {g[a]:>10.6f}, chi(a)*sqrt(p) = {expected:>10.6f}, "
                  f"match={abs(g[a]-expected) < 0.001}")

        # Check D^{2n} / (g(a)*g(b)) ratios
        print(f"\n    D^{{2n}} / (g(a)*g(b)) ratios:")
        for n in range(1, 4):
            ratios = []
            print(f"    n={n} (D^{2*n}):")
            for a in range(1, m+1):
                for b in range(a+1, m+1):
                    val = D2n_deg2(a, b, p, n)
                    prod = g[a] * g[b]
                    if abs(prod) > 0.01:
                        ratio = val / prod
                        ratios.append(ratio)
                        print(f"      ({a},{b}): D^{2*n}={val:>10.4f}, "
                              f"g*g={prod:>10.4f}, ratio={ratio:>10.6f}")

            if ratios:
                # Check if constant
                r0 = ratios[0]
                is_const = all(abs(r - r0) < 0.001 for r in ratios)
                print(f"      Constant? {is_const} (value ~ {r0:.6f})")

    # PART 3: D^4 = -2*T_pure algebraic derivation
    print("\n" + "=" * 70)
    print("PART 3: D^4 ALGEBRAIC STRUCTURE")
    print("=" * 70)
    print("D^4 deg-2 = T_pure + T_mixed")
    print("T_pure = sum_k 4*(S_a*S_b^3 + S_a^3*S_b)")
    print("T_mixed = sum_k 12*S_a*S_b*sum_{c!=a,b} S_c^2")
    print()
    print("CLAIM: T_mixed = -3 * T_pure")
    print("PROOF: sum_{c=1}^m S_c^2 = (p-1)/4  [half-angle identity]")
    print("  sum_{c!=a,b} S_c^2 = (p-1)/4 - S_a^2 - S_b^2")
    print("  T_mixed = 12*sum_k S_a*S_b*[(p-1)/4 - S_a^2 - S_b^2]")
    print("         = 3(p-1)*[D^2 deg-2] - 12*sum_k S_a*S_b*(S_a^2+S_b^2)")
    print("         = 3(p-1)*0 - 3*T_pure  [since D^2 deg-2 = 0]")
    print("  => D^4 deg-2 = T_pure - 3*T_pure = -2*T_pure")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        print(f"\n  p={p}: Verification:")
        for a in range(1, m+1):
            for b in range(a+1, m+1):
                T_pure = 0
                T_mixed = 0
                for k in range(1, p):
                    Sa = math.sin(2*math.pi*k*a/p)
                    Sb = math.sin(2*math.pi*k*b/p)
                    T_pure += 4 * (Sa * Sb**3 + Sa**3 * Sb)
                    for c in range(1, m+1):
                        if c == a or c == b:
                            continue
                        Sc = math.sin(2*math.pi*k*c/p)
                        T_mixed += 12 * Sa * Sb * Sc**2

                ratio = T_mixed / T_pure if abs(T_pure) > 0.01 else float('inf')
                D4 = D2n_deg2(a, b, p, 2)
                print(f"    ({a},{b}): T_mix/T_pure = {ratio:>8.4f}, "
                      f"D^4 = {D4:>10.3f}, -2*T_pure = {-2*T_pure:>10.3f}, "
                      f"match = {abs(D4 - (-2*T_pure)) < 0.01}")

    # PART 4: T_pure via delta functions
    print("\n" + "=" * 70)
    print("PART 4: T_PURE DELTA-FUNCTION FORMULA")
    print("=" * 70)
    print("T_pure = sum_k 4*S_a*S_b*(S_a^2 + S_b^2)")
    print("= 2*sum_k [cos(k(a-b)) - cos(k(a+b))]*[1 - cos(k(a-b))*cos(k(a+b))]")
    print()
    print("Expanding and using delta_sum(n) = p-1 if p|n, else -1:")
    print("T_pure depends on delta_sum at: a-b, a+b, 3a-b, a-3b, 3a+b, a+3b")

    def ds(n, p):
        return p - 1 if n % p == 0 else -1

    for p in [7, 11, 13, 19]:
        m = (p - 1) // 2
        print(f"\n  p={p}:")

        for a in range(1, m+1):
            for b in range(a+1, m+1):
                # Formula terms
                t1 = ds(a-b, p)
                t2 = -ds(a+b, p)
                t3 = -(ds(a+b, p)/2 + (ds(3*a-b, p) + ds(a-3*b, p))/4)
                t4 = ds(a-b, p)/2 + (ds(3*a+b, p) + ds(-(a+3*b), p))/4

                T_pure_formula = 2 * (t1 + t2 + t3 + t4)
                D4_formula = -2 * T_pure_formula

                # Direct computation for verification
                D4_direct = D2n_deg2(a, b, p, 2)
                chi_ab = legendre(a*b, p)

                # 3-resonance flags
                res = []
                if (3*a - b) % p == 0: res.append("3a=b")
                if (3*a + b) % p == 0: res.append("3a=-b")
                if (a - 3*b) % p == 0: res.append("a=3b")
                if (a + 3*b) % p == 0: res.append("a=-3b")
                res_str = ",".join(res) if res else "none"

                print(f"    ({a},{b}): D^4_formula={D4_formula:>8.2f}, "
                      f"D^4_direct={D4_direct:>8.2f}, "
                      f"chi={chi_ab:+d}, resonance={res_str}")

    # PART 5: Why does the sign work?
    print("\n" + "=" * 70)
    print("PART 5: SIGN MECHANISM FOR D^{2n}")
    print("=" * 70)
    print()
    print("At p=3 mod 4: g(a) = chi(a)*sqrt(p), so g(a)*g(b) = chi(ab)*p")
    print("If D^{2n} deg-2 = R_n * g(a)*g(b) with R_n constant,")
    print("then sign(D^{2n} deg-2) = sign(R_n) * chi(ab)")
    print()
    print("For sign(h_hat) = chi(ab), we need sign(R_n) to be CONSISTENT")
    print("across all contributing D^{2n} terms (or the dominant one wins).")

    for p in [7, 11]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        print(f"\n  p={p}:")

        g = {}
        for a in range(1, m+1):
            g[a] = sum(
                (1 if k in QR else -1) * math.sin(2*math.pi*k*a/p)
                for k in range(1, p)
            )

        max_n = 5 if p == 7 else 4
        R_values = {}
        for n in range(1, max_n+1):
            ratios = []
            for a in range(1, m+1):
                for b in range(a+1, m+1):
                    val = D2n_deg2(a, b, p, n)
                    prod = g[a] * g[b]
                    if abs(prod) > 0.01:
                        ratios.append(val / prod)

            if ratios:
                r0 = ratios[0]
                is_const = all(abs(r - r0) < 0.001 for r in ratios)
                R_values[n] = r0 if is_const else None
                print(f"    D^{2*n}: R_{n} = {r0:>12.6f}, constant={is_const}, "
                      f"sign={'+' if r0 > 0 else '-'}")

        # The coefficient of D^{2n} in h_hat via trace expansion
        # h_hat = sum_{k even} (1/k) * coeff_D2n * D^{2n}_deg2
        # where coeff_D2n comes from Re((-1/2+iD)^k)
        print(f"\n    Combined sign analysis:")
        print(f"    h_hat = sum_n w_n * R_n * g(a)*g(b) = [sum_n w_n * R_n] * chi(ab)*p")
        print(f"    Sign = sign(sum_n w_n * R_n)")

        # Compute effective weights from trace expansion
        # tr(A^k) contributes to h_hat via cycle counts
        # The degree-2 Walsh of c_k = (1/k)*tr(A^k)
        # c_k contributes to H via OCF: H = I(Omega, 2)
        # But the Walsh coefficients come from ALL even k

    # PART 6: The overlap weight interpretation
    print("\n" + "=" * 70)
    print("PART 6: OVERLAP WEIGHT = FEJER KERNEL OVERLAP")
    print("=" * 70)
    print()
    print("The overlap matrix W[k1,k2] = sum_{l=1}^m sin(2pi k1 l/p)*sin(2pi k2 l/p)")
    print("For p=3 mod 4, writing sin(kl) in terms of characters:")
    print("sin(kl) = Im(omega^{kl}) where omega = e^{2pi i/p}")
    print()
    print("The Fejer kernel F_m(t) = (1/m)|sum_{k=0}^{m-1} e^{ikt}|^2")
    print("is the spectral weight. Cycle overlaps are measured by F_m convolutions.")

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\n  p={p}, m={m}:")

        # Compute W matrix
        W = [[0.0]*(p-1) for _ in range(p-1)]
        for k1 in range(1, p):
            for k2 in range(1, p):
                val = sum(
                    math.sin(2*math.pi*k1*l/p) * math.sin(2*math.pi*k2*l/p)
                    for l in range(1, m+1)
                )
                W[k1-1][k2-1] = val

        print(f"    W diagonal: {W[0][0]:.4f} (expect (p-1)/4 = {(p-1)/4:.4f})")
        off_vals = [W[i][j] for i in range(p-1) for j in range(p-1) if i != j]
        if off_vals:
            print(f"    W off-diag: min={min(off_vals):.4f}, max={max(off_vals):.4f}")

        # W[k1,k2] = (1/2)[delta_cos(k1-k2) - delta_cos(k1+k2)]
        # where delta_cos(n) = sum_{l=1}^m cos(2pi*n*l/p)
        # = (1/2)[sum_{l=1}^m e^{2pi i n l/p} + e^{-2pi i n l/p}]
        # For n != 0 mod p: delta_cos(n) = -1/2 (from geometric sum)
        # For n = 0 mod p: delta_cos(0) = m = (p-1)/2
        print(f"\n    W[k1,k2] = (1/2)[cos_sum(k1-k2) - cos_sum(k1+k2)]")
        print(f"    cos_sum(n) = sum_l cos(2pi*n*l/p) = m if p|n, else -1/2")
        print(f"    => W[k,k] = (1/2)[m - cos_sum(2k)] = (1/2)[(p-1)/2 + 1/2] = (2p-1)/4")
        predicted_diag = (2*p-1)/4
        print(f"    Predicted diagonal: {predicted_diag:.4f}")

        # Hmm, let me verify cos_sum
        for n_test in [0, 1, 2, m]:
            cs = sum(math.cos(2*math.pi*n_test*l/p) for l in range(1, m+1))
            print(f"    cos_sum({n_test}) = {cs:.6f}")

    # PART 7: p=1 mod 4 analysis
    print("\n" + "=" * 70)
    print("PART 7: p=1 mod 4 (chi(-1)=+1, NO Gauss sum factorization)")
    print("=" * 70)

    for p in [13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        print(f"\n  p={p} ({p%4} mod 4), m={m}:")
        print(f"  QR = {sorted(QR)}")

        g = {}
        for a in range(1, m+1):
            g[a] = sum(
                (1 if k in QR else -1) * math.sin(2*math.pi*k*a/p)
                for k in range(1, p)
            )
            expected = legendre(a, p) * math.sqrt(p)
            print(f"    g({a}) = {g[a]:>10.6f}, chi(a)*sqrt(p) = {expected:>10.6f}")

        print(f"\n    g(a) ~ 0 at p=1 mod 4! The sine Gauss sum vanishes.")
        print(f"    This means the product law sign(h_hat)=chi(ab) FAILS at p=1 mod 4.")
        print(f"    At p=1 mod 4, the degree-2 Walsh is controlled by COSINE sums instead.")

        # Verify: D^{2n} at p=13
        print(f"\n    D^{{2n}} degree-2 Walsh at p={p}:")
        hdr = f"    {'(a,b)':>8} {'chi':>5}"
        for n in range(1, 4):
            hdr += f"  {'D^'+str(2*n):>12}"
        print(hdr)

        for a in range(1, m+1):
            for b in range(a+1, m+1):
                chi_ab = legendre(a*b, p)
                line = f"    ({a},{b}){' ':>3} {chi_ab:>+5d}"
                for n in range(1, 4):
                    val = D2n_deg2(a, b, p, n)
                    line += f"  {val:>12.3f}"
                print(line)


if __name__ == '__main__':
    main()
