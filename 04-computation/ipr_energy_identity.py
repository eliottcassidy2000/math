#!/usr/bin/env python3
"""
ipr_energy_identity.py -- Deriving the CORRECT IPR-Energy identity

Opus claimed: IPR(S) = (p*E(S) - m^4) / (m*(p-m))^2
But direct computation shows this is WRONG. Let's find the right formula.

DEFINITIONS:
- S_hat(k) = sum_{s in S} omega^{ks} where omega = e^{2*pi*i/p}
- |S_hat(k)|^2 = Fourier power spectrum
- IPR = sum_k |S_hat(k)|^4 / (sum_k |S_hat(k)|^2)^2
- E(S) = |{(a,b,c,d) in S^4 : a+b = c+d mod p}| = additive energy

PARSEVAL:
- sum_{k=0}^{p-1} |S_hat(k)|^2 = p*m (Parseval identity)
- sum_{k=0}^{p-1} |S_hat(k)|^4 = p*E(S) (Parseval for 4th moment!)

So IPR = p*E(S) / (p*m)^2 = E(S) / (p*m^2)

Wait, but we included k=0 where S_hat(0) = m. So:
- |S_hat(0)|^2 = m^2
- |S_hat(0)|^4 = m^4
- sum_{k>0} |S_hat(k)|^2 = p*m - m^2 = m(p-m)
- sum_{k>0} |S_hat(k)|^4 = p*E(S) - m^4

If IPR counts only k>0 (nontrivial frequencies):
  IPR_{k>0} = (p*E(S) - m^4) / (m(p-m))^2

If IPR counts ALL k:
  IPR_{all} = p*E(S) / (p*m)^2 = E(S) / (p*m^2)

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math


def additive_energy(S, p):
    S_set = set(S)
    e = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    e += 1
    return e


def main():
    print("=" * 70)
    print("IPR-ENERGY IDENTITY: FINDING THE CORRECT FORMULA")
    print("=" * 70)

    # PART 1: Verify Parseval for 4th moment
    print("\n--- PART 1: PARSEVAL IDENTITIES ---")
    print("sum |S_hat(k)|^2 = p*m  (Parseval)")
    print("sum |S_hat(k)|^4 = p*E(S)  (4th moment = energy)")

    for p in [7, 11, 13, 17, 19, 23]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        for name, S in [("Interval", S_int), ("Paley", QR)]:
            omega = cmath.exp(2j * cmath.pi / p)
            spectrum = []
            for k in range(p):
                val = sum(omega ** (k * s) for s in S)
                spectrum.append(abs(val)**2)

            L2 = sum(spectrum)
            L4 = sum(f**2 for f in spectrum)
            E = additive_energy(S, p)

            # Check Parseval 2nd moment
            parseval_2 = p * m
            # Check Parseval 4th moment
            parseval_4 = p * E

            print(f"  p={p:>2}, {name:8s}: sum|S|^2 = {L2:.1f} vs p*m = {parseval_2}, "
                  f"sum|S|^4 = {L4:.1f} vs p*E = {parseval_4}, "
                  f"match_2: {abs(L2 - parseval_2) < 0.01}, match_4: {abs(L4 - parseval_4) < 0.01}")

    # PART 2: The three IPR definitions
    print("\n--- PART 2: THREE IPR FORMULAS ---")
    print("IPR_all = E(S) / (p*m^2)")
    print("IPR_k>0 = (p*E(S) - m^4) / (m(p-m))^2")
    print("IPR_direct = sum|S_hat|^4 / (sum|S_hat|^2)^2")

    for p in [7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        S_int = list(range(1, m + 1))
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p}:")
        for name, S in [("Interval", S_int), ("Paley", QR)]:
            omega = cmath.exp(2j * cmath.pi / p)
            spectrum = []
            for k in range(p):
                val = sum(omega ** (k * s) for s in S)
                spectrum.append(abs(val)**2)

            E = additive_energy(S, p)

            # IPR_all (including k=0)
            L4_all = sum(f**2 for f in spectrum)
            L2_all = sum(spectrum)
            IPR_all_direct = L4_all / L2_all**2
            IPR_all_formula = E / (p * m**2)

            # IPR_nontrivial (k > 0 only)
            L4_nontriv = sum(f**2 for f in spectrum[1:])
            L2_nontriv = sum(spectrum[1:])
            IPR_nontriv_direct = L4_nontriv / L2_nontriv**2
            IPR_nontriv_formula = (p * E - m**4) / (m * (p - m))**2

            print(f"    {name:8s}: E={E:>6}")
            print(f"      IPR_all:     direct={IPR_all_direct:.8f}, formula={IPR_all_formula:.8f}, "
                  f"match={abs(IPR_all_direct - IPR_all_formula) < 1e-8}")
            print(f"      IPR_nontriv: direct={IPR_nontriv_direct:.8f}, formula={IPR_nontriv_formula:.8f}, "
                  f"match={abs(IPR_nontriv_direct - IPR_nontriv_formula) < 1e-8}")

    # PART 3: Which IPR correlates with H?
    print("\n--- PART 3: WHICH IPR PREDICTS H? ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}:")

        data = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            E = additive_energy(S, p)

            # Compute H
            S_set = set(S)
            A = [[0]*p for _ in range(p)]
            for i in range(p):
                for s in S_set:
                    A[i][(i + s) % p] = 1

            # Held-Karp for H
            dp = [[0]*p for _ in range(1 << p)]
            for v in range(p):
                dp[1 << v][v] = 1
            for mask in range(1, 1 << p):
                for v in range(p):
                    if not (mask & (1 << v)):
                        continue
                    if dp[mask][v] == 0:
                        continue
                    for w in range(p):
                        if mask & (1 << w):
                            continue
                        if A[v][w]:
                            dp[mask | (1 << w)][w] += dp[mask][v]
            full = (1 << p) - 1
            H = sum(dp[full][v] for v in range(p))

            # IPR_all
            IPR_all = E / (p * m**2)

            # IPR_nontrivial
            IPR_nontriv = (p * E - m**4) / (m * (p - m))**2

            data.append({'E': E, 'H': H, 'IPR_all': IPR_all, 'IPR_nontriv': IPR_nontriv, 'S': S})

        # Correlations
        def pearson(x_vals, y_vals):
            n = len(x_vals)
            mx = sum(x_vals) / n
            my = sum(y_vals) / n
            cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x_vals, y_vals)) / n
            sx = (sum((xi - mx)**2 for xi in x_vals) / n) ** 0.5
            sy = (sum((yi - my)**2 for yi in y_vals) / n) ** 0.5
            return cov / (sx * sy) if sx > 0 and sy > 0 else 0

        H_vals = [d['H'] for d in data]
        E_vals = [d['E'] for d in data]
        IPR_all_vals = [d['IPR_all'] for d in data]
        IPR_nontriv_vals = [d['IPR_nontriv'] for d in data]

        print(f"    corr(E, H)          = {pearson(E_vals, H_vals):+.6f}")
        print(f"    corr(IPR_all, H)    = {pearson(IPR_all_vals, H_vals):+.6f}")
        print(f"    corr(IPR_nontriv, H)= {pearson(IPR_nontriv_vals, H_vals):+.6f}")

        # Top 3 by H
        data.sort(key=lambda x: -x['H'])
        print(f"    Top 3 by H:")
        for d in data[:3]:
            print(f"      H={d['H']:>8}, E={d['E']:>4}, IPR_all={d['IPR_all']:.6f}, S={d['S']}")
        print(f"    Bottom 3 by H:")
        for d in data[-3:]:
            print(f"      H={d['H']:>8}, E={d['E']:>4}, IPR_all={d['IPR_all']:.6f}, S={d['S']}")

    # PART 4: The ratio E(S)/m^2 and its relationship to H
    print("\n--- PART 4: DEEPER ADDITIVE ENERGY ANALYSIS ---")
    print("E(S)/m^2 = normalized additive energy")
    print("For m-element subsets of Z_p:")
    print("  E_min >= m^2 (trivial pairs a=c, b=d)")
    print("  E(Interval) is max by Ruzsa-Szemeredi / rearrangement")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}, m={m}:")

        # All E values
        E_vals_all = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            E_vals_all[bits] = E

        E_unique = sorted(set(E_vals_all.values()))
        print(f"    Distinct E values: {E_unique}")
        print(f"    E(Interval)/m^2 = {E_unique[-1]}/{m**2} = {E_unique[-1]/m**2:.4f}")
        print(f"    E_min/m^2 = {E_unique[0]}/{m**2} = {E_unique[0]/m**2:.4f}")
        print(f"    E_gap = {E_unique[-1] - E_unique[-2]}" if len(E_unique) > 1 else "")

    # PART 5: Walsh degree-2 coefficient sum as function of E
    print("\n--- PART 5: WALSH DEGREE-2 SUM vs E(S) ---")
    print("f_2(sigma) = sum_{|S|=2} h_hat[S] * chi_S(sigma)")
    print("Is f_2 a function of E(sigma)?")

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        # Full Walsh decomposition
        n_orient = 1 << m
        H_vals = {}
        E_all = {}
        for bits in range(n_orient):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            A = build_adj_simple(p, S)
            H = count_ham_paths_simple(A, p)
            E = additive_energy(S, p)
            H_vals[bits] = H
            E_all[bits] = E

        # Walsh transform
        h_hat = {}
        for S_bits in range(n_orient):
            S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
            val = 0
            for bits in range(n_orient):
                sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
                chi = 1
                for i in S_set:
                    chi *= sigma[i]
                val += H_vals[bits] * chi
            h_hat[S_set] = val / n_orient

        # Compute f_2 and f_4 for each orientation
        print(f"\n  p={p}:")
        print(f"    {'bits':>5} {'E':>5} {'H':>8} {'f_2':>10} {'f_4':>10}")

        for bits in range(n_orient):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            f2 = 0
            f4 = 0
            for S_set, val in h_hat.items():
                chi = 1
                for i in S_set:
                    chi *= sigma[i]
                if len(S_set) == 2:
                    f2 += val * chi
                elif len(S_set) == 4:
                    f4 += val * chi

            print(f"    {bits:>5} {E_all[bits]:>5} {H_vals[bits]:>8} {f2:>10.1f} {f4:>10.1f}")

    # PART 6: NEW — The sign pattern of degree-2 Walsh coefficients
    print("\n--- PART 6: SIGN PATTERN OF DEGREE-2 WALSH COEFFICIENTS ---")
    print("h_hat[{i,j}] sign depends on whether pair_sum falls in QR or NQR")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        # Quick Walsh decomposition for degree-2
        n_orient = 1 << m
        H_vals = {}
        for bits in range(n_orient):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            A = build_adj_simple(p, S)
            H = count_ham_paths_simple(A, p)
            H_vals[bits] = H

        h_hat = {}
        for S_bits in range(n_orient):
            S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
            val = sum(H_vals[bits] * prod_chi(bits, S_set, m) for bits in range(n_orient))
            h_hat[S_set] = val / n_orient

        print(f"\n  p={p}:")
        print(f"    {'pair':>10} {'pair_sum':>8} {'sum in QR':>10} {'h_hat':>10} {'sign':>6}")

        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            ps = (pairs[i][0] + pairs[j][0]) % p
            in_qr = ps in QR if ps != 0 else "N/A"
            sign = "+" if val > 0 else "-" if val < 0 else "0"
            print(f"    ({i+1},{j+1}){'':<5} {ps:>8} {str(in_qr):>10} {val:>10.1f} {sign:>6}")

        # Check: is sign(h_hat[{i,j}]) = chi(pair_sum)?
        correct = 0
        total = 0
        for S_set, val in h_hat.items():
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            ps = (pairs[i][0] + pairs[j][0]) % p
            if ps == 0 or abs(val) < 0.01:
                continue
            total += 1
            expected_sign = 1 if ps in QR else -1
            actual_sign = 1 if val > 0 else -1
            if expected_sign == actual_sign:
                correct += 1

        print(f"    sign(h_hat) = chi(pair_sum): {correct}/{total} correct")


def build_adj_simple(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths_simple(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def prod_chi(bits, S_set, m):
    """Compute chi_S(sigma) = prod_{i in S} sigma_i."""
    chi = 1
    for i in S_set:
        if bits & (1 << i):
            chi *= 1
        else:
            chi *= -1
    return chi


if __name__ == '__main__':
    main()
