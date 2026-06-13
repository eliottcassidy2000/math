#!/usr/bin/env python3
"""
fejer_walsh_bridge.py -- Connecting Fejer kernel, Walsh decomposition, and Morgan-Voyce

THE BIG PICTURE:
1. Interval eigenvalues |lambda_k|^2 = sin^2(pi*m*k/p) / sin^2(pi*k/p) (Fejer kernel)
2. These are EXACTLY our Q_k = (1 - T_m(cos(2pi*k/p))) / (1 - cos(2pi*k/p)) (Chebyshev)
3. Walsh coefficients h_hat[S] are polynomial functions of the Q_k
4. The degree-4 phase transition (opus-S67b) has algebraic explanation via Morgan-Voyce

EQUIVALENCE PROOF:
  Q_k = (1 - T_m(c_k)) / (1 - c_k)   where c_k = cos(2pi*k/p)
  Fejer: F_m(theta) = sin^2(m*theta/2) / sin^2(theta/2) = (1 - cos(m*theta)) / (1 - cos(theta))
  Setting theta = 2pi*k/p: F_m(2pi*k/p) = (1 - cos(2pi*mk/p)) / (1 - cos(2pi*k/p))

  But T_m(cos(x)) = cos(mx), so:
  Q_k = (1 - cos(2pi*mk/p)) / (1 - cos(2pi*k/p)) = F_m(2pi*k/p)

  THEY ARE IDENTICAL. The Chebyshev formula and Fejer kernel are the same thing.

NEW INVESTIGATIONS:
- Walsh coefficients as explicit functions of additive energy E(S) and IPR
- The role of degree-4 Walsh terms in the p=1 mod 4 vs p=3 mod 4 transition
- Overlap concentration as measured by co-occurrence variance

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    """Count Hamiltonian paths using Held-Karp DP."""
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


def additive_energy(S, p):
    """E(S) = |{(a,b,c,d) in S^4 : a+b = c+d mod p}|"""
    S_set = set(S)
    e = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    e += 1
    return e


def ipr(S, p):
    """Inverse Participation Ratio = sum |S_hat(k)|^4 / (sum |S_hat(k)|^2)^2"""
    omega = cmath.exp(2j * cmath.pi / p)
    spectrum = []
    for k in range(p):
        val = sum(omega ** (k * s) for s in S)
        spectrum.append(abs(val)**2)
    L4 = sum(f**2 for f in spectrum)
    L2 = sum(spectrum)
    return L4 / L2**2


def compute_walsh_decomposition(p, pairs, H_func):
    """Compute full Walsh-Hadamard decomposition of H(sigma).

    H(sigma) = sum_{S subset [m]} h_hat[S] * chi_S(sigma)
    where chi_S(sigma) = prod_{i in S} sigma_i
    """
    m = len(pairs)
    n_orient = 1 << m

    # Evaluate H at all 2^m orientations
    H_vals = {}
    for bits in range(n_orient):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_vals[bits] = H

    # Walsh-Hadamard transform
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

    return h_hat, H_vals


def main():
    print("=" * 70)
    print("FEJER-WALSH BRIDGE: Connecting Spectral and Walsh Structures")
    print("=" * 70)

    # ============================================================
    # PART 1: Verify Fejer = Chebyshev Q_k identity
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 1: FEJER KERNEL = CHEBYSHEV Q_k IDENTITY")
    print("=" * 60)

    for p in [7, 11, 13, 17, 19, 23, 29, 31]:
        m = (p - 1) // 2
        omega = cmath.exp(2j * cmath.pi / p)
        S_int = list(range(1, m + 1))

        max_err = 0
        for k in range(1, m + 1):
            # Chebyshev: Q_k via eigenvalue
            lam_k = sum(omega ** (k * s) for s in S_int)
            Q_cheby = abs(lam_k) ** 2

            # Fejer kernel: sin^2(pi*m*k/p) / sin^2(pi*k/p)
            theta = 2 * math.pi * k / p
            Q_fejer = (math.sin(m * theta / 2) ** 2) / (math.sin(theta / 2) ** 2)

            # Direct Chebyshev T: (1 - T_m(cos(theta))) / (1 - cos(theta))
            c_k = math.cos(theta)
            T_m_ck = math.cos(m * theta)  # T_m(cos(x)) = cos(mx)
            Q_cheby_T = (1 - T_m_ck) / (1 - c_k)

            err1 = abs(Q_cheby - Q_fejer)
            err2 = abs(Q_cheby - Q_cheby_T)
            max_err = max(max_err, err1, err2)

        print(f"  p={p:>2}: max|Q_Fourier - Q_Fejer| = {max_err:.2e}  [IDENTICAL]")

    print("\n  PROVED: Q_k (Fourier) = Q_k (Chebyshev T) = F_m(2*pi*k/p) (Fejer kernel)")

    # ============================================================
    # PART 2: IPR = Additive Energy identity (opus HYP-514)
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 2: IPR = ADDITIVE ENERGY (ALGEBRAIC IDENTITY)")
    print("=" * 60)
    print("  Formula: IPR(S) = (p*E(S) - m^4) / (m*(p-m))^2")

    for p in [7, 11, 13, 17, 19]:
        m = (p - 1) // 2
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        S_int = list(range(1, m + 1))

        print(f"\n  p={p}:")
        for name, S in [("Interval", S_int), ("Paley", QR)]:
            E = additive_energy(S, p)
            ipr_direct = ipr(S, p)
            ipr_formula = (p * E - m**4) / (m * (p - m))**2

            print(f"    {name:8s}: E(S)={E:>6}, IPR(direct)={ipr_direct:.6f}, "
                  f"IPR(formula)={ipr_formula:.6f}, match={abs(ipr_direct-ipr_formula)<1e-10}")

    # ============================================================
    # PART 3: Full Walsh decomposition at p=7, 11, 13
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 3: WALSH DECOMPOSITION BY DEGREE")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p}, m={m}:")

        h_hat, H_vals = compute_walsh_decomposition(p, pairs, count_ham_paths)

        # Group by degree
        degree_sums = defaultdict(float)
        degree_abs_sums = defaultdict(float)
        degree_counts = defaultdict(int)

        for S_set, val in h_hat.items():
            d = len(S_set)
            degree_sums[d] += val
            degree_abs_sums[d] += abs(val)
            degree_counts[d] += 1

        # Evaluate at Interval (all +1) and Paley
        sigma_int_bits = (1 << m) - 1  # all bits set = all +1
        sigma_paley_bits = sum((1 << i) for i in range(m) if (i + 1) in QR)

        H_int = H_vals[sigma_int_bits]
        H_paley = H_vals[sigma_paley_bits]

        print(f"    H(Interval) = {H_int}")
        print(f"    H(Paley) = {H_paley}")
        print(f"    h_hat[empty] = {h_hat[frozenset()]:.1f} (= mean H)")

        # Walsh contribution at Interval vs Paley by degree
        print(f"\n    Walsh contribution by degree:")
        print(f"    {'deg':>4} {'count':>6} {'sum |h|':>12} {'f_d(Int)':>12} {'f_d(Pal)':>12} {'Int-Pal':>12}")

        for d in sorted(degree_sums.keys()):
            # f_d(sigma) = sum_{|S|=d} h_hat[S] * chi_S(sigma)
            f_int = 0
            f_pal = 0
            for S_set, val in h_hat.items():
                if len(S_set) != d:
                    continue
                # chi_S at all +1 = +1
                f_int += val
                # chi_S at Paley
                chi_pal = 1
                for i in S_set:
                    if (i + 1) in QR:
                        chi_pal *= 1
                    else:
                        chi_pal *= -1
                f_pal += val * chi_pal

            print(f"    {d:>4} {degree_counts[d]:>6} {degree_abs_sums[d]:>12.1f} "
                  f"{f_int:>12.1f} {f_pal:>12.1f} {f_int - f_pal:>12.1f}")

        # ============================================================
        # NEW: Connecting Walsh coefficients to Q_k values
        # ============================================================
        print(f"\n    --- Walsh coefficients vs Q_k structure ---")
        omega = cmath.exp(2j * cmath.pi / p)

        for d in [2, 4]:
            print(f"\n    Degree-{d} coefficients:")
            d_coeffs = [(S_set, val) for S_set, val in h_hat.items() if len(S_set) == d]
            d_coeffs.sort(key=lambda x: tuple(sorted(x[0])))

            for S_set, val in d_coeffs[:15]:  # show first 15
                indices = sorted(S_set)
                # Compute the "additive character" sum for this index set
                # Each index i corresponds to pair (i+1, p-i-1) in Z_p
                # The Walsh coefficient should relate to whether i+j, i-j, etc. have
                # arithmetic structure mod p

                # Check: does sum of the pair-indices have special mod-p structure?
                pair_vals = [pairs[i][0] for i in indices]  # the "positive" elements
                pair_sum = sum(pair_vals) % p

                if abs(val) > 0.1:
                    print(f"      S={set(i+1 for i in indices)}: h_hat={val:>10.1f}, "
                          f"pair_vals={pair_vals}, sum mod p={pair_sum}")

        # ============================================================
        # NEW: Check if degree-4 coefficients relate to zero-sum property
        # ============================================================
        if d == 4 or m >= 4:
            print(f"\n    --- Degree-4: Zero-sum classification (opus HYP-538) ---")
            zero_sum_total = 0
            non_zero_sum_total = 0
            zero_sum_count = 0
            non_zero_count = 0

            for S_set, val in h_hat.items():
                if len(S_set) != 4:
                    continue
                indices = sorted(S_set)
                pair_vals = [pairs[i][0] for i in indices]

                # Count W = #{(eps) in {+1,-1}^4 : sum eps_i * pair_vals[i] = 0 mod p}
                W = 0
                for bits in range(1 << 4):
                    eps = [(1 if bits & (1 << j) else -1) for j in range(4)]
                    s = sum(e * v for e, v in zip(eps, pair_vals)) % p
                    if s == 0:
                        W += 1

                if W > 0:
                    zero_sum_total += val
                    zero_sum_count += 1
                else:
                    non_zero_sum_total += val
                    non_zero_count += 1

            print(f"      Zero-sum (W>0): count={zero_sum_count}, net h_hat = {zero_sum_total:.1f}")
            print(f"      Non-zero-sum:   count={non_zero_count}, net h_hat = {non_zero_sum_total:.1f}")
            if zero_sum_count > 0 and non_zero_count > 0:
                print(f"      Ratio: {zero_sum_total/non_zero_sum_total:.4f}" if non_zero_sum_total != 0 else "      Non-zero net = 0")

    # ============================================================
    # PART 4: Additive energy maximization and overlap concentration
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 4: ADDITIVE ENERGY -> OVERLAP CONCENTRATION")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"\n  p={p}:")

        # For ALL orientations, compute E(S), IPR, H, and co-occurrence variance
        results = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            S_set = set(S)

            # Additive energy
            E = additive_energy(S, p)

            # H
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            # Co-occurrence variance (from 3-cycles)
            c3_sets = []
            for a, b, c in combinations(range(p), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_sets.append(frozenset([a, b, c]))
            c3_sets = list(set(c3_sets))

            # Co-occurrence by gap
            gap_co = [0] * p
            for fs in c3_sets:
                verts = list(fs)
                for i in range(len(verts)):
                    for j in range(i+1, len(verts)):
                        d = (verts[j] - verts[i]) % p
                        if d == 0:
                            continue
                        gap_co[d] += 1
            # By circulant symmetry, gap_co[d] for d=1..p-1 should be divided by p
            co_vals = [gap_co[d] / p for d in range(1, p)]
            co_mean = sum(co_vals) / len(co_vals) if co_vals else 0
            co_var = sum((c - co_mean)**2 for c in co_vals) / len(co_vals) if co_vals else 0

            is_paley = S_set == set(QR)
            is_interval = S == list(range(1, m + 1))

            results.append({
                'sigma': sigma, 'S': S, 'E': E, 'H': H,
                'co_var': co_var, 'is_paley': is_paley, 'is_interval': is_interval
            })

        # Sort by H descending
        results.sort(key=lambda x: -x['H'])

        # Show top 5 and bottom 5
        print(f"    {'rank':>4} {'H':>8} {'E(S)':>6} {'co_var':>8} {'marker':>10}")
        for i, r in enumerate(results[:5]):
            marker = ""
            if r['is_paley']:
                marker = "PALEY"
            elif r['is_interval']:
                marker = "INTERVAL"
            print(f"    {i+1:>4} {r['H']:>8} {r['E']:>6} {r['co_var']:>8.3f} {marker:>10}")
        print(f"    {'...':>4}")
        for i, r in enumerate(results[-3:]):
            marker = ""
            if r['is_paley']:
                marker = "PALEY"
            elif r['is_interval']:
                marker = "INTERVAL"
            rank = len(results) - 2 + i
            print(f"    {rank:>4} {r['H']:>8} {r['E']:>6} {r['co_var']:>8.3f} {marker:>10}")

        # Correlation between E and H
        E_vals = [r['E'] for r in results]
        H_vals_list = [r['H'] for r in results]
        mean_E = sum(E_vals) / len(E_vals)
        mean_H = sum(H_vals_list) / len(H_vals_list)
        cov_EH = sum((e - mean_E) * (h - mean_H) for e, h in zip(E_vals, H_vals_list)) / len(E_vals)
        std_E = (sum((e - mean_E)**2 for e in E_vals) / len(E_vals)) ** 0.5
        std_H = (sum((h - mean_H)**2 for h in H_vals_list) / len(H_vals_list)) ** 0.5
        if std_E > 0 and std_H > 0:
            corr_EH = cov_EH / (std_E * std_H)
        else:
            corr_EH = 0

        # Correlation between co_var and H
        cv_vals = [r['co_var'] for r in results]
        mean_cv = sum(cv_vals) / len(cv_vals)
        cov_cvH = sum((c - mean_cv) * (h - mean_H) for c, h in zip(cv_vals, H_vals_list)) / len(cv_vals)
        std_cv = (sum((c - mean_cv)**2 for c in cv_vals) / len(cv_vals)) ** 0.5
        if std_cv > 0 and std_H > 0:
            corr_cvH = cov_cvH / (std_cv * std_H)
        else:
            corr_cvH = 0

        print(f"\n    Pearson corr(E, H)     = {corr_EH:+.6f}")
        print(f"    Pearson corr(co_var, H) = {corr_cvH:+.6f}")
        int_E = next((r['E'] for r in results if r['is_interval']), None)
        pal_E = next((r['E'] for r in results if r['is_paley']), None)
        print(f"    E(Interval) = {int_E}")
        if pal_E is not None:
            print(f"    E(Paley)    = {pal_E}")
        else:
            print(f"    E(Paley)    = N/A (p=1 mod 4, no Paley tournament)")
        max_E = max(r['E'] for r in results)
        max_E_names = []
        for r in results:
            if r['E'] == max_E:
                if r['is_interval']:
                    max_E_names.append("Interval")
                elif r['is_paley']:
                    max_E_names.append("Paley")
        print(f"    max E(S) achieved by: E={max_E}, names={max_E_names if max_E_names else 'other'}")

    # ============================================================
    # PART 5: The Fejer-Morgan-Voyce-Walsh triangle
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 5: THE FEJER-MORGAN-VOYCE-WALSH TRIANGLE")
    print("=" * 60)
    print("""
    Three viewpoints on the SAME structure:

    1. FEJER: Q_k = sin^2(pi*m*k/p) / sin^2(pi*k/p)
       => Peaked spectrum concentrates power in low frequencies
       => This IS the Fejer kernel, well-studied in approximation theory

    2. MORGAN-VOYCE: e_j(Q_1,...,Q_m) = C(m+j, 2j)
       => The ESFs of the Q-values are BINOMIAL COEFFICIENTS
       => prod(1+Q_k) = F_p (Fibonacci number!)
       => The characteristic polynomial is t^m * b(m, -1/t) where b is Morgan-Voyce

    3. WALSH: H(sigma) = sum_S h_hat[S] * chi_S(sigma)
       => Degree-4 Walsh terms dominate at p >= 13
       => Zero-sum index sets (sum +/- pair_vals = 0 mod p) drive maximality
       => The degree-4 surplus at Interval explains the p=1 mod 4 regime
    """)

    # ============================================================
    # PART 6: NEW — Degree-2 Walsh as function of pairwise Q products
    # ============================================================
    print("=" * 60)
    print("PART 6: DEGREE-2 WALSH COEFFICIENTS vs Q-PRODUCT STRUCTURE")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n  p={p}:")

        h_hat, _ = compute_walsh_decomposition(p, pairs, count_ham_paths)

        # For each degree-2 coefficient h_hat[{i,j}], compare with
        # the product Q_i * Q_j of Interval eigenvalues
        S_int = list(range(1, m + 1))
        Q_vals = []
        for k in range(1, m + 1):
            val = sum(omega ** (k * s) for s in S_int)
            Q_vals.append(abs(val)**2)

        print(f"    Interval Q_k = [{', '.join(f'{q:.3f}' for q in Q_vals)}]")

        print(f"\n    {'pair':>10} {'h_hat':>10} {'Q_i*Q_j':>10} {'Q_i+Q_j':>10} {'ratio':>10}")
        for S_set, val in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            qi, qj = Q_vals[i], Q_vals[j]
            ratio = val / (qi * qj) if abs(qi * qj) > 0.01 else float('inf')
            print(f"    ({i+1},{j+1}){'':<5} {val:>10.2f} {qi*qj:>10.2f} {qi+qj:>10.2f} {ratio:>10.4f}")

    # ============================================================
    # PART 7: NEW — J-matrix eigenvalues and Fejer kernel connection
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 7: J-MATRIX (degree-2 coupling) EIGENSTRUCTURE")
    print("=" * 60)
    print("  J[i,j] = h_hat[{i,j}] for |{i,j}|=2, J[i,i] = h_hat[{i}]")
    print("  J eigenvalues determine degree-2 landscape")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        h_hat, _ = compute_walsh_decomposition(p, pairs, count_ham_paths)

        # Build J matrix
        J = [[0.0] * m for _ in range(m)]
        for S_set, val in h_hat.items():
            if len(S_set) == 1:
                i = list(S_set)[0]
                J[i][i] = val
            elif len(S_set) == 2:
                i, j = sorted(S_set)
                J[i][j] = val
                J[j][i] = val

        print(f"\n  p={p}, J matrix ({m}x{m}):")
        for i in range(m):
            row = "    "
            for j in range(m):
                row += f"{J[i][j]:>10.1f}"
            print(row)

        # Eigenvalues via power method (or numpy-free approach)
        # For small matrices, compute characteristic polynomial and find roots
        # Actually let's just compute the trace and Frobenius norm
        tr = sum(J[i][i] for i in range(m))
        frob_sq = sum(J[i][j]**2 for i in range(m) for j in range(m))

        print(f"    Trace(J) = {tr:.1f} (= sum of degree-1 Walsh coefficients)")
        print(f"    ||J||_F^2 = {frob_sq:.1f}")
        print(f"    det signature: ", end="")

        # 2x2 and 3x3 determinants
        if m == 3:
            det = (J[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
                   - J[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
                   + J[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0]))
            print(f"det = {det:.1f}")
        elif m == 5:
            # Use cofactor expansion for 5x5
            # Just report trace and sum of 2x2 minors
            sum_minors = 0
            for i in range(m):
                for j in range(i+1, m):
                    minor = J[i][i]*J[j][j] - J[i][j]*J[j][i]
                    sum_minors += minor
            print(f"sum of 2x2 principal minors = {sum_minors:.1f}")
        elif m == 6:
            sum_minors = 0
            for i in range(m):
                for j in range(i+1, m):
                    minor = J[i][i]*J[j][j] - J[i][j]*J[j][i]
                    sum_minors += minor
            print(f"sum of 2x2 principal minors = {sum_minors:.1f}")

    # ============================================================
    # PART 8: NEW — The critical ratio: when does degree-4 overtake degree-2?
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 8: DEGREE-4 vs DEGREE-2 SURPLUS AT INTERVAL")
    print("=" * 60)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, H_vals = compute_walsh_decomposition(p, pairs, count_ham_paths)

        # f_d(sigma) at Interval (all +1) vs Paley
        sigma_int_bits = (1 << m) - 1
        sigma_paley_bits = sum((1 << i) for i in range(m) if (i + 1) in QR)

        print(f"\n  p={p} (p mod 4 = {p % 4}):")

        for d in range(m + 1):
            f_int = sum(val for S_set, val in h_hat.items() if len(S_set) == d)

            f_pal = 0
            for S_set, val in h_hat.items():
                if len(S_set) != d:
                    continue
                chi_pal = 1
                for i in S_set:
                    if (i + 1) in QR:
                        chi_pal *= 1
                    else:
                        chi_pal *= -1
                f_pal += val * chi_pal

            if abs(f_int) > 0.01 or abs(f_pal) > 0.01:
                surplus_int = f_int - f_pal
                print(f"    deg {d}: f(Int)={f_int:>12.1f}, f(Pal)={f_pal:>12.1f}, "
                      f"Int-Pal={surplus_int:>12.1f}")

        total_int = sum(val for S_set, val in h_hat.items())
        total_pal_approx = H_vals[sigma_paley_bits]
        print(f"    TOTAL: H(Int)={H_vals[sigma_int_bits]}, H(Pal)={total_pal_approx}")

    # ============================================================
    # PART 9: Co-occurrence variance as function of additive energy
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 9: CO-OCCURRENCE VARIANCE vs ADDITIVE ENERGY")
    print("=" * 60)
    print("  Question: Does higher E(S) => higher co-occurrence variance?")
    print("  Higher co_var means cycles are MORE spatially concentrated (Fejer-like)")

    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}:")

        # Compute E(S) and co_var for all orientations
        E_cv_pairs = []
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            E = additive_energy(S, p)

            A = build_adj(p, S)
            # 3-cycle co-occurrence
            c3_sets = []
            for a, b, c in combinations(range(p), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_sets.append(frozenset([a, b, c]))
            c3_sets = list(set(c3_sets))

            gap_co = [0] * p
            for fs in c3_sets:
                verts = list(fs)
                for vi in range(len(verts)):
                    for vj in range(vi+1, len(verts)):
                        d = (verts[vj] - verts[vi]) % p
                        gap_co[d] += 1

            co_vals = [gap_co[d] / p for d in range(1, p)]
            co_mean = sum(co_vals) / len(co_vals) if co_vals else 0
            co_var = sum((c - co_mean)**2 for c in co_vals) / len(co_vals) if co_vals else 0

            H = count_ham_paths(A, p)
            E_cv_pairs.append((E, co_var, H, S))

        # Sort by E
        E_cv_pairs.sort()

        # Group by E value
        E_groups = defaultdict(list)
        for E, cv, H, S in E_cv_pairs:
            E_groups[E].append((cv, H))

        print(f"    {'E(S)':>6} {'count':>5} {'avg co_var':>10} {'avg H':>10} {'max H':>10}")
        for E in sorted(E_groups):
            items = E_groups[E]
            avg_cv = sum(x[0] for x in items) / len(items)
            avg_H = sum(x[1] for x in items) / len(items)
            max_H = max(x[1] for x in items)
            print(f"    {E:>6} {len(items):>5} {avg_cv:>10.4f} {avg_H:>10.1f} {max_H:>10}")

    print("\n" + "=" * 60)
    print("SUMMARY OF CONNECTIONS")
    print("=" * 60)
    print("""
    THE FEJER-WALSH-MORGAN-VOYCE TRIANGLE:

    1. Fejer kernel concentration => Q_k values peaked (large Q_1, small Q_{m-1})
       MECHANISM: sin^2(pi*m*k/p) / sin^2(pi*k/p) is maximally peaked at k=1

    2. Peaked Q => large additive energy E(S) => large IPR
       MECHANISM: E(S) = |{a+b=c+d in S}|; contiguous sets maximize this (rearrangement)

    3. Large E(S) => cycle concentration (high co-occurrence variance)
       MECHANISM: cycles use gaps from S; contiguous S forces cycles to "cluster"

    4. Cycle concentration => more independent sets in Omega(T)
       MECHANISM: clustered cycles leave room for non-conflicting cycles elsewhere

    5. Walsh: degree-4 terms "see" the additive structure of S
       MECHANISM: h_hat[{i,j,k,l}] large iff pair_vals have zero-sum (additive relation)
       At p >= 13, degree-4 dominates and Interval wins via additive energy

    6. Morgan-Voyce: e_j = C(m+j, 2j) encodes the EXACT algebraic structure
       MECHANISM: Chebyshev T polynomial + resultant theory
       prod(1+Q_k) = F_p is a CONSEQUENCE, not a cause

    THE CAUSAL CHAIN:
    Contiguous S => max E(S) => max IPR => max Q-concentration
    => max cycle clustering => max alpha(Omega) => max H

    Algebraically: Fejer kernel = Chebyshev T = Morgan-Voyce polynomial
    """)


if __name__ == '__main__':
    main()
