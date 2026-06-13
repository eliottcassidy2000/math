#!/usr/bin/env python3
"""
ghat_factorization_deep.py -- Deep analysis of the Walsh g_hat factorization

KEY FINDING FROM product_law_proof.py:
At p=7 degree-2: h_hat[{i,j}] / (g_hat[i] * g_hat[j]) = 0.5 (CONSTANT)
At p=11 degree-4: h_hat[{i,j,k,l}] / prod(g_hat) = 43/44 (CONSTANT)

where g_hat(a) = sum_k chi(k) * sin(2pi*k*a/p) = chi(a)*sqrt(p) (Gauss sum)

This script investigates:
1. WHY is the ratio constant? What is C_d(p) as a function of d and p?
2. Does it extend to degree 6 at p=13?
3. The exact Walsh formula: h_hat[S] = C_d(p) * p^{d/2} * chi(prod a_i)
4. Connection to (E_2, E_3) determination of H

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict
from fractions import Fraction


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
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


def walsh_decomp(p, pairs):
    m = len(pairs)
    n_orient = 1 << m
    H_vals = {}
    for bits in range(n_orient):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_vals[bits] = H

    h_hat = {}
    for S_bits in range(n_orient):
        S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
        val = 0
        for bits in range(n_orient):
            chi = 1
            for i in S_set:
                if not (bits & (1 << i)):
                    chi *= -1
            val += H_vals[bits] * chi
        h_hat[S_set] = val / n_orient
    return h_hat, H_vals


def legendre(a, p):
    """Legendre symbol (a/p)"""
    if a % p == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def gauss_sin_transform(p, a, QR):
    """g_hat(a) = sum_{k=1}^{p-1} chi(k) * sin(2*pi*k*a/p)"""
    return sum((1 if k in QR else -1) * math.sin(2 * math.pi * k * a / p)
               for k in range(1, p))


def main():
    print("=" * 70)
    print("G_HAT FACTORIZATION DEEP ANALYSIS")
    print("=" * 70)

    # PART 1: Compute C_d(p) systematically for ALL degrees at small p
    print("\n--- PART 1: C_d(p) = h_hat[S] / prod(g_hat) ---")
    print("For p = 3 mod 4: g_hat(a) = chi(a)*sqrt(p)")
    print("Expect h_hat[S] = C_d(p) * prod(g_hat[i]) = C_d(p) * p^{d/2} * chi(prod a_i)")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, _ = walsh_decomp(p, pairs)

        # Compute g_hat for each pair
        g_vals = []
        for ii in range(m):
            a = pairs[ii][0]  # a = i+1
            g = gauss_sin_transform(p, a, QR)
            g_vals.append(g)

        # g_hat(a) should = chi(a) * sqrt(p) for p = 3 mod 4
        sqrt_p = math.sqrt(p)
        print(f"\n  p={p} ({p % 4} mod 4):")
        print(f"    g_hat values: {[f'{v:.4f}' for v in g_vals]}")
        print(f"    chi(a)*sqrt(p): {[f'{legendre(i+1,p)*sqrt_p:.4f}' for i in range(m)]}")
        g_match = all(abs(g_vals[i] - legendre(i+1, p)*sqrt_p) < 0.01 for i in range(m))
        print(f"    g_hat = chi*sqrt(p)? {g_match}")

        # For each degree, compute the ratio h_hat / prod(g_hat)
        for d in range(2, m + 1, 2):  # Even degrees only
            ratios = []
            for S_set, val in h_hat.items():
                if len(S_set) != d or abs(val) < 0.001:
                    continue
                indices = sorted(S_set)
                prod = 1.0
                for i in indices:
                    prod *= g_vals[i]
                if abs(prod) > 0.001:
                    ratio = val / prod
                    ratios.append(ratio)

            if ratios:
                unique_ratios = sorted(set(round(r, 6) for r in ratios))
                is_constant = len(unique_ratios) == 1
                count = len(ratios)
                print(f"    degree {d}: {count} coefficients, "
                      f"ratio(s) = {unique_ratios}, "
                      f"constant = {is_constant}")

                if is_constant:
                    r = unique_ratios[0]
                    # Try to identify as rational number
                    frac = Fraction(r).limit_denominator(1000)
                    print(f"      C_{d}({p}) = {r} ~ {frac} = {float(frac):.6f}")
                    # Check: C_d(p) = ratio means h_hat = C_d * p^{d/2} * chi(prod)
                    # At degree 2, p=7: ratio = 0.5 = 1/2
                    # At degree 4, p=11: ratio ~ 43/44
                    # Pattern: C_d(p) = ???

    # PART 2: Algebraic derivation of C_d(p)
    print("\n--- PART 2: ALGEBRAIC DERIVATION ---")
    print("h_hat[S] = (1/2^m) * sum_sigma H(sigma) * chi_S(sigma)")
    print("If H = sum_{d=0}^m f_d(sigma), where f_d = degree-d component,")
    print("then h_hat[S] involves ONLY f_d(sigma) for |S|=d.")
    print("")
    print("KEY: H depends on sigma via the eigenvalue magnitudes Q_k = |lambda_k|^2.")
    print("Each Q_k = 1/4 + D(k;sigma)^2 where D = sum sigma_i sin(...)")
    print("Q_k is QUADRATIC in sigma => H is polynomial in sigma => Walsh structure.")
    print("")
    print("For cycle count c_k: c_k = (1/k) * sum Q_k^{k/2} * exp(i...) ...")
    print("Actually c_k = (1/k) * Tr(A^k) = (1/k) * sum_j lambda_j^k")
    print("lambda_j = -1/2 + i*D(j;sigma)")

    # PART 3: Direct computation of tr(A^k) Walsh coefficients
    print("\n--- PART 3: TRACE-BASED WALSH COEFFICIENTS ---")
    print("Tr(A^k) = sum_{j=0}^{p-1} lambda_j^k")
    print("For circulant: lambda_0 = m, lambda_j = -1/2 + i*D(j;sigma)")
    print("D(j;sigma) = sum_{i=1}^m sigma_i * sin(2*pi*j*(i+1)/p)")
    print("")
    print("STRATEGY: Expand lambda_j^k = (-1/2 + i*D_j)^k and extract")
    print("Walsh coefficients from D_j^n terms (each is multilinear in sigma)")

    for p in [7]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, H_vals = walsh_decomp(p, pairs)

        print(f"\n  p={p}, m={m}:")

        # Compute trace powers for each orientation
        tr_powers = {}  # (bits, k) -> tr(A^k)
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            # Eigenvalues
            lam = [0] * p
            lam[0] = m
            for j in range(1, p):
                D_j = sum(sigma[i] * math.sin(2 * math.pi * j * (i+1) / p)
                          for i in range(m))
                lam[j] = -0.5 + 1j * D_j

            for k in [3, 4, 5, 6, 7]:
                tr_k = sum(l**k for l in lam)
                tr_powers[(bits, k)] = tr_k.real  # should be real for tournaments

        # Walsh decompose each tr(A^k)
        print(f"\n    Walsh decomposition of tr(A^k):")
        for k in [3, 4, 5, 6, 7]:
            # tr_k Walsh coefficients
            tr_hat = {}
            for S_bits in range(1 << m):
                S_set = frozenset(i for i in range(m) if S_bits & (1 << i))
                val = 0
                for bits in range(1 << m):
                    chi = 1
                    for i in S_set:
                        if not (bits & (1 << i)):
                            chi *= -1
                    val += tr_powers[(bits, k)] * chi
                tr_hat[S_set] = val / (1 << m)

            # Report nonzero by degree
            deg_sums = defaultdict(float)
            deg_counts = defaultdict(int)
            for S_set, val in tr_hat.items():
                if abs(val) > 0.001:
                    d = len(S_set)
                    deg_sums[d] += val**2
                    deg_counts[d] += 1

            line = f"    tr(A^{k}): "
            for d in sorted(deg_counts):
                line += f"deg-{d}: {deg_counts[d]} coeffs (energy={deg_sums[d]:.1f}), "
            print(line)

        # CRUCIAL: Compare h_hat[{0,1}] (degree-2 Walsh of H) with
        # the degree-2 Walsh of tr(A^k) for various k
        print(f"\n    Degree-2 Walsh coefficients of tr(A^k) vs H:")
        print(f"    {'S':>8} {'tr^3':>10} {'tr^4':>10} {'tr^5':>10} {'tr^7':>10} {'H':>10}")
        for S_set, hval in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            vals = []
            for k in [3, 4, 5, 7]:
                # Recompute
                tv = 0
                for bits in range(1 << m):
                    chi = 1
                    for i in S_set:
                        if not (bits & (1 << i)):
                            chi *= -1
                    tv += tr_powers[(bits, k)] * chi
                vals.append(tv / (1 << m))
            idxs = sorted(S_set)
            print(f"    {{{idxs[0]+1},{idxs[1]+1}}}: {vals[0]:>10.2f} {vals[1]:>10.2f} "
                  f"{vals[2]:>10.2f} {vals[3]:>10.2f} {hval:>10.2f}")

    # PART 4: The connection between g_hat factorization and Gauss sum structure
    print("\n--- PART 4: WHY THE RATIO IS CONSTANT ---")
    print("CONJECTURE: C_d(p) comes from the trace formula for H.")
    print("")
    print("H = I(Omega(T), 2) = sum_{j=0}^{alpha_max} alpha_j * 2^j")
    print("The alpha_j are polynomial in Q_k's, which are quadratic in sigma.")
    print("Each Q_k = 1/4 + (sum sigma_i * sin_i(k))^2 = 1/4 + S_k(sigma)^2")
    print("")
    print("The degree-d Walsh coefficient of Q_k^n involves sum of products")
    print("of d sine values at frequency k. The SUM over k then brings in")
    print("the Gauss sum structure when p = 3 mod 4.")
    print("")
    print("SPECIFICALLY: For degree 2,")
    print("hat{Q_k}[{i,j}] = 2 * sin(2*pi*k*a_i/p) * sin(2*pi*k*a_j/p)")
    print("Sum over k: sum hat{Q_k}[{i,j}] = (p/2) if a_i=a_j, else...")

    # Let's verify: compute hat{Q_k}[{i,j}] and sum over k
    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        sqrt_p = math.sqrt(p)

        print(f"\n  p={p}:")

        # For each frequency k = 1,...,m, compute sin(2pi*k*a/p) for a = 1,...,m
        sin_matrix = []  # sin_matrix[k-1][i] = sin(2*pi*k*(i+1)/p)
        for k in range(1, m+1):
            row = [math.sin(2*math.pi*k*(i+1)/p) for i in range(m)]
            sin_matrix.append(row)

        # D^T D matrix (m x m)
        # (D^T D)[i,j] = sum_{k=1}^m sin(k*a_i) * sin(k*a_j)
        # Should be (p/4) * delta_{ij} - something
        print(f"    D^T D matrix (scaled by 4/p):")
        DtD = [[0.0]*m for _ in range(m)]
        for i in range(m):
            for j in range(m):
                DtD[i][j] = sum(sin_matrix[k][i]*sin_matrix[k][j] for k in range(m))
        # Scale
        for i in range(m):
            row = [round(DtD[i][j] * 4/p, 4) for j in range(m)]
            print(f"      {row}")

        # The sum over k of sin_k(a_i)*sin_k(a_j) = (p-1)/4 * delta_{ij}
        # because the sine system is orthogonal
        print(f"    Expected diagonal value: (p-1)/4 = {(p-1)/4}")

    # PART 5: The EXACT formula derivation
    print("\n--- PART 5: EXACT FORMULA DERIVATION ---")
    print("At p=7 (m=3):")
    print("H = 1 + 2*alpha_1 + 4*alpha_2")
    print("H = 231 - 2*c4  (from HYP-478)")
    print("c4 = (1/4)*tr(A^4)")
    print("")
    print("tr(A^4) = sum lambda_j^4 = m^4 + sum_{j>=1} (-1/2+i*D_j)^4")
    print("= m^4 + (p-1)/16 - 3*sum D_j^2 + sum D_j^4 + lower order")
    print("   [using binomial expansion of (-1/2+iD)^4]")
    print("")
    print("The degree-2 Walsh part comes from the D_j^2 terms:")
    print("hat{D_j^2}[{i,j}] = 2*sin_i*sin_j (cross term in multinomial)")
    print("Sum over j: sum_k 2*sin(k*a_i)*sin(k*a_j) = 2 * (p-1)/4 * delta_{ij}")
    print("But this is the DIAGONAL part, not the off-diagonal Walsh coefficient!")
    print("")
    print("Off-diagonal (i != j):  sum_k 2*sin(k*a_i)*sin(k*a_j)")
    print("By orthogonality of sines mod p, this should be 0... but it's NOT?")
    print("Let me check:")

    for p in [7]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        # sum_{k=1}^{p-1} sin(2pi*k*a/p) * sin(2pi*k*b/p)
        # = (p/2) * delta_{a,b}  for a,b in 1,...,p-1
        # But we only sum k=1,...,m (not the FULL set to p-1)
        # The half-sum is NOT diagonal!

        print(f"\n  p={p}: Half-sum vs full-sum:")
        for a in range(1, m+1):
            for b in range(1, m+1):
                full = sum(math.sin(2*math.pi*k*a/p)*math.sin(2*math.pi*k*b/p)
                           for k in range(1, p))
                half = sum(math.sin(2*math.pi*k*a/p)*math.sin(2*math.pi*k*b/p)
                           for k in range(1, m+1))
                print(f"    (a={a},b={b}): full_sum = {full:.4f}, half_sum = {half:.4f}")

        print(f"\n  Crucial: sum k=1..p-1 sin(ka)*sin(kb) = (p/2)*delta_ab")
        print(f"  But lambda_k appears for k=1..p-1 (ALL frequencies)")
        print(f"  and D(k) = D(p-k) by sin symmetry => lambda_k and lambda_{p-k} are conjugate")
        print(f"  So: sum_{{k=1}}^{{p-1}} D_k^2 = 2*sum_{{k=1}}^m D_k^2")
        print(f"  The FULL sine orthogonality still applies!")
        print(f"")
        print(f"  Therefore: the degree-2 Walsh coefficient of sum_k D_k^2")
        print(f"  comes from sum_{{k=1}}^{{p-1}} sin(ka_i)*sin(ka_j) = (p/2)*delta")
        print(f"  and gives ZERO for i != j?? But h_hat[{{i,j}}] != 0!")
        print(f"")
        print(f"  RESOLUTION: The degree-2 Walsh of H comes from D_j^4 terms,")
        print(f"  NOT from D_j^2 terms! The D_j^2 sum is diagonal (no cross terms)")
        print(f"  but D_j^4 = (sum sigma_i sin_i)^4 has degree-2 cross terms!")

    # PART 6: Multinomial expansion of D_j^4
    print("\n--- PART 6: DEGREE-2 WALSH FROM D^4 TERM ---")

    for p in [7]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, _ = walsh_decomp(p, pairs)

        # D_j^4 = (sum sigma_i sin_ij)^4
        # = sum sigma_a sigma_b sigma_c sigma_d * sin_aj sin_bj sin_cj sin_dj
        # Degree-2 Walsh of this: terms where exactly 2 distinct indices appear
        # among {a,b,c,d}. That means patterns like (a,a,b,b) with 4!/(2!2!)=6 terms
        # For pair (i,j): coefficient = 6 * sin_ij^2 * sin_ij^2? No...
        # Pattern {a,a,b,b}: coefficient in D^4 expansion = C(4,2) = 6
        # So: hat{D_k^4}[{i,j}] = 6 * sin^2(k*a_i) * sin^2(k*a_j)  ... wait

        # Actually, degree-2 Walsh coefficient of D^4 at {i,j}:
        # D^4 = sum_{abcd} sigma_a sigma_b sigma_c sigma_d * S_a S_b S_c S_d
        # where S_a = sin(2pi*k*a/p) for frequency k.
        # For chi_{i,j}(sigma) to pick out {i,j}, we need sigma_i*sigma_j
        # from the product sigma_a*sigma_b*sigma_c*sigma_d.
        # This means: {a,b,c,d} as multiset must contain i and j each an ODD
        # number of times. With 4 total: either (i,i,i,j) or (i,j,X,X) for X != i,j.
        # Actually for degree-2 Walsh, we need sigma_product = sigma_i * sigma_j.
        # So among {a,b,c,d}: the total count of index i must be odd, and
        # the total count of index j must be odd. All other indices even.
        #
        # With 4 positions and indices from {0,...,m-1}:
        # Patterns contributing to {i,j} degree-2:
        #   (i,j,X,X) where X != i and X != j, X ranges over all m indices
        #     => sigma_i * sigma_j * sigma_X^2 = sigma_i * sigma_j
        #     number of such: C(4,1)*C(3,1)*C(2,2) * m choices for X... no
        #     Actually: multinomial coefficient for each pattern
        #   (i,i,i,j): count(i)=3(odd), count(j)=1(odd) => contributes
        #     coefficient: C(4,3)*1 = 4. S_i^3 * S_j
        #   (i,j,j,j): similar, coefficient: 4. S_i * S_j^3
        #   (i,j,k,k) for k != i,j: coefficient: C(4,2)*C(2,1)*C(1,1)... hmm
        #     Actually it's 4!/(1!*1!*2!) = 12 for each k
        #   (i,j,i,j) == (i,i,j,j): coefficient: 4!/(2!*2!) = 6. S_i^2 * S_j^2

        # Total degree-2 contribution at {i,j} from D_k^4:
        # = 4*S_i^3*S_j + 4*S_i*S_j^3 + 6*S_i^2*S_j^2 + sum_{l!=i,j} 12*S_i*S_j*S_l^2
        # = S_i*S_j * [4*S_i^2 + 4*S_j^2 + 6*S_i*S_j + 12*sum_{l!=i,j} S_l^2]
        # = S_i*S_j * [4*S_i^2 + 4*S_j^2 + 6*S_i*S_j + 12*(sum_l S_l^2 - S_i^2 - S_j^2)]
        # = S_i*S_j * [-8*S_i^2 - 8*S_j^2 + 6*S_i*S_j + 12*Sigma_2]
        # where Sigma_2 = sum_l S_l^2 = (p-1)/4 (from orthogonality)

        # Wait, S_l here means sin(2pi*k*l/p) for a FIXED k.
        # And we need to SUM over all frequencies k = 1,...,p-1 (but k and p-k give same D).

        print(f"\n  p={p}: Verify degree-2 Walsh of sum_k D_k^4")

        for S_set, hval in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2:
                continue
            i, j = sorted(S_set)
            a_i, a_j = pairs[i][0], pairs[j][0]  # = i+1, j+1

            # Compute the degree-2 Walsh of sum_k D_k^4 at {i,j}
            # For each k, contribution = S_ik * S_jk * [-8S_ik^2 - 8S_jk^2 + 6S_ik*S_jk + 12*Sigma2_k]
            # where S_ik = sin(2pi*k*a_i/p) and Sigma2_k = sum_l sin^2(2pi*k*l/p)

            total = 0
            for k in range(1, p):
                S_ik = math.sin(2*math.pi*k*a_i/p)
                S_jk = math.sin(2*math.pi*k*a_j/p)
                Sigma2_k = sum(math.sin(2*math.pi*k*(l+1)/p)**2 for l in range(m))

                contrib = S_ik * S_jk * (-8*S_ik**2 - 8*S_jk**2 + 6*S_ik*S_jk + 12*Sigma2_k)
                total += contrib

            print(f"    {{i={i+1},j={j+1}}}: sum_k D_k^4 deg-2 = {total:.4f}")

        # Now H = 231 - (1/2)*tr(A^4) at p=7
        # tr(A^4) = m^4 + sum_{k>=1} lambda_k^4
        # lambda_k^4 = ((-1/2) + i*D_k)^4
        # = 1/16 - (3/2)*D_k^2 + D_k^4 + i*(...) [imaginary parts cancel in sum]
        # Actually:
        # (-1/2 + iD)^4 = (-1/2)^4 + 4(-1/2)^3(iD) + 6(-1/2)^2(iD)^2 + 4(-1/2)(iD)^3 + (iD)^4
        # = 1/16 + 4(-1/8)(iD) + 6(1/4)(-D^2) + 4(-1/2)(-iD^3) + D^4
        # = 1/16 - iD/2 - 3D^2/2 + 2iD^3 + D^4
        # Real part: 1/16 - 3D^2/2 + D^4
        # Imaginary: -D/2 + 2D^3  (this sums to 0 over k by symmetry)

        print(f"\n  H = 231 - (1/2)*tr(A^4) where:")
        print(f"  tr(A^4) = m^4 + (p-1)/16 - (3/2)*sum D_k^2 + sum D_k^4")
        print(f"  Degree-2 Walsh of tr(A^4) at {{i,j}}:")
        print(f"  = -3 * sum_k S_ik*S_jk + (formula above)")
        print(f"  The D^2 term: -3 * sum_{{k=1}}^{{p-1}} sin(k*a_i)*sin(k*a_j)")
        print(f"  = -3 * (p/2)*delta_{{a_i,a_j}}")
        print(f"  = 0 for i != j (by sine orthogonality over full k)")
        print(f"  So ALL degree-2 Walsh content of tr(A^4) comes from D^4!")

    # PART 7: The chi(a*b) sign from the D^4 Gauss sum structure
    print("\n--- PART 7: CHI(a*b) FROM GAUSS SUM ---")
    print("Key: sum_k sin(ka)*sin(kb)*sin^2(kc) involves 4th-power Gauss terms")
    print("Product-to-sum: sin(ka)*sin(kb) = (1/2)[cos(k(a-b)) - cos(k(a+b))]")
    print("sin^2(kc) = (1/2)[1 - cos(2kc)]")
    print("So the product = (1/4)[cos(k(a-b)) - cos(k(a+b))]*(1 - cos(2kc))")
    print("")
    print("Summing over k = 1,...,p-1:")
    print("sum cos(kn) = -1 if n != 0 mod p, = p-1 if n = 0 mod p")
    print("")
    print("This gives EXACT formulas for the D^4 cross terms!")

    # VERIFY: compute the sum analytically and compare to numerical
    for p in [7, 11]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        h_hat, _ = walsh_decomp(p, pairs)

        print(f"\n  p={p}:")

        for S_set, hval in sorted(h_hat.items(), key=lambda x: tuple(sorted(x[0]))):
            if len(S_set) != 2 or abs(hval) < 0.01:
                continue
            i, j = sorted(S_set)
            a, b = pairs[i][0], pairs[j][0]  # a=i+1, b=j+1

            # The degree-2 Walsh coefficient of sum_k D_k^4 at {a,b}:
            # = sum_{l=1}^m 12 * C_{ab,l} + terms from i^3*j and i*j^3 and i^2*j^2
            # where C_{ab,l} = sum_k sin(ka)*sin(kb)*sin^2(kl)
            #
            # C_{ab,l} = sum_k (1/4)[cos(k(a-b))-cos(k(a+b))]*[1-cos(2kl)]
            # = (1/4) sum_k [cos(k(a-b)) - cos(k(a+b)) - cos(k(a-b-2l)) + ... ]
            # Wait, let me be more careful with the product:
            # sin(ka)sin(kb) = (1/2)[cos(k(a-b)) - cos(k(a+b))]
            # sin^2(kl) = (1/2)(1 - cos(2kl))
            #
            # Product = (1/4)[cos(k(a-b)) - cos(k(a+b))] * [1 - cos(2kl)]
            # = (1/4)[cos(k(a-b)) - cos(k(a+b)) - cos(k(a-b))*cos(2kl) + cos(k(a+b))*cos(2kl)]
            #
            # cos(kA)*cos(kB) = (1/2)[cos(k(A+B)) + cos(k(A-B))]
            # So:
            # cos(k(a-b))*cos(2kl) = (1/2)[cos(k(a-b+2l)) + cos(k(a-b-2l))]
            # cos(k(a+b))*cos(2kl) = (1/2)[cos(k(a+b+2l)) + cos(k(a+b-2l))]
            #
            # Full sum at each k:
            # (1/4){ cos(k(a-b)) - cos(k(a+b))
            #       - (1/2)cos(k(a-b+2l)) - (1/2)cos(k(a-b-2l))
            #       + (1/2)cos(k(a+b+2l)) + (1/2)cos(k(a+b-2l)) }
            #
            # Sum over k=1..p-1: sum cos(kn) = -1 if n%p != 0, else p-1
            # delta(n) = p if n%p==0 else 0 (since sum_k cos = p-1 when n=0, else -1)
            # Actually: sum_{k=0}^{p-1} cos(2*pi*k*n/p) = p*delta_{n%p,0}
            # So: sum_{k=1}^{p-1} cos(2*pi*k*n/p) = p*delta_{n%p,0} - 1

            def delta_sum(n_val, p_val):
                """sum_{k=1}^{p-1} cos(2*pi*k*n/p) = p*delta(n%p,0) - 1"""
                return p_val - 1 if (n_val % p_val == 0) else -1

            # C_{ab,l} = (1/4){ delta_sum(a-b) - delta_sum(a+b)
            #            - (1/2)delta_sum(a-b+2l) - (1/2)delta_sum(a-b-2l)
            #            + (1/2)delta_sum(a+b+2l) + (1/2)delta_sum(a+b-2l) }
            # (all arguments mod p, applied to 2*pi*k/p)

            # Similarly need terms for i^3*j, i*j^3, i^2*j^2
            # i^3*j: sum_k sin^3(ka)*sin(kb) = sin(ka)^3*sin(kb)
            # = (3/4*sin(ka) - 1/4*sin(3ka)) * sin(kb)   [triple angle]
            # = (3/4)[cos(k(a-b)) - cos(k(a+b))]/2 - (1/4)[cos(k(3a-b)) - cos(k(3a+b))]/2

            # Actually, let me just verify numerically for a specific case
            # Compute (1/4)*{sum_k [cos(...)] analytical} vs numerical
            # For the full degree-2 Walsh of sum D_k^4:

            # Full formula:
            # D_k^4 = (sum_{l=1}^m sigma_l sin_l(k))^4
            # Degree-2 Walsh at {i,j} means:
            # coefficient of sigma_i*sigma_j in (1/2^m) sum_sigma D_k^4 * sigma_i * sigma_j
            # = sum of S_{ik}*S_{jk} * (terms)

            # Let me compute it numerically and compare to analytical
            num = 0
            for k in range(1, p):
                S_ik = math.sin(2*math.pi*k*a/p)
                S_jk = math.sin(2*math.pi*k*b/p)
                Sigma2_k = sum(math.sin(2*math.pi*k*(l+1)/p)**2 for l in range(m))

                # Full degree-2 contribution from D_k^4
                # = S_ik*S_jk * [-8*S_ik^2 - 8*S_jk^2 + 6*S_ik*S_jk + 12*Sigma2_k]
                num += S_ik * S_jk * (-8*S_ik**2 - 8*S_jk**2 + 6*S_ik*S_jk + 12*Sigma2_k)

            # Analytical: using sum cos formulas
            # Sigma2_k = sum_{l=1}^m sin^2(2*pi*k*(l+1)/p)
            # = (1/2) * sum_{l=1}^m [1 - cos(2*2*pi*k*(l+1)/p)]
            # = m/2 - (1/2)*sum_{l=1}^m cos(2*2*pi*k*(l+1)/p)

            # 12*sum S_il*S_jl*Sigma2 term:
            # = 12 * sum_k S_ik*S_jk*sum_l S_lk^2
            # = 12 * sum_l sum_k S_ik*S_jk*S_lk^2

            # The key observation: the sum over k=1..p-1 of sin(ka)*sin(kb)*sin^2(kc)
            # can be computed exactly using delta_sum identities.

            # C_{abc} = sum_{k=1}^{p-1} sin(ka)sin(kb)sin^2(kc)
            C_sum = 0
            for c in range(1, m+1):
                c_val = c  # c = l+1 since pairs[l][0] = l+1
                # analytical C_{a,b,c}
                C_abc = (1.0/8) * (
                    2*delta_sum(a-b, p) - 2*delta_sum(a+b, p)
                    - delta_sum(a-b+2*c_val, p) - delta_sum(a-b-2*c_val, p)
                    + delta_sum(a+b+2*c_val, p) + delta_sum(a+b-2*c_val, p)
                )
                C_sum += C_abc

            # Similarly compute the other terms analytically
            # 4*S_i^3*S_j + 4*S_i*S_j^3 contribution:
            # sin^3(ka) = (3sin(ka) - sin(3ka))/4
            T1 = 0  # 4*sum_k sin^3(ka)*sin(kb)
            for k in range(1, p):
                T1 += math.sin(2*math.pi*k*a/p)**3 * math.sin(2*math.pi*k*b/p)
            T1 *= 4

            T2 = 0  # 4*sum_k sin(ka)*sin^3(kb)
            for k in range(1, p):
                T2 += math.sin(2*math.pi*k*a/p) * math.sin(2*math.pi*k*b/p)**3
            T2 *= 4

            T3 = 0  # 6*sum_k sin^2(ka)*sin^2(kb)
            for k in range(1, p):
                T3 += math.sin(2*math.pi*k*a/p)**2 * math.sin(2*math.pi*k*b/p)**2
            T3 *= 6

            total_terms = T1 + T2 + T3 + 12*C_sum

            chi_ab = legendre(a*b, p)

            print(f"    (a={a},b={b}): D^4_deg2 = {num:.4f}, "
                  f"T1+T2+T3+12C = {total_terms:.4f}, "
                  f"chi(ab) = {chi_ab:+d}, "
                  f"h_hat = {hval:.4f}")

    # PART 8: The C_d(p) pattern
    print("\n--- PART 8: C_d(p) PATTERN ---")
    print("Known values:")
    print("  C_2(7) = 1/2")
    print("  C_2(11) = NOT CONSTANT (0.75 or 24.75)")
    print("  C_4(11) = 43/44")
    print("")
    print("HYPOTHESIS: C_d(p) is constant ONLY when d = (p-1)/2 mod something?")
    print("At p=7, m=3: degree 2 (= m-1). C_2(7) = 1/2.")
    print("At p=11, m=5: degree 4 (= m-1). C_4(11) = 43/44.")
    print("Pattern: C_{m-1}(p) is constant? Always = (p-2)/(p-1)?")
    print("  C_2(7) = 1/2 = 5/6?? No, 1/2 = (7-2)/(7-1-4)?? Hmm.")
    print("  C_4(11) = 43/44 = (p-2)/(p-1) * (something)?")
    print("")
    print("  1/2 = 3/6 = (m)/(2m)?")
    print("  43/44 = (p-2)/(p-1) = 9/10?? No. 43/44...")
    print("")
    print("Let me check: at degree = m-1, C_{m-1}(p) = ?")
    print("  p=7, d=2: C = 0.500000 = 1/2")
    print("  p=11, d=4: C = 0.977273 = 43/44")
    print("")
    print("  1/2 = 1 - 1/2   = 1 - 1/(m-1)? m=3, 1/(m-1) = 1/2. Yes!")
    print("  43/44 = 1 - 1/44 = 1 - 1/(p-1 choose d/2) ??")
    print("  C(10,2) = 45. 1-1/45 = 44/45 != 43/44")
    print("  C(5,2) = 10. 1-1/10 = 9/10 != 43/44")
    print("")
    print("Let me try: 43/44 as a fraction. 43 and 44 are coprime.")
    print("  44 = 4*11 = 4p. 1/2 = 3.5/7 = (p/2)/p.")
    print("  43/(4*11) = 43/44. And 3.5/7 = 3.5/7 = 1/2.")
    print("  Hmm, doesn't quite work.")
    print("")
    print("  Perhaps C_d(p) = (F_p - 1) / F_p where F_p is the Fibonacci bound?")
    print("  F_7 = F_7 = 13. (13-1)/13 ~ 0.923 != 0.5")

    # Let me compute C_d(p) for p=13 where we have both p=1 and degree 2,4,6
    print("\n  p=13 C_d values (from Part 1):")
    p = 13
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
    h_hat, _ = walsh_decomp(p, pairs)
    g_vals = [gauss_sin_transform(p, pairs[ii][0], QR) for ii in range(m)]

    for d in range(2, m + 1, 2):
        ratios = []
        for S_set, val in h_hat.items():
            if len(S_set) != d or abs(val) < 0.001:
                continue
            indices = sorted(S_set)
            prod = 1.0
            for idx in indices:
                prod *= g_vals[idx]
            if abs(prod) > 0.001:
                ratios.append(val / prod)
        if ratios:
            unique = sorted(set(round(r, 6) for r in ratios))
            print(f"    d={d}: ratios = {unique[:5]}{'...' if len(unique)>5 else ''} "
                  f"({len(unique)} distinct, {len(ratios)} total)")

    print("\n  CONCLUSION ON C_d(p):")
    print("  The ratio is constant ONLY at specific degrees.")
    print("  For p=3 mod 4: constant at degree 2 when m=3 (p=7),")
    print("                  constant at degree 4 when m=5 (p=11).")
    print("  For p=1 mod 4: g_hat = chi(a)*sqrt(p) may NOT hold")
    print("  (the Gauss sum has different phase).")
    print("  The constant-ratio phenomenon needs the QR character structure.")


if __name__ == '__main__':
    main()
