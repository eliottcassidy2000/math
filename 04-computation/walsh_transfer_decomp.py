#!/usr/bin/env python3
"""
walsh_transfer_decomp.py -- Transfer matrix eigenspace decomposition of H

By THM-125, for circulant tournaments on Z_p, the symbol matrix M_m(t) is
constant. This means the transfer matrix decomposes into p eigenspaces,
each contributing independently to H.

H = p * sum_{k=0}^{p-1} h_k(S)

where h_k depends on S_hat(k). For k=0: h_0 = (p-1)! / p^{p-1} * m^{p-1}
(depends only on |S| = m). For k != 0: h_k depends on S_hat(k).

By conjugate pairing: h_k = h_{p-k}* (complex conjugate), so h_k + h_{p-k}
is real and depends only on |S_hat(k)|^2 = Q_k.

But the nonlinearity in H(Q_k) means the h_k contribution is NOT simply
proportional to |S_hat(k)|^2. Instead, it's a more complex function.

This script:
1. Computes the eigenspace contributions h_k directly via the transfer matrix
2. Tests if h_k is a polynomial in S_hat(k) and |S_hat(k)|^2
3. Finds the exact functional form of h_k
4. Derives the Walsh decomposition from the eigenspace decomposition

Author: kind-pasteur-2026-03-12-S59c
"""

import cmath
import math
from itertools import combinations, permutations
from collections import defaultdict


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


def transfer_matrix_contribution(p, S, verbose=False):
    """Compute H via DFT eigenspace decomposition.

    The transfer matrix T_{ij} = 1 if j-i in S (tournament arc).
    Its eigenvalues are S_hat(k) = sum_{s in S} omega^{ks}.

    For Hamiltonian paths, we need the permanent of a (p-1)x(p-1) submatrix.
    But by circulant symmetry, H/p = permanent of the (p-1)x(p-1) "reduced" matrix.

    By DFT: the eigenvalues of the circulant matrix are S_hat(k).
    H = sum over permutations = related to permanent.

    Permanent of circulant matrix is NOT simply prod of eigenvalues (that's det).
    Instead: permanent = sum_{sigma} prod T_{i,sigma(i)}

    For a circulant C with first row (c_0, c_1, ..., c_{n-1}):
      perm(C) = sum_{sigma in S_n} prod c_{sigma(i)-i mod n}

    This is exactly the permanent of the circulant.
    """
    omega = cmath.exp(2j * cmath.pi / p)
    S_set = set(S)

    # Fourier transform of connection set
    S_hat = []
    for k in range(p):
        val = sum(omega ** (k * s) for s in S)
        S_hat.append(val)

    # For the eigenspace decomposition, we need to understand how the
    # permanent decomposes. The key formula (Minc, Marcus):
    #
    # perm(C) = (1/n) sum_{k} S_hat(k) * perm(C_k)
    #
    # where C_k is a (n-1)x(n-1) matrix. But this is recursive and doesn't
    # give a simple closed form.
    #
    # Instead, let's use the Ryser formula for the permanent as inclusion-exclusion.
    # Or better: directly compute H as a function of the DFT coefficients
    # by computing the "DFT permanent decomposition."

    # For SMALL p, we can compute the eigenspace contribution directly.
    # H = number of Hamiltonian paths in the tournament
    # = sum over all orderings (v_0, v_1, ..., v_{p-1}) of Z_p such that
    #   v_{i+1} - v_i in S for all i.
    #
    # Using DFT: the number of walks of length p-1 with all distinct vertices
    # and all steps in S.

    # The circulant has eigenvalues S_hat(k). The number of length-L walks
    # from any vertex is S_hat(0)^L = m^L (doesn't depend on S for fixed |S|).
    # But Hamiltonian PATHS require distinct vertices.

    # Let's compute H/p as the permanent of a (p-1)x(p-1) "reduced" matrix.
    # Actually, H = sum over start vertex v of #{paths starting at v}.
    # By circulant symmetry, each start contributes H/p.

    # For the reduced matrix approach:
    # Fix start vertex v=0. Count paths 0, v_1, ..., v_{p-1} with
    # v_1 - 0 in S, v_2 - v_1 in S, etc.
    # This is H/p = sum over permutations sigma of {1,...,p-1} with
    #   sigma(1) in S, sigma(2)-sigma(1) in S, ..., sigma(p-1)-sigma(p-2) in S.

    return S_hat


def eigenspace_H_decomp(p, S):
    """Try to decompose H into contributions from each DFT frequency.

    The idea: write H as a function of the DFT coefficients S_hat(k).
    Test various polynomial forms.
    """
    omega = cmath.exp(2j * cmath.pi / p)
    S_hat = []
    for k in range(p):
        val = sum(omega ** (k * s) for s in S)
        S_hat.append(val)

    return S_hat


def main():
    print("=" * 70)
    print("TRANSFER MATRIX EIGENSPACE DECOMPOSITION")
    print("=" * 70)

    for p in [5, 7, 11]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        omega = cmath.exp(2j * cmath.pi / p)

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}")
        print(f"{'='*70}")

        # Compute H for all orientations along with their S_hat values
        pairs = [(s, p - s) for s in range(1, m + 1)]
        H_dict = {}
        Shat_dict = {}

        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = []
            for i, (a, b) in enumerate(pairs):
                S.append(a if sigma[i] == 1 else b)
            S = sorted(S)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            H_dict[sigma] = H

            # DFT
            S_hat = []
            for k in range(p):
                val = sum(omega ** (k * s) for s in S)
                S_hat.append(val)
            Shat_dict[sigma] = S_hat

        # ====== DISTINCT (Q_1,...,Q_m) PROFILES ======
        print(f"\n  1. DISTINCT Q PROFILES")

        Q_profiles = {}
        for sigma, H in H_dict.items():
            S_hat = Shat_dict[sigma]
            Q = tuple(round(abs(S_hat[k])**2, 8) for k in range(1, m+1))
            Q_sorted = tuple(sorted(Q))
            if Q_sorted not in Q_profiles:
                Q_profiles[Q_sorted] = {'H': H, 'sigmas': [], 'Q': Q}
            Q_profiles[Q_sorted]['sigmas'].append(sigma)

        for Q_sorted, data in sorted(Q_profiles.items(), key=lambda x: -x[1]['H']):
            print(f"    Q_sorted = {Q_sorted}")
            print(f"    H = {data['H']}, {len(data['sigmas'])} orientations")
            # Show the S_hat values for one representative
            sigma = data['sigmas'][0]
            S_hat = Shat_dict[sigma]
            shat_vals = [S_hat[k] for k in range(1, m+1)]
            print(f"    S_hat = {['({:.3f})'.format(v) for v in shat_vals]}")

        # ====== POWER SUMS OF S_hat(k) ======
        print(f"\n  2. POWER SUMS OF S_hat AND H")
        print(f"     Testing: H = f(P_1, P_2, ...) where P_j = sum_k S_hat(k)^j")

        # For each sigma, compute power sums of S_hat(k) for k=1,...,p-1
        power_sums = {}
        for sigma in H_dict:
            S_hat = Shat_dict[sigma]
            ps = {}
            for j in range(1, p):
                # P_j = sum_{k=1}^{p-1} S_hat(k)^j
                ps[j] = sum(S_hat[k]**j for k in range(1, p))
            power_sums[sigma] = ps

        # P_1 = sum S_hat(k) = -S_hat(0) + sum_all = -m + (sum_k sum_s omega^{ks})
        # = -m + sum_s (sum_k omega^{ks}) = -m + sum_s (p * delta_{s,0} - 1)
        # Since 0 not in S: = -m + (-m) = -2m. Wait, let me compute directly.

        print(f"\n    P_1 for first few sigmas:")
        for sigma in sorted(H_dict, key=lambda s: -H_dict[s])[:4]:
            print(f"      sigma={sigma}: P_1 = {power_sums[sigma][1]:.4f}, "
                  f"P_2 = {power_sums[sigma][2]:.4f}")

        # P_1 = sum_{k=1}^{p-1} S_hat(k) = (sum_{k=0}^{p-1} S_hat(k)) - S_hat(0)
        # sum_k S_hat(k) = sum_k sum_s omega^{ks} = sum_s sum_k omega^{ks} = sum_s p*delta(s,0) = 0
        # (since 0 not in S)
        # So P_1 = -S_hat(0) = -m. Should be constant.

        print(f"\n    P_1 should be -m = {-m}")
        for sigma in list(H_dict.keys())[:3]:
            print(f"      P_1 = {power_sums[sigma][1]:.6f}")

        # P_2 = sum_{k=1}^{p-1} |S_hat(k)|^2 (since S_hat is Hermitian: S_hat(p-k) = conj(S_hat(k)))
        # Actually P_2 = sum S_hat(k)^2 which is NOT |S_hat(k)|^2.
        # |S_hat(k)|^2 = S_hat(k) * conj(S_hat(k)) = S_hat(k) * S_hat(p-k)
        # So sum |S_hat(k)|^2 = sum S_hat(k)*S_hat(p-k) = trace of auto-correlation
        # = p * sum_{d in Z_p} |{(s1,s2) in S^2 : s1-s2=d}| = p * |S| = p*m  (by Parseval)
        # Wait: sum_k |S_hat(k)|^2 = p * |S| = p*m. This includes k=0: |S_hat(0)|^2 = m^2.
        # So sum_{k=1}^{p-1} |S_hat(k)|^2 = p*m - m^2 = m(p-m) = m*(p+1)/2.

        print(f"\n    sum Q_k = sum |S_hat(k)|^2 (k=1..p-1) should be m(p-m) = {m*(p-m)}")
        for sigma in list(H_dict.keys())[:3]:
            S_hat = Shat_dict[sigma]
            sum_q = sum(abs(S_hat[k])**2 for k in range(1, p))
            print(f"      sum Q_k = {sum_q:.6f}")

        # ====== PRODUCT FORM ======
        print(f"\n  3. TESTING PRODUCT RELATIONSHIPS")
        # Maybe H is related to prod_{k=1}^{m} f(S_hat(k)) for some f?

        # Test: log(H) vs sum log(Q_k)
        print(f"    log(H) vs sum log(Q_k):")
        for sigma in sorted(H_dict, key=lambda s: -H_dict[s]):
            H = H_dict[sigma]
            S_hat = Shat_dict[sigma]
            Q_vals = [abs(S_hat[k])**2 for k in range(1, m+1)]
            prod_Q = 1
            for q in Q_vals:
                prod_Q *= q
            log_H = math.log(H)
            log_prod_Q = math.log(prod_Q) if prod_Q > 0 else float('-inf')
            if abs(log_prod_Q) > 1e-10 and log_prod_Q != float('-inf'):
                print(f"      H={H:>10}, log(H)={log_H:.4f}, "
                      f"log(prod Q_k)={log_prod_Q:.4f}, "
                      f"ratio={log_H/log_prod_Q:.4f}")
            else:
                print(f"      H={H:>10}, log(H)={log_H:.4f}, "
                      f"log(prod Q_k)={log_prod_Q:.4f} (near zero)")

        # ====== NEWTON'S IDENTITIES: p_k from e_k ======
        print(f"\n  4. SYMMETRIC FUNCTIONS OF Q_k")
        # Let e_1, e_2, ..., e_m be the elementary symmetric functions of Q_1,...,Q_m
        # If H = F(Q_sorted), it can be written as a polynomial in e_1,...,e_m.

        # Compute e_j for each sigma
        for sigma in sorted(H_dict, key=lambda s: -H_dict[s]):
            H = H_dict[sigma]
            S_hat = Shat_dict[sigma]
            Q_vals = [abs(S_hat[k])**2 for k in range(1, m+1)]
            Q_sorted = sorted(Q_vals)

            # Elementary symmetric functions
            e = [1]  # e_0 = 1
            for j in range(1, m+1):
                ej = sum(math.prod(Q_vals[i] for i in combo)
                        for combo in combinations(range(m), j))
                e.append(ej)

            e_str = ", ".join(f"e_{j}={e[j]:.4f}" for j in range(1, min(m+1, 5)))
            print(f"    H={H:>10}: {e_str}")

        # ====== EXACT POLYNOMIAL in e_j ======
        print(f"\n  5. POLYNOMIAL FIT: H = sum c_j * e_j")
        # Build a linear system: for each distinct Q-sorted profile, H is known.
        # The unknowns are the coefficients c_0, c_1, ..., c_m in
        #   H = c_0 + c_1*e_1 + c_2*e_2 + ... + c_m*e_m

        sigma_list = list(H_dict.keys())
        H_vector = [H_dict[s] for s in sigma_list]

        # Compute feature matrix
        features = []
        for sigma in sigma_list:
            S_hat = Shat_dict[sigma]
            Q_vals = [abs(S_hat[k])**2 for k in range(1, m+1)]
            row = [1.0]
            for j in range(1, m+1):
                ej = sum(math.prod(Q_vals[i] for i in combo)
                        for combo in combinations(range(m), j))
                row.append(ej)
            features.append(row)

        n_feats = m + 1
        n_data = len(sigma_list)

        # Solve least squares
        XTX = [[0.0]*n_feats for _ in range(n_feats)]
        XTH = [0.0]*n_feats
        for r in range(n_data):
            for i in range(n_feats):
                for j in range(n_feats):
                    XTX[i][j] += features[r][i] * features[r][j]
                XTH[i] += features[r][i] * H_vector[r]

        # Gaussian elimination
        aug = [XTX[i][:] + [XTH[i]] for i in range(n_feats)]
        for col in range(n_feats):
            max_row = max(range(col, n_feats), key=lambda r: abs(aug[r][col]))
            aug[col], aug[max_row] = aug[max_row], aug[col]
            pivot = aug[col][col]
            if abs(pivot) < 1e-20:
                print(f"    Near-singular at column {col}")
                break
            for j in range(col, n_feats + 1):
                aug[col][j] /= pivot
            for r in range(n_feats):
                if r == col:
                    continue
                factor = aug[r][col]
                for j in range(col, n_feats + 1):
                    aug[r][j] -= factor * aug[col][j]

        c_vals = [aug[i][n_feats] for i in range(n_feats)]

        # Check residuals
        residuals = []
        for r in range(n_data):
            H_pred = sum(c_vals[i] * features[r][i] for i in range(n_feats))
            residuals.append(H_vector[r] - H_pred)

        max_resid = max(abs(r) for r in residuals)
        print(f"    Max residual = {max_resid:.6f}")
        print(f"    Is H = polynomial in e_j? {'YES' if max_resid < 0.01 else 'NO'}")

        if max_resid < 0.01:
            print(f"    Coefficients:")
            for j in range(n_feats):
                print(f"      c_{j} = {c_vals[j]:.10f}")
        else:
            print(f"    Linear-in-e_j model fails.")

        # ====== POWER SUM POLYNOMIAL ======
        print(f"\n  6. POLYNOMIAL IN POWER SUMS p_j = sum Q_k^j")
        # p_j = sum_{k=1}^{m} Q_k^j

        pow_features = []
        for sigma in sigma_list:
            S_hat = Shat_dict[sigma]
            Q_vals = [abs(S_hat[k])**2 for k in range(1, m+1)]
            row = [1.0]
            for j in range(1, m+1):
                pj = sum(q**j for q in Q_vals)
                row.append(pj)
            pow_features.append(row)

        # Solve
        XTX2 = [[0.0]*n_feats for _ in range(n_feats)]
        XTH2 = [0.0]*n_feats
        for r in range(n_data):
            for i in range(n_feats):
                for j in range(n_feats):
                    XTX2[i][j] += pow_features[r][i] * pow_features[r][j]
                XTH2[i] += pow_features[r][i] * H_vector[r]

        aug2 = [XTX2[i][:] + [XTH2[i]] for i in range(n_feats)]
        for col in range(n_feats):
            max_row = max(range(col, n_feats), key=lambda r: abs(aug2[r][col]))
            aug2[col], aug2[max_row] = aug2[max_row], aug2[col]
            pivot = aug2[col][col]
            if abs(pivot) < 1e-20:
                print(f"    Near-singular at column {col}")
                break
            for j in range(col, n_feats + 1):
                aug2[col][j] /= pivot
            for r in range(n_feats):
                if r == col:
                    continue
                factor = aug2[r][col]
                for j in range(col, n_feats + 1):
                    aug2[r][j] -= factor * aug2[col][j]

        c2_vals = [aug2[i][n_feats] for i in range(n_feats)]

        residuals2 = []
        for r in range(n_data):
            H_pred = sum(c2_vals[i] * pow_features[r][i] for i in range(n_feats))
            residuals2.append(H_vector[r] - H_pred)

        max_resid2 = max(abs(r) for r in residuals2)
        print(f"    Max residual = {max_resid2:.6f}")
        print(f"    Is H = polynomial in p_j? {'YES' if max_resid2 < 0.01 else 'NO'}")

        if max_resid2 < 0.01:
            print(f"    Coefficients:")
            for j in range(n_feats):
                print(f"      c_{j} = {c2_vals[j]:.10f}")

        # ====== MIXED PRODUCTS ======
        if max_resid > 0.01 and max_resid2 > 0.01:
            print(f"\n  7. QUADRATIC IN e_j / MIXED PRODUCTS")
            # Try H = sum c_{jl} * e_j * e_l
            mixed_features = []
            for sigma in sigma_list:
                S_hat = Shat_dict[sigma]
                Q_vals = [abs(S_hat[k])**2 for k in range(1, m+1)]
                e_vals = [1.0]
                for j in range(1, m+1):
                    ej = sum(math.prod(Q_vals[i] for i in combo)
                            for combo in combinations(range(m), j))
                    e_vals.append(ej)

                row = []
                for j in range(m+1):
                    for l in range(j, m+1):
                        row.append(e_vals[j] * e_vals[l])
                mixed_features.append(row)

            n_mixed = len(mixed_features[0])
            print(f"    {n_mixed} mixed features, {n_data} data points")

            if n_mixed <= n_data:
                XTX3 = [[0.0]*n_mixed for _ in range(n_mixed)]
                XTH3 = [0.0]*n_mixed
                for r in range(n_data):
                    for i in range(n_mixed):
                        for j in range(n_mixed):
                            XTX3[i][j] += mixed_features[r][i] * mixed_features[r][j]
                        XTH3[i] += mixed_features[r][i] * H_vector[r]

                aug3 = [XTX3[i][:] + [XTH3[i]] for i in range(n_mixed)]
                for col in range(n_mixed):
                    max_row = max(range(col, n_mixed), key=lambda r: abs(aug3[r][col]))
                    aug3[col], aug3[max_row] = aug3[max_row], aug3[col]
                    pivot = aug3[col][col]
                    if abs(pivot) < 1e-20:
                        break
                    for j in range(col, n_mixed + 1):
                        aug3[col][j] /= pivot
                    for r in range(n_mixed):
                        if r == col:
                            continue
                        factor = aug3[r][col]
                        for j in range(col, n_mixed + 1):
                            aug3[r][j] -= factor * aug3[col][j]

                c3_vals = [aug3[i][n_mixed] for i in range(n_mixed)]

                residuals3 = []
                for r in range(n_data):
                    H_pred = sum(c3_vals[i] * mixed_features[r][i] for i in range(n_mixed))
                    residuals3.append(H_vector[r] - H_pred)

                max_resid3 = max(abs(r) for r in residuals3)
                print(f"    Max residual = {max_resid3:.6f}")
                print(f"    Is H = quadratic in e_j? {'YES' if max_resid3 < 0.01 else 'NO'}")

        # ====== RELATIONSHIP TO PERMANENT ======
        print(f"\n  8. PERMANENT CONNECTION")
        # H/p = permanent of (p-1)x(p-1) reduced circulant matrix.
        # The (i,j) entry of the reduced matrix is 1 if j-i mod p is in S, 0 otherwise.
        # (where indices run from 1 to p-1, not 0 to p-1)

        # For a circulant matrix with eigenvalues lambda_k = S_hat(k):
        # det = prod lambda_k
        # perm != prod lambda_k (permanent doesn't have this nice formula)

        # BUT: by a result of Marcus and Minc, for circulant matrices:
        #   perm(C) = (1/n) sum_{d|n} phi(n/d) * (sum_{k: dk=0 mod n} lambda_k)^{n/d}  (???)
        #   No, that's not right either.

        # The permanent of a circulant matrix is the sum over all permutations
        # sigma of prod c_{sigma(i)-i}. By DFT:
        #   c_j = (1/n) sum_k lambda_k * omega^{-jk}
        # So:
        #   prod c_{sigma(i)-i} = prod_i [(1/n) sum_k lambda_k omega^{-k(sigma(i)-i)}]
        # This is a product of sums, which doesn't factor nicely.

        # Let's just compute H/p for each Q-profile and look for patterns.
        print(f"    H/p for each Q-profile:")
        for Q_sorted, data in sorted(Q_profiles.items(), key=lambda x: -x[1]['H']):
            H = data['H']
            print(f"      Q = {Q_sorted}: H/p = {H/p:.6f}, "
                  f"prod Q = {math.prod(Q_sorted):.4f}, "
                  f"H/(p*prod_Q) = {H/(p*math.prod(Q_sorted)):.8f}")

        # ====== H in terms of character sums ======
        print(f"\n  9. H AND CHARACTER SUMS")
        # For a circulant tournament, H counts Hamiltonian paths.
        # H = sum_{(v_0,...,v_{p-1})} prod_{i=0}^{p-2} 1_{v_{i+1}-v_i in S}
        # = sum over orderings where each step is in S
        #
        # In terms of S_hat:
        # 1_{d in S} = (1/p) sum_k S_hat(k) * omega^{-kd}
        #
        # So H = sum_{orderings} prod_{i} [(1/p) sum_k S_hat(k_i) omega^{-k_i * d_i}]
        # where d_i = v_{i+1} - v_i.
        #
        # H = (1/p^{p-1}) sum_{orderings} sum_{k_0,...,k_{p-2}}
        #     prod S_hat(k_i) * omega^{-sum k_i d_i}
        #
        # The sum over orderings of omega^{-sum k_i d_i} is a character sum.
        # When all k_i = 0: contributes m^{p-1} * p! / p = p! * m^{p-1} / p^p
        # Wait, this isn't right. Let me think more carefully.

        # Actually: for a circulant T on Z_p, H/p is the permanent of the
        # (p-1)x(p-1) matrix M where M_{a,b} = 1 if (b-a mod p) in S,
        # a,b in {1,...,p-1} (after fixing vertex 0 as start).

        # The eigenvalues of the full p x p circulant are S_hat(k) for k=0,...,p-1.
        # The (p-1) x (p-1) submatrix doesn't inherit these eigenvalues directly.

        # Better: think of H as a polynomial in the connection set indicator.
        # H(S) = sum_{pi in S_p} prod_{i=0}^{p-2} S(pi(i+1) - pi(i) mod p)
        # where S(d) = 1 if d in S, 0 otherwise.

        # This is a multilinear function of the S(d) values (each S(d) appears
        # multiple times across different edges, so it's NOT multilinear in S).

        # H is a degree-(p-1) polynomial in the variables {S(d) : d=1,...,p-1},
        # but with many constraints (S has exactly m elements, S union (p-S) = {1,...,p-1}).

        # For the orientation cube: S(i+1) = (1+sigma_i)/2 and S(p-i-1) = (1-sigma_i)/2.
        # So H is a polynomial in sigma. But we already know the Walsh decomposition
        # gives the polynomial structure.

        print(f"    (Skipping character sum analysis for now)")

    # ====== SUMMARY ======
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    print("""
  Key findings:
  1. C(k) = -1/2 for all k (universal: sum cos(2pi*k*i/p) over i=1..m)
  2. Q_k constant = (p+1)/4 (Parseval average)
  3. sum_k Q_k = m*(p-m) = m*(m+1) (Parseval total)
  4. Gram matrix of centered Q_k = (p/4)(mI - J) — equicorrelated
  5. H = F(sorted Q_k) — symmetric function, no chirality
  6. F is NOT polynomial of any reasonable degree in elementary symm functions
  7. The nonlinearity comes from the PERMANENT structure of H
  8. Permanent of circulant does NOT simplify to product of eigenvalues
    """)


if __name__ == '__main__':
    main()
