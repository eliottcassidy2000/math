"""
hessian_sign_origin.py — Why does the Hessian sign flip at p mod 4?

THE QUESTION: H is an exact polynomial in e_k(y^2) on the Parseval simplex.
The coefficients {c_k} determine the Hessian eigenvalue lambda_H.
lambda_H < 0 for p=3 mod 4 (center=max), lambda_H > 0 for p=1 mod 4 (center=min).

WHY? What determines the sign of the c_k coefficients?

APPROACH:
  1. The e_k(y^2) are elementary symmetric functions of y_k^2 = Im(lambda_k)^2
  2. The eigenvalues lambda_k = sum_{s in S} omega^{ks} depend on the connection set S
  3. H = I(Omega(T), 2) where Omega is the conflict graph of all odd directed cycles
  4. The bridge is: cycles are counted by eigenvalue traces (c_k = tr(A^k)/k for k<=5)
     and the cycle ARRANGEMENT determines Omega

  Key observation: for circulant tournaments, the eigenvalues parametrize the
  conflict graph, and the e_k decomposition captures how the OCF (independence
  polynomial) depends on the spectral parameters.

  SPECIFIC HYPOTHESIS: The coefficient c_k in H = sum c_k * e_k(y^2) is related
  to the k-th alpha (independent k-sets) of Omega via the OCF formula.

  At p=7: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
          H = 170.625 + 2*e_2(y^2)
  The question: is there a formula relating alpha_k to e_k?

  ALTERNATIVE APPROACH: Use the permanent formula.
  H = perm(A) for the adjacency matrix of the tournament.
  For circulant matrices, Minc (1978) showed:
    perm(circ(a_0, a_1, ..., a_{n-1})) = sum_{S subset {0,...,n-1}, |S|=k}
                                           product_{i in S} lambda_i
  Wait, that's the determinant, not permanent!

  Actually, the permanent of a circulant matrix does NOT simplify to a product
  of eigenvalues. The permanent is #P-hard in general.

  BUT: for a 0-1 circulant matrix (which is what a circulant tournament gives),
  perm(A) counts the number of permutation matrices dominated by A.
  This is exactly the number of Hamiltonian cycles in the bipartite double cover...
  No, for tournaments perm(A) = H(T) (number of Hamiltonian paths,
  using the fact that A has 0 diagonal).

  Actually, perm(A) for a tournament adjacency matrix is NOT H(T).
  perm(A) counts perfect matchings in the bipartite double cover.
  H(T) = number of Hamiltonian PATHS = sum over permutations sigma of
  product_{i=0}^{n-2} A[sigma(i), sigma(i+1)] ... no, this is not the permanent.

  H(T) = sum over permutations sigma of product_{i=0}^{n-2} A[P(i), P(i+1)]
  where P is a path ordering. This is the PERMANENT of the TRANSFER matrix, not A.

  Actually, H(T) relates to the permanent differently. Let me think...

  For a digraph D, the number of Hamiltonian paths from any start to any end
  is the permanent of (J - I) restricted to the adjacency... no.

  H(T) = sum_sigma product A[sigma(i), sigma(i+1)] for i=0,...,n-2
       = sum_sigma product A[sigma(i), sigma(i+1)]

  This is the sum over all permutations of a product of n-1 adjacency entries
  along a path, not n entries around a cycle. So it's NOT the permanent.

  Let's just compute what we can and look for patterns.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import cmath
import math
from itertools import combinations
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    S_set = set(S)
    adj = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in S_set: adj[i][j] = True
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if dp[mask][v]==0 or not(mask&(1<<v)): continue
            for w in range(n):
                if mask&(1<<w): continue
                if adj[v][w]: dp[mask|(1<<w)][w] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j*cmath.pi/n)
    return [sum(omega**(k*s) for s in S) for k in range(n)]


def elementary_symmetric(xs, k):
    if k == 0: return 1.0
    if k > len(xs): return 0.0
    total = 0.0
    for subset in combinations(range(len(xs)), k):
        prod = 1.0
        for i in subset:
            prod *= xs[i]
        total += prod
    return total


def solve_linear_system(A, b):
    n = len(b)
    M = [row[:] + [bi] for row, bi in zip(A, b)]
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[max_row][col]) < 1e-12: return None
        M[col], M[max_row] = M[max_row], M[col]
        for row in range(col+1, n):
            f = M[row][col] / M[col][col]
            for j in range(col, n+1): M[row][j] -= f * M[col][j]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j]*x[j] for j in range(i+1,n))) / M[i][i]
    return x


def main():
    print("=" * 70)
    print("HESSIAN SIGN ORIGIN — WHY DOES IT FLIP AT p mod 4?")
    print("=" * 70)

    # ================================================================
    # SECTION 1: Coefficient analysis across primes
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 1: e_k COEFFICIENT TABLE ACROSS PRIMES")
    print(f"{'=' * 60}")

    all_coeffs = {}
    all_lambdaH = {}

    for p in [3, 5, 7, 11, 13]:
        m = (p-1)//2
        all_S = all_circulant_tournaments(p)

        data = []
        for S in all_S:
            H = ham_count_dp(p, S)
            eigs = circulant_eigenvalues(p, S)
            y2 = [eigs[k].imag**2 for k in range(1, m+1)]
            esyms = [elementary_symmetric(y2, k) for k in range(m+1)]
            data.append({'S': S, 'H': H, 'y2': y2, 'esyms': esyms})

        data.sort(key=lambda d: d['H'], reverse=True)
        H_vals = sorted(set(d['H'] for d in data), reverse=True)
        n_H = len(H_vals)

        # Get representatives
        reps = [next(d for d in data if d['H']==H_val) for H_val in H_vals]

        # Fit H = c_0 + c_2*e_2 + ... + c_m*e_m
        # Try increasing number of terms until exact
        best_coeffs = None
        for n_terms in range(1, min(n_H+1, m)):
            k_start = 2
            k_end = k_start + n_terms
            n_vars = 1 + n_terms

            if n_H < n_vars:
                continue

            X = [[1.0] + [d['esyms'][k] for k in range(k_start, k_end)] for d in reps]
            y = [float(d['H']) for d in reps]

            XtX = [[sum(X[i][r]*X[i][c] for i in range(n_H)) for c in range(n_vars)] for r in range(n_vars)]
            Xty = [sum(X[i][r]*y[i] for i in range(n_H)) for r in range(n_vars)]

            coeffs = solve_linear_system(XtX, Xty)
            if coeffs is None:
                continue

            residuals = []
            for d in reps:
                pred = coeffs[0]
                for j, k in enumerate(range(k_start, k_end)):
                    pred += coeffs[j+1] * d['esyms'][k]
                residuals.append(d['H'] - pred)

            max_res = max(abs(r) for r in residuals)
            if max_res < 0.01:
                best_coeffs = {0: coeffs[0]}
                for j, k in enumerate(range(k_start, k_end)):
                    best_coeffs[k] = coeffs[j+1]
                break

        if best_coeffs:
            all_coeffs[p] = best_coeffs

            # Compute Hessian eigenvalue
            v = p / 4.0
            lambda_H = 0
            for k, c in best_coeffs.items():
                if k >= 2:
                    binom = math.comb(m-2, k-2)
                    lambda_H -= c * binom * v**(k-2)

            all_lambdaH[p] = lambda_H

            print(f"\n  p = {p} (mod 4 = {p%4}), m = {m}")
            print(f"    H = {best_coeffs[0]:.4f}", end="")
            for k in sorted(k2 for k2 in best_coeffs if k2 >= 2):
                sign = "+" if best_coeffs[k] >= 0 else ""
                print(f" {sign}{best_coeffs[k]:.4f}*e_{k}", end="")
            print()
            print(f"    lambda_H = {lambda_H:.4f}")
            print(f"    Sign pattern: {['+'if best_coeffs.get(k,0)>=0 else '-' for k in range(2, max(best_coeffs.keys())+1)]}")

    # ================================================================
    # SECTION 2: Rationality of coefficients
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 2: ARE COEFFICIENTS RATIONAL?")
    print(f"{'=' * 60}")

    for p in [7, 11, 13]:
        if p not in all_coeffs:
            continue
        print(f"\n  p = {p}:")
        for k in sorted(all_coeffs[p]):
            c = all_coeffs[p][k]
            # Try to express as fraction
            frac = Fraction(c).limit_denominator(10000)
            err = abs(float(frac) - c)
            print(f"    c_{k} = {c:.10f} ~ {frac} (error = {err:.2e})")

    # ================================================================
    # SECTION 3: Coefficient growth analysis
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 3: COEFFICIENT STRUCTURE")
    print(f"{'=' * 60}")
    print("""
    Observation: The coefficients grow rapidly with k.
    At p=13: c_2 = 6596, c_3 = -3730, c_4 = 1835, c_5 = -401, c_6 = -700

    The alternating sign pattern (+,-,+,-,-) breaks at k=6.
    The Hessian formula is:
      lambda_H = -sum_k c_k * C(m-2, k-2) * (p/4)^{k-2}

    The term with largest MAGNITUDE dominates.
    At p=7: only one term (k=2), trivially negative.
    At p=11: three terms. The k=3 term (positive, from c_3 < 0) is large
             but k=2 and k=4 terms together outweigh it.
    At p=13: five terms. The k=4 and k=5 terms change the balance.
    """)

    for p in [7, 11, 13]:
        if p not in all_coeffs:
            continue
        m = (p-1)//2
        v = p/4.0
        print(f"  p = {p}, v = p/4 = {v:.2f}:")
        total = 0
        for k in sorted(k2 for k2 in all_coeffs[p] if k2 >= 2):
            c = all_coeffs[p][k]
            binom = math.comb(m-2, k-2)
            power = v**(k-2)
            term = -c * binom * power
            total += term
            print(f"    k={k}: c_k={c:>12.4f}, C({m-2},{k-2})={binom:>4}, "
                  f"v^{k-2}={power:>10.4f}, term={term:>12.4f} "
                  f"({'+' if term > 0 else '-'})")
        print(f"    SUM = lambda_H = {total:.4f} ({'+' if total > 0 else '-'})")

    # ================================================================
    # SECTION 4: Connection to Gauss sum structure
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 4: GAUSS SUM AND COEFFICIENT ORIGIN")
    print(f"{'=' * 60}")
    print("""
    For Paley T_p (p = 3 mod 4):
      y_k^2 = |g_k|^2 where g_k = sum_{a in QR} omega^{ak} (Gauss sum)
      g_k = chi(k) * g_1, so |g_k| = |g_1| = sqrt(p)/2
      => y_k^2 = p/4 for all k (spectral flatness)

    For p = 1 mod 4: -1 in QR, so S and -S overlap.
    The QR set is NOT a valid tournament connection set.
    No circulant tournament achieves flat spectrum.

    The coefficient signs seem to encode how H's landscape changes
    as we move away from the flat point on the spectral simplex.

    KEY INSIGHT: At the flat point (y_k^2 = p/4 for all k),
    e_k(y^2) = C(m,k) * (p/4)^k. This is the maximum of e_k (Schur).
    Moving away from flat DECREASES all e_k simultaneously.
    If c_k > 0: this decreases H (good for Paley being max).
    If c_k < 0: this INCREASES H (bad for Paley).

    The Hessian sign depends on the NET effect at second order.
    """)

    # Compute e_k at uniform point
    for p in [7, 11, 13]:
        m = (p-1)//2
        v = p/4.0
        print(f"\n  p = {p}: e_k at uniform point (v = {v:.2f}):")
        for k in range(1, m+1):
            ek_uniform = math.comb(m, k) * v**k
            print(f"    e_{k} = C({m},{k}) * {v}^{k} = {ek_uniform:.4f}")

    # ================================================================
    # SECTION 5: The dichotomy mechanism
    # ================================================================
    print(f"\n{'=' * 60}")
    print(f"SECTION 5: DICHOTOMY MECHANISM")
    print(f"{'=' * 60}")
    print("""
    CLAIM: The sign of lambda_H is determined by whether
    the ACHIEVABLE spectral configurations (circulant tournaments on Z_p)
    favor "expansion from center" or "contraction toward center".

    For p = 3 mod 4:
      - Paley (QR) achieves the CENTER (flat spectrum)
      - ALL other circulants are FARTHER from center
      - H DECREASES as we move away => lambda_H < 0

    For p = 1 mod 4:
      - No circulant achieves the center
      - The CLOSEST to center (Satake/NDR) has INTERMEDIATE H
      - The FARTHEST from center (cyclic interval) has MAXIMUM H
      - H INCREASES as we move away => lambda_H > 0

    BUT THIS IS THE EMPIRICAL OBSERVATION, NOT THE CAUSE.
    The CAUSE is in the coefficient signs c_k.
    And the coefficients are determined by the POLYNOMIAL structure
    of H(e_2, ..., e_m) which comes from the OCF and the conflict
    graph Omega's dependence on eigenvalues.

    The most promising route to understanding:
    H = I(Omega, 2) where Omega has alpha_k independent k-sets.
    The alpha_k are complicated functions of the cycle arrangement.
    For circulant tournaments, alpha_k depends on the eigenvalues
    through the cycle counts (which ARE trace functions of eigenvalues).

    The polynomial H(e_k) encodes how DISJOINT CYCLE COLLECTIONS
    vary with the spectral parameters. The sign of each c_k
    tells us whether adding more "spectral concentration" (lower e_k)
    helps or hurts the total count of independent cycle collections
    of each size.
    """)


if __name__ == '__main__':
    main()
