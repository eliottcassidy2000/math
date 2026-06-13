#!/usr/bin/env python3
"""
ihara_to_IP.py -- Connecting Ihara zeta function to Independence Polynomial

The Ihara zeta function Z_T(u) of a tournament T satisfies:
  Z_T(u)^{-1} = det(I - u*A + u^2*Q)  for undirected graphs
  For directed graphs: Z_T(u)^{-1} = det(I - u*A)  (simpler!)

For tournaments: all out-degrees = m, so Q = (m-1)*I.

The key: Z_T(u) encodes ALL cycle information.
The independence polynomial I(Omega, x) also encodes cycle information (via Omega).

QUESTION: Is there a closed-form relationship between Z_T(u) and I(Omega, 2)?

Approach: Expand log Z_T(u) in terms of prime cycles:
  log Z_T(u) = sum_{p prime cycle} sum_{k>=1} (u^{k*|p|}/k)
  = sum_{n>=1} (c_n/n) * u^n  where c_n = #{directed n-cycles}

Then: Z_T(u) = exp(sum c_n * u^n / n)
  = prod_{n odd} exp(c_n * u^n / n) * prod_{n even} exp(c_n * u^n / n)

For tournaments (no 2-cycles): c_2 = 0, and c_n for even n > 2 exists.

The OCF says H = I(Omega, 2) = sum_{S independent} 2^|S|.
And the conflict graph Omega has vertices = {directed odd cycles in T}.

Can we express I(Omega, 2) in terms of the c_n and their conflict structure?

Actually, I(Omega, 2) = sum_S 2^|S| where S ranges over vertex-disjoint
collections of odd cycles. This is:
  I(Omega, 2) = prod_{C odd cycle} (1 + 2*x_C)  summed over independent sets
  where x_C = 1 if C in S, 0 otherwise.

This is NOT simply expressible in terms of c_n alone -- it depends on the
DISJOINTNESS structure.

But there's a graph-theoretic identity:
  I(G, x) = 1 + x*alpha_1 + x^2*alpha_2 + ...
  = det(I + x*D - x*A_G)  if G is a claw-free graph (Chudnovsky-Seymour)

Wait, that's not right. The independence polynomial doesn't have a determinantal
formula in general. But for LINE GRAPHS:
  I(L(G), x) = matching polynomial of G
  = det(x*I - A_G) up to some scaling

For our Omega: Omega is claw-free at n <= 8 (from MEMORY), so we CAN use
the Chudnovsky-Seymour theory for small cases.

Let me just compute and compare.

Author: kind-pasteur-2026-03-12-S57
"""

import math
import cmath


def compute_ihara_zeta_coeffs(p, S, max_k=None):
    """Compute c_k (directed cycle counts) from Ihara expansion.

    For a digraph D: Z(u)^{-1} = det(I - u*A)
    = prod_j (1 - u*lambda_j) where lambda_j are eigenvalues.

    log Z(u) = -sum_j log(1 - u*lambda_j) = sum_j sum_{k>=1} lambda_j^k u^k / k
    = sum_{k>=1} (tr(A^k)/k) * u^k

    So the k-th coefficient of log Z is tr(A^k)/k, which counts all CLOSED WALKS
    of length k, not just prime cycles. For prime cycles:
    c_k^{prime} = (1/k) sum_{d|k} mu(k/d) tr(A^d)

    For k prime: c_k^{prime} = (tr(A^k) - tr(A)) / k = tr(A^k) / k
    (since tr(A) = 0 for tournaments).

    But these count DIRECTED cycles, not undirected. A directed cycle
    (v_0,...,v_{k-1},v_0) is counted once per starting vertex... no, tr(A^k)
    counts each directed cycle k times (once per start). So the number of
    distinct directed k-cycles is tr(A^k)/k (for prime k with no shorter dividing).
    """
    if max_k is None:
        max_k = p

    omega = cmath.exp(2j * cmath.pi / p)
    eigs = [sum(omega ** (j * s) for s in S) for j in range(p)]

    result = {}
    for k in range(3, max_k + 1):
        try:
            trk = round(sum(e ** k for e in eigs).real)
        except OverflowError:
            break
        # For ODD k: the "directed k-cycle" count is related to tr(A^k)
        # but includes walks that revisit vertices.
        # For k=3: tr(A^3) = 3 * c_3 (exact, no shorter cycles possible)
        # For k=5: tr(A^5) = 5*c_5 + ??? (some contribution from 3-cycle walks)
        result[k] = trk

    return result


def main():
    print("=" * 70)
    print("IHARA ZETA AND CYCLE STRUCTURE")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        S_int = list(range(1, m + 1))

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            print(f"\n--- p={p}, {name} ---")
            coeffs = compute_ihara_zeta_coeffs(p, S)

            print(f"  k  tr(A^k)  tr/k")
            for k, trk in sorted(coeffs.items()):
                print(f"  {k:2d}  {trk:12d}  {trk/k:12.2f}")

    # Now compute the Ihara zeta determinant explicitly for small p
    print("\n" + "=" * 70)
    print("IHARA ZETA DETERMINANT: Z(u)^{-1} = det(I - u*A)")
    print("=" * 70)

    for p in [7]:
        m = (p - 1) // 2
        S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
        S_int = list(range(1, m + 1))

        omega = cmath.exp(2j * cmath.pi / p)

        for name, S in [("Paley", S_qr), ("Interval", S_int)]:
            eigs = [sum(omega ** (j * s) for s in S) for j in range(p)]

            # det(I - u*A) = prod_j (1 - u*lambda_j)
            # Expand as polynomial in u
            # This is the characteristic polynomial evaluated at 1/u, times u^p
            print(f"\n  p={p}, {name}:")
            print(f"  Eigenvalues: {[f'{e.real:.2f}{e.imag:+.2f}j' for e in eigs]}")

            # Evaluate Z(u)^{-1} at u = 1/(2*m) (near convergence radius)
            # More interesting: evaluate at the "critical" u = 1/m
            for u_val in [0.1, 0.2, 0.3]:
                det_val = 1
                for e in eigs:
                    det_val *= (1 - u_val * e)
                print(f"  Z({u_val})^{{-1}} = det(I - {u_val}*A) = {det_val.real:.6f}")

    # Key insight: H = I(Omega, 2) cannot be extracted from Z alone
    # because Z encodes c_k (all cycles) but not the DISJOINTNESS structure.
    # The independence polynomial of Omega is a FINER invariant than Z.
    print("\n" + "=" * 70)
    print("KEY INSIGHT")
    print("=" * 70)
    print("""
  The Ihara zeta function Z_T(u) encodes COUNTS of directed k-cycles:
    log Z = sum c_k u^k / k

  But H = I(Omega, 2) depends on INDEPENDENT SETS of cycles (disjointness).

  Two tournaments can have identical c_k for all k but different H values,
  if their cycles have different conflict structures.

  At p=7:
    - Paley: c_3=14, c_5=42, c_7=24, alpha_1=80, alpha_2=7, H=189
    - Interval: c_3=14, c_5=28, c_7=17, alpha_1=59, alpha_2=14, H=175

  Same c_3, different c_5 and c_7, but crucially different alpha_2!
  Interval has DOUBLE the disjoint cycle pairs (14 vs 7).

  This means H-maximization is NOT purely a spectral problem.
  It requires understanding the GEOMETRIC structure of cycle overlaps,
  which is where the additive energy / flow picture becomes essential.
""")


if __name__ == '__main__':
    main()
