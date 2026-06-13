"""
euler_identity_tournaments.py -- kind-pasteur-2026-03-14-S96
Thinking about e, pi, and i TOGETHER in tournament theory.

EULER'S IDENTITY: e^{i*pi} + 1 = 0
This ties 0, 1, e, pi, i into one equation.

In tournaments:
  0 = H of the empty tournament
  1 = H of the transitive tournament (= the UNIT)
  2 = the OCF fugacity (= the GENERATOR)
  e = Szele limit (= the EFFICIENCY)
  pi = geometric normalization
  i = sqrt(-1) appears in eigenvalues and roots of unity

THE QUESTION: Is there a "tournament Euler identity"?
Something like: H evaluated at i gives a formula involving e and pi?

DEEPER: i appears when we ROTATE in the complex plane.
Tournament theory has a natural ROTATION: the arc-flip involution.
Flipping an arc is like multiplying by -1 = e^{i*pi}.
The GS involution is a REFLECTION.
Can we understand these as ROTATIONS in a complex tournament space?
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def C(n, k):
    if k < 0 or k > n: return 0
    return math.comb(n, k)

maxH = {1:1, 2:1, 3:3, 4:5, 5:15, 6:45, 7:189, 8:661, 9:3357, 10:15745, 11:95095}

def main():
    e = math.e
    pi = math.pi
    I = 1j  # the imaginary unit

    print("=" * 70)
    print("THE IMAGINARY UNIT i IN TOURNAMENT THEORY")
    print("kind-pasteur-2026-03-14-S96")
    print("=" * 70)

    # ============================================================
    # PART 1: EULER'S IDENTITY AND TOURNAMENTS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 1: EULER'S IDENTITY e^{{i*pi}} + 1 = 0")
    print(f"{'='*70}")

    print(f"\n  Euler: e^(i*pi) = -1")
    print(f"  So: e^(i*pi) + 1 = 0")
    print(f"")
    print(f"  THE FIVE CONSTANTS AND THEIR TOURNAMENT ROLES:")
    print(f"    0 = H(empty tournament) = the zero element")
    print(f"    1 = H(transitive) = the unit / identity")
    print(f"    e = max_H/mean_H limit = the efficiency")
    print(f"    pi = normalization in sqrt(2*pi*n) = the geometry")
    print(f"    i = sqrt(-1) = the ROTATION / PHASE")
    print(f"")
    print(f"  But tournaments also have a SIXTH fundamental constant:")
    print(f"    2 = OCF fugacity = the binary choice")
    print(f"")
    print(f"  TOURNAMENT EULER IDENTITY?")
    print(f"  e^(i*pi) = -1 relates e, pi, i")
    print(f"  H = I(Omega, 2) relates H, Omega, 2")
    print(f"  Can we find: f(e, pi, i, 2) = H for some f?")

    # The evaluation of F(T, x) at x = e^{i*theta}:
    # F(T, e^{i*theta}) = sum over paths: e^{i*theta*fwd(P)}
    # At theta = pi: F(T, e^{i*pi}) = F(T, -1) = signed Ham path count
    # F(T, -1) is always ODD (proved in S70)

    print(f"\n  F(T, e^{{i*pi}}) = F(T, -1) = signed Ham path count")
    print(f"  F(T, -1) is always ODD (proved)")
    print(f"  So: Euler's e^{{i*pi}} = -1 applied to F gives the SIGNED count.")
    print(f"")
    print(f"  At theta = pi/2: F(T, i) = sum e^{{i*pi*fwd/2}}")
    print(f"  = sum i^{{fwd}} = sum over paths: i^{{fwd(P)}}")
    print(f"  This decomposes paths by fwd mod 4:")
    print(f"    fwd ≡ 0 (mod 4): contributes +1 (real)")
    print(f"    fwd ≡ 1 (mod 4): contributes +i (imaginary)")
    print(f"    fwd ≡ 2 (mod 4): contributes -1 (real)")
    print(f"    fwd ≡ 3 (mod 4): contributes -i (imaginary)")

    # ============================================================
    # PART 2: F(T, i) — THE COMPLEX EVALUATION
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 2: F(T, i) — EVALUATING AT THE IMAGINARY UNIT")
    print(f"{'='*70}")

    # Compute F(T, i) for all tournaments at n=5
    n = 5
    m = C(n, 2)
    N = 2**m

    F_i_real_sum = 0
    F_i_imag_sum = 0
    F_i_data = []

    for bits in range(N):
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for ii in range(n):
            for jj in range(ii+1, n):
                if bits & (1 << idx): A[jj][ii] = 1
                else: A[ii][jj] = 1
                idx += 1

        # Compute F(T, i) = sum over Ham paths: i^{fwd(P)}
        F_i = complex(0, 0)
        H = 0
        for perm in __import__('itertools').permutations(range(n)):
            valid = all(A[perm[k]][perm[k+1]] for k in range(n-1))
            if valid:
                fwd = sum(1 for k in range(n-1) if perm[k] < perm[k+1])
                F_i += I**fwd
                H += 1

        F_i_real_sum += F_i.real
        F_i_imag_sum += F_i.imag
        F_i_data.append((bits, H, F_i))

    mean_F_i = complex(F_i_real_sum / N, F_i_imag_sum / N)
    print(f"\n  n={n}: Mean F(T, i) over all tournaments = {mean_F_i.real:.4f} + {mean_F_i.imag:.4f}i")
    print(f"  |Mean F(T, i)| = {abs(mean_F_i):.4f}")

    # Is the mean purely real? Purely imaginary?
    print(f"  Purely real? {abs(mean_F_i.imag) < 0.001}")
    print(f"  Purely imaginary? {abs(mean_F_i.real) < 0.001}")

    # Sum over ALL tournaments: sum F(T, i) = ?
    total_F_i = complex(F_i_real_sum, F_i_imag_sum)
    print(f"\n  Sum F(T, i) over all {N} tournaments = {total_F_i.real:.4f} + {total_F_i.imag:.4f}i")
    print(f"  = {int(total_F_i.real)} + {int(total_F_i.imag)}i (should be integers!)")

    # This should be: sum_T sum_P i^{fwd(P)} = sum_P (sum_T that contain P) * i^{fwd(P)}
    # For each permutation P: how many tournaments contain P as a Ham path?
    # A tournament contains P iff all n-1 consecutive arcs are present.
    # The remaining C(n,2) - (n-1) arcs are free.
    # So #tournaments containing P = 2^{C(n,2)-(n-1)} = 2^{m-n+1}

    free_arcs = m - (n - 1)
    factor = 2**free_arcs
    print(f"\n  Each permutation appears in 2^{{m-(n-1)}} = 2^{free_arcs} = {factor} tournaments")
    print(f"  Sum F(T, i) = {factor} * sum_P i^{{fwd(P)}} = {factor} * F(complete_graph, i)")

    # sum_P i^{fwd(P)} over ALL n! permutations:
    # This is the EULERIAN polynomial A_n(i) where A_n(x) = sum_k A(n,k) * x^k
    # and A(n,k) = Eulerian number

    eulerian_poly_at_i = sum(
        sum(1 for perm in __import__('itertools').permutations(range(n))
            if sum(1 for k in range(n-1) if perm[k] < perm[k+1]) == fwd_val) * I**fwd_val
        for fwd_val in range(n)
    )

    print(f"\n  Eulerian polynomial A_{n}(i) = {eulerian_poly_at_i.real:.0f} + {eulerian_poly_at_i.imag:.0f}i")
    print(f"  Predicted: Sum F(T,i) = {factor} * A_{n}(i) = "
          f"{factor * eulerian_poly_at_i.real:.0f} + {factor * eulerian_poly_at_i.imag:.0f}i")
    print(f"  Actual: {total_F_i.real:.0f} + {total_F_i.imag:.0f}i")
    print(f"  Match: {abs(total_F_i - factor * eulerian_poly_at_i) < 0.01}")

    # ============================================================
    # PART 3: THE EULERIAN POLYNOMIAL AT ROOTS OF UNITY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 3: EULERIAN POLYNOMIAL AT ROOTS OF UNITY")
    print(f"{'='*70}")

    # A_n(x) = sum_k A(n,k) x^k where A(n,k) = Eulerian numbers
    # A_n(1) = n! (all permutations)
    # A_n(-1) = 0 for n >= 2 (equal ascending/descending)
    # A_n(i) = ?

    for nn in range(2, 8):
        # Compute Eulerian numbers
        from itertools import permutations as perms
        fwd_dist = [0] * nn
        for perm in perms(range(nn)):
            fwd = sum(1 for k in range(nn-1) if perm[k] < perm[k+1])
            fwd_dist[fwd] += 1

        print(f"\n  n={nn}: Eulerian numbers A({nn},k) = {fwd_dist}")

        # Evaluate at special points
        A_1 = sum(fwd_dist)
        A_neg1 = sum(fwd_dist[k] * (-1)**k for k in range(nn))
        A_i = sum(fwd_dist[k] * I**k for k in range(nn))
        A_omega = sum(fwd_dist[k] * np.exp(2j*pi/3)**k for k in range(nn))

        print(f"    A({nn}, 1) = {A_1} = {nn}!")
        print(f"    A({nn}, -1) = {A_neg1}")
        print(f"    A({nn}, i) = {A_i.real:.0f} + {A_i.imag:.0f}i, |A| = {abs(A_i):.4f}")
        print(f"    A({nn}, omega_3) = {A_omega.real:.2f} + {A_omega.imag:.2f}i, |A| = {abs(A_omega):.4f}")

    # KEY: A_n(i) determines the sum of F(T, i) over all tournaments!
    # And A_n(i) involves the mod-4 structure of the Eulerian numbers.

    # ============================================================
    # PART 4: i AND THE PALEY EIGENVALUES
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 4: i AND THE PALEY EIGENVALUES")
    print(f"{'='*70}")

    # Paley T_p eigenvalue: mu = (-1 +/- i*sqrt(p))/2
    # This can be written as: mu = (1/2) * (-1 + i*sqrt(p))
    # = (1/2) * sqrt(p+1) * e^{i*theta}
    # where theta = pi - arctan(sqrt(p))

    for p in [3, 7, 11, 19, 23]:
        mu = complex(-0.5, math.sqrt(p)/2)
        mod = abs(mu)
        arg = math.atan2(mu.imag, mu.real)

        # Express mu in terms of e, pi, i:
        # mu = mod * e^{i*arg}
        print(f"\n  Paley T_{p}:")
        print(f"    mu = {mu.real:.4f} + {mu.imag:.4f}i")
        print(f"    = {mod:.4f} * e^(i * {arg:.4f})")
        print(f"    = {mod:.4f} * e^(i * {arg/pi:.4f} * pi)")
        print(f"    |mu|^2 = (p+1)/4 = {(p+1)/4}")

        # The MAGIC: |mu|^2 = (p+1)/4 means:
        # |mu| = sqrt((p+1)/4) = sqrt(p+1)/2
        # At p=3: |mu| = 1 (unit circle!)
        # At p=7: |mu| = sqrt(2) (the first algebraic irrational!)
        # At p=11: |mu| = sqrt(3) (the variance irrational!)

    print(f"\n  THE PALEY EIGENVALUE MODULUS PROGRESSION:")
    print(f"    p=3: |mu| = sqrt(1) = 1")
    print(f"    p=7: |mu| = sqrt(2) = {math.sqrt(2):.6f}")
    print(f"    p=11: |mu| = sqrt(3) = {math.sqrt(3):.6f}")
    print(f"    p=19: |mu| = sqrt(5) = {math.sqrt(5):.6f}")
    print(f"    p=23: |mu| = sqrt(6) = {math.sqrt(6):.6f}")
    print(f"")
    print(f"  (p+1)/4: 1, 2, 3, 5, 6 for p = 3, 7, 11, 19, 23")
    print(f"  These are the FIRST FEW non-square integers!")
    print(f"  The Paley eigenvalue modulus walks through sqrt(1), sqrt(2), sqrt(3), sqrt(5), sqrt(6)")
    print(f"  SKIPPING sqrt(4) = 2 because p=15 is not prime!")

    # ============================================================
    # PART 5: e^{i*theta} AS TOURNAMENT ROTATION
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 5: e^{{i*theta}} AS TOURNAMENT ROTATION")
    print(f"{'='*70}")

    print(f"""
  THE ARC-FLIP AS ROTATION:
    Flipping arc (i,j) maps x_{{ij}} -> 1 - x_{{ij}}
    In the {{-1, +1}} basis: y_{{ij}} -> -y_{{ij}}
    Multiplying by -1 = e^{{i*pi}} is a ROTATION BY pi.

  THE T^op INVOLUTION:
    T -> T^op flips ALL arcs simultaneously.
    In Fourier: chi_S -> (-1)^{{|S|}} chi_S
    This is multiplication by e^{{i*pi*|S|}} for each basis function.

  THE GS INVOLUTION:
    GS maps paired positions: (r,c) -> (r, n-r-c).
    In the tiling space: this is a REFLECTION.
    Reflection = rotation by pi along one axis.

  THE FULL SYMMETRY GROUP:
    T^op (order 2) x GS (order 2) = Z/2 x Z/2 (Klein four-group)
    In terms of e^{{i*theta}}:
      identity = e^(i*0) = 1
      T^op = e^(i*pi) = -1
      GS = reflection (not a pure rotation)
      T^op * GS = glide reflection

  DEEPER: Can we extend to a CONTINUOUS rotation?
    Define T(theta) where T(0) = T, T(pi) = T^op.
    At intermediate theta: T(theta) is a "quantum tournament"
    with COMPLEX amplitudes e^(i*theta*x_{ij}) per arc.
    The "quantum H" = F(T, e^{{i*theta}}).
    At theta = 0: F = H (real, positive)
    At theta = pi: F = F(T, -1) (real, signed)
    At theta = pi/2: F = F(T, i) (complex)
""")

    # ============================================================
    # PART 6: THE TOURNAMENT EULER IDENTITY
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 6: THE TOURNAMENT EULER IDENTITY")
    print(f"{'='*70}")

    print(f"""
  CLASSICAL: e^(i*pi) + 1 = 0

  TOURNAMENT VERSION:
  Consider F(T, e^{{i*theta}}) as theta varies:
    F(T, e^{{i*0}}) = F(T, 1) = n! (total permutations)
    F(T, e^{{i*pi}}) = F(T, -1) (signed count, always odd)
    F(T, e^{{i*pi/2}}) = F(T, i) (complex, involves Eulerian at i)

  THE IDENTITY:
    F(T, 1) + F(T, -1) = 2 * (even-fwd paths)
    F(T, 1) - F(T, -1) = 2 * (odd-fwd paths)
    F(T, 1) = n! for ALL T (universal)
    F(T, -1) = H(T) mod 2 is 1 (always odd, from Redei)

  COMBINING:
    F(T, e^{{i*pi}}) = F(T, -1) is always odd
    F(T, e^{{i*0}}) = n! is always the same
    F(T, e^{{i*pi}}) + F(T, e^{{i*0}}) = n! + F(T, -1) (always even + odd = odd)

  THE DEEPER IDENTITY — at the OCF point x = 2:
    F(T, 2) = H(T) = I(Omega(T), 2)
    And 2 = |e^{{i*theta}}|^2 when theta is the right angle... no, |e^{{i*theta}}| = 1 always.
    But: 2 = 1 + 1 = e^0 + e^0.
    Or: 2 = e^{{i*0}} + e^{{i*pi}} ... no, that's 1 + (-1) = 0.
    Actually: 2 = e^{{i*pi/3}} + e^{{-i*pi/3}} = 2*cos(pi/3) = 2*(1/2) = 1. No!
    2*cos(0) = 2. So 2 = 2*cos(0) = e^{{i*0}} + e^{{-i*0}}.

  HMMMM. The number 2 = e^0 + e^0 is trivially true.
  But: 2 = e^{{ln 2}} and ln(2) ≈ 0.693.
  The 2 in OCF is NOT directly e^{{i*anything}} on the unit circle.
  It's OUTSIDE the unit circle: |2| > 1.

  THIS IS THE LEE-YANG CONNECTION:
  The OCF evaluates at |z| = 2 > 1 (outside the unit circle).
  The Lee-Yang zeros approach |z| = 1 (the unit circle).
  So H = I(Omega, 2) evaluates the partition function
  in the "ordered phase" (outside the critical circle).

  THE TOURNAMENT EULER IDENTITY IS:
    e^{{i*pi}} = -1, and F(T, -1) is always odd.
    So: F(T, e^{{i*pi}}) ≡ 1 (mod 2) for all tournaments T.
    This is a MOD-2 Euler identity:
      "Evaluating F at Euler's number gives a Redei-like constraint."
""")

    # ============================================================
    # PART 7: F(T, i) DATA AND PATTERNS
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 7: F(T, i) DATA — PATTERNS IN THE COMPLEX PLANE")
    print(f"{'='*70}")

    # Group F(T, i) by H value
    Fi_by_H = {}
    for bits, H, Fi in F_i_data:
        if H not in Fi_by_H:
            Fi_by_H[H] = []
        Fi_by_H[H].append(Fi)

    print(f"\n  n=5: F(T, i) grouped by H:")
    for H in sorted(Fi_by_H.keys()):
        vals = Fi_by_H[H]
        moduli = sorted(set(round(abs(v), 4) for v in vals))
        args = sorted(set(round(math.atan2(v.imag, v.real) / pi, 4) for v in vals if abs(v) > 0.01))
        print(f"    H={H:3d}: {len(vals)} tournaments, |F(i)| in {moduli[:5]}, "
              f"arg/pi in {args[:5]}")

    # Is |F(T,i)|^2 related to H?
    print(f"\n  Is |F(T,i)|^2 determined by H?")
    for H in sorted(Fi_by_H.keys()):
        vals = Fi_by_H[H]
        mod_sq = sorted(set(round(abs(v)**2) for v in vals))
        print(f"    H={H:3d}: |F(i)|^2 in {mod_sq}")

    # ============================================================
    # PART 8: THE GAUSSIAN INTEGER STRUCTURE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 8: F(T, i) AND GAUSSIAN INTEGERS")
    print(f"{'='*70}")

    # Is F(T, i) always a Gaussian integer (a + bi with a, b integers)?
    all_gaussian = True
    for bits, H, Fi in F_i_data:
        if abs(Fi.real - round(Fi.real)) > 0.01 or abs(Fi.imag - round(Fi.imag)) > 0.01:
            all_gaussian = False
            break

    print(f"\n  F(T, i) is always a Gaussian integer? {all_gaussian}")

    if all_gaussian:
        # The Gaussian integer norms: N(a+bi) = a^2 + b^2
        norms = set()
        for bits, H, Fi in F_i_data:
            a, b = round(Fi.real), round(Fi.imag)
            norm = a*a + b*b
            norms.add(norm)

        print(f"  Gaussian norms N(F(T,i)) = a^2 + b^2: {sorted(norms)}")

        # Are these norms related to H?
        norm_by_H = {}
        for bits, H, Fi in F_i_data:
            a, b = round(Fi.real), round(Fi.imag)
            norm = a*a + b*b
            if H not in norm_by_H:
                norm_by_H[H] = set()
            norm_by_H[H].add(norm)

        print(f"\n  H -> Gaussian norms N(F(T,i)):")
        for H in sorted(norm_by_H.keys()):
            print(f"    H={H:3d}: N in {sorted(norm_by_H[H])}")

    # ============================================================
    # PART 9: SYNTHESIS — i AS THE BRIDGE
    # ============================================================
    print(f"\n{'='*70}")
    print("PART 9: SYNTHESIS — i AS THE BRIDGE BETWEEN e AND pi")
    print(f"{'='*70}")

    print(f"""
  i = sqrt(-1) is the BRIDGE that connects e and pi.

  Without i: e and pi are unrelated transcendental numbers.
  With i: e^(i*pi) = -1, tying them together.

  IN TOURNAMENT THEORY:

  Without i:
    e governs EFFICIENCY (max_H/mean_H -> e)
    pi governs GEOMETRY (sqrt(2*pi*n) normalization)
    They seem unrelated.

  With i:
    F(T, e^(i*pi)) = F(T, -1) = signed count (always odd) [CONNECTS e AND pi TO H]
    Paley eigenvalue mu = (-1 + i*sqrt(p))/2 [CONNECTS i TO spectral theory]
    Lee-Yang zeros on |z| = 1 = {{e^(i*theta)}} [CONNECTS i TO phase transitions]
    F(T, i) is a Gaussian integer [CONNECTS i TO number theory]

  i IS THE ROTATION THAT CREATES PHASE:
    At x = 1: F = n! (trivial, all permutations)
    At x = 2: F = H (the tournament invariant)
    At x = -1: F = signed count (Redei-like)
    At x = i: F = Gaussian integer (number-theoretic)
    At x = e^(i*theta): F traces a curve in C as theta varies

  THE TOURNAMENT LIVES ON A LINE IN C:
    F(T, x) for real x >= 0 gives a real positive function.
    Evaluating at x = 2 gives H.
    Extending to complex x via e^(i*theta) reveals:
    - The signed structure (at theta = pi)
    - The 4-fold symmetry (at theta = pi/2)
    - The Lee-Yang phase boundary (at |x| = 1)

  THE DEEPEST INSIGHT:
  The tournament IS a path in the complex plane.
  H = F(T, 2) is a SINGLE POINT on this path.
  The full path F(T, x) for x in C contains ALL the information.
  The values at x = 1, 2, -1, i are four "projections" of this path.
""")

    print(f"\n{'='*70}")
    print("DONE — i AS THE BRIDGE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
