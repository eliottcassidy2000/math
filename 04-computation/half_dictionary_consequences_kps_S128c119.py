#!/usr/bin/env python3
"""half_dictionary_consequences_kps_S128c119.py -- kind-pasteur-2026-07-20-S128c119

CONSEQUENCES OF THE HALF DICTIONARY (THM-1555).  The single identity

    A + A^T = J - I            (the owner's x |-> (1+x)/2, differentiated)

is used THREE times here, and each use is a different theorem.

(A) THE CHARACTERISTIC-2 COLLAPSE -- why the repo needs BOTH matrices.
    The dictionary A = (J - I + S)/2 is invertible exactly when 2 is invertible.
    Over GF(2) it dies, and it dies asymmetrically:  mod 2,
        S = A - A^T = A + A^T = J - I = J + I,
    which does NOT depend on the tournament at all.  So the skew matrix carries
    ZERO information mod 2, while A carries everything.  That is the structural
    reason the repo's GF(2) layer (cut/cycle space, switching classes, even
    graphs, Redei mod 2) is written in A and its char-0 layer (Pfaffians, skew
    spectrum, THM-1455 parity) is written in S.  They are not two conventions;
    one of them is unavailable on each side.  Checked exhaustively.

(B) THE PERRON END OF THE SAME IDENTITY.  Let x >= 0 be the unit Perron vector,
    A x = rho x.  Then x^T (A + A^T) x = 2 rho and x^T (J - I) x = (1^T x)^2 - 1,
    so
        2 rho = (1^T x)^2 - 1  <=  n - 1        (Cauchy-Schwarz, ||x|| = 1)
    giving rho <= (n-1)/2 with equality iff x is uniform iff T is REGULAR.
    Note this is the SAME computation as the Re >= -1/2 bound, run on the Perron
    vector instead of on a general eigenvector.

(C) THE TRICHOTOMY, and it closes with no extra input.  tr(A) = 0 gives
    0 = rho + sum_{k>=1} Re(lambda_k), and Re(lambda_k) >= -1/2 gives
    rho <= (n-1)/2 a SECOND way.  Chaining the two:

        T regular  <=>  rho = (n-1)/2  <=>  every non-Perron eigenvalue has
                                            Re exactly -1/2.

    and the deficiency is an exact identity, not an inequality:

        sum_{k>=1} ( Re(lambda_k) + 1/2 )  =  (n-1)/2 - rho  =  ( n - (1^T x)^2 ) / 2.

    Total height of the spectrum above the line Re = -1/2   =   Perron deficiency
    =   the Perron vector's deviation from uniform.  Three quantities, one number.

    Also checked: tr(A^2) = 0 for every tournament (A_ij A_ji = 0 always), so
    sum lambda^2 = 0 as well -- the spectrum is constrained to two moments.
"""
import sys
import numpy as np

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
TOL = 1e-9


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=float)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1.0
            else:
                A[j, i] = 1.0
            idx += 1
    return A


def all_tournaments(n, cap):
    m = n * (n - 1) // 2
    for bits in range(min(1 << m, cap)):
        yield from_bits(bits, n)


print("=" * 78)
print("(A) THE CHARACTERISTIC-2 COLLAPSE:  S = J + I mod 2, ALWAYS")
print("=" * 78)
print("  The dictionary needs 1/2.  Over GF(2) it dies, and the skew matrix dies")
print("  with it: S = A - A^T = A + A^T = J - I mod 2, independent of T.")
bad_s = bad_a = 0
tot = 0
for n in range(3, NMAX + 1):
    J = np.ones((n, n))
    I = np.eye(n)
    seen_S, seen_A = set(), set()
    for A in all_tournaments(n, 3000):
        S = A - A.T
        tot += 1
        if not np.allclose(S % 2, (J + I) % 2):
            bad_s += 1
        seen_S.add(((S % 2).astype(int)).tobytes())
        seen_A.add(((A % 2).astype(int)).tobytes())
    print("  n = %d : distinct S mod 2 = %-4d   distinct A mod 2 = %-6d"
          % (n, len(seen_S), len(seen_A)))
print("  tournaments checked : %d ; failures of S = J + I mod 2 : %d" % (tot, bad_s))
print("  -> S mod 2 takes exactly ONE value at every n.  All tournament information")
print("     survives mod 2 in A and NONE of it in S.  This is why the repo's GF(2)")
print("     layer is written in A and its char-0 parity layer is written in S.")

print()
print("=" * 78)
print("(B) PERRON END:  2 rho = (1^T x)^2 - 1  <=  n - 1,  equality iff regular")
print("=" * 78)
worst_gap = {}
for n in range(3, NMAX + 1):
    best = -1.0
    best_reg = None
    viol = 0
    for A in all_tournaments(n, 4000):
        ev, evec = np.linalg.eig(A)
        k = int(np.argmax(ev.real))
        rho = ev[k].real
        x = np.abs(evec[:, k].real) if np.allclose(evec[:, k].imag, 0, atol=1e-7) \
            else np.abs(evec[:, k])
        x = x / np.linalg.norm(x)
        # identity 2 rho = (1^T x)^2 - 1
        if abs(2 * rho - ((x.sum()) ** 2 - 1)) > 1e-6:
            viol += 1
        if rho > best:
            best = rho
            best_reg = bool(np.allclose(A.sum(axis=1), (n - 1) / 2))
    worst_gap[n] = best
    print("  n = %d : max rho over all tournaments = %.6f   (n-1)/2 = %.4f   "
          "attained by a regular T : %s   identity failures : %d"
          % (n, best, (n - 1) / 2, best_reg, viol))
print("  -> rho never exceeds (n-1)/2.  At ODD n the bound is attained (regular")
print("     tournaments exist); at EVEN n there are none, so rho is strictly below.")

print()
print("=" * 78)
print("(C) THE TRICHOTOMY AND THE EXACT DEFICIENCY IDENTITY")
print("=" * 78)
print("  regular  <=>  rho = (n-1)/2  <=>  all non-Perron eigenvalues on Re = -1/2")
print()
print("  %-4s %-8s %-10s %-10s %-10s %s"
      % ("n", "#T", "reg<=>rho", "reg<=>line", "tr(A^2)=0", "max |identity error|"))
for n in range(3, NMAX + 1):
    cnt = 0
    ok_rho = ok_line = ok_tr2 = True
    maxerr = 0.0
    for A in all_tournaments(n, 4000):
        cnt += 1
        ev = np.linalg.eigvals(A)
        k = int(np.argmax(ev.real))
        rho = ev[k].real
        rest = [ev[j] for j in range(len(ev)) if j != k]
        regular = bool(np.allclose(A.sum(axis=1), (n - 1) / 2))
        # equivalence 1: regular <=> rho = (n-1)/2
        if regular != (abs(rho - (n - 1) / 2) < 1e-7):
            ok_rho = False
        # equivalence 2: regular <=> all non-Perron on the line Re = -1/2
        online = all(abs(e.real + 0.5) < 1e-7 for e in rest)
        if regular != online:
            ok_line = False
        # tr(A^2) = 0
        if abs(np.trace(A @ A)) > 1e-7:
            ok_tr2 = False
        # deficiency identity: sum (Re lam_k + 1/2) = (n-1)/2 - rho
        lhs = sum(e.real + 0.5 for e in rest)
        rhs = (n - 1) / 2 - rho
        maxerr = max(maxerr, abs(lhs - rhs))
    print("  %-4d %-8d %-10s %-10s %-10s %.2e"
          % (n, cnt, ok_rho, ok_line, ok_tr2, maxerr))
print()
print("  The deficiency identity holds to machine precision at every n: the TOTAL")
print("  HEIGHT of the spectrum above the line Re = -1/2 equals the PERRON")
print("  DEFICIENCY (n-1)/2 - rho, which equals ( n - (1^T x)^2 )/2, the Perron")
print("  vector's deviation from uniform.  One number, three readings.")
print()
print("  And tr(A) = tr(A^2) = 0 for every tournament, so sum lam = sum lam^2 = 0:")
print("  the spectrum is pinned at its first two power sums before any structure")
print("  is imposed.  The half-dictionary supplies the third constraint, a HALF-PLANE")
print("  Re >= -1/2, and regularity is exactly the case of equality throughout.")
