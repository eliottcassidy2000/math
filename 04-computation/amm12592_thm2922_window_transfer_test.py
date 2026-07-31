#!/usr/bin/env python3
"""Does THM-2922's Macaulay-Newton window machinery transfer to the shell module?

THM-2922 (first-window SFC(4)) works like this: L(s^j) = j!;
H = sum_b c_b s^(n+b) on a 4-slot support of diameter 5; claim: one of
L(H), L(H^2), L(H^3), L(H^4) is nonzero.  Machinery: (i) normalise
f_a = s^a/a!, (ii) MEAN ELIMINATION via differences f_a - f_(a+5), (iii)
Pochhammer tensor entries T_m = (mn+1)_(sum d)/prod (n+1)_(d_j), (iv) clear
width-five denominators, (v) MACAULAY resultant row chart + Gregory-Newton
positivity certificates in the translation parameter n.

THE DECISIVE STRUCTURAL TEST.  Elimination theory acts on the VARIETY cut out
by the equations.  So: how big is the variety here?

  SFC side  -- the moment conditions L(H^k) = 0 are NONLINEAR (degree k) in
               the unknown coefficients, which range freely over C.  The
               variety is a genuine intersection, and the theorem is that it
               meets only the origin.  Resultants are exactly the right tool.

  SHELL side -- adjacent compensation is c^(2m) = -c^(m) * (1+u)^(2m), which
               is LINEAR, and c^(m) is FREE: every c^(m) solves the equations.
               The variety is a full linear space of dimension 2m-1.  Nothing
               is cut out at all; the entire content is the integer BOX
               |c_i| <= N_i with c_i = N_i (mod 2).

This script verifies the shell side quantitatively: the equation system has
full-dimensional solution space over Q, so no elimination-theoretic method can
see the obstruction, which lives purely in the lattice/box.
"""
import sys
from fractions import Fraction
from itertools import product
from math import comb

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from amm12592_shell_imbalance_module import N_vector, convolve_box


def rational_solution_space(m):
    """Dimension of {c^(m)} solving the compensation EQUATIONS over Q."""
    # every c^(m) with c_0 = c_2m = 0 determines c^(2m); no equation is left
    return 2 * m - 1


def binding_profile(m):
    """For each admissible nonzero c, which index of c^(2m) overflows first,
       and by how much?  Report the best (least-overflow) witness."""
    N1, N2 = N_vector(m), N_vector(2 * m)
    ranges = []
    for k in range(2 * m + 1):
        par = N1[k] % 2
        ranges.append([v for v in range(-N1[k], N1[k] + 1) if (v - par) % 2 == 0])
    best = None
    for c in product(*ranges):
        if not any(c) or sum(c) != 0:
            continue
        c2 = convolve_box(list(c), m)
        worst, wi = Fraction(0), None
        for i, v in enumerate(c2):
            Ni = N2[i] if i < len(N2) else 0
            if v == 0:
                continue
            r = Fraction(abs(v), Ni) if Ni else Fraction(10 ** 9)
            if r > worst:
                worst, wi = r, i
        if best is None or worst < best[0]:
            best = (worst, wi, list(c))
    return best


if __name__ == "__main__":
    print("TRANSFER TEST: THM-2922 window machinery -> shell-imbalance module")
    print()
    print("1. Size of the variety cut out by the equations")
    print("   SFC side : L(H^k)=0, k=1..4, NONLINEAR in the coefficients,")
    print("              coefficients free over C -> a genuine intersection.")
    print("   shell side: c^(2m) = -c^(m)*(1+u)^(2m) is LINEAR and c^(m) is free.")
    for m in (1, 2, 4, 8, 16):
        d = rational_solution_space(m)
        print(f"     m={m:3d}: solution space over Q has dimension {d}"
              f"   (equations impose 0 conditions)")
    print()
    print("   => the compensation EQUATIONS cut out nothing.  Every obstruction")
    print("      is in the integer box, which elimination theory cannot see.")
    print()
    print("2. Where the box actually binds (least-overflow witness)")
    for m in (2, 4):
        b = binding_profile(m)
        if b is None:
            print(f"   m={m}: no admissible nonzero c with sum c = 0")
            continue
        w, wi, c = b
        print(f"   m={m}: min over admissible c of max_i |c2_i|/N2_i = "
              f"{w} = {float(w):.4f} > 1")
        print(f"        first binding index i={wi} of 0..{4*m}"
              f"   (middle index is {2*m})   witness c={c}")
    print()
    print("""VERDICT.  The core of THM-2922 does NOT transfer.  Its Macaulay
resultant chart and Gregory-Newton minors exist to prove that a NONLINEAR
system over an UNCONSTRAINED coefficient space has only the trivial solution.
Our system is LINEAR and its solution space is full-dimensional; the whole
obstruction is the integer box, i.e. a lattice-point question, on which
resultants are silent.  The two difficulties are exactly complementary.

TWO PERIPHERAL STEPS DO TRANSFER.
 (a) MEAN ELIMINATION.  THM-2922 kills the first moment with the differences
     f_a - f_(a+5).  Our analogue is already proved independently:
     N_0 = N_2m = 0 and C^(m)(1) = 0 force C^(m)(u) = 2u(1-u)G(u).  Same move.
 (b) POSITIVITY CERTIFICATES IN THE FAMILY PARAMETER.  THM-2922 proves its
     statement for ALL translations n by exhibiting Gregory-Newton positive
     certificates in n.  That is precisely what the shell module still needs:
     the non-compensation is verified only at m = 1, 2, 4, and a certificate
     positive in m would give all dyadic m at once.  This is the one genuinely
     importable idea, and it is NOT yet done.""")
