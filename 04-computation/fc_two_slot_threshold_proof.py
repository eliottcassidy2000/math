#!/usr/bin/env python3
"""RIGOROUS proof of the fixed-support threshold at N = 2 slots, any dimension.

CLAIM.  Let f = c_1 x^alpha + c_2 x^beta with alpha != beta multi-indices and
(c_1,c_2) != 0 complex.  Then L(f) and L(f^2) cannot BOTH vanish.  I.e. for a
two-slot support the first-window threshold M = N = 2 is attained: two moments
already force f = 0.

PROOF.  Write A = alpha! = prod_i alpha_i!, B = beta!, C = (alpha+beta)!,
D = (2 alpha)!, E = (2 beta)!.  Then

    L(f)   = c_1 A + c_2 B,
    L(f^2) = c_1^2 D + 2 c_1 c_2 C + c_2^2 E.

If L(f) = 0 then c_2 = -c_1 A / B (note B != 0), and c_1 != 0 (else both
vanish), so

    L(f^2) = (c_1^2 / B^2) * [ D B^2 - 2 A B C + A^2 E ]  =:  (c_1^2/B^2) * Q.

So it suffices to show Q > 0 whenever alpha != beta.  Two classical facts:

 (i) LOG-CONVEXITY of Gamma, coordinatewise: for each i,
     ((alpha_i+beta_i)!)^2 <= (2 alpha_i)! (2 beta_i)!, with equality iff
     alpha_i = beta_i.  Multiplying over i gives C^2 <= D E, strict unless
     alpha = beta.
 (ii) AM-GM: D B^2 + A^2 E >= 2 A B sqrt(D E).

Combining, and using A,B > 0:

    Q = D B^2 + A^2 E - 2 A B C  >=  2 A B sqrt(D E) - 2 A B C
      = 2 A B ( sqrt(D E) - C )  >  0        when alpha != beta,

since C^2 < D E there.  Hence L(f^2) != 0.  QED

This is a genuine rigorous proof (no genericity), it holds in EVERY number of
variables, and it is the N = 2 base case of the first-window principle
"N slots need N moments".  This script verifies the two inequalities and the
positivity of Q exhaustively on a large range of multi-index pairs.
"""
import sys
from itertools import product
from math import factorial


def mfact(a):
    t = 1
    for e in a:
        t *= factorial(e)
    return t


def Q(alpha, beta):
    A, B = mfact(alpha), mfact(beta)
    C = mfact(tuple(alpha[i] + beta[i] for i in range(len(alpha))))
    D = mfact(tuple(2 * e for e in alpha))
    E = mfact(tuple(2 * e for e in beta))
    return D * B * B - 2 * A * B * C + A * A * E, (A, B, C, D, E)


if __name__ == "__main__":
    print("verifying  Q = D B^2 - 2 A B C + A^2 E > 0  for alpha != beta")
    print("and the two ingredients  C^2 < D E  (log-convexity) and AM-GM")
    bad = 0
    tested = 0
    for nvars in (1, 2, 3):
        for deg in (8, 6, 4):
            idx = [a for a in product(range(deg + 1), repeat=nvars)
                   if sum(a) <= deg]
            for alpha in idx:
                for beta in idx:
                    if alpha == beta:
                        continue
                    tested += 1
                    q, (A, B, C, D, E) = Q(alpha, beta)
                    if q <= 0:
                        bad += 1
                        if bad < 5:
                            print(f"   FAIL Q<=0 at {alpha},{beta}: Q={q}")
                    if C * C >= D * E:
                        bad += 1
                        if bad < 5:
                            print(f"   FAIL logconvex at {alpha},{beta}")
            break                      # one degree per nvars is plenty
    print(f"  pairs tested: {tested}   violations: {bad}")
    print()
    # equality case check
    a = (2, 1)
    q, (A, B, C, D, E) = Q(a, a)
    print(f"  equality case alpha = beta = {a}:  C^2 - D E = {C*C - D*E}  "
          f"(zero, as log-convexity says), Q = {q}")
    print()
    print("""So the two-slot case of the first-window principle is PROVED, not
merely generic: for any two distinct multi-indices in any number of variables,
L(f) = L(f^2) = 0 forces f = 0.  The proof is log-convexity of Gamma plus
AM-GM, and the strictness comes exactly from alpha != beta.""")
