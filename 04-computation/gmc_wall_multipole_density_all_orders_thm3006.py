#!/usr/bin/env python3
"""THM-3006 -- the first-gap wall is a four-band charge density, and its multipole
moments are known to ALL orders.

THM-2997 (24) states four wall-moment limits (w_0..w_3).  In the multipole
language of THM-3003 section 3 the wall W_M of THM-2997 eq (9) is a charge
distribution on [0,M]; its band profile is piecewise constant, so its moments are
the moments of an explicit density and are available at every order.

BAND PROFILE from eq (9), as a density in x = s/M:
    multiplicity 28 on (0, 1/3],  26 on (1/3, 1/2],  24 on (1/2, 2/3],
    23 on (2/3, 1),  plus O(1) endpoint atoms (the (n+1/2)^6, (n+1)^6,
    (n+2)^24 head and the (n+M)^20 tail, and the mod-10 root), which do not
    affect the leading term.

THEOREM (all orders).
    lim_{M->inf} w_k / M^(k+1)
      = [ 23 + (2/3)^(k+1) + 2 (1/2)^(k+1) + 2 (1/3)^(k+1) ] / (k+1).

PROOF.  w_k/M^(k+1) -> int_0^1 rho(x) x^k dx with rho the step density above:
    28/(k+1) (1/3)^(k+1)
  + 26/(k+1) [ (1/2)^(k+1) - (1/3)^(k+1) ]
  + 24/(k+1) [ (2/3)^(k+1) - (1/2)^(k+1) ]
  + 23/(k+1) [ 1          - (2/3)^(k+1) ]
  = [ 23 + (24-23)(2/3)^(k+1) + (26-24)(1/2)^(k+1) + (28-26)(1/3)^(k+1) ]/(k+1),
which is the stated formula: the coefficients 1, 2, 2 are exactly the DOWNWARD
JUMPS of the multiplicity profile at x = 2/3, 1/2, 1/3.  QED.

CONTROLS.  k = 0,1,2,3 reproduce THM-2997 (24)'s 76/3, 145/12, 2551/324,
1681/288 exactly.  Direct summation over eq (9) at M = 600, 1200, 2400 converges
to the predicted limits at every k <= 6.

NEW (extends (24) by three indices):
    w_4/M^5 -> 90211/19440,  w_5/M^6 -> 179795/46656,
    w_6/M^7 -> 3229771/979776.

SCOPE.  This is a statement about the WALL only, which is explicit.  It supplies
one half of THM-3000 section 7 / THM-3003 section 3's third-edge invoice; the
other half (the resultant jet P_4, or equivalently the core's own fourth power
sum) is NOT supplied here and remains the open obligation.  Independently, the
root-modulus route is now known to be unreachable from the three currently
published jets: a Hurwitz positive-coefficient polynomial exists with the same
degree and same first three power sums whose largest root modulus is Theta(M^2).

Reproduce: python3 04-computation/gmc_wall_multipole_density_all_orders_thm3006.py
"""

from fractions import Fraction as Fr


def wall_roots(M):
    a, b, c = M // 3, M // 2, (2 * M) // 3
    R = [Fr(1, 2)] * 6 + [Fr(1)] * 6 + [Fr(2)] * 24
    for s in range(3, a + 1):
        R += [Fr(s)] * 28
    for s in range(a + 1, b + 1):
        R += [Fr(s)] * 26
    for s in range(b + 1, c + 1):
        R += [Fr(s)] * 24
    for s in range(c + 1, M):
        R += [Fr(s)] * 23
    R += [Fr(M)] * 20
    if M % 10 == 1:
        R += [Fr(4 * M + 1, 5)]
    return R


def law(k):
    return Fr(23 + Fr(2, 3) ** (k + 1) + 2 * Fr(1, 2) ** (k + 1)
              + 2 * Fr(1, 3) ** (k + 1), k + 1)


KNOWN = {0: Fr(76, 3), 1: Fr(145, 12), 2: Fr(2551, 324), 3: Fr(1681, 288)}


def main():
    print("=" * 74)
    print("THM-3006  wall multipole moments to all orders")
    print("=" * 74)
    print("  lim w_k/M^(k+1) = [23 + (2/3)^(k+1) + 2(1/2)^(k+1) + 2(1/3)^(k+1)]/(k+1)")
    print("  (coefficients 1,2,2 are the downward jumps of the band profile at 2/3,1/2,1/3)")
    print()
    print("  k   law                    THM-2997 (24)          match    direct sum w_k/M^(k+1)")
    ok = True
    Ms = (600, 1200, 2400)
    walls = {M: wall_roots(M) for M in Ms}
    for k in range(0, 7):
        L = law(k)
        vals = []
        for M in Ms:
            wk = sum(r ** k for r in walls[M])
            vals.append(float(Fr(wk, 1) / Fr(M) ** (k + 1)))
        kn = KNOWN.get(k)
        m = (kn is None) or (L == kn)
        ok &= m
        # monotone convergence toward L
        conv = abs(vals[-1] - float(L)) < abs(vals[0] - float(L))
        ok &= conv
        print(f"  {k}   {str(L):22s} {str(kn) if kn else '(new)':22s} {str(m):7s}"
              f" {vals[0]:.6f} {vals[1]:.6f} {vals[2]:.6f} -> {float(L):.6f}  conv={conv}")
    print()
    print("  NEW DELIVERABLES beyond THM-2997 (24):")
    for k in (4, 5, 6):
        print(f"    w_{k}/M^{k + 1} -> {law(k)} = {float(law(k)):.9f}")
    print()
    print("=" * 74)
    print(f"SUMMARY  all four (24) constants reproduced and k<=6 converging: {ok}")
    print("=" * 74)


if __name__ == "__main__":
    main()
