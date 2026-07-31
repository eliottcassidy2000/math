#!/usr/bin/env python3
"""THM-3005 -- global no-return is not multiplicative, and the first-gap wall
already satisfies it.

Setting is THM-3000 section 1.  Call N RATIO-MONOTONE (equivalently: N has global
no-return) when its circuit c_k = log(R_k/R_(k-1)), k=2..d-1, is single-signed.

WHAT IS ESTABLISHED
  A. (PROVED) The ratio-monotone class is NOT closed under multiplication, and
     the reason is structural rather than accidental:
        * THM-3001: reversal negates the circuit, so N ratio-monotone implies
          N* ratio-monotone in the opposite direction;
        * the root multiset of N*N* is {r_i} union {1/r_i}, whose log-root
          measure is SYMMETRIC about 0 by construction;
        * THM-3003 section 1: a symmetric log-root measure forces
          R_k = R_(2d-k), i.e. the circuit is an ANTIPALINDROME, so it changes
          sign unless it vanishes identically.
     Hence N -> (N, N*) is a machine that turns any single ratio-monotone
     polynomial into a pair of ratio-monotone factors with a non-monotone
     product.  Minimal integer witness found by exhaustive search:
        A = (n+1)^2 (n+2),  B = A* = (n+1)(n+2)^2,  A*B = (n+1)^3 (n+2)^3,
     with circuits  '-',  '+',  and '++--'.
  B. (VERIFIED-EXACT) The first-gap wall W_M of THM-2997 eq (9) is ITSELF
     ratio-monotone -- its circuit is strictly positive with ZERO sign changes --
     for M = 6,8,10,12,14, at degrees up to 330.  Its two curvatures satisfy
     THM-3001's necessary condition C(mu) >= 0 >= C(mu*) with a growing margin:
     C(mu_W) rises from +0.137 to +0.194 while C(mu_W*) falls from -2.97 to
     -10.93.  So the wall is not the obstruction.
  C. (CONSEQUENCE) Because of A, wall-stripping is NECESSARY and not merely
     convenient: monotonicity of W_M and of N_M would NOT imply monotonicity of
     R_M = W_M N_M, and monotonicity of R_M would not imply it for N_M.
  D. (DICHOTOMY) The wall's root measure is a NEAR-CONTINUUM (roots 1/2,1,2,...,M
     with slowly varying multiplicities), the opposite extreme from THM-3004's
     well-separated clusters.  Sign-change count interpolates between 0 for a
     near-continuum and 2m-3 for m separated clusters.  Separation, not root
     count, is what creates reversals.

Reproduce: python3 04-computation/gmc_no_return_is_not_multiplicative_and_the_wall_is_monotone_thm3005.py
"""

from fractions import Fraction as Fr
from math import comb, gcd
import itertools


def circuit_signs(roots):
    """Exact sign of c_k = log(R_k/R_(k-1)) for k=2..d-1, roots given as Fractions."""
    d = len(roots)
    den = 1
    for r in roots:
        den = den * r.denominator // gcd(den, r.denominator)
    ri = [int(r * den) for r in roots]           # R_k is scale-invariant
    e = [0] * (d + 1)
    e[0] = 1
    for t in ri:
        for k in range(d, 0, -1):
            e[k] = e[k] + t * e[k - 1]
    out = []
    for k in range(2, d):
        L = e[k] ** 3 * e[k - 2] * comb(d, k - 1) ** 3 * comb(d, k + 1)
        R = e[k - 1] ** 3 * e[k + 1] * comb(d, k) ** 3 * comb(d, k - 2)
        out.append(1 if L > R else (-1 if L < R else 0))
    return out


def show(s):
    return ''.join('+' if x > 0 else ('-' if x < 0 else '0') for x in s)


def monotone(roots):
    s = [x for x in circuit_signs(roots) if x != 0]
    return len(s) == 0 or all(x > 0 for x in s) or all(x < 0 for x in s)


def sign_changes(s):
    t = [x for x in s if x != 0]
    return sum(1 for i in range(len(t) - 1) if t[i] != t[i + 1])


def curvature(roots):
    d = len(roots)
    m1 = sum(roots) / d
    m2 = sum(r * r for r in roots) / d
    m3 = sum(r ** 3 for r in roots) / d
    return float(3 * m2 ** 2 / m1 ** 4 - 2 * m3 / m1 ** 3 - 1)


def wall_roots(M):
    """THM-2997 eq (9): W_M = (n+1/2)^6 (n+1)^6 (n+2)^24 prod ... (n+M)^20."""
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


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


# --------------------------------------------------------------------------
def partA():
    rule("A. NO-RETURN IS NOT MULTIPLICATIVE -- the N times N* machine")
    print("  THM-3001: reversal negates the circuit.  THM-3003 sec 1: a SYMMETRIC")
    print("  log-root measure forces R_k = R_(D-k), an antipalindromic circuit.")
    print("  The root multiset of N*N* is {r} union {1/r}, whose log-root measure is")
    print("  symmetric about 0 BY CONSTRUCTION.  So for every N, the product of the")
    print("  two ratio-monotone polynomials N and N* is ratio-PALINDROMIC, hence not")
    print("  monotone unless its circuit vanishes identically.")
    print()
    ok = True
    print("  The ANTIPALINDROME is universal (it needs no hypothesis on N).  Whether")
    print("  BOTH FACTORS are monotone is an extra property of the chosen N, and is")
    print("  what turns the identity into a counterexample to multiplicativity.")
    print("  machine applied to several N (A = N, B = N* up to scaling):")
    for name, r in [("(n+1)^2(n+2)", [Fr(1), Fr(1), Fr(2)]),
                    ("(n+1)^3(n+5)", [Fr(1)] * 3 + [Fr(5)]),
                    ("(n+1)(n+2)(n+7)", [Fr(1), Fr(2), Fr(7)]),
                    ("(n+2)^2(n+9)^3", [Fr(2)] * 2 + [Fr(9)] * 3)]:
        rs = [Fr(1) / x for x in r]
        prod = r + rs
        sa, sb, sp = circuit_signs(r), circuit_signs(rs), circuit_signs(prod)
        anti = all(sp[i] == -sp[len(sp) - 1 - i] for i in range(len(sp)))
        # universal claim: antipalindromic product, hence never monotone
        ok &= anti and (not monotone(prod))
        both = monotone(r) and monotone(rs)
        tag = "counterexample" if both else "factors not both monotone (identity still holds)"
        print(f"    N = {name:18s} c(N)={show(sa):8s} c(N*)={show(sb):8s}"
              f" c(N N*)={show(sp):14s} antipalindrome={anti}"
              f" product monotone={monotone(prod)}  [{tag}]")
    print()
    print("  MINIMAL INTEGER WITNESS (exhaustive over integer roots <=6, degree <=5 each):")
    pool = [1, 2, 3, 4, 5, 6]
    cands = []
    for d in (3, 4, 5):
        for r in itertools.combinations_with_replacement(pool, d):
            if len(set(r)) < 2:
                continue
            rr = [Fr(x) for x in r]
            if monotone(rr):
                cands.append(rr)
    best = None
    for A in cands:
        for B in cands:
            if len(A) + len(B) > 8:
                continue
            if not monotone(A + B):
                if best is None or len(A) + len(B) < best[0]:
                    best = (len(A) + len(B), A, B)
    tot, A, B = best
    print(f"    total degree {tot}:  A roots {[str(x) for x in A]},  B roots {[str(x) for x in B]}")
    print(f"      c(A)   = {show(circuit_signs(A))}   monotone={monotone(A)}")
    print(f"      c(B)   = {show(circuit_signs(B))}   monotone={monotone(B)}")
    print(f"      c(A*B) = {show(circuit_signs(A + B))}   monotone={monotone(A + B)}")
    print("    B is exactly A's reversal up to scaling (A roots {1,1,2}, reciprocals")
    print("    {1,1,1/2} = (1/2){2,2,1} = B roots), so this is the machine at its")
    print("    smallest, not a sporadic coincidence.")
    ok &= (not monotone(A + B))
    print(f"  VERDICT A: {'NOT CLOSED UNDER MULTIPLICATION (proved + minimal witness)' if ok else 'FAILED'}")
    return ok


def partB():
    print()
    rule("B. THE FIRST-GAP WALL IS ITSELF RATIO-MONOTONE")
    ok = True
    print("   M   deg(W)  distinct scales  sign changes  monotone   C(mu_W)     C(mu_W*)   C>=0>=C*")
    for M in (6, 8, 10, 12, 14):
        W = wall_roots(M)
        Wi = [Fr(1) / r for r in W]
        s = circuit_signs(W)
        nc = sign_changes(s)
        mo = monotone(W)
        c1, c2 = curvature(W), curvature(Wi)
        nec = (c1 >= 0 >= c2)
        ok &= mo and (nc == 0) and nec
        print(f"  {M:3d}  {len(W):6d}  {len(set(W)):^15d}  {nc:^12d}  {str(mo):8s}"
              f" {c1:+9.6f}  {c2:+10.6f}   {nec}")
    print("  The wall's circuit is strictly POSITIVE throughout: zero reversals at")
    print("  degree 330.  THM-3001's necessary condition holds with a margin that")
    print("  GROWS in M on both sides.  The wall is not the obstruction.")
    print(f"  VERDICT B: {'WALL IS MONOTONE AND CLEARS THE NECESSARY CONDITION' if ok else 'FAILED'}")
    return ok


def partC():
    print()
    rule("C. DICHOTOMY: separation, not root count, creates reversals")
    print("  same number of distinct scales, two extremes:")
    for M in (10,):
        W = wall_roots(M)
        m = len(set(W))
        print(f"    near-continuum (the wall, M={M}): {m} distinct scales spread over")
        print(f"      [1/2, {M}] with ratios -> 1 :  sign changes = {sign_changes(circuit_signs(W))}")
    roots = []
    for i in range(11):
        roots += [Fr(10) ** (6 * i)] * 2
    print(f"    well separated: {len(set(roots))} distinct scales with ratio 10^6 :"
          f"  sign changes = {sign_changes(circuit_signs(roots))}"
          f"  (2m-3 = {2 * len(set(roots)) - 3})")
    print("  Same scale COUNT, opposite behaviour.  THM-3004's cluster law needs")
    print("  separation; a near-continuum root measure sits at the other extreme with")
    print("  no reversals at all.")
    return True


def main():
    a = partA()
    b = partB()
    c = partC()
    print()
    rule(f"SUMMARY  non-multiplicative={a}  wall-monotone={b}  dichotomy={c}")


if __name__ == "__main__":
    main()
