#!/usr/bin/env python3
"""HYP-8070 referee — the Jacobian-conjecture counterexample, verified exactly
(mac-mini-2026-07-19-S129; owner-supplied map, provenance unconfirmed by two
web searches at session time — the mathematics below stands on its own).

THE MAP  F: C^3 -> C^3,
  F1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
  F2 = y + 3x(1+xy)^2 z + 3xy^2 (4+3xy)
  F3 = 2x - 3x^2 y - x^3 z          (total degree 7)

CLAIMS VERIFIED (all exact, sympy):
  (1) det JF = -2 identically (constant, nonzero: a Keller map);
  (2) the fiber over (-1/4, 0, 0) is EXACTLY
      {(0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2)}  (Groebner-complete);
  hence F is a non-injective Keller map: THE JACOBIAN CONJECTURE IS FALSE
  in dimension 3 (and every n >= 3 by adding identity coordinates), and via
  DC_n => JC_n the DIXMIER CONJECTURE IS FALSE for n >= 3.  Survivors:
  JC_2, DC_1 (Dixmier's original A_1 question), DC_2.

ANATOMY (verified below): F intertwines the involution
  sigma(x,y,z) = (-x,-y,z)   with   tau(u,v,w) = (u,-v,-w):
  F o sigma = tau o F.  Fibers over tau-fixed targets (v=w=0) are
  sigma-stable, so their cardinality is congruent mod 2 to the number of
  preimages on the sigma-fixed axis x=y=0.  On that axis F(0,0,z)=(z,0,0)
  is linear: exactly ONE fixed-axis preimage => every tau-fixed fiber has
  ODD cardinality.  The witnessed fiber = 1 fixed point + 1 two-orbit —
  the non-injectivity certificate is organized by an involution with a
  one-point fixed layer: the Redei pairing mechanism, in the wild.
"""

import sympy as sp

x, y, z = sp.symbols('x y z')

F1 = (1 + x*y)**3 * z + y**2 * (1 + x*y) * (4 + 3*x*y)
F2 = y + 3*x*(1 + x*y)**2 * z + 3*x*y**2 * (4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
F = (F1, F2, F3)


def main():
    J = sp.Matrix([[sp.diff(f, v) for v in (x, y, z)] for f in F])
    det = sp.expand(J.det())
    assert det == -2, det
    print("(1) det JF = -2 identically  VERIFIED")

    deg = max(sp.total_degree(sp.expand(f)) for f in F)
    print(f"    total degree of the map: {deg}")

    eqs = [sp.expand(F1 + sp.Rational(1, 4)), sp.expand(F2), sp.expand(F3)]
    sols = sp.solve(eqs, [x, y, z], dict=True)
    got = sorted(tuple(sp.nsimplify(s[v]) for v in (x, y, z)) for s in sols)
    want = sorted([(0, 0, sp.Rational(-1, 4)),
                   (1, sp.Rational(-3, 2), sp.Rational(13, 2)),
                   (-1, sp.Rational(3, 2), sp.Rational(13, 2))])
    assert got == want, got
    print("(2) fiber over (-1/4,0,0) = exactly the three claimed points "
          "(Groebner-complete)  VERIFIED")
    print("    => non-injective Keller map: JC(3) FALSE; JC(n) false for all "
          "n >= 3; Dixmier DC(n) false for n >= 3 (DC_n => JC_n).")

    # equivariance
    sub = {x: -x, y: -y, z: z}
    assert sp.expand(F1.subs(sub) - F1) == 0
    assert sp.expand(F2.subs(sub) + F2) == 0
    assert sp.expand(F3.subs(sub) + F3) == 0
    print("(3) equivariance F(-x,-y,z) = (F1, -F2, -F3)  VERIFIED")
    print("    fixed axis maps linearly: F(0,0,z) =",
          tuple(sp.simplify(f.subs({x: 0, y: 0})) for f in F))
    print("    => every fiber over (u,0,0) has ODD cardinality (Redei-style")
    print("    parity: sigma-stable fiber, one fixed-axis preimage).")


if __name__ == "__main__":
    main()
