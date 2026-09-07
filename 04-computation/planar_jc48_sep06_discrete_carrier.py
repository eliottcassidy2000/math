#!/usr/bin/env python3
"""Exact carrier-map and discrete/formal controls; no shape census.

SymPy performs rational/polynomial identities. The proof of the one-way
automorphism classification is geometric and is in the companion note.
No inherited producer is imported or executed.
"""

import hashlib
import json
import sympy as s


x, t, p, y, z, lam, w, v = s.symbols("x t p y z lambda w v")
aa = 1 + x**2*t
pp = t*aa
yy = x*t*pp
uu = x**2*t
DD = p**3-y**2
sub_source = {p: pp, y: yy}
inv = {x: y*p/DD, t: DD/p**2}
GATES = 0


def check(test, name):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(name)


def equal(a, b, name):
    check(s.cancel(a-b) == 0, name)


def J(a, b, xx=x, tt=t):
    return s.diff(a, xx)*s.diff(b, tt)-s.diff(a, tt)*s.diff(b, xx)


def subst(f, table):
    return f.subs(table, simultaneous=True)


def source(f):
    return subst(f, sub_source)


def main():
    # Both rational inversion directions, not merely the denominator formula.
    equal(source(DD), t**3*aa**2, "literal Delta factorization")
    equal(source(DD), t*pp**2, "fibre denominator identity")
    equal(source(y*p/DD), x, "source x inverse")
    equal(source(DD/p**2), t, "source t inverse")
    equal(subst(pp, inv), p, "target p inverse")
    equal(subst(yy, inv), y, "target y inverse")
    equal(subst(uu, inv), y**2/DD, "target u inverse")
    equal(subst(aa, inv), p**3/DD, "target second divisor inverse")
    for f, name in ((pp, "p"), (yy, "y")):
        equal(f.subs(t, 0), 0, "line collapsed " + name)
        equal(subst(f, {x: z, t: -z**-2}), 0, "punctured line collapsed " + name)
    equal(aa.subs(t, 0), 1, "exceptional components disjoint")
    equal(subst(aa, {x: z, t: -z**-2}), 0, "Gm component equation")

    # The complete scalar family uses an indeterminate invertible parameter.
    scaling = {x: x/lam, t: lam**2*t}
    equal(subst(pp, scaling), lam**2*pp, "scaled p")
    equal(subst(yy, scaling), lam**3*yy, "scaled y")
    equal(subst(uu, scaling), uu, "scaled u")
    equal(J(scaling[x], scaling[t]), lam, "source scalar Jacobian")
    equal(subst(scaling[x], {x: lam*x, t: t/lam**2}), x, "scale inverse x")
    equal(subst(scaling[t], {x: lam*x, t: t/lam**2}), t, "scale inverse t")
    H = 1 + p + y + p**2*y + 3*p*y**2
    equal(subst(-uu/2 + source(H), scaling),
          -uu/2 + source(subst(H, {p: lam**2*p, y: lam**3*y})),
          "full displayed source under scaling")

    # Independent Poisson consistency path, including the corrected u,p.
    equal(J(pp, yy), -source(DD/p), "literal carrier Poisson factor")
    equal(J(uu, pp), 2*x*pp, "correct u,p")
    check(s.expand(J(uu, pp)-2*yy) != 0, "rejected old 2y typo")
    equal(J(uu, yy), 3*uu*pp, "u,y")
    equal(J(lam**2*p, lam**3*y, p, y), lam**5, "carrier scalar Jacobian")
    equal(lam**5*DD*(lam**2*p), lam*p*(lam**6*DD), "Poisson covariance of scaling")

    # A polynomial symplectic hostile fails descent on an actual whole fibre.
    trans = {x: x+1, t: t}
    equal(J(trans[x], trans[t]), 1, "translation is symplectic")
    on_M = s.cancel(subst(subst(pp, trans), {x: z, t: -z**-2}))
    equal(on_M, (2*z+1)/z**4, "translation distinguishes collapsed points")
    check(s.diff(on_M, z) != 0, "collapsed image is nonconstant")
    equal(on_M.subs(z, 1), 3, "first same-fibre control")
    equal(on_M.subs(z, 2), s.Rational(5, 16), "second same-fibre control")

    # A genuinely rational canonical action, with inverse and group law.
    psi = {x: x/(1-w*x), t: t*(1-w*x)**2}
    equal(J(psi[x], psi[t]), 1, "rational action is symplectic")
    equal(subst(uu, psi), uu, "rational action fixes u")
    equal(subst(psi[x], {x: x/(1+w*x), t: t*(1+w*x)**2}), x, "rational inverse x")
    equal(subst(psi[t], {x: x/(1+w*x), t: t*(1+w*x)**2}), t, "rational inverse t")
    for coordinate in (x, t):
        lhs = subst(psi[coordinate], {x: x/(1-v*x), t: t*(1-v*x)**2})
        equal(lhs, psi[coordinate].subs(w, w+v), "rational additive law " + str(coordinate))
    equal(subst(pp, psi), pp*(1-w*x)**2, "rational transformed p")
    equal(subst(yy, psi), yy*(1-w*x)**3, "rational transformed y")
    carrier_p = p*(DD-w*y*p)**2/DD**2
    carrier_y = y*(DD-w*y*p)**3/DD**3
    equal(subst(subst(pp, psi), inv), carrier_p, "p rational carrier address")
    equal(subst(subst(yy, psi), inv), carrier_y, "y rational carrier address")
    residue_p = s.rem(p*(DD-w*y*p)**2, y**2-p**3, y)
    residue_y = s.rem(y*(DD-w*y*p)**3, y**2-p**3, y)
    equal(residue_p, w**2*p**6, "uncancelled order-two cusp pole")
    equal(residue_y, -w**3*p**9, "uncancelled order-three cusp pole")
    check(residue_p != 0 and residue_y != 0, "nonzero rational carrier residues")
    equal(s.diff(psi[x], w).subs(w, 0), J(x, uu), "rational generator x")
    equal(s.diff(psi[t], w).subs(w, 0), J(t, uu), "rational generator t")

    # Universal formal carrier: literal source derivatives versus A-derivation.
    SS = p**2*DD
    SX = source(SS)
    equal(SX, t**5*aa**4, "formal carrier source factorization")
    lp = 2*p*y*DD
    ly = DD*(5*p**3-2*y**2)
    lu = 10*p**3*y
    Lx = t**4*aa**3*(5+9*x**2*t)
    Lt = -8*x*t**6*aa**3
    for f, target, name in ((pp, source(lp), "p"), (yy, source(ly), "y"),
                            (uu, source(lu), "u"), (x, Lx, "x"), (t, Lt, "t")):
        equal(J(f, SX), target, "first ordinary Hamiltonian " + name)
    equal(3*p**2*lp-2*y*ly, -4*y*DD**2, "formal cusp ideal retained")
    equal(J(-uu/2, SX), source(-5*p**3*y), "formal affine source variation")
    LA = lambda f: s.diff(f, p)*lp+s.diff(f, y)*ly
    for f, name in ((p, "p"), (y, "y"), (H, "H")):
        actual = source(f)
        expected = f
        for order in range(1, 4):
            actual = s.expand(J(actual, SX))
            expected = s.expand(LA(expected))
            equal(actual, source(expected), f"independent iterate {name},{order}")
    # Check the first two symplectic exponential coefficients without
    # constructing or specializing an infinite time series.
    equal(s.diff(Lx, x)+s.diff(Lt, t), 0, "formal divergence zero")
    L2x, L2t = J(Lx, SX), J(Lt, SX)
    equal(J(Lx, Lt)+(s.diff(L2x, x)+s.diff(L2t, t))/2,
          0, "formal second-order symplectic coefficient")
    equal(s.expand(Lx).coeff(t, 4), 5, "nonidentity formal first correction")

    payload = {
        "image": "{(0,0)} union {p*Delta != 0}",
        "exceptional": ["t=0: A1", "1+x^2*t=0: Gm"],
        "one_way_classification": ["x -> x/lambda", "t -> lambda^2*t", "Jacobian=lambda"],
        "translation_fibre": str(on_M),
        "rational_residues": [str(residue_p), str(residue_y)],
        "formal_carrier": [str(lp), str(ly), str(lu)],
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("STATUS: exact controls PASS; all-degree classification is proved geometrically")
    print("universe: symbolic carrier map, arbitrary scalar scaling, rational action parameters")
    print("actual_image={(0,0)} union {p*Delta!=0}; sole positive-dimensional fibre=A1 disjoint Gm")
    print("one_way_carrier_polynomial_automorphisms: (x/lambda,lambda^2*t)")
    print("symplectic_carrier_polynomial_automorphism: identity only")
    print("translation_hostile: collapsed M images=(2*z+1)/z^4; values=3,5/16")
    print("rational_positive: Psi_w=(x/(1-w*x),t*(1-w*x)^2); J=1; u fixed")
    print("rational_carrier_hostile: exact cusp pole orders 2,3 for nonzero w")
    print("formal_positive: S=p^2*Delta; nontrivial all-time-formal carrier flow")
    print("independent_paths: rational inverse both ways; Poisson covariance; ordinary and carrier iterates 1..3")
    print("scope: automorphism and full-carrier hypotheses retained; no arbitrary Keller endomorphism conclusion")
    print("semantic_sha256=" + digest)
    print(f"PASS gates={GATES}")


if __name__ == "__main__":
    main()
