#!/usr/bin/env python3
"""Exact controls for THM-3131.

The theorem is local and its proof is valuation/group theoretic.  This
companion checks the two polynomial identities used by the proof, the S4
intersection ledger, both sides of the 3|v(q) boundary, and the inherited
THM-3059 hostile.
"""

from itertools import permutations

import sympy as sp


def compose(a, b):
    """Permutation a after b, represented by image tuples."""
    return tuple(a[b[i]] for i in range(4))


def power(a, n):
    out = tuple(range(4))
    for _ in range(n):
        out = compose(a, out)
    return out


def generated_by(a):
    return {power(a, n) for n in range(12)}


def ord0(expr, t):
    """t-adic order of a nonzero rational function over Q(other symbols)."""
    num, den = sp.fraction(sp.cancel(expr))
    pn = sp.Poly(num, t)
    pd = sp.Poly(den, t)
    vn = min(m[0] for m, c in pn.terms() if c != 0)
    vd = min(m[0] for m, c in pd.terms() if c != 0)
    return vn - vd


def quartic_and_resolvent_discriminants():
    T, U, p, q, r = sp.symbols("T U p q r")
    f = T**4 + p * T**2 + q * T + r
    S = U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2
    df = sp.factor(sp.discriminant(f, T))
    dS = sp.factor(sp.discriminant(S, U))
    assert sp.expand(df - dS) == 0
    return df


def group_intersections():
    s4 = set(permutations(range(4)))
    identity = tuple(range(4))
    cycle = (1, 2, 0, 3)  # (0 1 2), fixing sheet 3
    inertia = generated_by(cycle)
    v4 = {
        identity,
        (1, 0, 3, 2),
        (2, 3, 0, 1),
        (3, 2, 1, 0),
    }
    moved_stabilizer = {g for g in s4 if g[0] == 0}
    fixed_stabilizer = {g for g in s4 if g[3] == 3}
    assert len(inertia) == 3
    assert len(inertia & v4) == 1
    assert len(inertia & moved_stabilizer) == 1
    assert inertia <= fixed_stabilizer
    return (
        len(inertia & v4),
        len(inertia & moved_stabilizer),
        len(inertia & fixed_stabilizer),
    )


def nondivisible_control():
    # f=T^4-t^n T has C3 orbit {s^n,zeta*s^n,zeta^2*s^n} plus 0,
    # with t=s^3.  Use n=1 for the exact transcript.
    T, U, t = sp.symbols("T U t")
    f = T**4 - t * T
    S = U**3 - t**2
    df = sp.factor(sp.discriminant(f, T))
    dS = sp.factor(sp.discriminant(S, U))
    assert df == dS == -27 * t**4
    return (1, sp.oo, sp.oo, ord0(df, t))


def divisible_boundary_control(m=1):
    # Pair sums beta_i=s^(3m)(1+zeta^i s).  Their squared orbit has
    # q=-t^(3m)(1+t), and the discriminant gains the strict +2 collision.
    T, U, t = sp.symbols("T U t")
    p = -sp.Rational(3, 2) * t ** (2 * m)
    q = -(t ** (3 * m)) * (1 + t)
    r = t ** (4 * m) * (-sp.Rational(3, 16) + sp.Rational(3, 2) * t)
    f = sp.expand(T**4 + p * T**2 + q * T + r)
    S = sp.expand(U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2)
    df = sp.factor(sp.discriminant(f, T))
    dS = sp.factor(sp.discriminant(S, U))
    expected = -27 * t ** (12 * m + 2) * (t - 8) ** 2
    assert df == dS == expected
    n = ord0(q, t)
    vp = ord0(p, t)
    ve2 = ord0(p**2 - 4 * r, t)
    vd = ord0(df, t)
    assert n == 3 * m
    assert 3 * vp == 2 * n
    assert 3 * ve2 == 4 * n
    assert vd > 4 * n
    return (n, vp, ve2, vd)


def thm3059_hostile():
    # THM-3059: p=-v/u, q=-1/u^2, r=w/u^2 and Disc=u^-8*H,
    # with v,w,H units at the generic u=0 point.
    u, v, w, H = sp.symbols("u v w H")
    p = -v / u
    q = -1 / u**2
    r = w / u**2
    n = ord0(q, u)
    vp = ord0(p, u)
    ve2 = ord0(p**2 - 4 * r, u)
    vd = ord0(H / u**8, u)
    assert n == -2 and 3 * vp > 2 * n
    assert 3 * ve2 > 4 * n and vd == 4 * n
    return (n, vp, ve2, vd)


def main():
    df = quartic_and_resolvent_discriminants()
    intersections = group_intersections()
    strict = nondivisible_control()
    boundary = divisible_boundary_control()
    hostile = thm3059_hostile()

    print("THM-3131 exact controls")
    print("disc_identity=PASS")
    print(f"disc_formula={df}")
    print("intersection_sizes I_V4,I_Hmoved,I_Hfixed=" + ",".join(map(str, intersections)))
    print("nondivisible_n,vp,ve2,vdisc=" + ",".join(map(str, strict)))
    print("divisible_n,vp,ve2,vdisc=" + ",".join(map(str, boundary)))
    print("thm3059_n,vp,ve2,vdisc=" + ",".join(map(str, hostile)))
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
