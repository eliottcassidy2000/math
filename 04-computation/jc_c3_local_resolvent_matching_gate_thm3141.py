#!/usr/bin/env python3
"""Exact controls for THM-3141 and its graph-quartic saturation addendum.

The theorem is local and its proof is valuation/group theoretic.  This
companion checks the two polynomial identities used by the proof, the S4
intersection ledger, both sides of the 3|v(q) boundary, and the inherited
THM-3059 hostile.  It also checks the exact depression formula, both leading
graph-root lanes, normal-parameter gauge invariance, and the terminal-scale
cancellation which identifies the missing prefactor cubeclass.
"""

from itertools import permutations

import sympy as sp


def require(condition, label):
    """Keep every exact check live under both normal Python and python -O."""
    if not bool(condition):
        raise AssertionError(label)


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
    require(sp.expand(df - dS) == 0, "quartic/resolvent discriminants")
    return df


def depressed_graph_coefficient():
    """Check q after depression, including resultant scaling and a hostile."""
    X, T = sp.symbols("X T")
    a3, a2, a1, a0 = sp.symbols("a3 a2 a1 a0")
    f = X**4 + a3 * X**3 + a2 * X**2 + a1 * X + a0
    depressed = sp.expand(f.subs(X, T - a3 / 4))
    q_dep = sp.expand(depressed.coeff(T, 1))
    expected = a1 - a2 * a3 / 2 + a3**3 / 8
    require(sp.expand(q_dep - expected) == 0, "depressed linear coefficient")

    c, r3, r2, r1 = sp.symbols("c r3 r2 r1", nonzero=True)
    scaled = sp.factor(expected.subs({a3: r3 / c, a2: r2 / c, a1: r1 / c}))
    expected_scaled = (8 * c**2 * r1 - 4 * c * r2 * r3 + r3**3) / (8 * c**3)
    require(sp.cancel(scaled - expected_scaled) == 0, "resultant-scaled q")

    u, t = sp.symbols("u t", nonzero=True)
    hostile = sp.expand((X - u) * (X**3 - t**-1))
    hostile_a3 = hostile.coeff(X, 3)
    hostile_a2 = hostile.coeff(X, 2)
    hostile_a1 = hostile.coeff(X, 1)
    hostile_q = sp.factor(
        hostile_a1 - hostile_a2 * hostile_a3 / 2 + hostile_a3**3 / 8
    )
    require(
        sp.cancel(hostile_q - (-1 / t - u**3 / 8)) == 0,
        "fixed-plus-cubic depressed q",
    )
    return expected, expected_scaled, hostile_q


def graph_leading_lanes():
    """Depress the four leading roots and recover c_m=-1 or 1/8 exactly."""
    A = sp.symbols("A", nonzero=True)
    zeta = -sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2

    constants = {}
    for residue in (1, 2):
        moved = [sp.expand(A * zeta ** (-i * residue)) for i in range(3)]
        trace = sp.simplify(sum(moved))
        require(trace == 0, f"moved trace cancellation p={residue} mod 3")
        beta = [sp.simplify(root - trace / 2) for root in moved]
        q_lead = sp.simplify(-sp.prod(beta))
        require(sp.simplify(q_lead / A**3) == -1, f"q leading p={residue}")
        constants[residue] = sp.simplify(q_lead / A**3)

    moved = [A, A, A]
    trace = sum(moved)
    beta = [sp.simplify(root - trace / 2) for root in moved]
    require(all(entry == -A / 2 for entry in beta), "divisible depressed pair sums")
    q_lead = sp.simplify(-sp.prod(beta))
    require(q_lead == A**3 / 8, "q leading p=0 mod 3")
    constants[0] = sp.simplify(q_lead / A**3)
    return constants


def graph_gauge_and_terminal_scale():
    """Check A^3 tau^p gauge invariance and rho cancellation for controls."""
    A, tau, h, rho, H, K = sp.symbols("A tau h rho H K", nonzero=True)
    for pole in (1, 2, 3, 4, 6):
        c_p = sp.Rational(1, 8) if pole % 3 == 0 else -1
        q0 = c_p * A**3 * tau**pole
        transformed = c_p * (A * h**pole) ** 3 * (tau * h**-3) ** pole
        require(sp.cancel(q0 - transformed) == 0, f"normal gauge p={pole}")

        terminal = sp.cancel(
            q0.subs({A: rho**(-pole) * H, tau: rho**3 * K})
        )
        require(
            terminal == c_p * H**3 * K**pole,
            f"terminal rho cancellation p={pole}",
        )
    return "m=1,2,3,4,6"


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
    require(len(inertia) == 3, "C3 inertia order")
    require(len(inertia & v4) == 1, "C3 intersects V4 trivially")
    require(
        len(inertia & moved_stabilizer) == 1,
        "C3 intersects moved stabilizer trivially",
    )
    require(inertia <= fixed_stabilizer, "C3 fixes the fourth sheet")
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
    require(df == dS == -27 * t**4, "nondivisible discriminant control")
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
    require(df == dS == expected, "divisible discriminant control")
    n = ord0(q, t)
    vp = ord0(p, t)
    ve2 = ord0(p**2 - 4 * r, t)
    vd = ord0(df, t)
    require(n == 3 * m, "divisible q value")
    require(3 * vp == 2 * n, "divisible p value")
    require(3 * ve2 == 4 * n, "divisible second resolvent value")
    require(vd > 4 * n, "divisible strict discriminant value")
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
    require(n == -2 and 3 * vp > 2 * n, "THM-3059 first Newton gate")
    require(
        3 * ve2 > 4 * n and vd == 4 * n,
        "THM-3059 second Newton/discriminant gate",
    )
    return (n, vp, ve2, vd)


def main():
    df = quartic_and_resolvent_discriminants()
    q_dep, q_scaled, hostile_q = depressed_graph_coefficient()
    graph_constants = graph_leading_lanes()
    gauge_controls = graph_gauge_and_terminal_scale()
    intersections = group_intersections()
    strict = nondivisible_control()
    boundary = divisible_boundary_control()
    hostile = thm3059_hostile()

    print("THM-3141 exact controls")
    print("disc_identity=PASS")
    print(f"disc_formula={df}")
    print(f"depressed_q={q_dep}")
    print(f"resultant_scaled_q={q_scaled}")
    print(f"fixed_plus_cubic_q={hostile_q}")
    print(
        "graph_constants mmod3=1,2,0="
        + ",".join(str(graph_constants[i]) for i in (1, 2, 0))
    )
    print(f"graph_gauge_terminal_scale_controls={gauge_controls}:PASS")
    print("intersection_sizes I_V4,I_Hmoved,I_Hfixed=" + ",".join(map(str, intersections)))
    print("nondivisible_n,vp,ve2,vdisc=" + ",".join(map(str, strict)))
    print("divisible_n,vp,ve2,vdisc=" + ",".join(map(str, boundary)))
    print("thm3059_n,vp,ve2,vdisc=" + ",".join(map(str, hostile)))
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
