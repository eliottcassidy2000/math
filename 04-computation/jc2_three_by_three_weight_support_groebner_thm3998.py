#!/usr/bin/env python3
"""Independent finite direct-source replay for THM-3998.

Unlike the all-degree theorem certificate, this script does not use the
weight-row reductions.  It substitutes generic coefficient polynomials
directly into k[x,t], expands J(A,C)-1, and asks for the exact characteristic
zero Groebner basis.  These finite slices are hostile controls for the
all-degree proof, not its logical basis.
"""

from __future__ import annotations

import sympy as sp


x, t = sp.symbols("x t")
u = x**2 * t


def jac(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, t)
                     - sp.diff(P, t) * sp.diff(Q, x))


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


def generic_poly(prefix: str, constant, degree: int):
    tail = sp.symbols(" ".join(f"{prefix}{j}" for j in range(1, degree + 1)))
    if degree == 1:
        tail = (tail,)
    return sp.expand(constant + sum(c * u**j for j, c in enumerate(tail, 1))), list(tail)


def replay(a_value: int, gamma_value: int, degree: int) -> None:
    av = sp.Rational(a_value)
    gv = sp.Rational(gamma_value)
    tag = f"a{a_value}_g{gamma_value}_D{degree}"

    f, fvars = generic_poly(f"f_{tag}_", gv**2, degree)
    F, Fvars = generic_poly(f"F_{tag}_", av, degree)
    g, gvars = generic_poly(f"g_{tag}_", gv**3, degree)
    k, kvars = generic_poly(f"k_{tag}_", sp.Rational(3, 2) * av * gv, degree)
    M, mvars = generic_poly(
        f"m_{tag}_", -sp.Rational(2, 3) / (gv * av), degree
    )
    m = sp.expand(u * (1 + u) * M)

    # The x^-2 row is polynomial because m has the exact u(1+u) factor.
    A = sp.cancel(x**2 * f + F + x**-2 * m)
    C = sp.expand(x**3 * g + x * k)
    gate(f"{tag}: A lies in k[x,t]", sp.denom(A) == 1)
    gate(f"{tag}: node A boundary", sp.expand(A.subs(t, 0) - (gv**2 * x**2 + av)) == 0)
    gate(
        f"{tag}: node C boundary",
        sp.expand(C.subs(t, 0) - (gv**3 * x**3 + sp.Rational(3, 2) * av * gv * x)) == 0,
    )
    gate(
        f"{tag}: forced A_t jet",
        sp.expand(sp.diff(A, t).subs({x: 0, t: 0})
                  + sp.Rational(2, 3) / (gv * av)) == 0,
    )

    equations = list(dict.fromkeys(sp.Poly(jac(A, C) - 1, x, t).coeffs()))
    variables = fvars + Fvars + gvars + kvars + mvars
    basis = sp.groebner(equations, *variables, order="grevlex")
    gate(f"{tag}: direct Groebner ideal is [1]", [p.as_expr() for p in basis.polys] == [1])
    print(
        f"REPLAY {tag}: variables={len(variables)} equations={len(equations)} "
        "survivors=0"
    )


print("STATUS=THM-3998 INDEPENDENT FINITE DIRECT-SOURCE REPLAY")
print("UNIVERSE=A weights {2,0,-2}; C weights {3,1}; coefficient caps shown per slice")
print("METHOD=expand J(A,C)-1 in x,t over QQ; exact Groebner basis")
gate("ambient orientation control J(x,t)=1", jac(x, t) == 1)

# One deeper slice and three parameter hostiles.  Fixing a,gamma is not a
# gauge claim; the all-parameter conclusion comes only from the companion
# typed proof.  These are deliberately independent specializations.
for case in [(1, 1, 4), (2, 1, 3), (1, 2, 3), (-1, 1, 3)]:
    replay(*case)

print("LIMITATION=finite coefficient caps and four (a,gamma) slices only")
print("THEOREM_ID=THM-3998")
print("ALL THM-3998 DIRECT GROEBNER REPLAYS PASSED")
