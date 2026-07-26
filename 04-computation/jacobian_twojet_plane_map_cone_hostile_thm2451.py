#!/usr/bin/env python3
"""Exact symbolic companion for THM-2451."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bracket(u: sp.Matrix, v: sp.Matrix, w: sp.Matrix) -> sp.Expr:
    return sp.expand(sp.Matrix.hstack(u, v, w).det())


def main() -> None:
    x, y, z = sp.symbols("x y z")
    a = sp.Matrix(sp.symbols("a0:3"))
    b = sp.Matrix(sp.symbols("b0:3"))
    n = sp.Matrix(sp.symbols("n0:3"))
    zero = sp.zeros(3, 1)

    A = a + x * n
    B = b + y * n
    C = zero

    Ax, Ay = A.diff(x), A.diff(y)
    Bx, By = B.diff(x), B.diff(y)
    Cx, Cy = C.diff(x), C.diff(y)

    D5 = 2 * bracket(Ax, Ay, A)
    D4 = (
        bracket(Ax, Ay, B)
        + 2 * bracket(Ax, By, A)
        + 2 * bracket(Bx, Ay, A)
    )
    D3 = (
        bracket(Ax, By, B)
        + bracket(Bx, Ay, B)
        + 2
        * (
            bracket(Ax, Cy, A)
            + bracket(Cx, Ay, A)
            + bracket(Bx, By, A)
        )
    )
    D2 = (
        bracket(Bx, By, B)
        + bracket(Ax, Cy, B)
        + bracket(Cx, Ay, B)
        + 2 * (bracket(Bx, Cy, A) + bracket(Cx, By, A))
    )
    D1 = (
        bracket(Bx, Cy, B)
        + bracket(Cx, By, B)
        + 2 * bracket(Cx, Cy, A)
    )
    D0 = bracket(Cx, Cy, B)
    graded = (D5, D4, D3, D2, D1, D0)
    require(all(sp.expand(value) == 0 for value in graded), "graded zero failed")

    F = A * z**2 + B * z + C
    direct_jacobian = sp.factor(F.jacobian((x, y, z)).det())
    require(direct_jacobian == 0, "direct Jacobian is not zero")

    w = A.cross(B)
    w_volume = sp.factor(bracket(w.diff(x), w.diff(y), w))
    delta = bracket(a, b, n)
    require(sp.factor(w_volume - delta**2) == 0, "cross-volume identity failed")

    witness = {
        a[0]: 1,
        a[1]: 0,
        a[2]: 0,
        b[0]: 0,
        b[1]: 1,
        b[2]: 0,
        n[0]: 0,
        n[1]: 0,
        n[2]: 1,
    }
    witness_A = tuple(sp.expand(value.subs(witness)) for value in A)
    witness_B = tuple(sp.expand(value.subs(witness)) for value in B)
    witness_w = tuple(sp.expand(value.subs(witness)) for value in w)
    witness_volume = sp.expand(w_volume.subs(witness))
    require(witness_A == (1, 0, x), "witness A changed")
    require(witness_B == (0, 1, y), "witness B changed")
    require(witness_w == (-x, -y, 1), "witness w changed")
    require(witness_volume == 1, "witness volume changed")

    print("THM-2451 exact symbolic companion")
    print("graded_brackets_D5_to_D0=0,0,0,0,0,0")
    print("direct_jacobian=0")
    print("plane_map_volume=[a,b,n]^2")
    print("integral_witness_A=(1,0,x)")
    print("integral_witness_B=(0,1,y)")
    print("integral_witness_w=(-x,-y,1)")
    print("integral_witness_volume=1")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
