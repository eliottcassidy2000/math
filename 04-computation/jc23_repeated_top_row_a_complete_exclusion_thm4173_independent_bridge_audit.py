#!/usr/bin/env python3
"""Clean-room local-intersection audit for THM-4173.

This referee proves the differential identity in a generic coefficient ring
and exhibits the two possible meanings of a repeated resultant root.  It is
independent of the large row-A resultant computation.
"""

import sympy as sp


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    X, T = sp.symbols("X T")
    coefficients = sp.symbols("a0:9")
    J = (
        coefficients[0]
        + coefficients[1] * X
        + coefficients[2] * T
        + coefficients[3] * X**2
        + coefficients[4] * X * T
        + coefficients[5] * T**2
        + coefficients[6] * X**3
        + coefficients[7] * X**2 * T
        + coefficients[8] * X * T**2
    )
    G = sp.expand(T * J)
    f = sp.diff(G, X) / T
    h = sp.diff(G, T)
    hessian = sp.det(sp.hessian(G, (X, T)))
    jacobian = sp.det(sp.Matrix([[sp.diff(f, X), sp.diff(f, T)],
                                 [sp.diff(h, X), sp.diff(h, T)]]))
    need(sp.factor(T * jacobian - hessian - f * sp.diff(G, X, T)) == 0,
         "generic differential identity")

    # Two reduced points with one projected coordinate: a repeated
    # resultant root need not mean a nonreduced source intersection.
    f_split = X**2 - 1
    h_split = T
    res_split = sp.resultant(f_split, h_split, X)
    need(sp.factor(res_split - T**2) == 0, "split resultant")
    split_jacobian = sp.det(sp.Matrix([[sp.diff(f_split, X), sp.diff(f_split, T)],
                                       [sp.diff(h_split, X), sp.diff(h_split, T)]]))
    need(split_jacobian.subs({X: 1, T: 0}) != 0
         and split_jacobian.subs({X: -1, T: 0}) != 0, "split points reduced")

    # One doubled point has the same projected resultant divisor, but its
    # local Jacobian vanishes.  Keller-Morse excludes precisely this branch.
    f_double = X**2
    h_double = T
    res_double = sp.resultant(f_double, h_double, X)
    need(sp.factor(res_double - T**2) == 0, "double resultant")
    double_jacobian = sp.det(sp.Matrix([[sp.diff(f_double, X), sp.diff(f_double, T)],
                                        [sp.diff(h_double, X), sp.diff(h_double, T)]]))
    need(double_jacobian.subs({X: 0, T: 0}) == 0, "doubled point singular")

    residual_length = 19 + 2
    restored_length = residual_length + 2
    need(restored_length == 23, "row-A length arithmetic")
    need(2 * 19 - restored_length - 1 + 3 == 17 < 18,
         "finite response arithmetic")
    need(2 * (25 - restored_length) == 4 < 18,
         "full response arithmetic")

    print("generic_bridge identity=0")
    print("projected_double_root split_points=2 reduced=yes")
    print("projected_double_root doubled_point=1 reduced=no")
    print("row_a_length 19+2+2=23")
    print("budgets finite=17<18 full=4<18")
    print("THM4173_INDEPENDENT_BRIDGE_ACCEPT")


if __name__ == "__main__":
    main()
