#!/usr/bin/env python3
"""Exact companion for THM-3063.

The calculation checks that a branch-cofactor vertex gauge reaches the three
quartic matching monomials only through the product of the four cofactors.
It then gives an inertia-equivariant product-one twist which preserves the
entire matching triple but destroys sheetwise Keller constancy.  A residue-7
C3 face makes the minimum face and all initial residues explicit.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def inv_mod(value: int, prime: int) -> int:
    value %= prime
    require(value != 0, "inverse requested at zero")
    return pow(value, -1, prime)


def main() -> None:
    z = sp.symbols("z0:4")
    c = sp.symbols("c0:4", nonzero=True)
    d = sp.symbols("d0:3", nonzero=True)

    def edge(i: int, j: int) -> sp.Expr:
        return z[i] - z[j]

    p1 = edge(0, 1) * edge(2, 3)
    p2 = edge(0, 2) * edge(1, 3)
    p3 = edge(0, 3) * edge(1, 2)
    require(sp.expand(p1 - p2 + p3) == 0, "oriented Pluecker identity")

    cofactor_product = sp.prod(c)
    gauged_edges = {
        (i, j): c[i] * c[j] * edge(i, j)
        for i in range(4)
        for j in range(i + 1, 4)
    }
    a1 = gauged_edges[(0, 1)] * gauged_edges[(2, 3)]
    a2 = -gauged_edges[(0, 2)] * gauged_edges[(1, 3)]
    a3 = gauged_edges[(0, 3)] * gauged_edges[(1, 2)]
    require(sp.expand(a1 - cofactor_product * p1) == 0, "first matching gauge")
    require(sp.expand(a2 + cofactor_product * p2) == 0, "second matching gauge")
    require(sp.expand(a3 - cofactor_product * p3) == 0, "third matching gauge")
    require(sp.expand(a1 + a2 + a3) == 0, "gauged Pluecker contraction")

    # The full product-one torus is invisible to all three matchings.
    twist = (d[0], d[1], d[2], 1 / (d[0] * d[1] * d[2]))
    require(sp.simplify(sp.prod(twist) - 1) == 0, "product-one twist")
    twisted_product = sp.prod(c[i] * twist[i] for i in range(4))
    require(sp.simplify(twisted_product - cofactor_product) == 0, "twist changed cofactor product")

    derivative_stars = tuple(sp.prod(edge(i, j) for j in range(4) if j != i) for i in range(4))
    keller_cofactors = tuple(1 / derivative_stars[i] for i in range(4))
    keller_values = tuple(sp.simplify(keller_cofactors[i] * derivative_stars[i]) for i in range(4))
    require(keller_values == (1, 1, 1, 1), "abstract Keller packet")
    twisted_values = tuple(
        sp.simplify(twist[i] * keller_cofactors[i] * derivative_stars[i]) for i in range(4)
    )
    require(twisted_values == twist, "twisted Jacobian packet")

    # Exact C3-compatible twist: constant t on the cubic orbit and t^-3 on
    # the fixed factor.  It has norm one but breaks orbit/fixed equality.
    t = Fraction(2, 1)
    c3_twist = (t, t, t, t ** -3)
    require(c3_twist == (Fraction(2), Fraction(2), Fraction(2), Fraction(1, 8)), "C3 twist values")
    require(sp.prod(c3_twist) == 1, "C3 twist norm")
    require(c3_twist[0] != c3_twist[3], "C3 twist must break sheetwise equality")

    # Residue-7 planar C3 face.  For M=-X^4+X^3+u with u=s^3, the three
    # small leading roots solve a^3=-1 and the fixed root has residue 1.
    q = 7
    small = (3, 5, 6)
    fixed = 1
    require(all(pow(value, 3, q) == q - 1 for value in small), "C3 residual roots")
    require(len(set(small)) == 3, "C3 residual roots are not distinct")

    signed_matching_residues = (
        (small[0] - small[1]) * (0 - fixed),
        -((small[0] - small[2]) * (0 - fixed)),
        (0 - fixed) * (small[1] - small[2]),
    )
    signed_matching_residues = tuple(value % q for value in signed_matching_residues)
    require(signed_matching_residues == (2, 4, 1), "C3 matching initial residues")
    require(sum(signed_matching_residues) % q == 0, "C3 initial face must cancel")

    edge_valuations = {
        (0, 1): 1,
        (0, 2): 1,
        (0, 3): 0,
        (1, 2): 1,
        (1, 3): 0,
        (2, 3): 0,
    }
    matching_valuations = (
        edge_valuations[(0, 1)] + edge_valuations[(2, 3)],
        edge_valuations[(0, 2)] + edge_valuations[(1, 3)],
        edge_valuations[(0, 3)] + edge_valuations[(1, 2)],
    )
    require(matching_valuations == (1, 1, 1), "C3 matching face valuations")

    # Derivative residues: small M'(sa)/s^2 -> 3a^2, fixed M'(1) -> -1.
    derivative_residues = tuple((3 * value * value) % q for value in small) + ((-1) % q,)
    require(derivative_residues == (6, 5, 3, 6), "C3 derivative residues")
    derivative_valuations = (2, 2, 2, 0)
    cofactor_valuations = tuple(-value for value in derivative_valuations)
    require(cofactor_valuations == (-2, -2, -2, 0), "Keller cofactor valuations")
    cofactor_residues = tuple(inv_mod(value, q) for value in derivative_residues)
    require(cofactor_residues == (6, 3, 5, 6), "Keller cofactor residues")
    require(
        tuple(a * b % q for a, b in zip(derivative_residues, cofactor_residues)) == (1, 1, 1, 1),
        "Keller residue equality",
    )
    require(sp.prod(cofactor_residues) % q == 1, "cofactor product residue")

    cofactor_product_valuation = sum(cofactor_valuations)
    gauged_matching_valuations = tuple(value + cofactor_product_valuation for value in matching_valuations)
    require(cofactor_product_valuation == -6, "cofactor product valuation")
    require(gauged_matching_valuations == (-5, -5, -5), "gauged C3 minimum face")

    twist_residues = tuple(int(value.numerator * inv_mod(value.denominator, q)) % q for value in c3_twist)
    require(twist_residues == (2, 2, 2, 1), "C3 twist residues")
    require(sp.prod(twist_residues) % q == 1, "C3 twist residue norm")
    broken_jacobian_residues = tuple(
        twist_residues[i] * cofactor_residues[i] * derivative_residues[i] % q
        for i in range(4)
    )
    require(broken_jacobian_residues == (2, 2, 2, 1), "twist did not break Keller residues")

    print("theorem=THM-3063")
    print("status=PROVED_VERIFIED_EXACT")
    print("contraction=(c0c1c2c3)*(P1,-P2,P3);sum=0_by_Pluecker")
    print("cofactor_quotient=matching_triple_depends_only_on_product_c0c1c2c3")
    print("Keller_predicate=c_i*product_(j!=i)(z_i-z_j)_constant_for_i=0..3")
    print("invisible_group={(d0,d1,d2,d3):product_di=1}")
    print("C3_invariant_twist=(2,2,2,1/8);product=1;orbit_vs_fixed_not_equal")
    print("C3_edge_order=(01,02,03,12,13,23);valuations=(1,1,0,1,0,0);matching_valuations=(1,1,1)")
    print("C3_signed_initial_residues_F7=(2,4,1);sigma=0")
    print("C3_derivative_valuations=(2,2,2,0);cofactor_valuations=(-2,-2,-2,0)")
    print("C3_Keller_residues=(1,1,1,1);twisted_residues=(2,2,2,1)")
    print("gauged_matching_face=(-5,-5,-5);twist_preserves_full_matching_triple")
    print("conclusion=THM3058_vertex_gauge_augmentation_cannot_test_sheetwise_Keller_constancy")


if __name__ == "__main__":
    main()
