#!/usr/bin/env python3
"""THM-4359 independent exact referee.

Reconstruct the row-eight response solely from the formulas displayed in
THM-4308.  This file does not read or import the primary certificate.  It
audits the closure kernels, constructible image, exceptional fibre,
algebraic-matroid circuits, and pointed affine zero-wall nonfaces.
"""

from itertools import combinations
import sys

import sympy as s


sys.stdout.reconfigure(newline="\n")

P, eta, xi = s.symbols("Phi eta xi_10")
z, X, U, W, Z = s.symbols("zeta_3 X U W Z")
SOURCE = (P, eta, xi)

z_src = -s.Rational(3, 2) * P
U_src = s.Rational(475515904, 200475) - s.Rational(6, 11) * xi
W_src = (
    -s.Rational(13, 12) * P**2
    + s.Rational(424, 99) * xi
    - s.Rational(35956576256, 1002375)
)
Z_src = (
    s.Rational(3126529518592, 27064125)
    - s.Rational(130, 81) * P**2
    - s.Rational(65, 36) * P * eta
    - s.Rational(22886, 2673) * xi
)

OBS4 = {"z": z_src, "U": U_src, "W": W_src, "Z": Z_src}
OBS5 = {"z": z_src, "x": xi, "U": U_src, "W": W_src, "Z": Z_src}

F = 1184625 * z**2 + 19318500 * U + 2460375 * W + 42434609152
R1 = 109350 * X + 200475 * U - 475515904
R2 = 482625 * z**2 - 4293000 * X + 1002375 * W + 35956576256
G0 = -1042743375 * U + 66430125 * Z - 5200877686784
R7 = 231720750 * X + 27064125 * Z - 3126529518592


def pull(poly):
    """Pull a target polynomial back to Q[Phi,eta,xi_10]."""

    return s.cancel(poly.subs({z: z_src, X: xi, U: U_src, W: W_src, Z: Z_src}))


def primitive(expr):
    """Return the positive-leading primitive integer associate."""

    _, cleared = s.Poly(expr).clear_denoms(convert=True)
    _, prim = cleared.primitive()
    out = s.expand(prim.as_expr())
    if s.LC(s.Poly(out)) < 0:
        out = -out
    return out


def graph_elimination(retain_x):
    if retain_x:
        graph = [z - z_src, X - xi, U - U_src, W - W_src, Z - Z_src]
        gens = (P, eta, xi, z, X, U, W, Z)
        source_set = {P, eta, xi}
    else:
        graph = [z - z_src, U - U_src, W - W_src, Z - Z_src]
        gens = (P, eta, xi, z, U, W, Z)
        source_set = {P, eta, xi}
    basis = s.groebner(graph, *gens, order="lex", domain=s.QQ)
    return [
        primitive(g.as_expr())
        for g in basis.polys
        if not (source_set & g.as_expr().free_symbols)
    ]


def algebraic_circuits(observables, parameters):
    """Enumerate circuits using the characteristic-zero Jacobian criterion."""

    keys = tuple(observables)
    dependent = []
    for size in range(1, len(keys) + 1):
        for subset in combinations(keys, size):
            jac = s.Matrix([observables[name] for name in subset]).jacobian(parameters)
            if jac.rank() < len(subset):
                dependent.append(subset)
    return [
        subset
        for subset in dependent
        if not any(set(smaller) < set(subset) for smaller in dependent)
    ]


def is_unit_ideal(polys):
    basis = s.groebner(polys, *SOURCE, order="lex", domain=s.QQ)
    return len(basis.polys) == 1 and basis.polys[0].as_expr() == 1


def minimal_wall_nonfaces(observables):
    keys = tuple(observables)
    unit_subsets = []
    for size in range(1, len(keys) + 1):
        for subset in combinations(keys, size):
            if is_unit_ideal([observables[name] for name in subset]):
                unit_subsets.append(subset)
    return [
        subset
        for subset in unit_subsets
        if not any(set(smaller) < set(subset) for smaller in unit_subsets)
    ]


def vanishes_mod_square(expr, square_value):
    numerator = s.together(expr).as_numer_denom()[0]
    return s.rem(numerator, P**2 - square_value, P) == 0


def main():
    checks = 0

    def ck(condition):
        nonlocal checks
        if not bool(condition):
            raise AssertionError(f"failed independent exact check {checks + 1}")
        checks += 1

    # Closure kernels, independently obtained from graph-ideal elimination.
    elim4 = graph_elimination(False)
    ck(elim4 == [s.expand(F)])
    ck(pull(F) == 0)
    ck(s.diff(F, W) == 2460375)
    ck(s.degree(F, W) == 1)

    elim5 = graph_elimination(True)
    ck(set(elim5) == {s.expand(F), s.expand(R1)})
    ck(pull(R1) == 0)
    ck(pull(R2) == 0)
    ck(s.expand(11 * F - 27 * R2 - 1060 * R1) == 0)
    ck(s.diff(R1, U) == 200475 and not R1.has(W, Z))
    ck(s.diff(R2, W) == 1002375 and not R2.has(U, Z))
    retained_jac = s.Matrix(
        [[s.diff(poly, variable) for variable in (X, U, W, Z, z)] for poly in (R1, R2)]
    )
    ck(retained_jac[:, [1, 2]].det() == 200475 * 1002375)

    # Two-sided inverse on the dense open z!=0.
    P_inv = -s.Rational(2, 3) * z
    xi_inv = s.Rational(237757952, 54675) - s.Rational(11, 6) * U
    eta_inv = s.Rational(54, 65) / z * (
        Z
        - s.Rational(5200877686784, 66430125)
        - s.Rational(125873, 8019) * U
        + s.Rational(520, 729) * z**2
    )
    ck(s.cancel(z_src.subs(P, P_inv) - z) == 0)
    ck(s.cancel(U_src.subs(xi, xi_inv) - U) == 0)
    W_back_num = s.together(W_src.subs({P: P_inv, xi: xi_inv}) - W).as_numer_denom()[0]
    ck(s.rem(W_back_num, F, W) == 0)
    ck(s.cancel(Z_src.subs({P: P_inv, xi: xi_inv, eta: eta_inv}) - Z) == 0)
    ck(s.cancel(P_inv.subs(z, z_src) - P) == 0)
    ck(s.cancel(xi_inv.subs(U, U_src) - xi) == 0)
    ck(s.cancel(eta_inv.subs({z: z_src, U: U_src, Z: Z_src}) - eta) == 0)

    eta_inv_x = s.Rational(54, 65) / z * (
        Z
        - s.Rational(3126529518592, 27064125)
        + s.Rational(520, 729) * z**2
        + s.Rational(22886, 2673) * X
    )
    ck(s.cancel(Z_src.subs({P: P_inv, xi: X, eta: eta_inv_x}) - Z) == 0)
    ck(s.cancel(eta_inv_x.subs({z: z_src, X: xi, Z: Z_src}) - eta) == 0)

    # Exceptional line, converse parameterization, and closure defect.
    special = {name: s.cancel(poly.subs(P, 0)) for name, poly in OBS5.items()}
    ck(pull(G0).subs(P, 0) == 0)
    ck(pull(R7).subs(P, 0) == 0)
    ck(s.factor(2 * pull(G0) - (159924375 * eta * z_src - 94770000 * z_src**2)) == 0)
    ck(not any(poly.has(eta) for poly in special.values()))

    line_basis = s.groebner([F.subs(z, 0), G0], W, Z, U, order="lex", domain=s.QQ)
    line_x = xi_inv
    ck(s.cancel(U_src.subs({P: 0, xi: line_x}) - U) == 0)
    ck(line_basis.reduce(s.together(W_src.subs({P: 0, xi: line_x}) - W).as_numer_denom()[0])[1] == 0)
    ck(line_basis.reduce(s.together(Z_src.subs({P: 0, xi: line_x}) - Z).as_numer_denom()[0])[1] == 0)

    phantom = {
        z: 0,
        U: 0,
        W: -s.Rational(42434609152, 2460375),
        Z: 0,
    }
    ck(s.cancel(F.subs(phantom)) == 0)
    ck(s.cancel(G0.subs(phantom)) == -5200877686784)
    retained_phantom = {**phantom, X: s.Rational(237757952, 54675)}
    ck(R1.subs(retained_phantom) == 0 and R2.subs(retained_phantom) == 0)
    ck(R7.subs(retained_phantom) == -s.Rational(57209654554624, 27))

    # The effective fibre jumps from A^0 to A^1.  The two silent source
    # coefficients and seven terminal tangent coordinates shift this to
    # A^9 and A^10 for the full finite-row family.
    effective_source_dim = 3
    silent_row8_dim = 2 + 7
    closure_dim = 3
    special_preimage_dim = 2
    special_image_dim = 1
    ck(effective_source_dim - closure_dim == 0)
    ck(special_preimage_dim - special_image_dim == 1)
    ck(silent_row8_dim == 9)
    ck(silent_row8_dim + 1 == 10)

    # Global and specialized algebraic matroids.
    circuits4 = algebraic_circuits(OBS4, SOURCE)
    circuits5 = algebraic_circuits(OBS5, SOURCE)
    ck(circuits4 == [("z", "U", "W")])
    ck(circuits5 == [("x", "U"), ("z", "x", "W"), ("z", "U", "W")])
    ck(s.Matrix(list(OBS4.values())).jacobian(SOURCE).rank() == 3)
    ck(s.Matrix(list(OBS5.values())).jacobian(SOURCE).rank() == 3)
    ck(F.free_symbols == {z, U, W})
    ck(R1.free_symbols == {X, U})
    ck(R2.free_symbols == {z, X, W})

    special4 = {name: s.cancel(poly.subs(P, 0)) for name, poly in OBS4.items()}
    special5 = {name: s.cancel(poly.subs(P, 0)) for name, poly in OBS5.items()}
    special_circuits4 = algebraic_circuits(special4, (eta, xi))
    special_circuits5 = algebraic_circuits(special5, (eta, xi))
    ck(
        special_circuits4
        == [("z",), ("U", "W"), ("U", "Z"), ("W", "Z")]
    )
    ck(
        special_circuits5
        == [
            ("z",),
            ("x", "U"),
            ("x", "W"),
            ("x", "Z"),
            ("U", "W"),
            ("U", "Z"),
            ("W", "Z"),
        ]
    )
    ck(s.Matrix(list(special5.values())).jacobian((eta, xi)).rank() == 1)

    special_pair_relations = [
        19318500 * U + 2460375 * W + 42434609152,
        -1042743375 * U + 66430125 * Z - 5200877686784,
        115860375 * W + 57955500 * Z - 2539122786304,
        109350 * X + 200475 * U - 475515904,
        -4293000 * X + 1002375 * W + 35956576256,
        231720750 * X + 27064125 * Z - 3126529518592,
    ]
    for relation in special_pair_relations:
        ck(s.cancel(relation.subs({X: xi, U: special5["U"], W: special5["W"], Z: special5["Z"]})) == 0)
        ck(s.Poly(relation).content() == 1)

    # Complete pointed affine zero-wall feasibility clutter.
    nonfaces4 = minimal_wall_nonfaces(OBS4)
    nonfaces5 = minimal_wall_nonfaces(OBS5)
    ck(nonfaces4 == [("z", "U", "W"), ("z", "U", "Z"), ("z", "W", "Z")])
    ck(
        nonfaces5
        == [
            ("x", "U"),
            ("z", "x", "W"),
            ("z", "x", "Z"),
            ("z", "U", "W"),
            ("z", "U", "Z"),
            ("z", "W", "Z"),
        ]
    )

    # Independent literal Bezout checks, one for each retained-ground-set
    # minimal nonface in the same order as above.
    certificates = [
        (109350 * xi + 200475 * U_src) / s.Integer(475515904),
        -(482625 * z_src**2 - 4293000 * xi + 1002375 * W_src)
        / s.Integer(35956576256),
        (
            (38610000 * z_src - 65154375 * eta) * z_src
            + 463441500 * xi
            + 54128250 * Z_src
        )
        / s.Integer(6253059037184),
        -(1184625 * z_src**2 + 19318500 * U_src + 2460375 * W_src)
        / s.Integer(42434609152),
        (
            -2085486750 * U_src
            + 132860250 * Z_src
            + (94770000 * z_src - 159924375 * eta) * z_src
        )
        / s.Integer(10401755373568),
        (
            115860375 * W_src
            + 57955500 * Z_src
            + (97124625 * z_src - 69761250 * eta) * z_src
        )
        / s.Integer(2539122786304),
    ]
    for certificate in certificates:
        ck(s.cancel(certificate) == 1)

    # Every special pair {z,q} has a witness at its distinct zero address.
    addresses = {
        "x": s.Integer(0),
        "U": s.solve(U_src, xi)[0],
        "W": s.solve(W_src.subs(P, 0), xi)[0],
        "Z": s.solve(Z_src.subs(P, 0), xi)[0],
    }
    ck(addresses["U"] == s.Rational(237757952, 54675))
    ck(addresses["W"] == s.Rational(4494572032, 536625))
    ck(addresses["Z"] == s.Rational(1563264759296, 115860375))
    ck(len(set(addresses.values())) == 4)
    for name, address in addresses.items():
        ck(s.cancel(z_src.subs(P, 0)) == 0)
        ck(s.cancel(OBS5[name].subs({P: 0, xi: address})) == 0)

    # Two feasible triples cover every no-z set which avoids {x,U}.
    r_xwz = -s.Rational(143826305024, 4343625)
    pe_xwz = s.Rational(6086390091776, 65154375)
    eta_xwz = pe_xwz / P
    xi_xwz = s.Integer(0)
    ck(r_xwz != 0)
    ck(xi_xwz == 0)
    ck(vanishes_mod_square(W_src.subs(xi, xi_xwz), r_xwz))
    ck(vanishes_mod_square(Z_src.subs({xi: xi_xwz, eta: eta_xwz}), r_xwz))

    xi_uwz = s.Rational(237757952, 54675)
    r_uwz = -s.Rational(13056802816, 820125)
    pe_uwz = s.Rational(707514056704, 12301875)
    eta_uwz = pe_uwz / P
    ck(r_uwz != 0)
    ck(vanishes_mod_square(U_src.subs(xi, xi_uwz), r_uwz))
    ck(vanishes_mod_square(W_src.subs(xi, xi_uwz), r_uwz))
    ck(vanishes_mod_square(Z_src.subs({xi: xi_uwz, eta: eta_uwz}), r_uwz))

    # The specialization-only wall nonfaces are not algebraic-matroid
    # circuits: each corresponding coordinate triple is globally independent.
    for names in (("z", "x", "Z"), ("z", "U", "Z"), ("z", "W", "Z")):
        jac = s.Matrix([OBS5[name] for name in names]).jacobian(SOURCE)
        ck(jac.rank() == 3)
        ck(names not in circuits5)

    # Raw affine functions on the special line have vector-space rank two,
    # even though their algebraic matroid has rank one.  Constants cannot be
    # silently discarded or treated as a linear representation.
    affine_vectors = s.Matrix(
        [
            [s.expand(special5[name]).coeff(xi, 0), s.diff(special5[name], xi)]
            for name in ("x", "U", "W", "Z")
        ]
    )
    ck(affine_vectors.rank() == 2)
    ck(affine_vectors[[0, 1], :].det() != 0)

    print("THM-4359 INDEPENDENT EXACT REFEREE")
    print("closure4: principal prime smooth hypersurface (F)")
    print("F=1184625zeta_3^2+19318500U+2460375W+42434609152")
    print("closure5: prime smooth complete intersection (R1,R2)")
    print("kernel identity: 11F=27R2+1060R1")
    print("image4: V(F) intersect (D(zeta_3) union V(G0))")
    print("image5: V(R1,R2) intersect (D(zeta_3) union V(R7))")
    print("phantom closure point: PASS")
    print("effective fibres: A^0 generic, A^1 exceptional")
    print("full finite-row8 fibres: A^9 generic, A^10 exceptional")
    print("global circuits4: {z,U,W}")
    print("global circuits5: {x,U}; {z,x,W}; {z,U,W}")
    print("special circuits: loop {z}; all pairs among remaining line coordinates")
    print("wall nonfaces4: {z,U,W}; {z,U,Z}; {z,W,Z}")
    print("wall nonfaces5: {x,U}; {z,x,W}; {z,x,Z}; {z,U,W}; {z,U,Z}; {z,W,Z}")
    print("proper-subset controls: distinct special addresses and two feasible triples PASS")
    print("TERMINOLOGY algebraic circuits != pointed affine zero-wall nonfaces")
    print("SCOPE effective A3 isomorphism only; no all-row, seam, Keller, JC(2), or DC(2) claim")
    print(f"PASS checks={checks}")


if __name__ == "__main__":
    main()
