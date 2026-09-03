#!/usr/bin/env python3
"""THM-4359 primary exact certificate.

Starting only from THM-4308's displayed source-normal row-eight response,
compute its scheme-theoretic closure, its actual constructible image, the
Phi=0 specialization, the affine-modification identity, and the complete
algebraic-matroid circuit and coordinate-wall nonface ledgers.

This proves only a finite row-eight response theorem.  It proves no all-row
lift, seam entry, Keller-pair existence, JC(2), or DC(2) consequence.
"""

from itertools import combinations
import sys

import sympy as s


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def check(condition, label):
    global CHECKS
    if not bool(condition):
        raise AssertionError(label)
    CHECKS += 1


def same_ideal(gens_a, gens_b, variables):
    """Compare reduced lexicographic Groebner bases over QQ."""

    ga = s.groebner(gens_a, *variables, order="lex", domain=s.QQ)
    gb = s.groebner(gens_b, *variables, order="lex", domain=s.QQ)
    return all(ga.reduce(s.expand(g))[1] == 0 for g in gens_b) and all(
        gb.reduce(s.expand(g))[1] == 0 for g in gens_a
    )


def pull(expr):
    return s.cancel(
        expr.subs(
            {
                z: zf,
                x: xi,
                U: Uf,
                W: Wf,
                Z: Zf,
            },
            simultaneous=True,
        )
    )


def matroid_circuits(functions, source_variables):
    """All circuits by the characteristic-zero Jacobian criterion."""

    names = list(functions)
    dependent = []
    for size in range(1, len(names) + 1):
        for subset in combinations(names, size):
            jac = s.Matrix(
                [[s.diff(functions[name], v) for v in source_variables] for name in subset]
            )
            if jac.rank() < size:
                dependent.append(subset)
    return [
        subset
        for subset in dependent
        if not any(set(smaller) < set(subset) for smaller in dependent)
    ]


def minimal_coordinate_nonfaces(functions):
    """All inclusion-minimal coordinate zero sets with unit source ideal."""

    names = list(functions)
    nonfaces = []
    for size in range(1, len(names) + 1):
        for subset in combinations(names, size):
            gb = s.groebner(
                [functions[name] for name in subset],
                eta,
                phi,
                xi,
                order="lex",
            )
            if any(p.as_expr() == 1 for p in gb.polys):
                nonfaces.append(subset)
    return [
        subset
        for subset in nonfaces
        if not any(set(smaller) < set(subset) for smaller in nonfaces)
    ]


phi, eta, xi = s.symbols("Phi eta xi")
z, x, U, W, Z = s.symbols("zeta xi_t U W Z")

DELTA = s.Rational(896, 15)
THETA = s.Rational(512, 75)
K = s.Rational(-32, 5)
UPSILON5 = s.Rational(-731648, 2025)

for name, value in (
    ("Delta", DELTA),
    ("Theta", THETA),
    ("K", K),
    ("upsilon5", UPSILON5),
):
    check(value != 0, f"fixed nonzero constant {name}")

# THM-4308 response, first in its displayed form and then simplified.
zf = -s.Rational(3, 2) * phi
U_displayed = (s.Integer(475515904) - s.Integer(109350) * xi) / s.Integer(200475)
W_displayed = -(
    s.Integer(4343625) * phi**2
    - s.Integer(17172000) * xi
    + s.Integer(143826305024)
) / s.Integer(4009500)
Z_displayed = (
    s.Integer(12506118074368)
    - s.Integer(173745000) * phi**2
    - s.Integer(195463125) * phi * eta
    - s.Integer(926883000) * xi
) / s.Integer(108256500)

Uf = s.Rational(475515904, 200475) - s.Rational(6, 11) * xi
Wf = (
    -s.Rational(13, 12) * phi**2
    + s.Rational(424, 99) * xi
    - s.Rational(35956576256, 1002375)
)
Zf = (
    s.Rational(3126529518592, 27064125)
    - s.Rational(130, 81) * phi**2
    - s.Rational(65, 36) * phi * eta
    - s.Rational(22886, 2673) * xi
)

check(s.cancel(U_displayed - Uf) == 0, "simplified U")
check(s.cancel(W_displayed - Wf) == 0, "simplified W")
check(s.cancel(Z_displayed - Zf) == 0, "simplified Z")

# Response-only closure kernel.
F = (
    s.Integer(1184625) * z**2
    + s.Integer(19318500) * U
    + s.Integer(2460375) * W
    + s.Integer(42434609152)
)
check(pull(F) == 0, "F pulls back to zero")

graph4 = [z - zf, U - Uf, W - Wf, Z - Zf]
gb4 = s.groebner(graph4, eta, phi, xi, z, U, W, Z, order="lex")
elim4 = [
    p.as_expr()
    for p in gb4.polys
    if not p.as_expr().has(eta, phi, xi)
]
check(same_ideal(elim4, [F], (z, U, W, Z)), "response kernel is (F)")

jac_det = s.factor(
    s.Matrix(
        [[s.diff(f, v) for v in (phi, xi, eta)] for f in (zf, Uf, Zf)]
    ).det()
)
check(jac_det == -s.Rational(65, 44) * phi, "z,U,Z algebraic independence")
check(s.Poly(F, W).degree() == 1, "F irreducible by W-linearity")
check(s.diff(F, W) == 2460375, "V(F) smooth")

# Explicit isomorphism on z!=0.
XI_U = s.Rational(237757952, 54675)
Z_U = s.Rational(5200877686784, 66430125)
phi_inverse = -s.Rational(2, 3) * z
xi_inverse = XI_U - s.Rational(11, 6) * U
eta_inverse = s.Rational(54, 65) / z * (
    Z
    - Z_U
    - s.Rational(125873, 8019) * U
    + s.Rational(520, 729) * z**2
)

inverse_sub = {z: zf, U: Uf, Z: Zf}
check(s.cancel(phi_inverse.subs(inverse_sub) - phi) == 0, "inverse recovers Phi")
check(s.cancel(xi_inverse.subs(inverse_sub) - xi) == 0, "inverse recovers xi")
check(s.cancel(eta_inverse.subs(inverse_sub) - eta) == 0, "inverse recovers eta")

source_from_target = {
    phi: phi_inverse,
    xi: xi_inverse,
    eta: eta_inverse,
}
U_back = s.cancel(Uf.subs(source_from_target, simultaneous=True) - U)
W_back = s.cancel(Wf.subs(source_from_target, simultaneous=True) - W)
Z_back = s.cancel(Zf.subs(source_from_target, simultaneous=True) - Z)
W_on_F = -(
    s.Integer(1184625) * z**2
    + s.Integer(19318500) * U
    + s.Integer(42434609152)
) / s.Integer(2460375)
check(U_back == 0, "target reconstruction U")
check(s.cancel(W_back.subs(W, W_on_F)) == 0, "target reconstruction W modulo F")
check(Z_back == 0, "target reconstruction Z")

# Exact Phi=0 special image.
F0 = s.expand(F.subs(z, 0))
G0 = (
    -s.Integer(1042743375) * U
    + s.Integer(66430125) * Z
    - s.Integer(5200877686784)
)
check(s.cancel(pull(F0).subs(phi, 0)) == 0, "special F0 relation")
check(s.cancel(pull(G0).subs(phi, 0)) == 0, "special G0 relation")

special_graph = [
    z,
    U - Uf.subs(phi, 0),
    W - Wf.subs(phi, 0),
    Z - Zf.subs(phi, 0),
]
special_gb = s.groebner(special_graph, eta, xi, z, U, W, Z, order="lex")
special_elim = [
    p.as_expr() for p in special_gb.polys if not p.as_expr().has(eta, xi)
]
check(
    same_ideal(special_elim, [z, F, G0], (z, U, W, Z)),
    "specialized image ideal is (z,F,G0)",
)

# The converse parametrization of the special line.
W_special_target = s.solve(F0, W)[0]
Z_special_target = s.solve(G0, Z)[0]
check(
    s.cancel(Wf.subs({phi: 0, xi: xi_inverse}) - W_special_target) == 0,
    "special line converse for W",
)
check(
    s.cancel(Zf.subs({phi: 0, xi: xi_inverse}) - Z_special_target) == 0,
    "special line converse for Z",
)

affine_modification_identity = s.cancel(
    2 * pull(G0)
    - zf * (s.Integer(159924375) * eta - s.Integer(94770000) * zf)
)
check(affine_modification_identity == 0, "affine-modification identity")

W_U = -s.Rational(42434609152, 2460375)
phantom = {z: 0, U: 0, W: W_U, Z: 0}
check(F.subs(phantom) == 0, "phantom point lies in closure")
check(G0.subs(phantom) == -5200877686784, "phantom point misses image")

# Extended closure retaining xi_10 as an observable.
L_xU = s.Integer(109350) * x + s.Integer(200475) * U - s.Integer(475515904)
L_zxW = (
    s.Integer(482625) * z**2
    - s.Integer(4293000) * x
    + s.Integer(1002375) * W
    + s.Integer(35956576256)
)
check(pull(L_xU) == 0, "extended relation x,U")
check(pull(L_zxW) == 0, "extended relation z,x,W")

graph5 = [z - zf, x - xi, U - Uf, W - Wf, Z - Zf]
gb5 = s.groebner(graph5, eta, phi, xi, z, x, U, W, Z, order="lex")
elim5 = [
    p.as_expr()
    for p in gb5.polys
    if not p.as_expr().has(eta, phi, xi)
]
check(
    same_ideal(elim5, [L_xU, L_zxW], (z, x, U, W, Z)),
    "extended closure kernel",
)

response_functions = {"zeta": zf, "U": Uf, "W": Wf, "Z": Zf}
extended_functions = {
    "zeta": zf,
    "xi": xi,
    "U": Uf,
    "W": Wf,
    "Z": Zf,
}
special_functions = {
    name: s.expand(value.subs(phi, 0)) for name, value in extended_functions.items()
}

response_circuits = matroid_circuits(response_functions, (phi, eta, xi))
extended_circuits = matroid_circuits(extended_functions, (phi, eta, xi))
special_circuits = matroid_circuits(special_functions, (eta, xi))
check(
    response_circuits == [("zeta", "U", "W")],
    "complete response matroid circuit ledger",
)
check(
    extended_circuits
    == [("xi", "U"), ("zeta", "xi", "W"), ("zeta", "U", "W")],
    "complete extended matroid circuit ledger",
)
check(
    special_circuits
    == [
        ("zeta",),
        ("xi", "U"),
        ("xi", "W"),
        ("xi", "Z"),
        ("U", "W"),
        ("U", "Z"),
        ("W", "Z"),
    ],
    "complete special matroid circuit ledger",
)

# Primitive circuit equations on the Phi=0 line.
special_pair_circuits = [
    s.Integer(109350) * x + s.Integer(200475) * U - s.Integer(475515904),
    -s.Integer(4293000) * x + s.Integer(1002375) * W + s.Integer(35956576256),
    s.Integer(231720750) * x + s.Integer(27064125) * Z - s.Integer(3126529518592),
    s.Integer(19318500) * U + s.Integer(2460375) * W + s.Integer(42434609152),
    -s.Integer(1042743375) * U + s.Integer(66430125) * Z - s.Integer(5200877686784),
    s.Integer(115860375) * W + s.Integer(57955500) * Z - s.Integer(2539122786304),
]
for index, circuit in enumerate(special_pair_circuits, 1):
    check(
        s.cancel(pull(circuit).subs(phi, 0)) == 0,
        f"special pair circuit {index}",
    )

# Exhaustive coordinate-wall nonfaces.
extended_nonfaces = minimal_coordinate_nonfaces(extended_functions)
response_nonfaces = minimal_coordinate_nonfaces(response_functions)
expected_extended_nonfaces = [
    ("xi", "U"),
    ("zeta", "xi", "W"),
    ("zeta", "xi", "Z"),
    ("zeta", "U", "W"),
    ("zeta", "U", "Z"),
    ("zeta", "W", "Z"),
]
expected_response_nonfaces = [
    ("zeta", "U", "W"),
    ("zeta", "U", "Z"),
    ("zeta", "W", "Z"),
]
check(
    extended_nonfaces == expected_extended_nonfaces,
    "complete extended nonface ledger",
)
check(
    response_nonfaces == expected_response_nonfaces,
    "complete response-only nonface ledger",
)

# Literal source-ring Bezout certificates, in the same order as the six
# extended nonfaces above.
certificates = [
    (s.Integer(109350) * xi + s.Integer(200475) * Uf) / s.Integer(475515904),
    -(
        s.Integer(482625) * zf**2
        - s.Integer(4293000) * xi
        + s.Integer(1002375) * Wf
    )
    / s.Integer(35956576256),
    (
        (s.Integer(38610000) * zf - s.Integer(65154375) * eta) * zf
        + s.Integer(463441500) * xi
        + s.Integer(54128250) * Zf
    )
    / s.Integer(6253059037184),
    -(
        s.Integer(1184625) * zf**2
        + s.Integer(19318500) * Uf
        + s.Integer(2460375) * Wf
    )
    / s.Integer(42434609152),
    (
        -s.Integer(2085486750) * Uf
        + s.Integer(132860250) * Zf
        + (s.Integer(94770000) * zf - s.Integer(159924375) * eta) * zf
    )
    / s.Integer(10401755373568),
    (
        s.Integer(115860375) * Wf
        + s.Integer(57955500) * Zf
        + (s.Integer(97124625) * zf - s.Integer(69761250) * eta) * zf
    )
    / s.Integer(2539122786304),
]
for index, certificate in enumerate(certificates, 1):
    check(s.cancel(certificate - 1) == 0, f"unit certificate {index}")


def source_values(point):
    return {
        "zeta": s.simplify(zf.subs(point)),
        "xi": s.simplify(xi.subs(point)),
        "U": s.simplify(Uf.subs(point)),
        "W": s.simplify(Wf.subs(point)),
        "Z": s.simplify(Zf.subs(point)),
    }


XI_W = s.Rational(4494572032, 536625)
XI_Z = s.Rational(1563264759296, 115860375)
rho_xW = s.sqrt(-s.Rational(143826305024, 4343625))
rho_UW = s.sqrt(-s.Rational(13056802816, 820125))

# Every admissible pair in the five-coordinate ground set.  The omitted
# tenth pair {xi,U} is the minimal two-wall nonface.
pair_witnesses = {
    ("zeta", "xi"): {phi: 0, eta: 0, xi: 0},
    ("zeta", "U"): {phi: 0, eta: 0, xi: XI_U},
    ("zeta", "W"): {phi: 0, eta: 0, xi: XI_W},
    ("zeta", "Z"): {phi: 0, eta: 0, xi: XI_Z},
    ("xi", "W"): {phi: rho_xW, eta: 0, xi: 0},
    ("xi", "Z"): {
        phi: 1,
        eta: s.Rational(12505944329368, 195463125),
        xi: 0,
    },
    ("U", "W"): {phi: rho_UW, eta: 0, xi: XI_U},
    ("U", "Z"): {
        phi: 1,
        eta: s.Rational(1600237252472, 36905625),
        xi: XI_U,
    },
    ("W", "Z"): {
        phi: 1,
        eta: s.Rational(781201309507, 32197500),
        xi: s.Rational(143830648649, 17172000),
    },
}
for pair, point in pair_witnesses.items():
    vals = source_values(point)
    for wall in pair:
        check(vals[wall] == 0, f"pair witness {pair}: {wall}")

# The only two triples that neither contain {xi,U} nor are declared
# nonfaces.  They prove the subset ledger is exhaustive without relying
# solely on the Groebner enumeration.
triple_UWZ = {
    phi: rho_UW,
    eta: s.Rational(707514056704, 12301875) / rho_UW,
    xi: XI_U,
}
triple_xWZ = {
    phi: rho_xW,
    eta: s.Rational(6086390091776, 65154375) / rho_xW,
    xi: 0,
}
for label, point, walls in (
    ("UWZ", triple_UWZ, ("U", "W", "Z")),
    ("xWZ", triple_xWZ, ("xi", "W", "Z")),
):
    vals = source_values(point)
    for wall in walls:
        check(vals[wall] == 0, f"live triple {label}: {wall}")
    check(vals["zeta"] != 0, f"live triple {label}: zeta nonzero")

# Hostiles to the two most tempting overclaims.
check(
    s.expand(F.subs(W, W + 1) - F) == 2460375,
    "perturbing W exits the closure",
)
check(
    pull(G0).subs({phi: 1, eta: 0, xi: 0}) == -106616250,
    "G0 is not a global kernel equation",
)

zero_addresses = [s.Integer(0), XI_U, XI_W, XI_Z]
for index, address in enumerate(zero_addresses):
    check(address != 0 or index == 0, f"zero address {index}")
for i, j in combinations(range(len(zero_addresses)), 2):
    check(zero_addresses[i] != zero_addresses[j], f"distinct zero addresses {i},{j}")

print("THM-4359 PRIMARY: ROW-EIGHT CONSTRUCTIBLE RESPONSE / AFFINE MODIFICATION")
print("fixed Delta=896/15 Theta=512/75 K=-32/5 upsilon5=-731648/2025")
print("closure kernel: (F)")
print("F=1184625*zeta^2+19318500*U+2460375*W+42434609152")
print("smoothness: dF/dW=2460375")
print("generic Jacobian det d(zeta,U,Z)/d(Phi,xi,eta)=-65*Phi/44")
print("on zeta!=0: response is an isomorphism onto V(F)")
print("inverse Phi=-2*zeta/3")
print("inverse xi=237757952/54675-11*U/6")
print(
    "inverse eta=54/(65*zeta)*(Z-5200877686784/66430125"
    "-125873*U/8019+520*zeta^2/729)"
)
print("special image: zeta=0, F=0, G0=0")
print("G0=-1042743375*U+66430125*Z-5200877686784")
print("affine modification: 2*G0=zeta*(159924375*eta-94770000*zeta)")
print("effective source ring=A[G0/zeta], A=Q[zeta,U,Z]")
print("actual image=V(F) intersect (D(zeta) union V(G0))")
print("missing closure locus=V(F,zeta) intersect D(G0)")
print(
    "phantom=(0,0,-42434609152/2460375,0): "
    "F=0 G0=-5200877686784"
)
print(
    "extended kernel: "
    "109350*xi+200475*U-475515904, "
    "482625*zeta^2-4293000*xi+1002375*W+35956576256"
)
print("response circuits:", response_circuits)
print("extended circuits:", extended_circuits)
print("special circuits:", special_circuits)
print("response-only minimal nonfaces:", response_nonfaces)
print("extended minimal nonfaces:", extended_nonfaces)
print("special zero addresses: 0,xi_U,xi_W,xi_Z pairwise distinct")
print(
    "live hostile U=W=Z=0: xi=237757952/54675, "
    "Phi^2=-13056802816/820125"
)
print(
    "live hostile xi=W=Z=0: "
    "Phi^2=-143826305024/4343625"
)
print(
    "specialization firewall: (zeta,F) is strictly contained in "
    "(zeta,F,G0)"
)
print("SCOPE finite-row-eight-only; no all-row lift, seam entry, Keller pair, JC(2), or DC(2)")
print(f"PASS checks={CHECKS}")
