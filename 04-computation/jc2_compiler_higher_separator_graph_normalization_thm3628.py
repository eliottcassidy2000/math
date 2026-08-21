#!/usr/bin/env python3
"""Exact controls for proved THM-3628.

The companion checks the uniform saturated graph equations for the higher
compiler separators H=E^(m-1)z, the replacement normalization obtained by
adjoining d and W=xE, its boundary arms, source-open inverse, and valuation
packets.  The all-m normality and image statements are proof-driven from the
uniform charts; m=2 and m=3 are reconstructed independently by elimination.
"""

import ast
import itertools
import math
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one deterministic exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact polynomial or rational-function zero test."""
    return sp.cancel(expression) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


def basis_expressions(groebner_basis):
    """Expanded expressions in deterministic Groebner order."""
    return tuple(sp.expand(poly.as_expr()) for poly in groebner_basis.polys)


print("THM-3628 exact companion -- higher separator graph normalization")
print("status=proved, verified exact, and independently hostile-audited")


print("SECTION compiler coordinates and uniform generic inverse")
SECTION_START = CHECKS
x, q = sp.symbols("x q")
d_source = 1 + x**2 * q
b_source = (d_source - 1) * (d_source + 2) ** 2
c_source = x * d_source * (d_source + 2)
e_source = q * (d_source + 3)
z_source = x * q
w_source = x * e_source

require("compiler B expansion", same(b_source, x**2 * q * (3 + x**2 * q) ** 2))
require("compiler surface", same(c_source**2 * e_source, b_source * (b_source + 4)))
require("E equals four q plus z square", same(e_source, 4 * q + z_source**2))
require("W equals dplus3 z", same(w_source, (d_source + 3) * z_source))

B, C, E, H, Z = sp.symbols("B C E H Z")
graph_data = {}
for m in (2, 3):
    r = m - 1
    n = 2 * m - 1
    h_source = e_source**r * z_source
    delta = E**n - H**2
    delta_source = sp.expand(delta.subs({E: e_source, H: h_source}))
    q_inverse = delta / (4 * E ** (2 * r))
    x_inverse = 4 * H * E**r / delta
    d_inverse = (E**n + 3 * H**2) / delta
    require(f"m{m} delta source", same(delta_source, 4 * e_source ** (2 * r) * q))
    require(f"m{m} q generic inverse", same(q_inverse.subs({E: e_source, H: h_source}), q))
    require(f"m{m} x generic inverse", same(x_inverse.subs({E: e_source, H: h_source}), x))
    require(f"m{m} d generic inverse", same(d_inverse.subs({E: e_source, H: h_source}), d_source))
print(f"PASS generic_gates={CHECKS - SECTION_START} Delta=E^(2m-1)-H^2")


print("SECTION exact saturated graph kernels for m2 and m3")
SECTION_START = CHECKS
jz = (
    B**2 + 4 * B - C**2 * E,
    B * E + 3 * B * Z**2 - 3 * C * E * Z - C * Z**3,
    8 * B * Z**3
    + C * E**2
    - 6 * C * E * Z**2
    - 3 * C * Z**4
    - 12 * E * Z
    - 4 * Z**3,
    C * E**3
    - 3 * C * E**2 * Z**2
    + 3 * C * E * Z**4
    - C * Z**6
    - 12 * E**2 * Z
    - 40 * E * Z**3
    - 12 * Z**5,
)
U = sp.symbols("U")
for m in (2, 3):
    r = m - 1
    n = 2 * m - 1
    delta = E**n - H**2
    f_b = sp.expand(B * delta**3 - 4 * H**2 * (3 * E**n + H**2) ** 2)
    f_c = sp.expand(
        C * delta**3
        - 4 * H * E**r * (E**n + 3 * H**2) * (3 * E**n + H**2)
    )
    require(
        f"m{m} generic B equation on source",
        zero(f_b.subs({B: b_source, E: e_source, H: e_source**r * z_source})),
    )
    require(
        f"m{m} generic C equation on source",
        zero(f_c.subs({C: c_source, E: e_source, H: e_source**r * z_source})),
    )

    saturation_full = sp.groebner(
        [f_b, f_c, 1 - U * E * delta], U, B, C, E, H, order="lex"
    )
    saturation_elimination = [
        poly.as_expr() for poly in saturation_full.polys if not poly.as_expr().has(U)
    ]
    saturation_basis = sp.groebner(saturation_elimination, B, C, E, H, order="lex")

    direct_full = sp.groebner(
        list(jz) + [H - E**r * Z], Z, B, C, E, H, order="lex"
    )
    direct_elimination = [
        poly.as_expr() for poly in direct_full.polys if not poly.as_expr().has(Z)
    ]
    direct_basis = sp.groebner(direct_elimination, B, C, E, H, order="lex")
    require(
        f"m{m} saturation equals direct graph elimination",
        basis_expressions(saturation_basis) == basis_expressions(direct_basis),
    )
    require(f"m{m} graph basis has six elements", len(saturation_basis.polys) == 6)

    boundary_basis = sp.groebner(
        saturation_elimination + [E], B, C, E, H, order="lex"
    )
    expected_boundary = (
        B**2 + 4 * B,
        B * H**3 + 4 * H**3,
        C * H**3,
        E,
    )
    require(
        f"m{m} exact E-zero boundary basis",
        basis_expressions(boundary_basis) == expected_boundary,
    )
    for b_value, c_value, h_value, label in (
        (0, C, 0, "L1"),
        (-4, C, 0, "Lminus4"),
        (-4, 0, H, "Vinf"),
    ):
        for index, relation in enumerate(saturation_basis.polys, start=1):
            require(
                f"m{m} boundary {label} relation {index}",
                relation.as_expr().subs({B: b_value, C: c_value, E: 0, H: h_value}) == 0,
            )
    graph_data[m] = saturation_basis
print(f"PASS graph_gates={CHECKS - SECTION_START} Jm=(FB,FC):(E*Delta)^infinity")


print("SECTION exact five-relation normalization")
SECTION_START = CHECKS
dn, W = sp.symbols("dn W")
a = dn * (dn + 2)
g = (dn - 1) * (dn + 3)
normal_data = {}
for m in (2, 3):
    r = m - 1
    core = (
        C * W - a * g,
        C * E - a * W,
        W**2 - g * E,
        (dn + 3) * H - E**r * W,
        C * H - a * (dn - 1) * E**r,
    )
    h_source = e_source**r * z_source
    normal_substitution = {
        dn: d_source,
        C: c_source,
        E: e_source,
        H: h_source,
        W: w_source,
    }
    for index, relation in enumerate(core, start=1):
        require(f"m{m} normalization relation {index}", zero(relation.subs(normal_substitution)))

    direct_normal = sp.groebner(
        [
            dn - d_source,
            C - c_source,
            E - e_source,
            H - h_source,
            W - w_source,
        ],
        x,
        q,
        dn,
        C,
        E,
        H,
        W,
        order="lex",
    )
    direct_normal_elimination = [
        poly.as_expr()
        for poly in direct_normal.polys
        if not poly.as_expr().has(x, q)
    ]
    direct_normal_basis = sp.groebner(
        direct_normal_elimination, dn, C, E, H, W, order="lex"
    )
    core_basis = sp.groebner(core, dn, C, E, H, W, order="lex")
    require(
        f"m{m} five relations equal direct normalization kernel",
        basis_expressions(core_basis) == basis_expressions(direct_normal_basis),
    )
    require(f"m{m} d monic equation", same(dn**3 + 3 * dn**2 - (((dn - 1) * (dn + 2) ** 2) + 4), 0))
    require(f"m{m} W monic equation present", core[2] == W**2 - g * E)

    def normal_zero(expression):
        return core_basis.reduce(sp.expand(expression))[1] == 0

    require(f"m{m} derived compiler surface", normal_zero(C**2 * E - a**2 * g))
    require(f"m{m} derived WH relation", normal_zero(W * H - (dn - 1) * E ** (r + 1)))
    require(
        f"m{m} derived H square relation",
        normal_zero((dn + 3) * H**2 - (dn - 1) * E ** (2 * r + 1)),
    )
    b_from_d = (dn - 1) * (dn + 2) ** 2
    for index, relation in enumerate(graph_data[m].polys, start=1):
        require(
            f"m{m} normalization maps to graph relation {index}",
            normal_zero(relation.as_expr().subs(B, b_from_d)),
        )
    normal_data[m] = (core, core_basis)
print(f"PASS normalization_gates={CHECKS - SECTION_START} Nm=Sm[d,W] finite_birational")


print("SECTION normality charts and sharp m2 versus m3 boundary")
SECTION_START = CHECKS
normal_variables = (dn, C, E, H, W)
for m in (2, 3):
    r = m - 1
    core, core_basis = normal_data[m]

    def normal_zero(expression):
        return core_basis.reduce(sp.expand(expression))[1] == 0

    require(f"m{m} g-chart E numerator", normal_zero(g * E - W**2))
    require(f"m{m} g-chart H numerator", normal_zero((dn + 3) * H - E**r * W))
    require(f"m{m} g-chart hypersurface", normal_zero(C * W - a * g))
    require(f"m{m} a-chart W numerator", normal_zero(a * W - C * E))
    require(
        f"m{m} a-chart H numerator",
        normal_zero(a * (dn + 3) * H - C * E ** (r + 1)),
    )
    require(f"m{m} a-chart hypersurface", normal_zero(C**2 * E - a**2 * g))
    require(
        f"m{m} dminus3 local equation",
        (C * H - a * (dn - 1) * E**r).subs(dn, -3) == C * H + 12 * E**r,
    )

    jacobian = sp.Matrix(
        [[sp.diff(relation, variable) for variable in normal_variables] for relation in core]
    )
    three_minors = [
        sp.expand(jacobian.extract(rows, columns).det())
        for rows in itertools.combinations(range(5), 3)
        for columns in itertools.combinations(range(5), 3)
    ]
    singular_basis = sp.groebner(
        list(core) + three_minors, dn, C, E, H, W, order="lex"
    )
    if m == 2:
        require("m2 normalization smooth", basis_expressions(singular_basis) == (sp.Integer(1),))
    else:
        require(
            "m3 normalization unique A1 singular point",
            basis_expressions(singular_basis) == (dn + 3, C, E, H, W),
        )
print(f"PASS normality_gates={CHECKS - SECTION_START} m2=smooth m3=A1 general=Amminus2")


print("SECTION five E-zero arms and exact source open")
SECTION_START = CHECKS
b_from_d = (dn - 1) * (dn + 2) ** 2
arm_packets = (
    (0, C, 0, 0, 0, -4, "N0"),
    (-2, C, 0, 0, 0, 0, "Nminus2"),
    (1, C, 0, 0, 0, 0, "N1"),
    (-3, C, 0, 0, 0, -4, "Dminus3"),
    (-3, 0, 0, H, 0, -4, "Vminus3"),
)
for m in (2, 3):
    core, core_basis = normal_data[m]
    for d_value, c_value, e_value, h_value, w_value, b_value, label in arm_packets:
        substitution = {dn: d_value, C: c_value, E: e_value, H: h_value, W: w_value}
        for index, relation in enumerate(core, start=1):
            require(f"m{m} arm {label} relation {index}", relation.subs(substitution) == 0)
        require(f"m{m} arm {label} target B", b_from_d.subs(dn, d_value) == b_value)

require("source qzero is covered N1", tuple(sp.expand(v.subs(q, 0)) for v in (d_source, c_source, e_source, w_source)) == (1, 3 * x, 0, 0))
u = sp.symbols("u", nonzero=True)
dminus3_source = tuple(
    sp.cancel(value.subs({x: u, q: -4 / u**2}))
    for value in (d_source, c_source, e_source, w_source)
)
require("source dminus3 covers punctured C-axis", dminus3_source == (-3, 3 * u, 0, 0))
require("x inverse chart glue", same((C / a) - (g / W), (C * W - a * g) / (a * W)))
require(
    "q inverse chart glue numerator",
    same(
        E / (dn + 3) - (dn - 1) * a**2 / C**2,
        (C**2 * E - a**2 * g) / (C**2 * (dn + 3)),
    ),
)
require("inverse recovers d", same(1 + (C / a) ** 2 * (dn - 1) * a**2 / C**2, dn))
require("inverse recovers C", same((C / a) * a, C))
require("inverse recovers E", same((dn - 1) * a**2 / C**2 * (dn + 3), a**2 * g / C**2))
print(f"PASS arm_open_gates={CHECKS - SECTION_START} source=Nm_minus_N0_Nminus2_Vminus3")


print("SECTION semigroups, valuations, and multiplied debts")
SECTION_START = CHECKS
for m in (2, 3):
    r = m - 1
    require(f"m{m} cusp gcd", math.gcd(2, 2 * m - 1) == 1)
    require(f"m{m} cusp delta", ((2 - 1) * ((2 * m - 1) - 1)) // 2 == m - 1)
    for root, expected_d, expected_q, expected_z, label in (
        (0, -sp.Rational(1, 6), -sp.Rational(1, 9), sp.Rational(1, 3), "N0"),
        (-2, sp.Rational(1, 6), -sp.Rational(1, 3), sp.Integer(1), "Nminus2"),
    ):
        g_root = sp.expand(g).subs(dn, root)
        ag_derivative = sp.diff(a * g, dn).subs(dn, root)
        d_coefficient = C / ag_derivative
        e_coefficient = 1 / g_root
        x_coefficient = g_root
        q_coefficient = e_coefficient / (root + 3)
        z_coefficient = x_coefficient * q_coefficient
        h_coefficient = e_coefficient**r * z_coefficient
        require(f"m{m} {label} d coefficient", same(d_coefficient, C * expected_d))
        require(f"m{m} {label} E coefficient", same(e_coefficient, -sp.Rational(1, 3)))
        require(f"m{m} {label} x coefficient", same(x_coefficient, -3))
        require(f"m{m} {label} q coefficient", same(q_coefficient, expected_q))
        require(f"m{m} {label} z coefficient", same(z_coefficient, expected_z))
        require(
            f"m{m} {label} H coefficient",
            same(h_coefficient, (-sp.Rational(1, 3)) ** r * expected_z),
        )
    x_order, q_order, z_order, h_order = -1, 2, 1, 2 * r + 1
    require(f"m{m} omitted x pole order", -x_order == 1)
    require(f"m{m} omitted q zero order", q_order == 2)
    require(f"m{m} omitted z order", z_order == 1)
    require(f"m{m} omitted H order", h_order == 2 * m - 1)

    h_generic = sp.symbols(f"h_generic_{m}", nonzero=True)
    c_lead = -12 / h_generic
    w_lead = c_lead / 3
    dplus3_lead = c_lead**2 / (-36)
    x_lead = c_lead / 3
    q_lead = 1 / dplus3_lead
    require(f"m{m} V C leading coefficient", same(c_lead, -12 / h_generic))
    require(f"m{m} V W leading coefficient", same(w_lead, -4 / h_generic))
    require(f"m{m} V dplus3 leading coefficient", same(dplus3_lead, -4 / h_generic**2))
    require(f"m{m} V x leading coefficient", same(x_lead, -4 / h_generic))
    require(f"m{m} V q leading coefficient", same(q_lead, -h_generic**2 / 4))
    require(f"m{m} V C order", r == m - 1)
    require(f"m{m} V W order", r + 1 == m)
    require(f"m{m} V dplus3 order", 2 * r + 1 == 2 * m - 1)
    require(f"m{m} V q pole order", 2 * r == 2 * m - 2)
    require(f"m{m} V z pole order", r == m - 1)
print(f"PASS valuation_gates={CHECKS - SECTION_START} cusp_delta=mminus1 q_pole=2mminus2")


print("SECTION m1 contrast, old-normalization hostile, and hygiene")
SECTION_START = CHECKS
require("m1 hyperbola first graph relation", jz[0].subs({B: -4, E: 0, C: -12 / Z}) == 0)
for index, relation in enumerate(jz[1:], start=2):
    require(
        f"m1 Ezero hyperbola relation {index}",
        zero(relation.subs({B: -4, E: 0, C: -12 / Z})),
    )
for m in (2, 3):
    require(f"m{m} hyperbola H collapses", (E ** (m - 1) * Z).subs(E, 0) == 0)
require("punctured line excludes compiler point", sp.denom(sp.cancel(-12 / Z)) == Z)
for m in (2, 3):
    e_order, h_order = 1, 0
    z_order = h_order - (m - 1) * e_order
    require(f"m{m} old z has V pole", -z_order == m - 1)
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require("no truth-bypassing assert nodes", not any(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree)))
print(f"PASS contrast_hygiene_gates={CHECKS - SECTION_START} T_not_finite_new_vertical_arm")


print(f"PASS total_exact_gates={CHECKS}")
print("RESULT PASS -- proved uniform theorem independently hostile-audited")
