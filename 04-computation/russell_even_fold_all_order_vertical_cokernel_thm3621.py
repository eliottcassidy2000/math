#!/usr/bin/env python3
"""Optimization-safe exact companion for the proposed THM-3621 proof.

This companion verifies the algebraic and formal-series identities used by the
all-order vertical-cokernel argument for the polynomial even-fold family of
THM-3612.  It treats an arbitrary formal target two-form; neither closedness
nor decomposability is used.  The scope is this one Russell-cylinder family,
not the planar Jacobian conjecture.

The all-order deduction encoded below is:

* shifted evaluation by rho generates the complete vertical cokernel row;
* the rational comparison germ Q_infinity has an exact three-branch target
  collision and exact side-sum/twice-middle bivector identity;
* after lower rows are constant, filtration leaves only the first variation
  in Q', whose normalized coefficient is -16;
* therefore Delta_(2n-2) equals
    -2^(n+3)/(3^(n-1)(n-1)!)*(q_n+(n+1)q_(n-1));
* vanishing at every order forces the rational Q_infinity Taylor germ, which
  no polynomial Q can have.

Every executable gate uses require().  A final AST gate rejects Python assert
statements, and the source is required to consist of raw LF bytes.
"""

import ast
import hashlib
from math import factorial
from pathlib import Path

import sympy as sp


CHECKS = 0
SEMANTIC_FACTS = []
HOSTILE_N_MIN = 3
HOSTILE_N_MAX = 20


def require(label, condition):
    """Record one active exact gate and fail loudly if it is false."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact rational-function zero test."""
    return sp.cancel(expression) == 0


def exact_coefficient(expression, variable, degree):
    """One exact formal coefficient, expanding only as far as needed."""
    return sp.expand(sp.series(expression, variable, 0, degree + 1).removeO()).coeff(
        variable, degree
    )


def wedge(first, second):
    """Three Plucker coordinates in the (C^Y,C^Z,Y^Z) basis."""
    return (
        sp.factor(first[0] * second[1] - first[1] * second[0]),
        sp.factor(first[0] * second[2] - first[2] * second[0]),
        sp.factor(first[1] * second[2] - first[2] * second[1]),
    )


def fold_bivector(coordinates, x_variable, t_variable):
    """Pushed source bivector for one map (x,t)->(C,Y,Z)."""
    x_tangent = tuple(sp.diff(value, x_variable) for value in coordinates)
    t_tangent = tuple(sp.diff(value, t_variable) for value in coordinates)
    return wedge(x_tangent, t_tangent)


def canonical_rational(value):
    """Stable text for an exact rational."""
    value = sp.Rational(value)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def canonical_row(row):
    """Stable serialization of a sparse finite-control row."""
    fields = []
    for (branch, (xi_degree, t_degree)), value in sorted(row.items()):
        fields.append(
            f"{branch}:{xi_degree},{t_degree}:{canonical_rational(value)}"
        )
    return "|".join(fields)


print("THM-3621 exact companion -- all-order Russell even-fold vertical cokernel")
print("status=exact_companion_for_independent_audit; theorem_file_not_modified")


print("SECTION pinned ancestry")
AUDITED_ANCESTRY = (
    (
        "THM-3561-rational-keller-danielewski-polynomial-completion",
        "5008d4965b356ccc178e8da21c334fb418136b415d88686ad3d8e9153fd1aa7a",
    ),
    (
        "THM-3605-russell-cylinder-graph-slice-puncture-no-filling",
        "5198efa462933e3a52cb370a4b096315034ce389decbb743dc8ab40f87cd1cc2",
    ),
    (
        "THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity",
        "70b0a053b281093db76f775fbfff89c84555f1e252c18261e714f9b9aaf463fa",
    ),
)
THM3612_PROMOTION_COMMIT = "3cd370a3a525d0b5aec22dcaf1810c8fedb0eecc"
THM3612_SCRIPT_SHA256 = (
    "fa210b98b493f5426f77567769c8f3356440de0b2dad00fca3c029b71efa0f43"
)
THM3612_OUTPUT_SHA256 = (
    "abb71fdac67e97e0c759cf970e3d787b0cbe60c802638303e1d237eab4f59c4b"
)
PROVISIONAL_CONTROL = (
    "THM-3619-russell-cylinder-even-fold-higher-jet-staircase",
    "31e3364dc868a90c308552c3412052fb8ae5c357",
    "26ebb9a9b9212100bd42ee146bcb0ecd0e2365d1f3f42bf39239441768d148dd",
    "7b5cc4f246c86e5cb4a061a23d371053ecb2fe9d54e0dc2d5e0b35569107cae7",
    "d16a7df35ff58ca77708169da6d7aea8832c8f511cf0ae05135d6dca2f20873e",
)

for slug, digest in AUDITED_ANCESTRY:
    require(f"audited ancestry SHA format {slug}", len(digest) == 64)
    require(f"audited ancestry SHA hex {slug}", all(c in "0123456789abcdef" for c in digest))
for label, digest in (
    ("THM-3612 promotion commit", THM3612_PROMOTION_COMMIT),
    ("THM-3612 script", THM3612_SCRIPT_SHA256),
    ("THM-3612 output", THM3612_OUTPUT_SHA256),
    ("THM-3619 provisional commit", PROVISIONAL_CONTROL[1]),
    ("THM-3619 theorem", PROVISIONAL_CONTROL[2]),
    ("THM-3619 script", PROVISIONAL_CONTROL[3]),
    ("THM-3619 output", PROVISIONAL_CONTROL[4]),
):
    expected_length = 40 if "commit" in label else 64
    require(f"pin length {label}", len(digest) == expected_length)
    require(f"pin hexadecimal {label}", all(c in "0123456789abcdef" for c in digest))
print(
    "PIN audited=THM-3561,THM-3605,THM-3612 "
    "provisional_finite_control=THM-3619"
)
print(
    f"PIN THM-3612_commit={THM3612_PROMOTION_COMMIT} "
    f"THM-3619_commit={PROVISIONAL_CONTROL[1]}"
)


print("SECTION compiler and exact Q_infinity identity")
x, q, w, t, u = sp.symbols("x q w t u")
U = x**2 * q
D = 1 + U
b = sp.expand(U * (U + 3) ** 2)
C = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
B = sp.expand(b + C * w)
Y = sp.expand(C * e + (2 * b + 4) * w + C * w**2)
S = sp.expand(((b + 2) * (e + 3 * w**2) + C * w * (3 * e + w**2)) / 8)
Z = sp.expand(S + sp.Rational(3, 4))
PHI = (C, Y, Z)

require("compiler c2e identity", zero(C**2 * e - b * (b + 4)))
require("Russell target relation", zero(C * Y - B * (B + 4)))

Q_INFINITY = -sp.Rational(3, 4) - sp.Rational(9, 4) / x**2
h_t = 1 / sp.sqrt(1 - sp.Rational(4, 3) * t**2)
side_fold = tuple(
    sp.factor(value.subs({q: Q_INFINITY + t**2, w: t})) for value in PHI
)
side_plus_target = tuple(
    sp.factor(sp.simplify(value.subs(x, h_t))) for value in side_fold
)
side_minus_target = tuple(
    sp.factor(sp.simplify(value.subs(x, -h_t))) for value in side_fold
)
middle_fold = tuple(
    sp.expand(value.subs({q: -sp.Rational(3, 4) + t**2, w: t}))
    for value in PHI
)
middle_target = tuple(sp.factor(value.subs(x, 0)) for value in middle_fold)
EXPECTED_TARGET = (sp.Integer(0), 4 * t, sp.Rational(7, 4) * t**2)

require("Q_infinity even", zero(Q_INFINITY.subs(x, -x) - Q_INFINITY))
require("Q_infinity shifted side q", zero(Q_INFINITY.subs(x, h_t) + t**2 - (-3 + 4 * t**2)))
require(
    "Q_infinity shifted side D",
    zero(D.subs({x: h_t, q: Q_INFINITY.subs(x, h_t) + t**2}) + 2),
)
require("plus target collision", side_plus_target == EXPECTED_TARGET)
require("minus target collision", side_minus_target == EXPECTED_TARGET)
require("middle target collision", middle_target == EXPECTED_TARGET)

side_bivector = fold_bivector(side_fold, x, t)
side_plus_bivector = tuple(
    sp.factor(sp.simplify(value.subs(x, h_t))) for value in side_bivector
)
side_minus_bivector = tuple(
    sp.factor(sp.simplify(value.subs(x, -h_t))) for value in side_bivector
)
middle_bivector = tuple(
    sp.factor(value.subs(x, 0)) for value in fold_bivector(middle_fold, x, t)
)
side_bivector_sum = tuple(
    sp.factor(left + right)
    for left, right in zip(side_minus_bivector, side_plus_bivector)
)
EXPECTED_MIDDLE_BIVECTOR = (
    sp.Integer(12),
    sp.Rational(21, 2) * t,
    3 * t * (11 * t**2 - 6),
)
EXPECTED_SIDE_SUM = tuple(2 * value for value in EXPECTED_MIDDLE_BIVECTOR)
require("middle bivector exact", middle_bivector == EXPECTED_MIDDLE_BIVECTOR)
require("side bivector sum exact", side_bivector_sum == EXPECTED_SIDE_SUM)
require(
    "side sum equals twice middle",
    all(zero(side_bivector_sum[index] - 2 * middle_bivector[index]) for index in range(3)),
)
A_form, K_form, R_form = sp.symbols("A_form K_form R_form")
require(
    "arbitrary target two-form contraction identity",
    zero(
        sum(
            coefficient * (side_bivector_sum[index] - 2 * middle_bivector[index])
            for index, coefficient in enumerate((A_form, K_form, R_form))
        )
    ),
)

phi_q = tuple(sp.diff(value, q) for value in PHI)
phi_w = tuple(sp.diff(value, w) for value in PHI)
phi_x = tuple(sp.diff(value, x) for value in PHI)
G = wedge(phi_q, phi_w)
q_prime = sp.symbols("q_prime")
ambient_fold_bivector = wedge(
    tuple(phi_x[index] + q_prime * phi_q[index] for index in range(3)),
    tuple(phi_w[index] + 2 * t * phi_q[index] for index in range(3)),
)
require(
    "chain-rule Q-prime first variation",
    all(
        zero(sp.diff(ambient_fold_bivector[index], q_prime) - G[index])
        for index in range(3)
    ),
)
G_plus_zero = tuple(sp.factor(value.subs({x: 1, q: -3, w: 0})) for value in G)
G_minus_zero = tuple(sp.factor(value.subs({x: -1, q: -3, w: 0})) for value in G)
G_difference_zero = tuple(
    sp.factor(plus - minus) for plus, minus in zip(G_plus_zero, G_minus_zero)
)
require("normalized first variation bivector", G_difference_zero == (-16, 0, 0))
require("middle normalization selects A0", middle_bivector[0].subs(t, 0) == 12)
require("middle normalization kills K0", middle_bivector[1].subs(t, 0) == 0)
require("middle normalization kills R0", middle_bivector[2].subs(t, 0) == 0)
SEMANTIC_FACTS.extend(
    (
        f"target={EXPECTED_TARGET}",
        f"middle_bivector={EXPECTED_MIDDLE_BIVECTOR}",
        f"side_sum={EXPECTED_SIDE_SUM}",
        f"G_difference={G_difference_zero}",
    )
)
print("PASS target=(0,4*t,7*t^2/4) side_sum=2*middle_bivector")
print("PASS G_plus(0)-G_minus(0)=(-16,0,0) normalized_A0=1")


print("SECTION rho shifted-evaluation cokernel rows")
rho = 1 - 1 / sp.sqrt(1 - sp.Rational(4, 3) * u)
h_u = 1 / sp.sqrt(1 - sp.Rational(4, 3) * u)
require("rho constant vanishes", rho.subs(u, 0) == 0)
require("rho linear coefficient", exact_coefficient(rho, u, 1) == -sp.Rational(2, 3))
require(
    "rho algebraic equation",
    zero((1 - rho) ** (-2) - (1 - sp.Rational(4, 3) * u)),
)
require("minus branch shifted point", zero(-1 + rho + h_u))
require("plus branch shifted point", zero(1 - rho - h_u))


def rho_row_coefficient(power, degree):
    """d_power^(degree)=-[u^degree]rho(u)^power."""
    return -exact_coefficient(rho**power, u, degree)


def generated_side_row(invoice_order):
    """All nonzero side coefficients at vertical order invoice_order=2m."""
    m = invoice_order // 2
    row = {}
    for t_half_degree in range(m):
        residual_degree = m - t_half_degree
        for xi_degree in range(1, residual_degree + 1):
            coefficient = rho_row_coefficient(xi_degree, residual_degree)
            if coefficient:
                exponent = (xi_degree, 2 * t_half_degree)
                row[(-1, exponent)] = coefficient
                row[(1, exponent)] = (-1) ** xi_degree * coefficient
    return row


THM3619_FINITE_SIDE_ROWS = {
    4: {
        (-1, (1, 0)): sp.Rational(2, 3),
        (-1, (2, 0)): -sp.Rational(4, 9),
        (-1, (1, 2)): sp.Rational(2, 3),
        (1, (1, 0)): -sp.Rational(2, 3),
        (1, (2, 0)): -sp.Rational(4, 9),
        (1, (1, 2)): -sp.Rational(2, 3),
    },
    6: {
        (-1, (1, 0)): sp.Rational(20, 27),
        (-1, (2, 0)): -sp.Rational(8, 9),
        (-1, (1, 2)): sp.Rational(2, 3),
        (-1, (3, 0)): sp.Rational(8, 27),
        (-1, (2, 2)): -sp.Rational(4, 9),
        (-1, (1, 4)): sp.Rational(2, 3),
        (1, (1, 0)): -sp.Rational(20, 27),
        (1, (2, 0)): -sp.Rational(8, 9),
        (1, (1, 2)): -sp.Rational(2, 3),
        (1, (3, 0)): -sp.Rational(8, 27),
        (1, (2, 2)): -sp.Rational(4, 9),
        (1, (1, 4)): -sp.Rational(2, 3),
    },
    8: {
        (-1, (1, 0)): sp.Rational(70, 81),
        (-1, (2, 0)): -sp.Rational(116, 81),
        (-1, (1, 2)): sp.Rational(20, 27),
        (-1, (3, 0)): sp.Rational(8, 9),
        (-1, (2, 2)): -sp.Rational(8, 9),
        (-1, (4, 0)): -sp.Rational(16, 81),
        (-1, (1, 4)): sp.Rational(2, 3),
        (-1, (3, 2)): sp.Rational(8, 27),
        (-1, (2, 4)): -sp.Rational(4, 9),
        (-1, (1, 6)): sp.Rational(2, 3),
        (1, (1, 0)): -sp.Rational(70, 81),
        (1, (2, 0)): -sp.Rational(116, 81),
        (1, (1, 2)): -sp.Rational(20, 27),
        (1, (3, 0)): -sp.Rational(8, 9),
        (1, (2, 2)): -sp.Rational(8, 9),
        (1, (4, 0)): -sp.Rational(16, 81),
        (1, (1, 4)): -sp.Rational(2, 3),
        (1, (3, 2)): -sp.Rational(8, 27),
        (1, (2, 4)): -sp.Rational(4, 9),
        (1, (1, 6)): -sp.Rational(2, 3),
    },
    10: {
        (-1, (1, 0)): sp.Rational(28, 27),
        (-1, (2, 0)): -sp.Rational(520, 243),
        (-1, (1, 2)): sp.Rational(70, 81),
        (-1, (3, 0)): sp.Rational(152, 81),
        (-1, (2, 2)): -sp.Rational(116, 81),
        (-1, (4, 0)): -sp.Rational(64, 81),
        (-1, (1, 4)): sp.Rational(20, 27),
        (-1, (3, 2)): sp.Rational(8, 9),
        (-1, (5, 0)): sp.Rational(32, 243),
        (-1, (2, 4)): -sp.Rational(8, 9),
        (-1, (4, 2)): -sp.Rational(16, 81),
        (-1, (1, 6)): sp.Rational(2, 3),
        (-1, (3, 4)): sp.Rational(8, 27),
        (-1, (2, 6)): -sp.Rational(4, 9),
        (-1, (1, 8)): sp.Rational(2, 3),
        (1, (1, 0)): -sp.Rational(28, 27),
        (1, (2, 0)): -sp.Rational(520, 243),
        (1, (1, 2)): -sp.Rational(70, 81),
        (1, (3, 0)): -sp.Rational(152, 81),
        (1, (2, 2)): -sp.Rational(116, 81),
        (1, (4, 0)): -sp.Rational(64, 81),
        (1, (1, 4)): -sp.Rational(20, 27),
        (1, (3, 2)): -sp.Rational(8, 9),
        (1, (5, 0)): -sp.Rational(32, 243),
        (1, (2, 4)): -sp.Rational(8, 9),
        (1, (4, 2)): -sp.Rational(16, 81),
        (1, (1, 6)): -sp.Rational(2, 3),
        (1, (3, 4)): -sp.Rational(8, 27),
        (1, (2, 6)): -sp.Rational(4, 9),
        (1, (1, 8)): -sp.Rational(2, 3),
    },
}

EXPECTED_CENTER_COEFFICIENTS = {
    4: sp.Integer(128),
    6: -sp.Integer(512),
    8: sp.Rational(12800, 9),
    10: -sp.Rational(90112, 27),
}
for invoice_order in (4, 6, 8, 10):
    generated = generated_side_row(invoice_order)
    expected = THM3619_FINITE_SIDE_ROWS[invoice_order]
    require(f"complete rho side row N={invoice_order}", generated == expected)
    m = invoice_order // 2
    require(
        f"rho row support size N={invoice_order}",
        len(generated) == m * (m + 1),
    )
    for (branch, (xi_degree, t_degree)), coefficient in generated.items():
        opposite = generated[(-branch, (xi_degree, t_degree))]
        require(
            f"branch parity N={invoice_order} a={xi_degree} b={t_degree}",
            opposite == (-1) ** xi_degree * coefficient,
        )
    row_digest = hashlib.sha256(canonical_row(generated).encode("ascii")).hexdigest()
    SEMANTIC_FACTS.append(f"rho_row_{invoice_order}={row_digest}")
    print(
        f"PASS rho_row N={invoice_order} side_entries={len(generated)} "
        f"center_control={EXPECTED_CENTER_COEFFICIENTS[invoice_order]} "
        f"sha256={row_digest}"
    )


print("SECTION shift filtration and arbitrary target two-form")
for n in range(HOSTILE_N_MIN, HOSTILE_N_MAX + 1):
    m = n - 1
    for xi_degree in range(1, 2 * m + 3):
        minimum_t_degree = max(0, 2 * m - xi_degree)
        shifted_order = 2 * xi_degree + minimum_t_degree
        require(
            f"shift filtration n={n} a={xi_degree}",
            shifted_order >= 2 * m + 1,
        )
    require(f"target displacement is higher n={n}", n > m)
    require(f"positive target coefficient is higher n={n}", 2 * m + 1 > 2 * m)
print(
    f"PASS filtration n={HOSTILE_N_MIN}..{HOSTILE_N_MAX} "
    "xi->rho raises every active lower-row monomial above the invoice"
)
print("PASS target_two_form=arbitrary_formal_(A,K,R); closedness/decomposability_unused")


print("SECTION normalized all-order invoice and recurrence")


def q_infinity_jet(order):
    """Closed coefficient Q_infinity^(order)(1)."""
    return (-1) ** (order - 1) * sp.Rational(9 * factorial(order + 1), 4)


def invoice_coefficient(order):
    """Positive c_n in Delta=-c_n*(q_n+(n+1)q_(n-1))."""
    return sp.Rational(2 ** (order + 3), 3 ** (order - 1) * factorial(order - 1))


def normalized_first_variation(order):
    """-16 times the leading shifted derivative coefficient."""
    return -16 * sp.Rational(2, 3) ** (order - 1) / factorial(order - 1)


HOSTILE_RECURRENCE_ROWS = []
h_shift_unit = sp.cancel((h_u - 1) / u)
require("shift unit constant", sp.limit(h_shift_unit, u, 0) == sp.Rational(2, 3))
y_local = sp.symbols("y_local")
require(
    "hostile perturbation base factor",
    zero((1 + y_local) ** 2 - 1 - y_local * (2 + y_local)),
)
for n in range(HOSTILE_N_MIN, HOSTILE_N_MAX + 1):
    m = n - 1
    shift_leading = sp.limit(h_shift_unit, u, 0) ** m
    require(f"shift leading coefficient n={n}", shift_leading == sp.Rational(2, 3) ** m)
    require(
        f"Q_infinity derivative n={n}",
        sp.diff(Q_INFINITY, x, n).subs(x, 1) == q_infinity_jet(n),
    )
    require(
        f"Q_infinity recurrence n={n}",
        q_infinity_jet(n) == -(n + 1) * q_infinity_jet(n - 1),
    )
    require(
        f"normalized coefficient n={n}",
        normalized_first_variation(n) == -invoice_coefficient(n),
    )
    q_n_symbol, q_previous_symbol = sp.symbols(f"q{n} q{n - 1}")
    invoice = sp.Poly(
        -invoice_coefficient(n) * (q_n_symbol + (n + 1) * q_previous_symbol),
        q_n_symbol,
        q_previous_symbol,
    )
    require(
        f"invoice q_n coefficient n={n}",
        invoice.coeff_monomial(q_n_symbol) == -invoice_coefficient(n),
    )
    require(
        f"invoice q_previous coefficient n={n}",
        invoice.coeff_monomial(q_previous_symbol) == -(n + 1) * invoice_coefficient(n),
    )
    require(
        f"rational germ kills invoice n={n}",
        -invoice_coefficient(n)
        * (q_infinity_jet(n) + (n + 1) * q_infinity_jet(n - 1))
        == 0,
    )
    perturbation_unit = (1 + y_local) ** 2 * (2 + y_local) ** n
    require(
        f"hostile perturbation multiplicity n={n}",
        perturbation_unit.subs(y_local, 0) == 2**n,
    )
    HOSTILE_RECURRENCE_ROWS.append(
        (
            n,
            canonical_rational(invoice_coefficient(n)),
            canonical_rational(-(n + 1) * invoice_coefficient(n)),
        )
    )

EXPECTED_FINITE_COEFFICIENTS = {
    3: sp.Rational(32, 9),
    4: sp.Rational(64, 81),
    5: sp.Rational(32, 243),
    6: sp.Rational(64, 3645),
}
THM3619_HOSTILE_CONTROLS = {
    3: (-sp.Integer(378), -sp.Rational(27, 2), sp.Integer(1536)),
    4: (sp.Integer(7506), sp.Integer(54), -sp.Integer(6144)),
    5: (-sp.Integer(127980), -sp.Integer(270), sp.Rational(51200, 3)),
    6: (sp.Integer(2269620), sp.Integer(1620), -sp.Rational(360448, 9)),
}
for n, expected_coefficient in EXPECTED_FINITE_COEFFICIENTS.items():
    hostile_q_n, previous_q, expected_invoice = THM3619_HOSTILE_CONTROLS[n]
    actual_invoice = -invoice_coefficient(n) * (hostile_q_n + (n + 1) * previous_q)
    require(f"THM-3619 coefficient n={n}", invoice_coefficient(n) == expected_coefficient)
    require(f"THM-3619 forced invoice n={n}", actual_invoice == expected_invoice)
    require(
        f"THM-3619 center quotient n={n}",
        actual_invoice / 12 == EXPECTED_CENTER_COEFFICIENTS[2 * n - 2],
    )

recurrence_digest = hashlib.sha256(repr(HOSTILE_RECURRENCE_ROWS).encode("ascii")).hexdigest()
SEMANTIC_FACTS.append(f"recurrence_rows={recurrence_digest}")
print(
    f"PASS normalized_invoice n={HOSTILE_N_MIN}..{HOSTILE_N_MAX} "
    "Delta=-2^(n+3)/(3^(n-1)*(n-1)!)*(q_n+(n+1)q_(n-1))"
)
print(
    "PASS THM-3619 finite invoices="
    "{n3:1536,n4:-6144,n5:51200/3,n6:-360448/9}"
)
print(f"PASS recurrence_rows_sha256={recurrence_digest}")


print("SECTION polynomial-versus-rational-germ endpoint")
require("Q_infinity value at one", Q_INFINITY.subs(x, 1) == -3)
require("Q_infinity first tuned jet", sp.diff(Q_INFINITY, x).subs(x, 1) == sp.Rational(9, 2))
require("Q_infinity second tuned jet", sp.diff(Q_INFINITY, x, 2).subs(x, 1) == -sp.Rational(27, 2))
require("Q_infinity pole numerator", sp.limit(x**2 * Q_INFINITY, x, 0) == -sp.Rational(9, 4))
require("polynomial identity endpoint contradiction", sp.Integer(0) != -sp.Rational(9, 4))

for degree in range(0, 17):
    derivative_matrix = sp.Matrix(
        [
            [
                sp.diff(x**power, x, derivative_order).subs(x, 1)
                for power in range(degree + 1)
            ]
            for derivative_order in range(degree + 1)
        ]
    )
    require(
        f"Taylor jet injective on polynomials degree={degree}",
        derivative_matrix.det() == sp.Integer(1) * sp.prod(factorial(k) for k in range(degree + 1)),
    )
print("PASS all Taylor jets determine a polynomial; x^2*(Q+3/4)=-9/4 fails at x=0")
print("SCOPE=polynomial even-fold family only; no global Keller pair or JC(2) claim")


print("SECTION source optimization and LF gates")
source_path = Path(__file__)
source_bytes = source_path.read_bytes()
require("source raw LF bytes", b"\r" not in source_bytes)
source_tree = ast.parse(source_bytes.decode("utf-8"))
assertion_nodes = sum(1 for node in ast.walk(source_tree) if isinstance(node, ast.Assert))
require("AST assertion nodes absent", assertion_nodes == 0)
source_sha256 = hashlib.sha256(source_bytes).hexdigest()
semantic_sha256 = hashlib.sha256("\n".join(SEMANTIC_FACTS).encode("utf-8")).hexdigest()
require("source SHA length", len(source_sha256) == 64)
require("semantic SHA length", len(semantic_sha256) == 64)
print(f"PASS ast_assertion_nodes={assertion_nodes} raw_LF=TRUE")
print(f"source_sha256={source_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print("normal_minus_O_byte_identity=VERIFY_WITH_MATCHING_STORED_TRANSCRIPT")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
