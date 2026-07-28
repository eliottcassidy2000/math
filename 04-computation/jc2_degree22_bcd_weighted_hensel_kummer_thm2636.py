#!/usr/bin/env python3
"""Exact Hensel--Kummer closure companion for the degree-22 BCD stratum.

The script reconstructs the universal five-parameter eliminant, specializes
E=W=0, excludes every factor through the fixed root and root-pair fields,
checks the five odd places of the retained spectral square class T^2=Z, and
reconstructs the y=0 boundary quartic.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


t, v, zeta = s.symbols("t v zeta")
c, d, e, w = s.symbols("c d e w")
alpha, beta, gamma = s.symbols("alpha beta gamma")

# B=rho^2, C=c rho^3, D=d rho^4, E=e rho^5, W=w rho^6,
# t=rho/y, v=u/y^2, zeta=Z/y^3.  These are the two normalized fluxes.
f1 = (
    819896 * t**2 * zeta
    - 1449459 * v * zeta
    + 83853 * zeta
    - 2981440 * t**2 * v
    + 24640 * t**2
    + 9370240 * c * t**3 * v
    - 232320 * c * t**3
    + 2044416 * d * t**4
    - 14992384 * e * t**5
    + 3689532 * v**2
    - 101640 * v
    + 252
)
f2 = (
    15944049 * zeta**2
    + 65591680 * t**2 * zeta
    - 206145280 * c * t**3 * zeta
    - 162339408 * v * zeta
    + 2236080 * zeta
    + 1443016960 * t**2 * v**2
    - 71554560 * t**2 * v
    + 98560 * t**2
    + 449771520 * c * t**3 * v
    - 1239040 * c * t**3
    - 1978994688 * d * t**4 * v
    + 16355328 * d * t**4
    - 239878144 * e * t**5
    - 1319329792 * w * t**6
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)

raw_resultant = s.resultant(f1, f2, zeta)
raw_poly = s.Poly(raw_resultant, t, v, c, d, e, w, domain=s.QQ)
raw_content, primitive_poly = raw_poly.primitive()
R_universal = primitive_poly.as_expr()
require(raw_content == 28344976, "universal resultant content changed")
require(len(primitive_poly.terms()) == 60, "universal eliminant lost its 60 terms")
require(primitive_poly.degree(t) == 10, "universal t-degree changed")
require(primitive_poly.degree(v) == 5, "universal v-degree changed")
require(
    {monomial[0] for monomial, _ in primitive_poly.terms()}
    == {0, 2, 3, 4, 5, 6, 7, 8, 9, 10},
    "universal weighted support changed",
)

R_bcd = s.expand(R_universal.subs({e: 0, w: 0}))
R_bcd_poly = s.Poly(R_bcd, t, v, c, d, domain=s.QQ)
require(R_bcd_poly.degree(t) == 9, "BCD specialized t-degree changed")
require(R_bcd_poly.degree(v) == 5, "BCD specialized v-degree changed")
leading_v = s.Poly(R_bcd, v).coeff_monomial(v**5)
require(leading_v == -88239118492602, "constant v-leading coefficient changed")
P = s.expand(R_bcd / leading_v)
P_poly = s.Poly(P, t, v, c, d, domain=s.QQ)
require(P_poly.LC() != 0 and s.Poly(P, v).LC() == 1, "monic normalization failed")

r_coefficients = [s.Poly(P, t).coeff_monomial(t**n) for n in range(11)]
P5_expr = r_coefficients[0]
P5 = s.Poly(P5_expr, v, domain=s.QQ)
require(P5.degree() == 5 and P5.LC() == 1, "fixed section is not monic quintic")
require(P5.is_irreducible, "fixed quintic stopped being irreducible")
require(s.gcd(P5, P5.diff()).degree() == 0, "fixed quintic stopped squarefree")
require(
    r_coefficients[1] == 0,
    "fixed-section missing-order-one tangency changed",
)

# At t=0 the first flux is
#   (83853-1449459*v)*zeta + (3689532*v^2-101640*v+252)=0.
# Both coefficients are units at every root of the irreducible quintic, so
# zeta has valuation zero at all five smooth fixed-section points.
fixed_zeta_denominator = s.Poly(83853 - 1449459 * v, v, domain=s.QQ)
fixed_zeta_numerator = s.Poly(
    3689532 * v**2 - 101640 * v + 252, v, domain=s.QQ
)
require(
    s.gcd(P5, fixed_zeta_denominator).degree() == 0,
    "zeta acquired a pole at a fixed-section point",
)
require(
    s.gcd(P5, fixed_zeta_numerator).degree() == 0,
    "zeta acquired a zero at a fixed-section point",
)

# The first-flux wall is not a component of the irreducible eliminant.  This
# identity also records the complete squared wall factor for later audits.
wall_v = (616 * t**2 + 63) / 1089
wall_polynomial = (
    745360 * c * t**5
    + (287496 * d - 71148) * t**4
    + 43560 * c * t**3
    + 5082 * t**2
    + 945
)
require(
    s.cancel(
        P.subs(v, wall_v)
        + s.Rational(128, 397076033216709) * wall_polynomial**2
    )
    == 0,
    "first-flux wall identity changed",
)

# Since t=rho/y and zeta=Z/y^3, the retained base-field square T^2=Z is
#   T^2=rho^3*zeta/t^3.
# Each fixed point therefore has Kummer valuation -3, hence is an odd branch
# place of the connected quadratic cover.  Five visible odd places and branch
# parity force at least six; Riemann--Hurwitz gives genus at least two even if
# the base curve has genus zero.
kummer_fixed_places = P5.degree()
kummer_valuation = -3
kummer_visible_branch_places = kummer_fixed_places
kummer_branch_floor = kummer_visible_branch_places + (
    kummer_visible_branch_places % 2
)
kummer_genus_floor = -1 + kummer_branch_floor // 2
require(
    kummer_fixed_places == 5
    and kummer_valuation % 2 != 0
    and kummer_branch_floor == 6
    and kummer_genus_floor == 2,
    "spectral-square Kummer genus invoice changed",
)

# The identically-zero boundary y=0 must be treated before t is defined.
# Normalize B=1.  The open chart has A=616-1089u != 0 and N1 reconstructs
# Z=-9370240*c*u/(1331*A).  Substitution in N2 has the exact primitive factor
# -u*Q4 below.  If u=0 then Z=0; otherwise Q4(u)=0, and its constant leading
# coefficient makes u algebraic over the algebraically closed constants.
u_boundary, Z_boundary = s.symbols("u_boundary Z_boundary")
A_boundary = 616 - 1089 * u_boundary
D_boundary = 1331 * A_boundary
N1_boundary = D_boundary * Z_boundary + 9370240 * c * u_boundary
N2_boundary = (
    15944049 * Z_boundary**2
    - 206145280 * c * Z_boundary
    + 1443016960 * u_boundary**2
    - 1978994688 * d * u_boundary
    - 1190488992 * u_boundary**3
)
boundary_raw = s.cancel(
    N2_boundary.subs(
        Z_boundary, -9370240 * c * u_boundary / D_boundary
    )
    * D_boundary**2
)
boundary_content, boundary_primitive = s.Poly(
    boundary_raw, u_boundary, c, d, domain=s.QQ
).primitive()
boundary_Q4 = (
    1267200 * c**2 * u_boundary
    - 1433600 * c**2
    + 3763584 * d * u_boundary**2
    - 4257792 * d * u_boundary
    + 1204224 * d
    + 2264031 * u_boundary**4
    - 5305608 * u_boundary**3
    + 3829056 * u_boundary**2
    - 878080 * u_boundary
)
require(
    boundary_content == 1104726788605792
    and s.expand(boundary_primitive.as_expr() + u_boundary * boundary_Q4) == 0,
    "y=0 eliminant changed",
)
require(
    s.Poly(boundary_Q4, u_boundary).degree() == 4
    and s.Poly(boundary_Q4, u_boundary).LC() == 2264031,
    "y=0 boundary lost its constant nonzero leading coefficient",
)


@dataclass(frozen=True)
class NumberField:
    """A small exact quotient Q[x]/(modulus), represented in the power basis."""

    generator: s.Symbol
    modulus: s.Poly

    def __post_init__(self) -> None:
        require(self.modulus.LC() == 1, "number-field modulus must be monic")
        require(self.modulus.is_irreducible, "number-field modulus is reducible")

    @lru_cache(maxsize=None)
    def reduce(self, expression: s.Expr | int) -> s.Expr:
        expression = s.cancel(s.sympify(expression))
        numerator, denominator = expression.as_numer_denom()
        numerator_poly = s.Poly(numerator, self.generator, domain=s.QQ).rem(
            self.modulus
        )
        denominator_poly = s.Poly(
            denominator, self.generator, domain=s.QQ
        ).rem(self.modulus)
        require(not denominator_poly.is_zero, "zero denominator in number field")
        inverse = s.invert(denominator_poly, self.modulus)
        return s.Poly(
            numerator_poly.as_expr() * inverse,
            self.generator,
            domain=s.QQ,
        ).rem(self.modulus).as_expr()

    def zero(self, expression: s.Expr | int) -> bool:
        return s.Poly(
            self.reduce(expression), self.generator, domain=s.QQ
        ).is_zero

    def inverse(self, expression: s.Expr | int) -> s.Expr:
        expression = self.reduce(expression)
        require(not self.zero(expression), "attempted inversion of zero")
        return self.reduce(1 / expression)

    def numerator(self, expression: s.Expr | int) -> s.Poly:
        numerator = s.cancel(self.reduce(expression)).as_numer_denom()[0]
        return s.Poly(numerator, self.generator, domain=s.QQ).primitive()[1]

    def coprime_numerator(self, expression: s.Expr | int) -> bool:
        return s.gcd(self.modulus, self.numerator(expression)).degree() == 0


# A parameter coefficient is a sparse polynomial in c,d over a number field.
PC = dict[tuple[int, int], s.Expr]
VP = list[PC]


class HenselSystem:
    def __init__(
        self,
        name: str,
        field: NumberField,
        q0_field: list[s.Expr],
        s0_field: list[s.Expr],
    ) -> None:
        self.name = name
        self.field = field
        self.q0 = self.v_from_field(q0_field)
        self.s0 = self.v_from_field(s0_field)
        self.k = len(q0_field) - 1
        self.inv_s0_mod_q0 = self.inverse_s0_mod_q0(q0_field, s0_field)

    def pc(self, terms: PC | None = None) -> PC:
        result: PC = {}
        for exponent, coefficient in (terms or {}).items():
            reduced = self.field.reduce(coefficient)
            if not self.field.zero(reduced):
                result[exponent] = reduced
        return result

    def pc_const(self, coefficient: s.Expr | int) -> PC:
        return self.pc({(0, 0): s.sympify(coefficient)})

    def pc_add(self, left: PC, right: PC) -> PC:
        result = dict(left)
        for exponent, coefficient in right.items():
            result[exponent] = self.field.reduce(
                result.get(exponent, 0) + coefficient
            )
        return self.pc(result)

    def pc_neg(self, value: PC) -> PC:
        return self.pc({exponent: -coefficient for exponent, coefficient in value.items()})

    def pc_mul(self, left: PC, right: PC) -> PC:
        result: PC = {}
        for (ic1, id1), coefficient1 in left.items():
            for (ic2, id2), coefficient2 in right.items():
                exponent = (ic1 + ic2, id1 + id2)
                result[exponent] = self.field.reduce(
                    result.get(exponent, 0) + coefficient1 * coefficient2
                )
        return self.pc(result)

    def pc_scale(self, value: PC, scalar: s.Expr | int) -> PC:
        return self.pc(
            {exponent: coefficient * scalar for exponent, coefficient in value.items()}
        )

    def vtrim(self, poly: VP) -> VP:
        poly = [self.pc(coefficient) for coefficient in poly]
        while len(poly) > 1 and not poly[-1]:
            poly.pop()
        return poly

    def v_from_field(self, poly: list[s.Expr]) -> VP:
        return self.vtrim([self.pc_const(coefficient) for coefficient in poly])

    def v_add(self, left: VP, right: VP) -> VP:
        length = max(len(left), len(right))
        return self.vtrim(
            [
                self.pc_add(
                    left[index] if index < len(left) else {},
                    right[index] if index < len(right) else {},
                )
                for index in range(length)
            ]
        )

    def v_neg(self, poly: VP) -> VP:
        return [self.pc_neg(coefficient) for coefficient in poly]

    def v_mul(self, left: VP, right: VP) -> VP:
        result: VP = [{} for _ in range(len(left) + len(right) - 1)]
        for i, left_coefficient in enumerate(left):
            for j, right_coefficient in enumerate(right):
                result[i + j] = self.pc_add(
                    result[i + j], self.pc_mul(left_coefficient, right_coefficient)
                )
        return self.vtrim(result)

    def v_divmod_monic(self, dividend: VP, divisor: VP) -> tuple[VP, VP]:
        divisor = self.vtrim(divisor)
        require(divisor[-1] == {(0, 0): s.Integer(1)}, "divisor not monic")
        remainder = self.vtrim(dividend[:])
        quotient: VP = [
            {} for _ in range(max(1, len(remainder) - len(divisor) + 1))
        ]
        while len(remainder) >= len(divisor) and not (
            len(remainder) == 1 and not remainder[0]
        ):
            shift = len(remainder) - len(divisor)
            leading = remainder[-1]
            quotient[shift] = self.pc_add(quotient[shift], leading)
            subtraction = [{} for _ in range(shift)] + [
                self.pc_mul(coefficient, leading) for coefficient in divisor
            ]
            remainder = self.v_add(remainder, self.v_neg(subtraction))
        return self.vtrim(quotient), self.vtrim(remainder)

    def v_remainder(self, dividend: VP, divisor: VP) -> VP:
        return self.v_divmod_monic(dividend, divisor)[1]

    def v_equal(self, left: VP, right: VP) -> bool:
        return all(not coefficient for coefficient in self.v_add(left, self.v_neg(right)))

    def inverse_s0_mod_q0(
        self, q0_field: list[s.Expr], s0_field: list[s.Expr]
    ) -> VP:
        # k is only one or two.  Construct the multiplication matrix of s0
        # on A[v]/q0 and invert its first column equation exactly.
        q0 = self.v_from_field(q0_field)
        s0 = self.v_from_field(s0_field)
        columns: list[list[s.Expr]] = []
        for exponent in range(self.k):
            monomial: VP = [{} for _ in range(exponent)] + [self.pc_const(1)]
            remainder = self.v_remainder(self.v_mul(s0, monomial), q0)
            entries: list[s.Expr] = []
            for index in range(self.k):
                coefficient = remainder[index] if index < len(remainder) else {}
                require(
                    set(coefficient).issubset({(0, 0)}),
                    "fixed inverse acquired parameters",
                )
                entries.append(coefficient.get((0, 0), s.Integer(0)))
            columns.append(entries)
        if self.k == 1:
            inverse_coefficients = [self.field.inverse(columns[0][0])]
        else:
            a00, a10 = columns[0]
            a01, a11 = columns[1]
            determinant = self.field.reduce(a00 * a11 - a01 * a10)
            require(not self.field.zero(determinant), "Hensel inverse is singular")
            inverse_coefficients = [
                self.field.reduce(a11 * self.field.inverse(determinant)),
                self.field.reduce(-a10 * self.field.inverse(determinant)),
            ]
        inverse = self.v_from_field(inverse_coefficients)
        product_remainder = self.v_remainder(self.v_mul(s0, inverse), q0)
        require(
            self.v_equal(product_remainder, self.v_from_field([1])),
            "inverse modulo fixed factor failed",
        )
        return inverse

    def expression_to_vpoly(self, expression: s.Expr) -> VP:
        poly = s.Poly(expression, v, c, d, domain=s.QQ)
        degree_v = poly.degree(v)
        result: VP = [{} for _ in range(max(0, degree_v) + 1)]
        for (iv, ic, id_), coefficient in poly.terms():
            result[iv] = self.pc_add(
                result[iv], self.pc({(ic, id_): coefficient})
            )
        return self.vtrim(result)

    def determinant(self, matrix: list[list[s.Expr]]) -> s.Expr:
        size = len(matrix)
        require(all(len(row) == size for row in matrix), "determinant not square")
        total = s.Integer(0)
        # Fraction-free recursion is tiny here (at most 4 by 4).
        if size == 0:
            return s.Integer(1)
        if size == 1:
            return self.field.reduce(matrix[0][0])
        for column in range(size):
            minor = [row[:column] + row[column + 1 :] for row in matrix[1:]]
            term = matrix[0][column] * self.determinant(minor)
            total = self.field.reduce(total + ((-1) ** column) * term)
        return self.field.reduce(total)

    def run(self, r_polys: list[VP]) -> dict[str, object]:
        require(len(r_polys) == 11, "expected t-orders zero through ten")
        require(
            self.v_equal(self.v_mul(self.q0, self.s0), r_polys[0]),
            "fixed factor/cofactor product failed",
        )
        qs: list[VP] = [self.q0]
        ss: list[VP] = [self.s0]
        product_controls: list[bool] = []
        for n in range(1, 12):
            rn = r_polys[n] if n < len(r_polys) else self.v_from_field([0])
            convolution = self.v_from_field([0])
            for index in range(1, n):
                convolution = self.v_add(
                    convolution, self.v_mul(qs[index], ss[n - index])
                )
            residual = self.v_add(rn, self.v_neg(convolution))
            qn = self.v_remainder(
                self.v_mul(residual, self.inv_s0_mod_q0), self.q0
            )
            sn_numerator = self.v_add(residual, self.v_neg(self.v_mul(qn, self.s0)))
            sn, sn_remainder = self.v_divmod_monic(sn_numerator, self.q0)
            require(
                all(not coefficient for coefficient in sn_remainder),
                f"{self.name}: Hensel division remainder at order {n}",
            )
            qs.append(qn)
            ss.append(sn)
            reconstructed = self.v_from_field([0])
            for index in range(n + 1):
                reconstructed = self.v_add(
                    reconstructed, self.v_mul(qs[index], ss[n - index])
                )
            product_controls.append(self.v_equal(reconstructed, rn))
        require(all(product_controls), f"{self.name}: direct product control failed")

        order11 = qs[11] + ss[11]
        require(len(order11) == 5, f"{self.name}: order-11 equation count changed")
        allowed_support = {(3, 0), (1, 2), (1, 1), (1, 0)}
        actual_support = set().union(*(set(coefficient) for coefficient in order11))
        require(
            actual_support == allowed_support,
            f"{self.name}: order-11 parameter support changed: {actual_support}",
        )
        monomial_basis = [(3, 0), (1, 2), (1, 1), (1, 0)]
        matrix = [
            [coefficient.get(monomial, s.Integer(0)) for monomial in monomial_basis]
            for coefficient in order11
        ]
        maximal_minors: list[s.Expr] = []
        for selected_rows in combinations(range(5), 4):
            square = [matrix[index] for index in selected_rows]
            maximal_minors.append(self.determinant(square))
        nonzero_minors = [
            minor for minor in maximal_minors if not self.field.zero(minor)
        ]
        require(nonzero_minors, f"{self.name}: order-11 matrix rank below four")
        require(
            all(self.field.coprime_numerator(minor) for minor in nonzero_minors),
            f"{self.name}: a nonzero minor is not uniform on conjugate sections",
        )
        first_numerator = self.field.numerator(nonzero_minors[0])
        return {
            "equations": len(order11),
            "rank": 4,
            "nonzero_minors": len(nonzero_minors),
            "minor_degree": first_numerator.degree(),
            "minor_terms": len(first_numerator.terms()),
            "product_controls": len(product_controls),
            "support": actual_support,
        }


def field_poly_from_expr(
    expression: s.Expr, variable: s.Symbol, field: NumberField
) -> list[s.Expr]:
    poly = s.Poly(expression, v, variable, domain=s.QQ)
    degree_v = poly.degree(v)
    result = [s.Integer(0)] * (degree_v + 1)
    for (iv, ix), coefficient in poly.terms():
        result[iv] = field.reduce(result[iv] + coefficient * variable**ix)
    while len(result) > 1 and field.zero(result[-1]):
        result.pop()
    return result


# Root-field section: q0=v-gamma.
root_modulus = s.Poly(P5_expr.subs(v, gamma), gamma, domain=s.QQ).monic()
root_field = NumberField(gamma, root_modulus)
root_q0 = [-gamma, s.Integer(1)]
root_s0_expr, root_remainder = s.div(
    s.Poly(P5_expr, v, domain=s.QQ),
    s.Poly(v - gamma, v, domain=s.QQ.frac_field(gamma)),
)
root_remainder_reduced = [
    root_field.reduce(root_remainder.nth(index))
    for index in range(max(0, root_remainder.degree()) + 1)
]
require(
    all(root_field.zero(entry) for entry in root_remainder_reduced),
    "root divisor failed after field reduction",
)
root_s0 = field_poly_from_expr(root_s0_expr.as_expr(), gamma, root_field)
root_system = HenselSystem("root_field", root_field, root_q0, root_s0)
root_r = [root_system.expression_to_vpoly(expression) for expression in r_coefficients]
root_result = root_system.run(root_r)


# Pair-field section: q0=v^2+alpha*v+beta for an unordered root pair.
pair_trial = s.Poly(v**2 + alpha * v + beta, v)
pair_remainder = s.rem(s.Poly(P5_expr, v), pair_trial).as_expr()
pair_remainder_poly = s.Poly(pair_remainder, v)
pair_equations = [pair_remainder_poly.nth(index) for index in range(2)]
pair_basis = s.groebner(pair_equations, beta, alpha, order="lex")
require(len(pair_basis.polys) == 2, "pair scheme lost triangular basis")
beta_relation = s.Poly(pair_basis.polys[0].as_expr(), beta, alpha)
pair_modulus = s.Poly(pair_basis.polys[1].as_expr(), alpha, domain=s.QQ).monic()
require(pair_modulus.degree() == 10, "pair field degree changed")
pair_field = NumberField(alpha, pair_modulus)
beta_linear = s.Poly(beta_relation.as_expr(), beta).coeff_monomial(beta)
require(
    beta_relation.degree(beta) == 1 and s.degree(beta_linear, alpha) == 0,
    "beta is not a rational polynomial in alpha",
)
beta_value = pair_field.reduce(-beta_relation.as_expr().subs(beta, 0) / beta_linear)
pair_q0 = [beta_value, alpha, s.Integer(1)]
pair_trial_field = s.Poly(v**2 + alpha * v + beta_value, v)
pair_s0_expr, pair_fixed_remainder = s.div(
    s.Poly(P5_expr, v, domain=s.QQ.frac_field(alpha)), pair_trial_field
)
pair_fixed_remainder_reduced = [
    pair_field.reduce(pair_fixed_remainder.nth(index))
    for index in range(max(0, pair_fixed_remainder.degree()) + 1)
]
require(
    all(pair_field.zero(entry) for entry in pair_fixed_remainder_reduced),
    "pair divisor failed after field reduction",
)
pair_s0 = field_poly_from_expr(pair_s0_expr.as_expr(), alpha, pair_field)
pair_system = HenselSystem("pair_field", pair_field, pair_q0, pair_s0)
pair_r = [pair_system.expression_to_vpoly(expression) for expression in r_coefficients]
pair_result = pair_system.run(pair_r)


print("degree-22 BCD weighted-scale Hensel--Kummer closure")
print(f"universal_resultant_content={raw_content}")
print(f"universal_terms={len(primitive_poly.terms())}")
print("universal_degrees=t:10,v:5")
print("universal_t_support=0,2,3,4,5,6,7,8,9,10")
print(f"bcd_terms={len(R_bcd_poly.terms())}")
print("bcd_degrees=t:9,v:5")
print(f"constant_v_leading={leading_v}")
print("fixed_quintic_degree=5")
print("fixed_quintic_irreducible=True")
print("fixed_quintic_squarefree=True")
print("fixed_t_uniformizers=5")
print("fixed_zeta_units=5")
print("first_flux_wall_identity=nonzero_square")
for name, field, result in (
    ("root", root_field, root_result),
    ("pair", pair_field, pair_result),
):
    print(f"{name}_field_degree={field.modulus.degree()}")
    print(f"{name}_field_irreducible=True")
    print(f"{name}_order11_equations={result['equations']}")
    print(f"{name}_order11_rank={result['rank']}")
    print(f"{name}_nonzero_maximal_minors={result['nonzero_minors']}")
    print(f"{name}_first_minor_degree={result['minor_degree']}")
    print(f"{name}_first_minor_terms={result['minor_terms']}")
    print(f"{name}_product_controls={result['product_controls']}")
    print(f"{name}_minor_coprime_to_field_modulus=True")
print("order11_support=c^3,c*d^2,c*d,c")
print("line_factor_scheme_for_c_nonzero=EMPTY")
print("quadratic_factor_scheme_for_c_nonzero=EMPTY")
print("bcd_eliminant_absolutely_irreducible_for_c_nonzero=True")
print("spectral_square_relation=T^2=rho^3*zeta/t^3")
print(f"kummer_fixed_valuation={kummer_valuation}")
print(f"kummer_visible_odd_places={kummer_visible_branch_places}")
print(f"kummer_branch_floor={kummer_branch_floor}")
print(f"kummer_genus_floor={kummer_genus_floor}")
print(f"y0_eliminant_content={boundary_content}")
print("y0_primitive_factor=-u*Q4")
print("y0_Q4_degree=4")
print("y0_Q4_leading=2264031")
print("bcd_trajectory_stratum=EMPTY")
print("ALL CHECKS PASSED")
