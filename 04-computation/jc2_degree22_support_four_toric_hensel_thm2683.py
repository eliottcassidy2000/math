#!/usr/bin/env python3
"""Exact support-four Hensel/toric certificates downstream of THM-2671.

This companion reconstructs THM-2671's five-coefficient pre-scale eliminant,
but uses a genuine three-parameter sparse coefficient ring over each of the
fixed degree-five root and degree-ten unordered-root-pair fields.  It closes
the BDEW and CDEW charts by full terminal linear rank.  On BCDW the terminal
matrix has a one-dimensional kernel; its monomials are

    C * (1, D, D^2, C^2, W).

The kernel line violates the toric identity X0*X2=X1^2 uniformly in both
fields.  On BCEW a two-dimensional kernel misses the intersection of two
necessary toric quadrics by a nonzero Sylvester resultant.  On BCDE three
necessary toric quadrics have a degree-four Macaulay map of full row rank on
the three-dimensional terminal kernel.  No finite-field specialization or
encoded exponent is used.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import io
from itertools import combinations
from pathlib import Path

import sympy as s


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py"
BASE_SHA256 = "0866a29f665aedc6d2c226f35943852e56907ff821e705a0dbca2651e71fa15c"
PARENT_PATH = HERE / "jc2_degree22_complete_support_three_hensel_thm2671.py"
PARENT_SHA256 = "6245dd4cc85d0a70bdbc8e56a0511ffad7889b6274130aed759e5729f92472e6"
require(hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256,
        "audited THM-2636 dependency changed")
require(hashlib.sha256(PARENT_PATH.read_bytes()).hexdigest() == PARENT_SHA256,
        "audited THM-2671 parent changed")

spec = importlib.util.spec_from_file_location("thm2636_base_support4", BASE_PATH)
require(spec is not None and spec.loader is not None, "cannot load THM-2636")
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)

t, v, zeta = base.t, base.v, base.zeta
p, q, r = s.symbols("p q r")
B0, C0, D0, E0, W0 = s.symbols("B0 C0 D0 E0 W0")

# THM-2671's normalized fluxes before the coefficient-root scale is chosen.
F1 = (
    819896 * B0 * zeta
    - 1449459 * v * zeta
    + 83853 * zeta
    - 2981440 * B0 * v
    + 24640 * B0
    + 9370240 * C0 * v
    - 232320 * C0
    + 2044416 * D0
    - 14992384 * E0
    + 3689532 * v**2
    - 101640 * v
    + 252
)
F2 = (
    15944049 * zeta**2
    + 65591680 * B0 * zeta
    - 206145280 * C0 * zeta
    - 162339408 * v * zeta
    + 2236080 * zeta
    + 1443016960 * B0 * v**2
    - 71554560 * B0 * v
    + 98560 * B0
    + 449771520 * C0 * v
    - 1239040 * C0
    - 1978994688 * D0 * v
    + 16355328 * D0
    - 239878144 * E0
    - 1319329792 * W0
    - 1190488992 * v**3
    + 147581280 * v**2
    - 1219680 * v
    + 672
)
general_raw = s.resultant(F1, F2, zeta)
general_content, general_primitive = s.Poly(
    general_raw, B0, C0, D0, E0, W0, v, domain=s.QQ
).primitive()
R_GENERAL = general_primitive.as_expr()
require(general_content == 28344976, "general resultant content changed")
require(len(general_primitive.terms()) == 60,
        "general eliminant support changed")
require(
    s.expand(R_GENERAL.subs({
        B0: t**2,
        C0: base.c * t**3,
        D0: base.d * t**4,
        E0: base.e * t**5,
        W0: base.w * t**6,
    }) - base.R_universal) == 0,
    "pre-scale eliminant does not recover THM-2636",
)


class HenselSystem3(base.HenselSystem):
    """The audited Hensel engine over K[p,q,r], with sparse triple exponents."""

    ZERO = (0, 0, 0)

    def pc(self, terms=None):
        result = {}
        for exponent, coefficient in (terms or {}).items():
            require(len(exponent) == 3, "parameter exponent is not a triple")
            reduced = self.field.reduce(coefficient)
            if not self.field.zero(reduced):
                result[exponent] = reduced
        return result

    def pc_const(self, coefficient):
        return self.pc({self.ZERO: s.sympify(coefficient)})

    def pc_add(self, left, right):
        result = dict(left)
        for exponent, coefficient in right.items():
            result[exponent] = self.field.reduce(
                result.get(exponent, 0) + coefficient
            )
        return self.pc(result)

    def pc_neg(self, value):
        return self.pc({exponent: -coefficient
                        for exponent, coefficient in value.items()})

    def pc_mul(self, left, right):
        result = {}
        for exponent1, coefficient1 in left.items():
            for exponent2, coefficient2 in right.items():
                exponent = tuple(a + b for a, b in zip(exponent1, exponent2))
                result[exponent] = self.field.reduce(
                    result.get(exponent, 0) + coefficient1 * coefficient2
                )
        return self.pc(result)

    def pc_scale(self, value, scalar):
        return self.pc({exponent: coefficient * scalar
                        for exponent, coefficient in value.items()})

    def v_divmod_monic(self, dividend, divisor):
        divisor = self.vtrim(divisor)
        require(divisor[-1] == {self.ZERO: s.Integer(1)},
                "divisor is not monic")
        remainder = self.vtrim(dividend[:])
        quotient = [
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

    def inverse_s0_mod_q0(self, q0_field, s0_field):
        q0 = self.v_from_field(q0_field)
        s0 = self.v_from_field(s0_field)
        columns = []
        for exponent in range(self.k):
            monomial = [{} for _ in range(exponent)] + [self.pc_const(1)]
            remainder = self.v_remainder(self.v_mul(s0, monomial), q0)
            entries = []
            for index in range(self.k):
                coefficient = remainder[index] if index < len(remainder) else {}
                require(set(coefficient).issubset({self.ZERO}),
                        "fixed inverse acquired parameters")
                entries.append(coefficient.get(self.ZERO, s.Integer(0)))
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
        require(
            self.v_equal(
                self.v_remainder(self.v_mul(s0, inverse), q0),
                self.v_from_field([1]),
            ),
            "inverse modulo fixed factor failed",
        )
        return inverse

    def expression_to_vpoly(self, expression):
        poly = s.Poly(expression, v, p, q, r, domain=s.QQ)
        degree_v = poly.degree(v)
        result = [{} for _ in range(max(0, degree_v) + 1)]
        for (iv, ip, iq, ir), coefficient in poly.terms():
            result[iv] = self.pc_add(
                result[iv], self.pc({(ip, iq, ir): coefficient})
            )
        return self.vtrim(result)


def scaled_eliminant(substitutions):
    expression = s.expand(R_GENERAL.subs(substitutions))
    raw_poly = s.Poly(expression, t, v, p, q, r, domain=s.QQ)
    content, primitive = raw_poly.primitive()
    leading = s.Poly(primitive.as_expr(), v).coeff_monomial(v**5)
    require(leading != 0 and not leading.free_symbols,
            "v-leading term is not constant")
    monic = s.expand(primitive.as_expr() / leading)
    require(s.Poly(monic, v).LC() == 1, "monic normalization failed")
    degree_t = s.Poly(monic, t).degree()
    coefficients = [
        s.Poly(monic, t).coeff_monomial(t**n)
        for n in range(degree_t + 1)
    ]
    require(s.expand(coefficients[0] - base.P5_expr) == 0,
            "fixed quintic changed")
    return primitive, content, leading, degree_t, coefficients


def terminal_equations(system, coefficients, terminal):
    require(
        system.v_equal(
            system.v_mul(system.q0, system.s0),
            system.expression_to_vpoly(coefficients[0]),
        ),
        f"{system.name}: fixed factorization failed",
    )
    qs = [system.q0]
    ss = [system.s0]
    zero = system.v_from_field([0])
    for n in range(1, terminal + 1):
        rn = (system.expression_to_vpoly(coefficients[n])
              if n < len(coefficients) else zero)
        convolution = zero
        for index in range(1, n):
            convolution = system.v_add(
                convolution, system.v_mul(qs[index], ss[n - index])
            )
        residual = system.v_add(rn, system.v_neg(convolution))
        qn = system.v_remainder(
            system.v_mul(residual, system.inv_s0_mod_q0), system.q0
        )
        numerator = system.v_add(
            residual, system.v_neg(system.v_mul(qn, system.s0))
        )
        sn, remainder = system.v_divmod_monic(numerator, system.q0)
        require(all(not coefficient for coefficient in remainder),
                f"{system.name}: Hensel division failed at order {n}")
        qs.append(qn)
        ss.append(sn)
        reconstructed = zero
        for index in range(n + 1):
            reconstructed = system.v_add(
                reconstructed, system.v_mul(qs[index], ss[n - index])
            )
        require(system.v_equal(reconstructed, rn),
                f"{system.name}: product control failed at order {n}")
    equations = qs[terminal] + ss[terminal]
    require(len(equations) == 5,
            f"{system.name}: terminal equation count changed")
    return equations


def coefficient_matrix(equations, support):
    actual = set().union(*(set(coefficient) for coefficient in equations))
    require(actual == set(support), f"terminal support changed: {actual}")
    return [[coefficient.get(monomial, s.Integer(0)) for monomial in support]
            for coefficient in equations]


def full_rank_certificate(system, matrix):
    width = len(matrix[0])
    minors = []
    for rows in combinations(range(5), width):
        minors.append(system.determinant([matrix[index] for index in rows]))
    nonzero = [minor for minor in minors if not system.field.zero(minor)]
    require(nonzero, f"{system.name}: matrix lost full column rank")
    require(all(system.field.coprime_numerator(minor) for minor in nonzero),
            f"{system.name}: a maximal minor vanishes at a conjugate")
    numerator = system.field.numerator(nonzero[0])
    return len(nonzero), numerator.degree(), len(numerator.terms())


def bcdw_toric_certificate(system, matrix):
    require(len(matrix) == 5 and len(matrix[0]) == 5,
            "BCDW terminal matrix is not 5 by 5")
    selected_rows = None
    kernel = None
    for rows in combinations(range(5), 4):
        four = [matrix[index] for index in rows]
        candidate = []
        for column in range(5):
            minor = [row[:column] + row[column + 1:] for row in four]
            candidate.append(system.field.reduce(
                (-1) ** column * system.determinant(minor)
            ))
        if any(not system.field.zero(value) for value in candidate):
            selected_rows, kernel = rows, candidate
            break
    require(kernel is not None, f"{system.name}: terminal rank below four")
    require(all(
        system.field.zero(sum(row[index] * kernel[index]
                              for index in range(5)))
        for row in matrix
    ), f"{system.name}: cofactor vector is not in the full kernel")
    require(all(system.field.coprime_numerator(value) for value in kernel),
            f"{system.name}: a kernel coordinate vanishes at a conjugate")
    defect = system.field.reduce(kernel[0] * kernel[2] - kernel[1] ** 2)
    require(not system.field.zero(defect),
            f"{system.name}: kernel line lies in the toric quadric")
    require(system.field.coprime_numerator(defect),
            f"{system.name}: toric defect vanishes at a conjugate")
    numerator = system.field.numerator(defect)
    digest = hashlib.sha256(str(numerator.as_expr()).encode()).hexdigest()
    return {
        "rows": selected_rows,
        "degree": numerator.degree(),
        "terms": len(numerator.terms()),
        "digest": digest,
    }


def bcew_projective_toric_certificate(system, matrix):
    """Intersect the two-dimensional terminal kernel with two toric quadrics."""
    require(len(matrix) == 5 and len(matrix[0]) == 6,
            "BCEW terminal matrix is not 5 by 6")
    reduced = [[system.field.reduce(value) for value in row] for row in matrix]
    pivot_columns = []
    pivot_row = 0
    for column in range(6):
        pivot = next(
            (row for row in range(pivot_row, 5)
             if not system.field.zero(reduced[row][column])),
            None,
        )
        if pivot is None:
            continue
        reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
        pivot_value = reduced[pivot_row][column]
        require(system.field.coprime_numerator(pivot_value),
                f"{system.name}: pivot vanishes at a conjugate")
        inverse = system.field.inverse(pivot_value)
        reduced[pivot_row] = [
            system.field.reduce(value * inverse)
            for value in reduced[pivot_row]
        ]
        for row in range(5):
            if row == pivot_row or system.field.zero(reduced[row][column]):
                continue
            multiplier = reduced[row][column]
            reduced[row] = [
                system.field.reduce(
                    reduced[row][index]
                    - multiplier * reduced[pivot_row][index]
                )
                for index in range(6)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == 5:
            break
    require(pivot_columns == [0, 1, 2, 3],
            f"{system.name}: pivot pattern changed to {pivot_columns}")
    require(all(system.field.zero(value) for value in reduced[4]),
            f"{system.name}: terminal rank exceeds four")

    basis = []
    for free_column in (4, 5):
        vector = [s.Integer(0)] * 6
        vector[free_column] = s.Integer(1)
        for row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = system.field.reduce(
                -reduced[row][free_column]
            )
        require(all(
            system.field.zero(sum(matrix[row][column] * vector[column]
                                  for column in range(6)))
            for row in range(5)
        ), f"{system.name}: RREF basis is not in the full kernel")
        basis.append(vector)
    ku, kv = basis

    def binary_quadric(i, j, k, ell):
        coefficients = (
            system.field.reduce(ku[i] * ku[j] - ku[k] * ku[ell]),
            system.field.reduce(
                ku[i] * kv[j] + kv[i] * ku[j]
                - ku[k] * kv[ell] - kv[k] * ku[ell]
            ),
            system.field.reduce(kv[i] * kv[j] - kv[k] * kv[ell]),
        )
        require(all(not system.field.zero(value) for value in coefficients),
                f"{system.name}: a toric-quadric coefficient vanished")
        require(all(system.field.coprime_numerator(value)
                    for value in coefficients),
                f"{system.name}: a toric coefficient vanished at a conjugate")
        return coefficients

    # In support order (E,EW,C,CW,C^2E,C^3), every monomial point obeys
    # E*CW=EW*C and C*(C^2E)=E*C^3.
    f0, f1, f2 = binary_quadric(0, 3, 1, 2)
    g0, g1, g2 = binary_quadric(2, 4, 0, 5)
    sylvester = [
        [f0, f1, f2, s.Integer(0)],
        [s.Integer(0), f0, f1, f2],
        [g0, g1, g2, s.Integer(0)],
        [s.Integer(0), g0, g1, g2],
    ]
    resultant = system.field.reduce(system.determinant(sylvester))
    require(not system.field.zero(resultant),
            f"{system.name}: projective toric quadrics acquired a common root")
    require(system.field.coprime_numerator(resultant),
            f"{system.name}: Sylvester resultant vanishes at a conjugate")
    numerator = system.field.numerator(resultant)
    digest = hashlib.sha256(str(numerator.as_expr()).encode()).hexdigest()
    return {
        "rank": 4,
        "pivots": tuple(pivot_columns),
        "coefficient_uniform": True,
        "degree": numerator.degree(),
        "terms": len(numerator.terms()),
        "digest": digest,
    }


def homogeneous_exponents(total):
    return tuple(
        (i, j, total - i - j)
        for i in range(total, -1, -1)
        for j in range(total - i, -1, -1)
    )


def field_pivot_columns(system, matrix):
    """Exact Gaussian rank profile over the fixed number field."""
    reduced = [[system.field.reduce(value) for value in row] for row in matrix]
    rows = len(reduced)
    columns = len(reduced[0])
    pivot_columns = []
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows)
             if not system.field.zero(reduced[row][column])),
            None,
        )
        if pivot is None:
            continue
        reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
        inverse = system.field.inverse(reduced[pivot_row][column])
        reduced[pivot_row] = [
            system.field.reduce(value * inverse)
            for value in reduced[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row or system.field.zero(reduced[row][column]):
                continue
            multiplier = reduced[row][column]
            reduced[row] = [
                system.field.reduce(
                    reduced[row][index]
                    - multiplier * reduced[pivot_row][index]
                )
                for index in range(columns)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    return tuple(pivot_columns), reduced


def field_determinant(system, matrix):
    """Cubic-time exact determinant over the fixed number field."""
    size = len(matrix)
    require(all(len(row) == size for row in matrix),
            "field determinant is not square")
    work = [[system.field.reduce(value) for value in row] for row in matrix]
    determinant = s.Integer(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size)
             if not system.field.zero(work[row][column])),
            None,
        )
        if pivot is None:
            return s.Integer(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = system.field.reduce(-determinant)
        pivot_value = work[column][column]
        determinant = system.field.reduce(determinant * pivot_value)
        inverse = system.field.inverse(pivot_value)
        work[column] = [
            system.field.reduce(value * inverse) for value in work[column]
        ]
        for row in range(column + 1, size):
            if system.field.zero(work[row][column]):
                continue
            multiplier = work[row][column]
            work[row] = [
                system.field.reduce(
                    work[row][index] - multiplier * work[column][index]
                )
                for index in range(size)
            ]
    return system.field.reduce(determinant)


def bcde_macaulay_toric_certificate(system, matrix):
    """Show three kernel-plane toric quadrics generate every quartic."""
    require(len(matrix) == 5 and len(matrix[0]) == 7,
            "BCDE terminal matrix is not 5 by 7")
    pivots, reduced = field_pivot_columns(system, matrix)
    require(pivots == (0, 1, 2, 3),
            f"{system.name}: BCDE terminal pivots changed to {pivots}")
    require(all(system.field.zero(value) for value in reduced[4]),
            f"{system.name}: BCDE terminal rank exceeds four")

    basis = []
    for free_column in (4, 5, 6):
        vector = [s.Integer(0)] * 7
        vector[free_column] = s.Integer(1)
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = system.field.reduce(
                -reduced[row][free_column]
            )
        require(all(
            system.field.zero(sum(matrix[row][column] * vector[column]
                                  for column in range(7)))
            for row in range(5)
        ), f"{system.name}: BCDE kernel basis failed")
        basis.append(vector)

    def ternary_quadric(i, j, k, ell):
        coefficients = {}
        units = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        for left in range(3):
            for right in range(3):
                exponent = tuple(
                    units[left][index] + units[right][index]
                    for index in range(3)
                )
                value = (
                    basis[left][i] * basis[right][j]
                    - basis[left][k] * basis[right][ell]
                )
                coefficients[exponent] = system.field.reduce(
                    coefficients.get(exponent, 0) + value
                )
        coefficients = {
            exponent: value for exponent, value in coefficients.items()
            if not system.field.zero(value)
        }
        require(coefficients and set(coefficients).issubset(
                    set(homogeneous_exponents(2))
                ), f"{system.name}: BCDE quadric support is malformed")
        require(all(system.field.coprime_numerator(value)
                    for value in coefficients.values()),
                f"{system.name}: BCDE quadric lost conjugate uniformity")
        return coefficients

    # Support order (E,DE,C,CD,CD^2,C^2E,C^3).
    quadrics = (
        ternary_quadric(0, 3, 1, 2),
        ternary_quadric(2, 4, 3, 3),
        ternary_quadric(2, 5, 0, 6),
    )
    degree_two = homogeneous_exponents(2)
    degree_four = homogeneous_exponents(4)
    columns = []
    for quadric in quadrics:
        for multiplier in degree_two:
            columns.append({
                tuple(exponent[index] + multiplier[index]
                      for index in range(3)): coefficient
                for exponent, coefficient in quadric.items()
            })
    require(len(degree_four) == 15 and len(columns) == 18,
            "BCDE Macaulay dimensions changed")
    macaulay = [
        [column.get(monomial, s.Integer(0)) for column in columns]
        for monomial in degree_four
    ]
    macaulay_pivots, _ = field_pivot_columns(system, macaulay)
    expected_pivots = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15)
    require(macaulay_pivots == expected_pivots,
            f"{system.name}: Macaulay pivots changed to {macaulay_pivots}")
    square = [[row[column] for column in macaulay_pivots] for row in macaulay]
    minor = field_determinant(system, square)
    require(not system.field.zero(minor),
            f"{system.name}: degree-four Macaulay minor vanished")
    require(system.field.coprime_numerator(minor),
            f"{system.name}: Macaulay minor vanished at a conjugate")
    numerator = system.field.numerator(minor)
    digest = hashlib.sha256(str(numerator.as_expr()).encode()).hexdigest()
    return {
        "rank": 4,
        "pivots": pivots,
        "quadrics": 3,
        "macaulay_shape": (15, 18),
        "macaulay_pivots": macaulay_pivots,
        "degree": numerator.degree(),
        "terms": len(numerator.terms()),
        "digest": digest,
    }


CHARTS = (
    {
        "name": "BCDE",
        "substitution": {
            B0: t**2, C0: p * t**3, D0: q * t**4,
            E0: r * t**5, W0: 0,
        },
        "terms": 54,
        "support": (
            (0, 0, 1),  # E
            (0, 1, 1),  # D E
            (1, 0, 0),  # C
            (1, 1, 0),  # C D
            (1, 2, 0),  # C D^2
            (2, 0, 1),  # C^2 E
            (3, 0, 0),  # C^3
        ),
        "kind": "macaulay_toric",
    },
    {
        "name": "BCDW",
        "substitution": {
            B0: t**2, C0: p * t**3, D0: q * t**4,
            E0: 0, W0: r * t**6,
        },
        "terms": 50,
        "support": (
            (1, 0, 0),  # C
            (1, 1, 0),  # C D
            (1, 2, 0),  # C D^2
            (3, 0, 0),  # C^3
            (1, 0, 1),  # C W
        ),
        "kind": "toric",
    },
    {
        "name": "BCEW",
        "substitution": {
            B0: t**2, C0: p * t**3, D0: 0,
            E0: q * t**5, W0: r * t**6,
        },
        "terms": 47,
        "support": (
            (0, 1, 0),  # E
            (0, 1, 1),  # E W
            (1, 0, 0),  # C
            (1, 0, 1),  # C W
            (2, 1, 0),  # C^2 E
            (3, 0, 0),  # C^3
        ),
        "kind": "projective_toric",
    },
    {
        "name": "BDEW",
        "substitution": {
            B0: t**2, C0: 0, D0: p * t**4,
            E0: q * t**5, W0: r * t**6,
        },
        "terms": 42,
        "support": ((0, 1, 0), (0, 1, 1), (1, 1, 0)),
        "kind": "linear",
    },
    {
        "name": "CDEW",
        "substitution": {
            B0: 0, C0: t**3, D0: p * t**4,
            E0: q * t**5, W0: r * t**6,
        },
        "terms": 28,
        "support": ((0, 1, 0), (0, 1, 1), (2, 0, 0)),
        "kind": "linear",
    },
)

certificates = []
for chart in CHARTS:
    primitive, content, leading, degree_t, coefficients = scaled_eliminant(
        chart["substitution"]
    )
    if chart["terms"] is not None:
        require(len(primitive.terms()) == chart["terms"],
                f"{chart['name']}: term count changed")
    require(degree_t == 10, f"{chart['name']}: t-degree changed")
    terminal = degree_t + 1
    result = {
        "name": chart["name"],
        "terms": len(primitive.terms()),
        "content": content,
        "leading": leading,
        "degree_t": degree_t,
        "terminal": terminal,
        "support": chart["support"],
        "fields": {},
    }
    for field_name, field, q0, s0 in (
        ("root", base.root_field, base.root_q0, base.root_s0),
        ("pair", base.pair_field, base.pair_q0, base.pair_s0),
    ):
        system = HenselSystem3(
            f"{chart['name']}_{field_name}", field, q0, s0
        )
        equations = terminal_equations(system, coefficients, terminal)
        matrix = coefficient_matrix(equations, chart["support"])
        if chart["kind"] == "linear":
            count, degree, terms = full_rank_certificate(system, matrix)
            field_result = {
                "rank": len(chart["support"]),
                "nonzero_minors": count,
                "minor_degree": degree,
                "minor_terms": terms,
            }
        elif chart["kind"] == "toric":
            toric = bcdw_toric_certificate(system, matrix)
            field_result = {"rank": 4, **toric}
        elif chart["kind"] == "projective_toric":
            field_result = bcew_projective_toric_certificate(system, matrix)
        else:
            field_result = bcde_macaulay_toric_certificate(system, matrix)
        result["fields"][field_name] = field_result
    certificates.append(result)

print("degree-22 support-four terminal Hensel/toric closure")
print(f"base_thm2636_sha256={BASE_SHA256}")
print(f"parent_thm2671_sha256={PARENT_SHA256}")
print(f"general_resultant_content={general_content} terms={len(general_primitive.terms())}")
print("fixed_quintic_degree=5 root_field_degree=5 pair_field_degree=10")
for certificate in certificates:
    support = ";".join(",".join(map(str, exponent))
                       for exponent in certificate["support"])
    print(
        f"{certificate['name']}:terms={certificate['terms']}:"
        f"degrees={certificate['degree_t']},5:"
        f"terminal={certificate['terminal']}:support={support}:"
        f"content={certificate['content']}:v5={certificate['leading']}"
    )
    for field_name in ("root", "pair"):
        result = certificate["fields"][field_name]
        if certificate["name"] == "BCDE":
            print(
                f"BCDE_{field_name}:rank={result['rank']}:"
                f"pivots={result['pivots']}:toric_quadrics={result['quadrics']}:"
                f"macaulay_shape={result['macaulay_shape']}:"
                f"macaulay_pivots={result['macaulay_pivots']}:"
                f"macaulay_minor_degree={result['degree']}:"
                f"macaulay_minor_terms={result['terms']}:"
                f"macaulay_minor_digest={result['digest']}:uniform=True"
            )
        elif certificate["name"] == "BCDW":
            print(
                f"BCDW_{field_name}:rank={result['rank']}:"
                f"kernel_rows={result['rows']}:all_kernel_coords_uniform=True:"
                f"toric_defect_degree={result['degree']}:"
                f"toric_defect_terms={result['terms']}:"
                f"toric_defect_digest={result['digest']}:uniform=True"
            )
        elif certificate["name"] == "BCEW":
            print(
                f"BCEW_{field_name}:rank={result['rank']}:"
                f"pivots={result['pivots']}:"
                f"two_toric_quadrics_coefficients_uniform=True:"
                f"sylvester_degree={result['degree']}:"
                f"sylvester_terms={result['terms']}:"
                f"sylvester_digest={result['digest']}:uniform=True"
            )
        else:
            print(
                f"{certificate['name']}_{field_name}:rank={result['rank']}:"
                f"nonzero_minors={result['nonzero_minors']}:"
                f"first_minor_degree={result['minor_degree']}:"
                f"first_minor_terms={result['minor_terms']}:uniform=True"
            )
print("closed_support_four_charts=BCDE,BCDW,BCEW,BDEW,CDEW")
print("all_five_support_four_strata=EMPTY")
print("remaining_degree22_coefficient_chart=full_support_BCDEW")
print("ALL CHECKS PASSED")
