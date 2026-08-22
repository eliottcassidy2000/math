#!/usr/bin/env python3
"""Closed mixed moments and a projective support-four Macaulay certificate.

This is an exact companion to THM-3310's still-open support-four/five frontier.
It has three logically separate layers.

1.  A closed coefficient formula for ``mu(a,b)=<z^a zbar^b>`` is derived from

        product_{lambda^3=1} (1-X lambda-Y lambda^{-1})
            = 1-X^3-Y^3-3XY.

    An independent exact barycentric expansion checks every pair a+b <= 10.

2.  The tempting continuous torus ``z -> t z, zbar -> t^-1 zbar`` is audited.
    It is a coordinate automorphism of ``uv=r^3`` but not a symmetry of the
    simplex moment functional.  Thus only projective scaling can normalize a
    coefficient; no second coefficient is set to one.

3.  For each of the five coordinate hyperplanes in the degree-four cyclic
    eigenspace, the guarded reductions modulo p=101 of M_3,...,M_21 generate
    the entire degree-21 piece in four coefficient variables.  Hence the
    common affine cone is only the origin, including every projective boundary
    and point at infinity.  Full rank modulo p exhibits a nonzero rational
    maximal minor, so the projective exclusion lifts to characteristic zero.

The modular guard is load-bearing: m <= 21 and deg(g) <= 4 require
``p > 4*21+2 = 86``.  No floating point, randomness, assertion-sensitive
test, or affine-only resultant is used.
"""

from functools import lru_cache
from fractions import Fraction
from math import comb, factorial

import flint
from flint import nmod_mat


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compositions(total, parts, prefix=()):
    """Yield weak compositions in a fixed lexicographic order."""
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(total + 1):
        yield from compositions(total - first, parts - 1, prefix + (first,))


def multinomial(total, exponents):
    value = factorial(total)
    for exponent in exponents:
        value //= factorial(exponent)
    return value


# ------------------------------------------------ exact Q(omega) control layer


class Cyclo:
    """Exact p+q*omega arithmetic with omega^2=-1-omega."""

    __slots__ = ("p", "q")

    def __init__(self, p=0, q=0):
        self.p = Fraction(p)
        self.q = Fraction(q)

    def __add__(self, other):
        other = as_cyclo(other)
        return Cyclo(self.p + other.p, self.q + other.q)

    __radd__ = __add__

    def __neg__(self):
        return Cyclo(-self.p, -self.q)

    def __sub__(self, other):
        return self + (-as_cyclo(other))

    def __mul__(self, other):
        other = as_cyclo(other)
        return Cyclo(
            self.p * other.p - self.q * other.q,
            self.p * other.q + self.q * other.p - self.q * other.q,
        )

    __rmul__ = __mul__

    def __pow__(self, exponent):
        require(exponent >= 0, "nonnegative cyclotomic exponent")
        answer = Cyclo(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power >>= 1
        return answer

    def __eq__(self, other):
        other = as_cyclo(other)
        return self.p == other.p and self.q == other.q

    def is_zero(self):
        return self.p == 0 and self.q == 0


def as_cyclo(value):
    return value if isinstance(value, Cyclo) else Cyclo(value)


ONE = Cyclo(1)
OMEGA = Cyclo(0, 1)
OMEGA2 = OMEGA * OMEGA
ROOTS = (ONE, OMEGA, OMEGA2)
CONJ_ROOTS = (ONE, OMEGA2, OMEGA)
require((OMEGA * OMEGA + OMEGA + ONE).is_zero(), "omega relation")


def poly2_mul(left, right):
    """Multiply bivariate polynomials with Cyclo coefficients."""
    out = {}
    for (a, b), lc in left.items():
        for (c, d), rc in right.items():
            key = (a + c, b + d)
            out[key] = out.get(key, Cyclo()) + lc * rc
    return {key: value for key, value in out.items() if not value.is_zero()}


def denominator_cubic_control():
    product = {(0, 0): ONE}
    for root, conjugate in zip(ROOTS, CONJ_ROOTS):
        product = poly2_mul(
            product,
            {(0, 0): ONE, (1, 0): -root, (0, 1): -conjugate},
        )
    expected = {
        (0, 0): Cyclo(1),
        (3, 0): Cyclo(-1),
        (0, 3): Cyclo(-1),
        (1, 1): Cyclo(-3),
    }
    require(product == expected, "Hesse denominator cubic identity")


@lru_cache(maxsize=None)
def coefficient_c(a, b):
    """[X^a Y^b] 1/(1-X^3-Y^3-3XY), as a nonnegative integer."""
    require(a >= 0 and b >= 0, "nonnegative mixed exponents")
    answer = 0
    for k in range(min(a, b) + 1):
        if (a - k) % 3 != 0 or (b - k) % 3 != 0:
            continue
        i = (a - k) // 3
        j = (b - k) // 3
        answer += multinomial(i + j + k, (i, j, k)) * 3**k
    return answer


@lru_cache(maxsize=None)
def mu_exact(a, b):
    """Closed exact mixed simplex moment <z^a zbar^b>."""
    return Fraction(
        2 * factorial(a) * factorial(b) * coefficient_c(a, b),
        factorial(a + b + 2),
    )


def direct_barycentric_mu(a, b):
    """Independent six-index Dirichlet expansion, used only as a control."""
    answer = Cyclo()
    for left in compositions(a, 3):
        left_coefficient = Cyclo(multinomial(a, left))
        for index in range(3):
            left_coefficient = left_coefficient * ROOTS[index] ** left[index]
        for right in compositions(b, 3):
            coefficient = left_coefficient * multinomial(b, right)
            for index in range(3):
                coefficient = coefficient * CONJ_ROOTS[index] ** right[index]
            barycentric_moment = Fraction(
                2
                * factorial(left[0] + right[0])
                * factorial(left[1] + right[1])
                * factorial(left[2] + right[2]),
                factorial(a + b + 2),
            )
            answer = answer + coefficient * barycentric_moment
    return answer


def mixed_moment_controls():
    denominator_cubic_control()
    checked = 0
    for total in range(11):
        for a in range(total + 1):
            b = total - a
            direct = direct_barycentric_mu(a, b)
            require(direct.q == 0, f"mu({a},{b}) is rational")
            require(direct.p == mu_exact(a, b), f"closed mu({a},{b})")
            require(
                (mu_exact(a, b) == 0) == ((a - b) % 3 != 0),
                f"selection rule at ({a},{b})",
            )
            checked += 1
    named = {
        (0, 0): Fraction(1),
        (1, 1): Fraction(1, 4),
        (3, 0): Fraction(1, 10),
        (0, 3): Fraction(1, 10),
        (2, 2): Fraction(1, 10),
        (1, 4): Fraction(2, 35),
        (4, 1): Fraction(2, 35),
        (6, 0): Fraction(1, 28),
        (3, 3): Fraction(29, 560),
    }
    for key, value in named.items():
        require(mu_exact(*key) == value, f"named moment {key}")

    # Coefficient extraction from the rational denominator gives a constant-
    # work recurrence per lattice point.  This is also a hostile boundary test:
    # every negative index is interpreted as zero, not wrapped by Python.
    recurrence_checked = 0
    for a in range(85):
        for b in range(85):
            previous = 0
            if a >= 3:
                previous += coefficient_c(a - 3, b)
            if b >= 3:
                previous += coefficient_c(a, b - 3)
            if a >= 1 and b >= 1:
                previous += 3 * coefficient_c(a - 1, b - 1)
            if a == 0 and b == 0:
                previous += 1
            require(coefficient_c(a, b) == previous,
                    f"mixed-coefficient recurrence ({a},{b})")
            recurrence_checked += 1
    return checked, recurrence_checked


# ----------------------------------------------------- coefficient moment ring


LABELS = ("A", "B", "C", "D", "E")
BASIS = ((0, 1), (2, 0), (1, 2), (3, 1), (0, 4))
TORUS_WEIGHTS = tuple(a - b for a, b in BASIS)
MOMENT_ORDERS = tuple(range(3, 22, 3))
PRIME = 101
MAX_ORDER = max(MOMENT_ORDERS)
MAX_SIMPLEX_DEGREE = 4 * MAX_ORDER
require(PRIME > MAX_SIMPLEX_DEGREE + 2, "MISTAKE-363 modular guard")


@lru_cache(maxsize=None)
def mu_mod(a, b, prime=PRIME):
    require(prime > a + b + 2, f"moment denominator guard ({a},{b})")
    numerator = 2 * factorial(a) * factorial(b) * coefficient_c(a, b)
    denominator = factorial(a + b + 2)
    denominator_mod = denominator % prime
    require(denominator_mod != 0, f"invertible moment denominator ({a},{b})")
    return numerator % prime * pow(denominator_mod, -1, prime) % prime


def moment_polynomial(order, retained, prime=PRIME):
    """Reduction of M_order in the four retained coefficient variables."""
    require(len(retained) == 4, "support-four coordinate hyperplane")
    polynomial = {}
    for counts in compositions(order, 4):
        z_power = sum(counts[j] * BASIS[index][0] for j, index in enumerate(retained))
        zbar_power = sum(
            counts[j] * BASIS[index][1] for j, index in enumerate(retained)
        )
        require(z_power + zbar_power <= 4 * order, "degree-four moment bound")
        coefficient = (
            multinomial(order, counts) * mu_mod(z_power, zbar_power, prime)
        ) % prime
        if coefficient:
            polynomial[counts] = coefficient
    require(polynomial, f"nonzero M_{order} on retained chart")
    return polynomial


def macaulay_matrix(generators, target_degree, variable_count, prime=PRIME):
    """Build the homogeneous Macaulay matrix in one target degree."""
    columns = list(compositions(target_degree, variable_count))
    column_index = {monomial: index for index, monomial in enumerate(columns)}
    row_count = 0
    for degree, polynomial in generators:
        require(degree <= target_degree, "generator degree in Macaulay target")
        for exponent in polynomial:
            require(sum(exponent) == degree, "homogeneous generator")
            require(len(exponent) == variable_count, "generator arity")
        row_count += comb(target_degree - degree + variable_count - 1,
                          variable_count - 1)

    matrix = nmod_mat(row_count, len(columns), prime)
    row = 0
    for degree, polynomial in generators:
        for shift in compositions(target_degree - degree, variable_count):
            for exponent, coefficient in polynomial.items():
                monomial = tuple(
                    exponent[index] + shift[index]
                    for index in range(variable_count)
                )
                matrix[row, column_index[monomial]] = coefficient
            row += 1
    require(row == row_count, "Macaulay row count")
    return row_count, len(columns), matrix


def macaulay_rank(generators, target_degree, variable_count, prime=PRIME):
    """Exact rank of the homogeneous Macaulay map in one target degree."""
    rows, columns, matrix = macaulay_matrix(
        generators, target_degree, variable_count, prime
    )
    return rows, columns, matrix.rank()


def compress_ranges(indices):
    require(indices, "nonempty index set")
    pieces = []
    first = indices[0]
    last = first
    for index in indices[1:]:
        if index == last + 1:
            last = index
            continue
        pieces.append(str(first) if first == last else f"{first}-{last}")
        first = index
        last = index
    pieces.append(str(first) if first == last else f"{first}-{last}")
    return ",".join(pieces)


def maximal_minor_certificate(generators, target_degree, variable_count,
                              prime=PRIME):
    """Select and evaluate one explicit maximal minor modulo ``prime``.

    Rows and columns of the Macaulay matrix use the displayed generator order
    and ``compositions`` order.  Pivot rows are the leading columns of the
    unique RREF of the transpose, so the compressed row list identifies the
    minor reproducibly rather than merely asserting that one exists.
    """
    rows, columns, matrix = macaulay_matrix(
        generators, target_degree, variable_count, prime
    )
    reduced, rank = matrix.transpose().rref()
    require(rank == columns, "full-column-rank maximal-minor request")
    pivot_rows = []
    scan_column = 0
    for row in range(rank):
        while scan_column < reduced.ncols() and reduced[row, scan_column] == 0:
            scan_column += 1
        require(scan_column < reduced.ncols(), "RREF pivot exists")
        require(reduced[row, scan_column] == 1, "normalized RREF pivot")
        pivot_rows.append(scan_column)
        scan_column += 1
    require(len(pivot_rows) == columns, "maximal-minor row count")

    minor = nmod_mat(columns, columns, prime)
    for minor_row, source_row in enumerate(pivot_rows):
        for column in range(columns):
            minor[minor_row, column] = matrix[source_row, column]
    determinant = int(minor.det())
    require(determinant != 0, "displayed maximal minor is nonzero")
    return rows, columns, rank, compress_ranges(pivot_rows), determinant


def macaulay_engine_controls():
    # (x0,x1,x2,x3) has no projective point and spans every cubic.
    empty_generators = []
    for index in range(4):
        exponent = tuple(1 if position == index else 0 for position in range(4))
        empty_generators.append((1, {exponent: 1}))
    empty_shape = macaulay_rank(empty_generators, 3, 4)
    require(empty_shape == (40, 20, 20), "known-empty projective control")

    # (x0,x1,x2) leaves [0:0:0:1], hence exactly x3^3 is missing.
    hostile_shape = macaulay_rank(empty_generators[:3], 3, 4)
    require(hostile_shape == (30, 20, 19), "known-nonempty hostile control")
    return empty_shape, hostile_shape


def support_four_certificates():
    records = []
    for deleted in range(5):
        retained = tuple(index for index in range(5) if index != deleted)
        generators = [
            (order, moment_polynomial(order, retained))
            for order in MOMENT_ORDERS
        ]
        chart_records = []
        minor_record = None
        for target_degree in range(18, 23):
            active = [item for item in generators if item[0] <= target_degree]
            if target_degree == 21:
                rows, columns, rank, pivot_ranges, determinant = (
                    maximal_minor_certificate(active, target_degree, 4)
                )
                minor_record = (pivot_ranges, determinant)
            else:
                rows, columns, rank = macaulay_rank(active, target_degree, 4)
            chart_records.append((target_degree, rows, columns, rank))
        require(chart_records[-2] == (21, 2926, 2024, 2024),
                f"degree-21 full rank after deleting {LABELS[deleted]}")
        require(chart_records[-1] == (22, 3514, 2300, 2300),
                f"degree-22 constancy after deleting {LABELS[deleted]}")
        require(minor_record is not None, "degree-21 maximal minor recorded")
        records.append((LABELS[deleted], chart_records, minor_record))
    return records


def fraction_text(value):
    return str(value.numerator) if value.denominator == 1 else str(value)


def main():
    checked, recurrence_checked = mixed_moment_controls()
    empty_shape, hostile_shape = macaulay_engine_controls()

    require(TORUS_WEIGHTS == (-1, 2, -1, 2, -4), "torus weights")
    cube_weights = tuple(3 * weight for weight in TORUS_WEIGHTS)
    cube_coefficients = tuple(mu_exact(3 * a, 3 * b) for a, b in BASIS)
    require(all(value != 0 for value in cube_coefficients), "nonzero pure M3 terms")
    require(len(set(cube_weights)) == 3, "continuous torus is non-covariant")

    records = support_four_certificates()

    print("DEGREE-FOUR CYCLIC EIGENSPACE: CLOSED MOMENTS + SUPPORT-4")
    print(f"python_flint={flint.__version__}")
    print("status_closed_moment=PROVED algebraic identity; exact controls")
    print("status_support4=FINITE-EXACT projective Macaulay certificate")
    print("status_support5=OPEN")
    print()
    print("CLOSED MIXED MOMENT")
    print("denominator_product=1-X^3-Y^3-3XY")
    print("C(a,b)=sum multinomial(i+j+k;i,j,k)*3^k")
    print("constraints=3i+k=a,3j+k=b")
    print("mu(a,b)=2*a!*b!*C(a,b)/(a+b+2)!")
    print(f"independent_barycentric_pairs_checked={checked} (a+b<=10)")
    print("C(a,b)=C(a-3,b)+C(a,b-3)+3*C(a-1,b-1)+[a=b=0]")
    print(f"recurrence_lattice_points_checked={recurrence_checked} (0<=a,b<=84)")
    print("selection_rule=mu(a,b)=0 iff a-b is nonzero mod 3")
    print()
    print("INVARIANT / NORMALIZATION AUDIT")
    print("r=z*zbar, u=z^3, v=zbar^3, relation=u*v-r^3=0")
    print("z*g=A*r+B*u+C*r^2+D*r*u+E*r*v")
    print(f"formal_torus_weights={TORUS_WEIGHTS}")
    print(f"pure_M3_torus_weights={cube_weights}")
    print("pure_M3_coefficients=(" + ",".join(
        fraction_text(value) for value in cube_coefficients
    ) + ")")
    print("continuous_torus_covariance=FAILS (three distinct nonzero weights)")
    print("t^3=1_action=common coefficient scalar; projectively trivial")
    print("lawful_normalization=one projective coefficient only")
    print()
    print("MODULAR / PROJECTIVE GUARDS")
    print(f"prime={PRIME}")
    print(f"max_order={MAX_ORDER}")
    print(f"max_simplex_degree_plus_2={MAX_SIMPLEX_DEGREE + 2}")
    print(f"guard={PRIME}>{MAX_SIMPLEX_DEGREE + 2}: PASS")
    print(f"empty_control_rows_cols_rank={empty_shape}")
    print(f"hostile_control_rows_cols_rank={hostile_shape}")
    print()
    print("SUPPORT-4 MACAULAY RANKS")
    print("columns=all monomials of displayed target degree in four variables")
    print("delete degree rows columns rank deficiency")
    for deleted, chart_records, _minor_record in records:
        for degree, rows, columns, rank in chart_records:
            print(f"{deleted:>6} {degree:>6} {rows:>4} {columns:>7} "
                  f"{rank:>4} {columns-rank:>10}")
    print()
    print("DISPLAYED DEGREE-21 MAXIMAL MINORS MOD 101")
    print("ordering=generators M3..M21; shifts/columns use compositions order")
    print("rows=leading columns of RREF(Macaulay_transpose), zero-based")
    print("delete determinant pivot_rows")
    for deleted, _chart_records, (pivot_ranges, determinant) in records:
        print(f"{deleted:>6} {determinant:>11} {pivot_ranges}")
    print()
    print("degree21=FULL for every deletion; every degree-21 monomial is in ideal")
    print("projective_infinity=EXCLUDED, including each pure 21st power")
    print("degree22=FULL for every deletion (fresh-build constancy guard)")
    print("rank_lift=each displayed nonzero maximal minor mod 101 is nonzero over Q")
    print("conclusion=no nonzero common zero on any coefficient hyperplane")
    print("combined_with_THM3310=support<=4 excluded; support5 remains open")


if __name__ == "__main__":
    main()
