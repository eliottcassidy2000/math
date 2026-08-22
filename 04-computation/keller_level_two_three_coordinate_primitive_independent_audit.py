#!/usr/bin/env python3
"""Independent exact audit of the three source coordinates of F^2.

This is a clean-room regular-representation calculation.  It imports no
repository computation module and uses neither a precomputed eliminant nor a
candidate theorem script.  The only external algebra engine is python-flint:
FLINT computes exact rational matrix ranks, determinants, characteristic and
minimal polynomials, and univariate discriminants.

The specialization is the target (a,b,c)=(1,1,1).  Its outer cubic is
squarefree but reducible over Q, so the lawful object is a rank-three etale
Q-algebra, not a cubic number field.  Adjoining the inner x-coordinate gives
the rank-nine fibre algebra.  Nonzero specialization determinants are enough
for the generic primitivity argument; connectedness of this one fibre is not
needed.

Reproduce with

    python3 -B 04-computation/keller_level_two_three_coordinate_primitive_independent_audit.py
    python3 -B -O 04-computation/keller_level_two_three_coordinate_primitive_independent_audit.py
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

from flint import fmpq, fmpq_mat, fmpq_poly


def fail(message: str) -> None:
    raise SystemExit(f"FAIL: {message}")


def require(condition: bool, message: str) -> None:
    # Deliberately not an assert: python -O must execute every audit gate.
    if not condition:
        fail(message)


Q0 = fmpq(0)
Q1 = fmpq(1)


def matrix_from_columns(columns: list[list[fmpq]]) -> fmpq_mat:
    require(bool(columns), "matrix needs at least one column")
    nrows = len(columns[0])
    require(all(len(column) == nrows for column in columns), "ragged columns")
    return fmpq_mat([[columns[j][i] for j in range(len(columns))] for i in range(nrows)])


@dataclass(frozen=True)
class CubicEtale:
    """Q[t]/(t^3+t/25-2/25), in the basis 1,t,t^2."""

    c: tuple[fmpq, fmpq, fmpq]

    @staticmethod
    def coerce(value: object) -> "CubicEtale":
        if isinstance(value, CubicEtale):
            return value
        return CubicEtale((fmpq(value), Q0, Q0))

    @staticmethod
    def t() -> "CubicEtale":
        return CubicEtale((Q0, Q1, Q0))

    def __add__(self, other: object) -> "CubicEtale":
        rhs = CubicEtale.coerce(other)
        return CubicEtale(tuple(self.c[i] + rhs.c[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self) -> "CubicEtale":
        return CubicEtale(tuple(-entry for entry in self.c))

    def __sub__(self, other: object) -> "CubicEtale":
        return self + (-CubicEtale.coerce(other))

    def __rsub__(self, other: object) -> "CubicEtale":
        return CubicEtale.coerce(other) - self

    def __mul__(self, other: object) -> "CubicEtale":
        rhs = CubicEtale.coerce(other)
        raw = [Q0 for _ in range(5)]
        for i, left in enumerate(self.c):
            for j, right in enumerate(rhs.c):
                raw[i + j] += left * right
        # t^3 = 2/25 - t/25.  Descending reduction also handles t^4.
        for degree in range(4, 2, -1):
            coefficient = raw[degree]
            raw[degree] = Q0
            raw[degree - 3] += fmpq(2, 25) * coefficient
            raw[degree - 2] -= fmpq(1, 25) * coefficient
        return CubicEtale(tuple(raw[:3]))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "CubicEtale":
        require(exponent >= 0, "negative powers use inverse() explicitly")
        value = CubicEtale.coerce(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                value *= base
            base *= base
            power >>= 1
        return value

    def multiplication_matrix(self) -> fmpq_mat:
        basis = [CubicEtale.coerce(1), CubicEtale.t(), CubicEtale.t() ** 2]
        return matrix_from_columns([list((self * basis_element).c) for basis_element in basis])

    def inverse(self) -> "CubicEtale":
        matrix = self.multiplication_matrix()
        require(matrix.det() != 0, f"attempted to invert nonunit {self}")
        solution = matrix.solve(fmpq_mat([[Q1], [Q0], [Q0]]))
        inverse = CubicEtale(tuple(solution[i, 0] for i in range(3)))
        require(self * inverse == CubicEtale.coerce(1), "cubic-etale inverse check")
        return inverse

    def __truediv__(self, other: object) -> "CubicEtale":
        return self * CubicEtale.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "CubicEtale":
        return CubicEtale.coerce(other) / self

    def is_zero(self) -> bool:
        return all(entry == 0 for entry in self.c)

    def __str__(self) -> str:
        return f"({self.c[0]})+({self.c[1]})*t+({self.c[2]})*t^2"


# These coefficients are assigned after reconstructing the middle fibre.
INNER_LINEAR: CubicEtale | None = None
INNER_CONSTANT: CubicEtale | None = None


@dataclass(frozen=True)
class RankNineEtale:
    """B[x]/(x^3+INNER_LINEAR*x+INNER_CONSTANT), B=CubicEtale."""

    c: tuple[CubicEtale, CubicEtale, CubicEtale]

    @staticmethod
    def coerce(value: object) -> "RankNineEtale":
        if isinstance(value, RankNineEtale):
            return value
        if isinstance(value, CubicEtale):
            return RankNineEtale((value, CubicEtale.coerce(0), CubicEtale.coerce(0)))
        return RankNineEtale((CubicEtale.coerce(value), CubicEtale.coerce(0), CubicEtale.coerce(0)))

    @staticmethod
    def x() -> "RankNineEtale":
        return RankNineEtale((CubicEtale.coerce(0), CubicEtale.coerce(1), CubicEtale.coerce(0)))

    def __add__(self, other: object) -> "RankNineEtale":
        rhs = RankNineEtale.coerce(other)
        return RankNineEtale(tuple(self.c[i] + rhs.c[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self) -> "RankNineEtale":
        return RankNineEtale(tuple(-entry for entry in self.c))

    def __sub__(self, other: object) -> "RankNineEtale":
        return self + (-RankNineEtale.coerce(other))

    def __rsub__(self, other: object) -> "RankNineEtale":
        return RankNineEtale.coerce(other) - self

    def __mul__(self, other: object) -> "RankNineEtale":
        require(INNER_LINEAR is not None and INNER_CONSTANT is not None, "inner relation initialized")
        rhs = RankNineEtale.coerce(other)
        raw = [CubicEtale.coerce(0) for _ in range(5)]
        for i, left in enumerate(self.c):
            for j, right in enumerate(rhs.c):
                raw[i + j] += left * right
        # x^3 = -INNER_LINEAR*x - INNER_CONSTANT.
        for degree in range(4, 2, -1):
            coefficient = raw[degree]
            raw[degree] = CubicEtale.coerce(0)
            raw[degree - 2] -= coefficient * INNER_LINEAR
            raw[degree - 3] -= coefficient * INNER_CONSTANT
        return RankNineEtale(tuple(raw[:3]))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "RankNineEtale":
        require(exponent >= 0, "negative powers use inverse() explicitly")
        value = RankNineEtale.coerce(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                value *= base
            base *= base
            power >>= 1
        return value

    def vector(self) -> list[fmpq]:
        # Basis: 1,t,t^2,x,tx,t^2x,x^2,tx^2,t^2x^2.
        return [entry for x_coefficient in self.c for entry in x_coefficient.c]

    def multiplication_matrix(self) -> fmpq_mat:
        return matrix_from_columns([(self * basis_element).vector() for basis_element in rank_nine_basis()])

    def inverse(self) -> "RankNineEtale":
        matrix = self.multiplication_matrix()
        require(matrix.det() != 0, "attempted to invert a rank-nine nonunit")
        rhs = fmpq_mat([[Q1]] + [[Q0] for _ in range(8)])
        solution = matrix.solve(rhs)
        inverse = rank_nine_from_vector([solution[i, 0] for i in range(9)])
        require(self * inverse == RankNineEtale.coerce(1), "rank-nine inverse check")
        return inverse

    def __truediv__(self, other: object) -> "RankNineEtale":
        return self * RankNineEtale.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "RankNineEtale":
        return RankNineEtale.coerce(other) / self

    def is_zero(self) -> bool:
        return all(entry.is_zero() for entry in self.c)


def rank_nine_from_vector(vector: list[fmpq]) -> RankNineEtale:
    require(len(vector) == 9, "rank-nine vector length")
    return RankNineEtale(
        tuple(CubicEtale(tuple(vector[3 * j + i] for i in range(3))) for j in range(3))
    )


def rank_nine_basis() -> list[RankNineEtale]:
    basis = []
    for x_degree in range(3):
        for t_degree in range(3):
            entries = [Q0 for _ in range(9)]
            entries[3 * x_degree + t_degree] = Q1
            basis.append(rank_nine_from_vector(entries))
    return basis


def L_polynomial(a: object, b: object, c: object):
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


def F_map(x: object, y: object, z: object):
    u = 1 + x * y
    return (
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


def inverse_section(a: object, b: object, c: object, x: object):
    """Reconstruct s=xy,y,z from the degree-one conic subresultant."""
    D = -3 * x**2 * a * c + x**2 * b**2 * c - x**2 * b + 2 * x * b * c - 2 * x + c
    N = (
        -3 * x**3 * a * b * c
        + 4 * x**3 * a
        - 6 * x**2 * a * c
        + x**2 * b**2 * c
        - x**2 * b
        + 2 * x * b * c
        - 2 * x
        + c
    )
    s = -N / D
    y = s / x
    z = (x * (2 - 3 * s) - c) / x**3
    return D, N, s, y, z


def check_inverse_equations(a: object, b: object, c: object, x: object, D: object, N: object, s: object) -> None:
    core = L_polynomial(a, b, c) * x**3 + (4 - 3 * b * c) * x - 2 * c
    conic_1 = 3 * a * x**2 - (1 + s) * (b * x - s)
    conic_2 = a * x**3 + c * (1 + s) ** 3 - x * (1 + s) * (2 + s)
    require(core.is_zero(), "inverse x-core")
    require((D * s + N).is_zero(), "degree-one subresultant equation")
    require(conic_1.is_zero(), "first inverse conic")
    require(conic_2.is_zero(), "second inverse conic")


def require_triple_equal(actual: tuple[object, object, object], expected: tuple[object, object, object], label: str) -> None:
    for index, (left, right) in enumerate(zip(actual, expected, strict=True), start=1):
        require((left - right).is_zero(), f"{label}, coordinate {index}")


def power_basis_matrix(element: RankNineEtale) -> fmpq_mat:
    powers = [RankNineEtale.coerce(1)]
    for _ in range(1, 9):
        powers.append(powers[-1] * element)
    return matrix_from_columns([power.vector() for power in powers])


def polynomial_digest(polynomial: fmpq_poly) -> str:
    payload = ";".join(str(coefficient) for coefficient in polynomial.coeffs()).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def rational_square_root(value: fmpq, label: str) -> fmpq:
    require(value > 0, f"{label} must be positive to be a rational square")
    try:
        root = value.sqrt()
    except Exception as exc:  # FLINT raises DomainError for a nonsquare.
        fail(f"{label} is not a rational square: {exc}")
    require(root * root == value, f"{label} square-root replay")
    return root


def F_mod(point: tuple[int, int, int], prime: int) -> tuple[int, int, int]:
    """The original map evaluated directly in F_prime."""
    x, y, z = point
    u = (1 + x * y) % prime
    return (
        (u**3 * z + y**2 * u * (4 + 3 * x * y)) % prime,
        (y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) % prime,
        (2 * x - 3 * x**2 * y - x**3 * z) % prime,
    )


def direct_finite_field_control() -> tuple[int, tuple[int, int, int], list[tuple[tuple[int, int, int], tuple[int, int, int]]]]:
    """Find a split nine-point fibre directly from F^2 equations.

    This control is intentionally disjoint from inverse_section and from the
    rational tower classes.  In the split fibre algebra F_p^9, the power-basis
    matrix of a coordinate is a Vandermonde matrix.  It has rank nine exactly
    when the nine coordinate values are distinct.  The first coordinate of
    the middle point is the degree-three hostile and should take three values,
    each on one block of three source points.
    """
    for prime in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        fibres: dict[
            tuple[int, int, int],
            list[tuple[tuple[int, int, int], tuple[int, int, int]]],
        ] = {}
        for x in range(prime):
            for y in range(prime):
                for z in range(prime):
                    source = (x, y, z)
                    middle = F_mod(source, prime)
                    target = F_mod(middle, prime)
                    fibres.setdefault(target, []).append((source, middle))
        for target in sorted(fibres):
            fibre = fibres[target]
            if len(fibre) != 9:
                continue
            source_value_counts = [len({entry[0][axis] for entry in fibre}) for axis in range(3)]
            middle_points = {entry[1] for entry in fibre}
            middle_x_values = [entry[1][0] for entry in fibre]
            middle_x_set = set(middle_x_values)
            middle_x_multiplicities = sorted(middle_x_values.count(value) for value in middle_x_set)
            if source_value_counts == [9, 9, 9] and len(middle_points) == 3:
                if len(middle_x_set) == 3 and middle_x_multiplicities == [3, 3, 3]:
                    return prime, target, fibre
    fail("no direct finite-field split control found in the registered prime bank")


def vandermonde_determinant_mod(values: list[int], prime: int) -> int:
    determinant = 1
    for j in range(len(values)):
        for i in range(j):
            determinant = determinant * (values[j] - values[i]) % prime
    return determinant


def main() -> None:
    global INNER_LINEAR, INNER_CONSTANT

    print("engine=python-flint regular representations over exact Q")
    print("candidate_script_imported=False")
    print("target=(1,1,1)")

    # Outer inverse stage.
    t = CubicEtale.t()
    outer_core = 25 * t**3 + t - 2
    require(outer_core.is_zero(), "outer cubic relation")
    outer_poly = fmpq_poly([-2, 1, 0, 25])
    outer_factorization = outer_poly.factor()
    require(outer_poly.gcd(outer_poly.derivative()).degree() == 0, "outer cubic squarefree")
    print(f"outer_core={outer_poly}")
    print(f"outer_factorization={outer_factorization}")
    print(f"outer_discriminant={outer_poly.discriminant()}")

    outer_D, outer_N, outer_s, middle_y, middle_z = inverse_section(
        CubicEtale.coerce(1), CubicEtale.coerce(1), CubicEtale.coerce(1), t
    )
    check_inverse_equations(
        CubicEtale.coerce(1), CubicEtale.coerce(1), CubicEtale.coerce(1),
        t, outer_D, outer_N, outer_s,
    )
    middle = (t, middle_y, middle_z)
    require_triple_equal(
        F_map(*middle),
        (CubicEtale.coerce(1), CubicEtale.coerce(1), CubicEtale.coerce(1)),
        "outer inverse replay",
    )
    print(f"outer_D_norm={outer_D.multiplication_matrix().det()}")
    print(f"outer_x_norm={t.multiplication_matrix().det()}")
    print("outer_inverse_equations=PASS")

    # Inner inverse stage over the rank-three etale algebra.
    inner_L = L_polynomial(*middle)
    inner_p = 4 - 3 * middle_y * middle_z
    inner_L_norm = inner_L.multiplication_matrix().det()
    require(inner_L_norm != 0, "inner leading coefficient is a unit")
    INNER_LINEAR = inner_p / inner_L
    INNER_CONSTANT = -2 * middle_z / inner_L
    print(f"Norm_outer(inner_L)={inner_L_norm}")

    x = RankNineEtale.x()
    inner_D, inner_N, inner_s, y, z = inverse_section(
        RankNineEtale.coerce(t), RankNineEtale.coerce(middle_y), RankNineEtale.coerce(middle_z), x
    )
    check_inverse_equations(
        RankNineEtale.coerce(t), RankNineEtale.coerce(middle_y), RankNineEtale.coerce(middle_z),
        x, inner_D, inner_N, inner_s,
    )
    source = (x, y, z)
    embedded_middle = tuple(RankNineEtale.coerce(entry) for entry in middle)
    require_triple_equal(F_map(*source), embedded_middle, "inner inverse replay")
    require_triple_equal(
        F_map(*embedded_middle),
        (RankNineEtale.coerce(1), RankNineEtale.coerce(1), RankNineEtale.coerce(1)),
        "outer map replay in rank-nine algebra",
    )
    require_triple_equal(
        F_map(*F_map(*source)),
        (RankNineEtale.coerce(1), RankNineEtale.coerce(1), RankNineEtale.coerce(1)),
        "F^2 replay",
    )
    print(f"inner_D_norm={inner_D.multiplication_matrix().det()}")
    print(f"inner_x_norm={x.multiplication_matrix().det()}")
    print("inner_inverse_equations=PASS")
    print("F2_equations=PASS")

    # Exact power-basis gates and the intermediate-coordinate hostile.
    coordinate_data: dict[str, dict[str, object]] = {}
    for name, coordinate in (("x", x), ("y", y), ("z", z)):
        power_matrix = power_basis_matrix(coordinate)
        multiplication_matrix = coordinate.multiplication_matrix()
        characteristic = multiplication_matrix.charpoly()
        minimal = multiplication_matrix.minpoly()
        discriminant = characteristic.discriminant()
        coordinate_data[name] = {
            "power_matrix": power_matrix,
            "volume": power_matrix.det(),
            "characteristic": characteristic,
            "minimal": minimal,
            "discriminant": discriminant,
        }
        require(power_matrix.rank() == 9, f"{name} power-basis rank")
        require(power_matrix.det() != 0, f"{name} power-basis volume")
        require(characteristic.degree() == 9, f"{name} characteristic degree")
        require(minimal == characteristic, f"{name} is primitive in specialized algebra")
        require(discriminant != 0, f"{name} eliminant is squarefree")
        print(f"{name}_power_rank={power_matrix.rank()}")
        print(f"{name}_power_volume={power_matrix.det()}")
        print(f"{name}_minpoly_degree={minimal.degree()}")
        print(f"{name}_minpoly_sha256={polynomial_digest(minimal)}")
        print(f"{name}_discriminant={discriminant}")

    intermediate = RankNineEtale.coerce(t)
    hostile_power = power_basis_matrix(intermediate)
    hostile_multiplication = intermediate.multiplication_matrix()
    hostile_charpoly = hostile_multiplication.charpoly()
    hostile_minpoly = hostile_multiplication.minpoly()
    require(hostile_power.rank() == 3, "intermediate hostile has rank exactly three")
    require(hostile_power.det() == 0, "intermediate hostile determinant vanishes")
    require(hostile_minpoly.degree() == 3, "intermediate hostile minimal degree")
    require(hostile_charpoly == hostile_minpoly**3, "intermediate characteristic polynomial is cubic^3")
    print(f"intermediate_power_rank={hostile_power.rank()}")
    print(f"intermediate_power_volume={hostile_power.det()}")
    print(f"intermediate_minpoly={hostile_minpoly}")
    print(f"intermediate_charpoly_equals_minpoly_cubed={hostile_charpoly == hostile_minpoly**3}")

    # Basis-volume congruence and discriminant square classes.
    disc_x = coordinate_data["x"]["discriminant"]
    volume_x = coordinate_data["x"]["volume"]
    require(isinstance(disc_x, fmpq) and isinstance(volume_x, fmpq), "typed x invariants")
    H_value = fmpq(951326441195)
    print(f"H_at_target={H_value}")
    for name in ("x", "y", "z"):
        disc = coordinate_data[name]["discriminant"]
        volume = coordinate_data[name]["volume"]
        require(isinstance(disc, fmpq) and isinstance(volume, fmpq), f"typed {name} invariants")
        ratio_to_x = disc / disc_x
        basis_ratio = volume / volume_x
        require(ratio_to_x == basis_ratio**2, f"{name}/x discriminant-basis ratio")
        ratio_to_H = disc / H_value
        root_to_H = rational_square_root(ratio_to_H, f"Disc_{name}/H")
        print(f"Disc_{name}_over_Disc_x={ratio_to_x}")
        print(f"basis_volume_{name}_over_x={basis_ratio}")
        print(f"Disc_{name}_over_H={ratio_to_H}")
        print(f"sqrt_Disc_{name}_over_H={root_to_H}")

    # The old norm identity is recomputed from multiplication by L(q), not
    # imported from a precomputed H polynomial or resultant ledger.
    expected_L_norm = H_value / (64 * 25)
    require(inner_L_norm == expected_L_norm, "independent specialized norm-to-H replay")
    print(f"Norm_outer(inner_L)_over_H={inner_L_norm / H_value}")
    print("specialized_norm_to_H=PASS")

    prime, finite_target, finite_fibre = direct_finite_field_control()
    finite_source_counts = tuple(
        len({source[axis] for source, _middle in finite_fibre}) for axis in range(3)
    )
    finite_middle_x_values = [middle[0] for _source, middle in finite_fibre]
    finite_middle_counts = tuple(
        sorted(finite_middle_x_values.count(value) for value in set(finite_middle_x_values))
    )
    finite_middle_points = sorted({middle for _source, middle in finite_fibre})
    finite_L_values = (
        int(L_polynomial(*finite_target) % prime),
        tuple((middle, int(L_polynomial(*middle) % prime)) for middle in finite_middle_points),
    )
    sorted_finite_fibre = sorted(finite_fibre)
    finite_vandermonde = tuple(
        vandermonde_determinant_mod([source[axis] for source, _middle in sorted_finite_fibre], prime)
        for axis in range(3)
    )
    require(finite_source_counts == (9, 9, 9), "direct finite-field coordinate ranks")
    require(finite_middle_counts == (3, 3, 3), "direct finite-field hostile blocks")
    require(finite_L_values[0] != 0, "direct finite-field outer target is off L")
    require(all(value != 0 for _middle, value in finite_L_values[1]),
            "direct finite-field middle targets are off L")
    require(all(determinant != 0 for determinant in finite_vandermonde),
            "direct finite-field Vandermonde determinants")
    require(all(F_mod(F_mod(source, prime), prime) == finite_target for source, _middle in finite_fibre),
            "direct finite-field F^2 replay")
    print("finite_field_prime_bank=(11,13,17,19,23,29,31,37,41,43)")
    print(f"finite_field_prime={prime}")
    print(f"finite_field_target={finite_target}")
    print(f"finite_field_fibre_size={len(finite_fibre)}")
    print(f"finite_field_fibre={sorted_finite_fibre}")
    print(f"finite_field_L_values={finite_L_values}")
    print(f"finite_field_source_coordinate_distinct_counts={finite_source_counts}")
    print(f"finite_field_source_vandermonde_determinants={finite_vandermonde}")
    print(f"finite_field_intermediate_x_multiplicities={finite_middle_counts}")
    print("direct_finite_field_control=PASS")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
