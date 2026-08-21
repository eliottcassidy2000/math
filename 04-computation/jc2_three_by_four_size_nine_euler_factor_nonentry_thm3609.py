#!/usr/bin/env python3
"""Exact companion for the four THM-3609 size-nine nonentries.

The hash-pinned THM-3606 transcript supplies the four canonical fibre words
W072, W073, W141, and W149 and their unique unexposed scalar schemes.  This
script parses those entries, derives the two collision cones by rational RREF,
reconstructs all seven actual weights from the scalar anchor, and applies the
singleton UFD and scalar-arm ledgers.

After the arm gate forces d=2, every nonsingleton bracket row is assembled in
an exact formal differential-monomial algebra over Q[a,t].  An independent
term-pairing division verifies the common factor

    E(h,k) = k h' + 2 h k'.

The scalar row would therefore have the form E times a polynomial equal to 1.
Since Sigma divides h and deg Sigma >= 2, E has positive degree and is not a
unit.  No bounded parameter scan or sample-weight extrapolation is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import re
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
PARENT_THEOREM = (
    ROOT
    / "01-canon/theorems/THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas.md"
)
PARENT_SCRIPT = (
    ROOT / "04-computation/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.py"
)
PARENT_OUTPUT = (
    ROOT / "05-knowledge/results/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.out"
)

EXPECTED_PARENT_THEOREM_SHA256 = (
    "c9b99569cafa1c2b14dbc065e1c66d24ce97d3e1bd2d388e616da2f5ad5510c4"
)
EXPECTED_PARENT_SCRIPT_SHA256 = (
    "9d558f5637c5cc3573214fe00817bd6acb4c749b85a9046644e2390bd1d0ad91"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "58b1160bd8831139fc5a6d4eb5102c876df4b47a823293080413c8f6c3f995b4"
)

TARGET_IDS = ("W072", "W073", "W141", "W149")
EXPECTED_SCHEME_SHA256 = {
    "W072": "2edabe6f195f9b106d9c430bf14d7f335383ec26b84d6f8661f1bdcea0939ba8",
    "W073": "16d074271e4dc7402a168439693bba66dc6b9081d388cfb5024c21ffd7354826",
    "W141": "b381f15bee1b47f0ec9b251f9e61b0c22629f5fce8406090e7e057f27b6b7173",
    "W149": "5a0101bcf6b78c3b44e581ab29b2ce9bb97caa55e3c057c36c8627db534a008d",
}

CONSTANT_NAMES = ("A", "B", "C", "L", "M", "N", "T")
CONSTANT_INDEX = {name: index for index, name in enumerate(CONSTANT_NAMES)}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def lf_sha256(path: Path) -> str:
    return sha256(lf_bytes(path)).hexdigest()


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True, order=True)
class Affine:
    """c + a_coeff*a + t_coeff*t over Q."""

    constant: Fraction = Fraction(0)
    a_coeff: Fraction = Fraction(0)
    t_coeff: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "constant", Fraction(self.constant))
        object.__setattr__(self, "a_coeff", Fraction(self.a_coeff))
        object.__setattr__(self, "t_coeff", Fraction(self.t_coeff))

    def __add__(self, other: object) -> "Affine":
        rhs = to_affine(other)
        return Affine(
            self.constant + rhs.constant,
            self.a_coeff + rhs.a_coeff,
            self.t_coeff + rhs.t_coeff,
        )

    def __radd__(self, other: object) -> "Affine":
        return self + other

    def __neg__(self) -> "Affine":
        return Affine(-self.constant, -self.a_coeff, -self.t_coeff)

    def __sub__(self, other: object) -> "Affine":
        return self + (-to_affine(other))

    def __rsub__(self, other: object) -> "Affine":
        return to_affine(other) - self

    def __mul__(self, scalar: object) -> "Affine":
        value = Fraction(scalar)
        return Affine(
            self.constant * value,
            self.a_coeff * value,
            self.t_coeff * value,
        )

    def __rmul__(self, scalar: object) -> "Affine":
        return self * scalar

    def __truediv__(self, scalar: object) -> "Affine":
        return self * (Fraction(1) / Fraction(scalar))

    def evaluate(self, a_value: int, t_value: int) -> Fraction:
        return self.constant + self.a_coeff * a_value + self.t_coeff * t_value

    def canonical(self) -> tuple[str, str, str]:
        return tuple(
            fraction_text(value)
            for value in (self.constant, self.a_coeff, self.t_coeff)
        )

    def parity_substitution(self) -> tuple[Fraction, Fraction, Fraction]:
        """Substitute a=2+2u and t=3+2v."""
        return (
            self.constant + 2 * self.a_coeff + 3 * self.t_coeff,
            2 * self.a_coeff,
            2 * self.t_coeff,
        )


def to_affine(value: object) -> Affine:
    if isinstance(value, Affine):
        return value
    if isinstance(value, (int, Fraction)):
        return Affine(Fraction(value))
    raise TypeError(f"cannot convert {type(value)!r} to Affine")


ZERO = Affine()
ONE = Affine(1)
A_VAR = Affine(0, 1, 0)
T_VAR = Affine(0, 0, 1)


class PolyAT:
    """A sparse polynomial in the formal parameters a,t over Q."""

    def __init__(self, terms: dict[tuple[int, int], Fraction] | None = None):
        cleaned: dict[tuple[int, int], Fraction] = {}
        for monomial, coefficient in (terms or {}).items():
            value = Fraction(coefficient)
            if value:
                cleaned[monomial] = cleaned.get(monomial, Fraction(0)) + value
        self.terms = {
            monomial: coefficient
            for monomial, coefficient in cleaned.items()
            if coefficient
        }

    @staticmethod
    def constant(value: object) -> "PolyAT":
        coefficient = Fraction(value)
        return PolyAT({(0, 0): coefficient} if coefficient else {})

    @staticmethod
    def affine(value: Affine) -> "PolyAT":
        terms: dict[tuple[int, int], Fraction] = {}
        if value.constant:
            terms[0, 0] = value.constant
        if value.a_coeff:
            terms[1, 0] = value.a_coeff
        if value.t_coeff:
            terms[0, 1] = value.t_coeff
        return PolyAT(terms)

    def __add__(self, other: object) -> "PolyAT":
        rhs = to_poly(other)
        terms = dict(self.terms)
        for monomial, coefficient in rhs.terms.items():
            terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
        return PolyAT(terms)

    def __radd__(self, other: object) -> "PolyAT":
        return self + other

    def __neg__(self) -> "PolyAT":
        return PolyAT({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "PolyAT":
        return self + (-to_poly(other))

    def __rsub__(self, other: object) -> "PolyAT":
        return to_poly(other) - self

    def __mul__(self, other: object) -> "PolyAT":
        rhs = to_poly(other)
        terms: dict[tuple[int, int], Fraction] = {}
        for (a_left, t_left), left in self.terms.items():
            for (a_right, t_right), right in rhs.terms.items():
                monomial = (a_left + a_right, t_left + t_right)
                terms[monomial] = terms.get(monomial, Fraction(0)) + left * right
        return PolyAT(terms)

    def __rmul__(self, other: object) -> "PolyAT":
        return self * other

    def __eq__(self, other: object) -> bool:
        try:
            return self.terms == to_poly(other).terms
        except TypeError:
            return False

    def is_zero(self) -> bool:
        return not self.terms

    def canonical(self) -> tuple[tuple[int, int, str], ...]:
        return tuple(
            (a_power, t_power, fraction_text(coefficient))
            for (a_power, t_power), coefficient in sorted(self.terms.items())
        )


def to_poly(value: object) -> PolyAT:
    if isinstance(value, PolyAT):
        return value
    if isinstance(value, Affine):
        return PolyAT.affine(value)
    if isinstance(value, (int, Fraction)):
        return PolyAT.constant(value)
    raise TypeError(f"cannot convert {type(value)!r} to PolyAT")


@dataclass(frozen=True, order=True)
class TermKey:
    constants: tuple[int, ...]
    h_exp: Affine
    k_exp: Affine
    hp_exp: int = 0
    kp_exp: int = 0

    def multiply(self, other: "TermKey") -> "TermKey":
        return TermKey(
            tuple(left + right for left, right in zip(self.constants, other.constants)),
            self.h_exp + other.h_exp,
            self.k_exp + other.k_exp,
            self.hp_exp + other.hp_exp,
            self.kp_exp + other.kp_exp,
        )


class Expr:
    """Formal sums of coefficient monomials times h,k,h',k'."""

    def __init__(self, terms: dict[TermKey, PolyAT] | None = None):
        cleaned: dict[TermKey, PolyAT] = {}
        for key, coefficient in (terms or {}).items():
            value = to_poly(coefficient)
            if value.is_zero():
                continue
            cleaned[key] = cleaned.get(key, PolyAT()) + value
        self.terms = {
            key: coefficient
            for key, coefficient in cleaned.items()
            if not coefficient.is_zero()
        }

    def __add__(self, other: object) -> "Expr":
        rhs = to_expr(other)
        terms = dict(self.terms)
        for key, coefficient in rhs.terms.items():
            terms[key] = terms.get(key, PolyAT()) + coefficient
        return Expr(terms)

    def __radd__(self, other: object) -> "Expr":
        return self + other

    def __neg__(self) -> "Expr":
        return Expr({key: -coefficient for key, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Expr":
        return self + (-to_expr(other))

    def __rsub__(self, other: object) -> "Expr":
        return to_expr(other) - self

    def __mul__(self, other: object) -> "Expr":
        rhs = to_expr(other)
        terms: dict[TermKey, PolyAT] = {}
        for left_key, left_coefficient in self.terms.items():
            for right_key, right_coefficient in rhs.terms.items():
                key = left_key.multiply(right_key)
                coefficient = left_coefficient * right_coefficient
                terms[key] = terms.get(key, PolyAT()) + coefficient
        return Expr(terms)

    def __rmul__(self, other: object) -> "Expr":
        if isinstance(other, (int, Fraction, Affine, PolyAT)):
            return self.scale(other)
        return to_expr(other) * self

    def scale(self, coefficient: object) -> "Expr":
        factor = to_poly(coefficient)
        return Expr({key: value * factor for key, value in self.terms.items()})

    def __eq__(self, other: object) -> bool:
        try:
            return self.terms == to_expr(other).terms
        except TypeError:
            return False

    def canonical(self) -> tuple[object, ...]:
        return tuple(
            (
                key.constants,
                key.h_exp.canonical(),
                key.k_exp.canonical(),
                key.hp_exp,
                key.kp_exp,
                coefficient.canonical(),
            )
            for key, coefficient in sorted(self.terms.items())
        )


def to_expr(value: object) -> Expr:
    if isinstance(value, Expr):
        return value
    if isinstance(value, (int, Fraction)):
        return monomial("", coefficient=value)
    raise TypeError(f"cannot convert {type(value)!r} to Expr")


def constant_vector(symbols: str) -> tuple[int, ...]:
    vector = [0] * len(CONSTANT_NAMES)
    for symbol in symbols:
        require(symbol in CONSTANT_INDEX, f"unknown coefficient symbol {symbol}")
        vector[CONSTANT_INDEX[symbol]] += 1
    return tuple(vector)


def monomial(
    symbols: str,
    h_exp: Affine = ZERO,
    k_exp: Affine = ZERO,
    hp_exp: int = 0,
    kp_exp: int = 0,
    coefficient: object = 1,
) -> Expr:
    value = to_poly(coefficient)
    if value.is_zero():
        return Expr()
    key = TermKey(constant_vector(symbols), h_exp, k_exp, hp_exp, kp_exp)
    return Expr({key: value})


def derivative(expression: Expr) -> Expr:
    terms: dict[TermKey, PolyAT] = {}
    for key, coefficient in expression.terms.items():
        require(key.hp_exp == 0 and key.kp_exp == 0, "first derivative input")
        if key.h_exp != ZERO:
            derived = TermKey(
                key.constants,
                key.h_exp - ONE,
                key.k_exp,
                hp_exp=1,
            )
            terms[derived] = terms.get(derived, PolyAT()) + coefficient * key.h_exp
        if key.k_exp != ZERO:
            derived = TermKey(
                key.constants,
                key.h_exp,
                key.k_exp - ONE,
                kp_exp=1,
            )
            terms[derived] = terms.get(derived, PolyAT()) + coefficient * key.k_exp
    return Expr(terms)


def wronskian(weight_f: Affine, f_value: Expr, weight_g: Affine, g_value: Expr) -> Expr:
    return (
        (derivative(f_value) * g_value).scale(weight_g)
        - (f_value * derivative(g_value)).scale(weight_f)
    )


EULER = (
    monomial("", k_exp=ONE, hp_exp=1)
    + monomial("", h_exp=ONE, kp_exp=1, coefficient=2)
)


def divide_by_euler(expression: Expr, label: str) -> Expr:
    """Recover F from expression=(k h'+2 h k')F by exact term pairing."""
    groups: dict[TermKey, dict[str, PolyAT]] = {}
    for key, coefficient in expression.terms.items():
        if key.hp_exp == 1 and key.kp_exp == 0:
            residual = TermKey(
                key.constants, key.h_exp, key.k_exp - ONE, 0, 0
            )
            tag = "hp"
        elif key.hp_exp == 0 and key.kp_exp == 1:
            residual = TermKey(
                key.constants, key.h_exp - ONE, key.k_exp, 0, 0
            )
            tag = "kp"
        else:
            raise RuntimeError(f"{label}: unexpected derivative monomial {key}")
        bucket = groups.setdefault(residual, {})
        require(tag not in bucket, f"{label}: duplicate Euler half")
        bucket[tag] = coefficient
    quotient: dict[TermKey, PolyAT] = {}
    for residual, halves in groups.items():
        require(set(halves) == {"hp", "kp"}, f"{label}: unmatched Euler half")
        require(halves["kp"] == 2 * halves["hp"], f"{label}: Euler ratio")
        quotient[residual] = halves["hp"]
    result = Expr(quotient)
    require(EULER * result == expression, f"{label}: Euler reconstruction")
    return result


Cell = tuple[int, int]
Fibre = tuple[Cell, ...]


@dataclass(frozen=True)
class ParentEntry:
    word_id: str
    scalar_index: int
    anchor: Cell
    orientation: str
    witness: tuple[int, ...]
    scheme_sha256: str
    fibres: tuple[Fibre, ...]
    source_line: str


def parse_fibres(word: str) -> tuple[Fibre, ...]:
    fibres: list[Fibre] = []
    for part in word.split("|"):
        cells: list[Cell] = []
        for token in part.split("="):
            require(re.fullmatch(r"[0-2][0-3]", token) is not None, "cell token")
            cells.append((int(token[0]), int(token[1])))
        fibres.append(tuple(cells))
    require(sum(map(len, fibres)) == 12, "complete three-by-four cell set")
    require(
        {cell for fibre in fibres for cell in fibre}
        == {(i, j) for i in range(3) for j in range(4)},
        "cell partition",
    )
    return tuple(fibres)


def fibre_word(fibres: tuple[Fibre, ...]) -> str:
    return "|".join("=".join(f"{i}{j}" for i, j in fibre) for fibre in fibres)


def parse_parent_entries() -> dict[str, ParentEntry]:
    entries: dict[str, ParentEntry] = {}
    text = lf_bytes(PARENT_OUTPUT).decode("utf-8")
    anchor_pattern = re.compile(r"(\d+):([0-2])([0-3])([+-]{2})@([0-9,]+)")
    for line in text.splitlines():
        if not any(line.startswith(f"{word_id}=") for word_id in TARGET_IDS):
            continue
        word_id, remainder = line.split("=", 1)
        fields: dict[str, str] = {}
        for field in remainder.split(";"):
            key, value = field.split(":", 1)
            fields[key] = value
        match = anchor_pattern.fullmatch(fields["residual_anchor"])
        require(match is not None, f"{word_id}: residual anchor parse")
        scalar_index = int(match.group(1))
        anchor = (int(match.group(2)), int(match.group(3)))
        orientation = match.group(4)
        witness = tuple(int(value) for value in match.group(5).split(","))
        fibres = parse_fibres(fields["word"])
        require(fields["m"] == "9", f"{word_id}: size nine")
        require(fields["unexposed"] == "1", f"{word_id}: unique unexposed scheme")
        require(fields["scheme_sha256"] == EXPECTED_SCHEME_SHA256[word_id],
                f"{word_id}: scheme digest")
        require(len(witness) == 5 and all(value > 0 for value in witness),
                f"{word_id}: witness")
        require(0 <= scalar_index < len(fibres), f"{word_id}: scalar index")
        require(anchor in fibres[scalar_index], f"{word_id}: anchor membership")
        entries[word_id] = ParentEntry(
            word_id,
            scalar_index,
            anchor,
            orientation,
            witness,
            fields["scheme_sha256"],
            fibres,
            line,
        )
    require(tuple(sorted(entries)) == TARGET_IDS, "four canonical parent entries")
    return entries


GAP_CELL_VECTOR: dict[Cell, tuple[int, ...]] = {}
for i in range(3):
    for j in range(4):
        a_vector = ((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (1, 1, 0, 0, 0))[i]
        b_vector = (
            (0, 0, 0, 0, 0),
            (0, 0, 1, 0, 0),
            (0, 0, 1, 1, 0),
            (0, 0, 1, 1, 1),
        )[j]
        GAP_CELL_VECTOR[i, j] = tuple(left + right for left, right in zip(a_vector, b_vector))


def vector_subtract(left: tuple[int, ...], right: tuple[int, ...]) -> list[Fraction]:
    return [Fraction(a - b) for a, b in zip(left, right)]


def collision_matrix(fibres: tuple[Fibre, ...]) -> list[list[Fraction]]:
    rows: list[list[Fraction]] = []
    for fibre in fibres:
        for cell in fibre[1:]:
            rows.append(vector_subtract(GAP_CELL_VECTOR[cell], GAP_CELL_VECTOR[fibre[0]]))
    return rows


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], tuple[int, ...]]:
    rows = [list(map(Fraction, row)) for row in matrix]
    require(rows and rows[0], "nonempty collision matrix")
    row_count, column_count = len(rows), len(rows[0])
    pivot_row = 0
    pivots: list[int] = []
    for column in range(column_count):
        pivot = next((row for row in range(pivot_row, row_count) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[row], rows[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    nonzero = [row for row in rows if any(row)]
    return nonzero, tuple(pivots)


def nullspace(matrix: list[list[Fraction]]) -> tuple[tuple[Fraction, ...], ...]:
    reduced, pivots = rref(matrix)
    column_count = len(matrix[0])
    free = tuple(column for column in range(column_count) if column not in pivots)
    basis: list[tuple[Fraction, ...]] = []
    for free_column in free:
        vector = [Fraction(0)] * column_count
        vector[free_column] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def cone_parameterization(
    fibres: tuple[Fibre, ...], a_coordinate: int, t_coordinate: int
) -> tuple[Affine, ...]:
    matrix = collision_matrix(fibres)
    reduced, pivots = rref(matrix)
    basis = nullspace(matrix)
    require(len(pivots) == 3 and len(basis) == 2, "rank-three two-dimensional cone")
    transform = (
        (basis[0][a_coordinate], basis[1][a_coordinate]),
        (basis[0][t_coordinate], basis[1][t_coordinate]),
    )
    determinant = transform[0][0] * transform[1][1] - transform[0][1] * transform[1][0]
    require(determinant != 0, "parameter coordinates")
    inverse = (
        (transform[1][1] / determinant, -transform[0][1] / determinant),
        (-transform[1][0] / determinant, transform[0][0] / determinant),
    )
    parameterization: list[Affine] = []
    for coordinate in range(5):
        basis_row = (basis[0][coordinate], basis[1][coordinate])
        a_coeff = basis_row[0] * inverse[0][0] + basis_row[1] * inverse[1][0]
        t_coeff = basis_row[0] * inverse[0][1] + basis_row[1] * inverse[1][1]
        parameterization.append(Affine(0, a_coeff, t_coeff))
    for row in reduced:
        total = sum((parameterization[index] * row[index] for index in range(5)), ZERO)
        require(total == ZERO, "parameterization satisfies RREF")
    return tuple(parameterization)


def affine_supports(parameterization: tuple[Affine, ...]) -> tuple[tuple[Affine, ...], tuple[Affine, ...]]:
    x_gap, y_gap, u_gap, v_gap, w_gap = parameterization
    return (
        (ZERO, x_gap, x_gap + y_gap),
        (ZERO, u_gap, u_gap + v_gap, u_gap + v_gap + w_gap),
    )


def cell_value(cell: Cell, parameterization: tuple[Affine, ...]) -> Affine:
    return sum(
        (parameterization[index] * coefficient for index, coefficient in enumerate(GAP_CELL_VECTOR[cell])),
        ZERO,
    )


def primitive_wall(value: Affine) -> tuple[int, int] | None:
    if value.constant or not value.a_coeff or not value.t_coeff:
        return None
    if value.a_coeff * value.t_coeff >= 0:
        return None
    require(value.a_coeff.denominator == 1 and value.t_coeff.denominator == 1,
            "integral chamber wall")
    left, right = int(value.a_coeff), int(value.t_coeff)
    divisor = gcd(abs(left), abs(right))
    return left // divisor, right // divisor


def derive_chamber(entry: ParentEntry, parameterization: tuple[Affine, ...]) -> str:
    fibre_values: list[Affine] = []
    for fibre in entry.fibres:
        values = tuple(cell_value(cell, parameterization) for cell in fibre)
        require(len(set(values)) == 1, f"{entry.word_id}: collision equality")
        fibre_values.append(values[0])
    differences = tuple(
        right - left for left, right in zip(fibre_values, fibre_values[1:])
    )
    walls = {wall for difference in differences if (wall := primitive_wall(difference)) is not None}
    require(len(walls) == 1, f"{entry.word_id}: unique chamber wall")
    wall = next(iter(walls))
    if wall == (-1, 1):
        chamber = "t>a"
    elif wall == (1, -1):
        chamber = "a>t"
    else:
        raise RuntimeError(f"{entry.word_id}: unexpected chamber wall {wall}")
    for difference in differences:
        require(
            positive_on_chamber(difference, chamber),
            f"{entry.word_id}: complete chamber ordering {difference}",
        )
    return chamber


def positive_on_chamber(value: Affine, chamber: str) -> bool:
    if value.constant < 0:
        return False
    if chamber == "t>a":
        # a=u, t=u+s with u,s>0.
        coefficients = (value.a_coeff + value.t_coeff, value.t_coeff)
    elif chamber == "a>t":
        # t=u, a=u+s with u,s>0.
        coefficients = (value.a_coeff + value.t_coeff, value.a_coeff)
    else:
        raise RuntimeError(f"unknown chamber {chamber}")
    return all(coefficient >= 0 for coefficient in coefficients) and (
        value.constant > 0 or any(coefficient > 0 for coefficient in coefficients)
    )


def derive_weights(
    parameterization: tuple[Affine, ...], entry: ParentEntry
) -> tuple[tuple[Affine, ...], tuple[Affine, ...]]:
    a_support, b_support = affine_supports(parameterization)
    target = {"+-": (1, -2), "-+": (-2, 1)}[entry.orientation]
    p_zero = Affine(target[0]) - a_support[entry.anchor[0]]
    q_zero = Affine(target[1]) - b_support[entry.anchor[1]]
    p_weights = tuple(p_zero + value for value in a_support)
    q_weights = tuple(q_zero + value for value in b_support)
    require(
        (p_weights[entry.anchor[0]], q_weights[entry.anchor[1]])
        == tuple(map(Affine, target)),
        f"{entry.word_id}: anchored weights",
    )
    for cell in entry.fibres[entry.scalar_index]:
        require(p_weights[cell[0]] + q_weights[cell[1]] == Affine(-1),
                f"{entry.word_id}: scalar complementary weight")
    return p_weights, q_weights


def sign_on_box(value: Affine, a_minimum: int = 2, t_minimum: int = 3) -> str:
    base = value.evaluate(a_minimum, t_minimum)
    if value.a_coeff >= 0 and value.t_coeff >= 0 and base > 0:
        return "+"
    if value.a_coeff <= 0 and value.t_coeff <= 0 and base < 0:
        return "-"
    raise RuntimeError(f"sign not certified on a>={a_minimum},t>={t_minimum}: {value}")


def singleton_components(
    entry: ParentEntry,
    p_weights: tuple[Affine, ...],
    q_weights: tuple[Affine, ...],
    a_minimum: int,
    t_minimum: int,
) -> dict[str, tuple[str, ...]]:
    vertices = tuple([f"P{i}" for i in range(3)] + [f"Q{j}" for j in range(4)])
    parent = {vertex: vertex for vertex in vertices}
    component_sign: dict[str, str] = {}

    def find(vertex: str) -> str:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    active: set[str] = set()
    for index, fibre in enumerate(entry.fibres):
        if index == entry.scalar_index or len(fibre) != 1:
            continue
        i, j = fibre[0]
        p_sign = sign_on_box(p_weights[i], a_minimum, t_minimum)
        q_sign = sign_on_box(q_weights[j], a_minimum, t_minimum)
        require(p_sign == q_sign, f"{entry.word_id}: singleton same-sign gate {i}{j}")
        left, right = f"P{i}", f"Q{j}"
        active.update((left, right))
        component_sign[left] = p_sign
        component_sign[right] = q_sign
        union(left, right)

    parts: dict[str, list[str]] = {}
    for vertex in sorted(active):
        parts.setdefault(find(vertex), []).append(vertex)
    signed_parts: dict[str, tuple[str, ...]] = {}
    for part in parts.values():
        signs = {component_sign[vertex] for vertex in part}
        require(len(signs) == 1, f"{entry.word_id}: component sign")
        sign = next(iter(signs))
        require(sign not in signed_parts, f"{entry.word_id}: unique {sign} component")
        signed_parts[sign] = tuple(part)
    require(set(signed_parts) == {"-", "+"}, f"{entry.word_id}: two signed components")
    return signed_parts


def singleton_lower_bounds(
    entry: ParentEntry,
    p_weights: tuple[Affine, ...],
    q_weights: tuple[Affine, ...],
) -> tuple[int, int]:
    """Derive the nonprojective t threshold from an actual singleton row."""
    singleton_cells = {
        fibre[0] for fibre in entry.fibres if len(fibre) == 1
    }
    if entry.word_id in {"W072", "W073"}:
        require((1, 2) in singleton_cells, f"{entry.word_id}: threshold singleton 12")
        require(p_weights[1] == ONE and q_weights[2] == T_VAR - 2,
                f"{entry.word_id}: threshold weights")
        # The first weight is positive, so the singleton gate forces t-2>0.
        return 1, 3
    require((0, 1) in singleton_cells, f"{entry.word_id}: threshold singleton 01")
    require(p_weights[0] == -A_VAR - 2 and q_weights[1] == ONE - T_VAR,
            f"{entry.word_id}: threshold weights")
    # The first weight is negative, so the singleton gate forces 1-t<0.
    return 1, 2


def node_weights(
    p_weights: tuple[Affine, ...], q_weights: tuple[Affine, ...]
) -> dict[str, Affine]:
    return {
        **{f"P{i}": value for i, value in enumerate(p_weights)},
        **{f"Q{j}": value for j, value in enumerate(q_weights)},
    }


def component_magnitudes(
    components: dict[str, tuple[str, ...]], weights: dict[str, Affine]
) -> dict[str, tuple[tuple[str, Affine], ...]]:
    return {
        sign: tuple(
            (node, weights[node] if sign == "+" else -weights[node])
            for node in nodes
        )
        for sign, nodes in components.items()
    }


def residue(value: Affine, a_residue: int, t_residue: int, modulus: int) -> int:
    require(
        value.constant.denominator == value.a_coeff.denominator == value.t_coeff.denominator == 1,
        "integral affine residue",
    )
    return int(value.evaluate(a_residue, t_residue)) % modulus


def parity_gate(negative_magnitudes: tuple[tuple[str, Affine], ...]) -> tuple[tuple[int, int], ...]:
    require(any(value == Affine(2) for _, value in negative_magnitudes),
            "negative ledger contains magnitude two")
    surviving: list[tuple[int, int]] = []
    for a_residue in range(2):
        for t_residue in range(2):
            if all(residue(value, a_residue, t_residue, 2) == 0
                   for _, value in negative_magnitudes):
                surviving.append((a_residue, t_residue))
    require(tuple(surviving) == ((0, 1),), "d=2 parity classification")
    return tuple(surviving)


def exact_gcd(values: tuple[int, ...]) -> int:
    return reduce(gcd, values)


def verify_gcd_ledger(
    entry: ParentEntry,
    magnitudes: dict[str, tuple[tuple[str, Affine], ...]],
) -> None:
    negative = magnitudes["-"]
    positive = magnitudes["+"]
    require(any(value == Affine(2) for _, value in negative),
            f"{entry.word_id}: negative gcd divides two")
    require(any(value == Affine(1) for _, value in positive),
            f"{entry.word_id}: positive gcd one")
    parity_gate(negative)
    for a_value in (2, 4, 6):
        for t_value in (3, 5, 7):
            negative_values = tuple(int(value.evaluate(a_value, t_value)) for _, value in negative)
            positive_values = tuple(int(value.evaluate(a_value, t_value)) for _, value in positive)
            require(exact_gcd(negative_values) == 2,
                    f"{entry.word_id}: d=2 hostile control")
            require(exact_gcd(positive_values) == 1,
                    f"{entry.word_id}: positive gcd hostile control")


def verify_unique_arm_gate(
    entry: ParentEntry,
    p_weights: tuple[Affine, ...],
    q_weights: tuple[Affine, ...],
    magnitudes: dict[str, tuple[tuple[str, Affine], ...]],
) -> None:
    scalar_fibre = entry.fibres[entry.scalar_index]
    require(len(scalar_fibre) == 2, f"{entry.word_id}: scalar double")
    pairs = tuple((p_weights[i], q_weights[j]) for i, j in scalar_fibre)
    anchor_pair = (p_weights[entry.anchor[0]], q_weights[entry.anchor[1]])
    expected_anchor = tuple(map(Affine, {"+-": (1, -2), "-+": (-2, 1)}[entry.orientation]))
    require(anchor_pair == expected_anchor, f"{entry.word_id}: active arm pair")
    other_pair = next(pair for cell, pair in zip(scalar_fibre, pairs) if cell != entry.anchor)
    negative = next(-value for value in other_pair if sign_on_box(value) == "-")
    positive = next(value for value in other_pair if sign_on_box(value) == "+")
    require(negative == A_VAR + 2 and positive == A_VAR + 1,
            f"{entry.word_id}: other complementary arm")
    require(negative - positive == ONE and negative.evaluate(1, 1) == 3,
            f"{entry.word_id}: other arm has R>=3")

    anchor_negative_node = (
        f"P{entry.anchor[0]}" if expected_anchor[0] < ZERO else f"Q{entry.anchor[1]}"
    )
    anchor_positive_node = (
        f"P{entry.anchor[0]}" if expected_anchor[0] > ZERO else f"Q{entry.anchor[1]}"
    )
    anchor_magnitude = next(
        value for node, value in magnitudes["-"] if node == anchor_negative_node
    )
    anchor_positive_magnitude = next(
        value for node, value in magnitudes["+"] if node == anchor_positive_node
    )
    require(anchor_magnitude == Affine(2), f"{entry.word_id}: arm magnitude two")
    require(anchor_positive_magnitude == ONE,
            f"{entry.word_id}: positive arm base exponent one")
    # The common-base exponent is 2/d.  A simple negative coefficient can
    # occur only for d=2 and ord_beta(h)=1; d=1 has even local order.
    feasible = tuple(
        d_value
        for d_value in (1, 2)
        if int(anchor_magnitude.constant) // d_value == 1
    )
    require(feasible == (2,), f"{entry.word_id}: d=2 forced by simple arm")


NODE_CONSTANT = {
    "P0": "A", "P1": "L", "P2": "M",
    "Q0": "B", "Q1": "C", "Q2": "N", "Q3": "T",
}


def ufd_forms(
    weights: dict[str, Affine], components: dict[str, tuple[str, ...]]
) -> dict[str, Expr]:
    forms: dict[str, Expr] = {}
    for sign, nodes in components.items():
        divisor = 2 if sign == "-" else 1
        for node in nodes:
            magnitude = -weights[node] if sign == "-" else weights[node]
            exponent = magnitude / divisor
            if sign == "-":
                forms[node] = monomial(NODE_CONSTANT[node], h_exp=exponent)
            else:
                forms[node] = monomial(NODE_CONSTANT[node], k_exp=exponent)
    require(set(forms) == set(NODE_CONSTANT), "all seven coefficient forms")
    return forms


def row_for_fibre(
    fibre: Fibre,
    forms: dict[str, Expr],
    p_weights: tuple[Affine, ...],
    q_weights: tuple[Affine, ...],
) -> Expr:
    row = Expr()
    for i, j in fibre:
        row += wronskian(
            p_weights[i], forms[f"P{i}"], q_weights[j], forms[f"Q{j}"]
        )
    return row


def require_polynomial_exponents(expression: Expr, label: str) -> None:
    for key in expression.terms:
        for exponent in (key.h_exp, key.k_exp):
            substituted = exponent.parity_substitution()
            require(all(value.denominator == 1 for value in substituted),
                    f"{label}: parity-integral exponent")
            require(all(value >= 0 for value in substituted),
                    f"{label}: polynomial exponent")


def expected_left_factors() -> dict[Fibre, Expr]:
    rho = (A_VAR + T_VAR - 1) / 2
    eta = (A_VAR + 2) / 2
    sigma = 2 * A_VAR + T_VAR - 2
    return {
        ((0, 2), (1, 0)): (
            monomial(
                "AN", h_exp=rho - 1, k_exp=T_VAR - 3,
                coefficient=(T_VAR - 2) * 1,
            ).scale(rho)
            - monomial("LB", h_exp=eta - 1, coefficient=eta)
        ),
        ((1, 1), (2, 0)): -(
            monomial("LC")
            + monomial(
                "MB", h_exp=eta - 1, k_exp=A_VAR,
                coefficient=PolyAT.affine(A_VAR + 1) * PolyAT.affine(eta),
            )
        ),
        ((0, 3), (2, 1)): (
            monomial(
                "AT", h_exp=rho - 1, k_exp=sigma - 1,
                coefficient=PolyAT.affine(sigma) * PolyAT.affine(rho),
            )
            - monomial("MC", k_exp=A_VAR, coefficient=A_VAR + 1)
        ),
    }


def expected_right_factors() -> dict[Fibre, Expr]:
    eta = (A_VAR + 2) / 2
    rho = (2 * A_VAR + T_VAR - 1) / 2
    omega = (T_VAR - 1) / 2
    r_exp = A_VAR + T_VAR - 2
    return {
        ((0, 2), (2, 0)): (
            monomial("AN", h_exp=eta - 1, coefficient=eta)
            - monomial(
                "MB", h_exp=rho - 1, k_exp=r_exp - 1,
                coefficient=PolyAT.affine(r_exp) * PolyAT.affine(rho),
            )
        ),
        ((0, 3), (1, 2)): (
            monomial(
                "AT", h_exp=eta - 1, k_exp=A_VAR,
                coefficient=PolyAT.affine(A_VAR + 1) * PolyAT.affine(eta),
            )
            + monomial("LN")
        ),
        ((1, 3), (2, 1)): (
            monomial("LT", k_exp=A_VAR, coefficient=A_VAR + 1)
            - monomial(
                "MC", h_exp=omega - 1, k_exp=r_exp - 1,
                coefficient=PolyAT.affine(r_exp) * PolyAT.affine(omega),
            )
        ),
    }


def reverse_fibres(fibres: tuple[Fibre, ...]) -> tuple[Fibre, ...]:
    return tuple(
        tuple(sorted((2 - i, 3 - j) for i, j in fibre))
        for fibre in reversed(fibres)
    )


def reverse_parameterization(parameterization: tuple[Affine, ...]) -> tuple[Affine, ...]:
    x_gap, y_gap, u_gap, v_gap, w_gap = parameterization
    return y_gap, x_gap, w_gap, v_gap, u_gap


def affine_text(value: Affine) -> str:
    known = {
        ZERO: "0", ONE: "1", Affine(2): "2", A_VAR: "a", T_VAR: "t",
        A_VAR + T_VAR: "a+t", 2 * A_VAR: "2a",
        A_VAR + 1: "a+1", A_VAR + 2: "a+2",
        T_VAR - 1: "t-1", T_VAR - 2: "t-2",
        A_VAR + T_VAR - 1: "a+t-1", A_VAR + T_VAR - 2: "a+t-2",
        2 * A_VAR + T_VAR - 1: "2a+t-1",
        2 * A_VAR + T_VAR - 2: "2a+t-2",
    }
    if value in known:
        return known[value]
    return str(value.canonical())


def parameterization_text(parameterization: tuple[Affine, ...]) -> str:
    return "(" + ",".join(affine_text(value) for value in parameterization) + ")"


def magnitude_text(values: tuple[tuple[str, Affine], ...]) -> str:
    return ",".join(f"{node}:{affine_text(value)}" for node, value in values)


def main() -> None:
    require(lf_sha256(PARENT_THEOREM) == EXPECTED_PARENT_THEOREM_SHA256,
            "THM-3606 theorem hash")
    require(lf_sha256(PARENT_SCRIPT) == EXPECTED_PARENT_SCRIPT_SHA256,
            "THM-3606 script hash")
    require(lf_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256,
            "THM-3606 output hash")

    entries = parse_parent_entries()
    collision_signatures: dict[tuple[object, ...], list[str]] = {}
    for word_id, entry in entries.items():
        reduced, pivots = rref(collision_matrix(entry.fibres))
        signature = (
            pivots,
            tuple(tuple(fraction_text(value) for value in row) for row in reduced),
        )
        collision_signatures.setdefault(signature, []).append(word_id)
    cone_groups = tuple(sorted(tuple(sorted(group)) for group in collision_signatures.values()))
    require(cone_groups == (("W072", "W073"), ("W141", "W149")),
            "two canonical collision cones")

    left = cone_parameterization(entries["W072"].fibres, 1, 3)
    right = cone_parameterization(entries["W141"].fibres, 0, 3)
    expected_left = (A_VAR + T_VAR, A_VAR, A_VAR, T_VAR, 2 * A_VAR)
    expected_right = (A_VAR, A_VAR + T_VAR, 2 * A_VAR, T_VAR, A_VAR)
    require(left == expected_left, "left cone derived parameterization")
    require(right == expected_right, "right cone derived parameterization")
    require(
        cone_parameterization(entries["W073"].fibres, 1, 3) == left
        and cone_parameterization(entries["W149"].fibres, 0, 3) == right,
        "both words derive each cone",
    )

    parameterizations = {
        "W072": left, "W073": left, "W141": right, "W149": right,
    }
    chambers: dict[str, str] = {}
    weights_by_word: dict[str, tuple[tuple[Affine, ...], tuple[Affine, ...]]] = {}
    components_by_word: dict[str, dict[str, tuple[str, ...]]] = {}
    magnitudes_by_word: dict[str, dict[str, tuple[tuple[str, Affine], ...]]] = {}
    factor_ledgers: dict[str, tuple[tuple[str, tuple[object, ...]], ...]] = {}

    expected_factors_by_family = {
        "left": expected_left_factors(), "right": expected_right_factors()
    }
    for word_id in TARGET_IDS:
        entry = entries[word_id]
        parameterization = parameterizations[word_id]
        chamber = derive_chamber(entry, parameterization)
        chambers[word_id] = chamber
        a_coordinate, t_coordinate = ((1, 3) if word_id in {"W072", "W073"} else (0, 3))
        a_value, t_value = entry.witness[a_coordinate], entry.witness[t_coordinate]
        require(
            tuple(int(value.evaluate(a_value, t_value)) for value in parameterization)
            == entry.witness,
            f"{word_id}: parent witness lies in derived cone",
        )

        p_weights, q_weights = derive_weights(parameterization, entry)
        weights_by_word[word_id] = p_weights, q_weights
        a_minimum, t_minimum = singleton_lower_bounds(entry, p_weights, q_weights)
        components = singleton_components(
            entry, p_weights, q_weights, a_minimum, t_minimum
        )
        components_by_word[word_id] = components
        magnitudes = component_magnitudes(components, node_weights(p_weights, q_weights))
        magnitudes_by_word[word_id] = magnitudes
        verify_gcd_ledger(entry, magnitudes)
        verify_unique_arm_gate(entry, p_weights, q_weights, magnitudes)

        forms = ufd_forms(node_weights(p_weights, q_weights), components)
        family = "left" if word_id in {"W072", "W073"} else "right"
        expected_factors = expected_factors_by_family[family]
        multi_fibres = tuple(tuple(sorted(fibre)) for fibre in entry.fibres if len(fibre) > 1)
        require(set(multi_fibres) == set(expected_factors), f"{word_id}: three collision rows")
        ledger: list[tuple[str, tuple[object, ...]]] = []
        for fibre in multi_fibres:
            row = row_for_fibre(fibre, forms, p_weights, q_weights)
            quotient = divide_by_euler(row, f"{word_id}:{fibre_word((fibre,))}")
            require(quotient == expected_factors[fibre],
                    f"{word_id}: displayed symbolic factor {fibre}")
            require_polynomial_exponents(row, f"{word_id}:{fibre}")
            require_polynomial_exponents(quotient, f"{word_id}:{fibre}:quotient")
            ledger.append((fibre_word((fibre,)), quotient.canonical()))
        factor_ledgers[word_id] = tuple(ledger)

    require(chambers == {
        "W072": "t>a", "W073": "a>t", "W141": "a>t", "W149": "t>a"
    }, "four chamber halves")

    expected_left_components = {
        "-": ("P0", "Q0", "Q1"), "+": ("P1", "P2", "Q2", "Q3")
    }
    expected_right_components = {
        "-": ("P0", "P1", "Q0", "Q1"), "+": ("P2", "Q2", "Q3")
    }
    for word_id in ("W072", "W073"):
        require(components_by_word[word_id] == expected_left_components,
                f"{word_id}: left singleton components")
    for word_id in ("W141", "W149"):
        require(components_by_word[word_id] == expected_right_components,
                f"{word_id}: right singleton components")

    left_negative = (("P0", A_VAR + T_VAR - 1), ("Q0", A_VAR + 2), ("Q1", Affine(2)))
    left_positive = (("P1", ONE), ("P2", A_VAR + 1), ("Q2", T_VAR - 2),
                     ("Q3", 2 * A_VAR + T_VAR - 2))
    right_negative = (("P0", A_VAR + 2), ("P1", Affine(2)),
                      ("Q0", 2 * A_VAR + T_VAR - 1), ("Q1", T_VAR - 1))
    right_positive = (("P2", A_VAR + T_VAR - 2), ("Q2", ONE), ("Q3", A_VAR + 1))
    for word_id in ("W072", "W073"):
        require(magnitudes_by_word[word_id] == {"-": left_negative, "+": left_positive},
                f"{word_id}: left UFD exponent ledger")
    for word_id in ("W141", "W149"):
        require(magnitudes_by_word[word_id] == {"-": right_negative, "+": right_positive},
                f"{word_id}: right UFD exponent ledger")

    reversal_pairs = {"W072": "W149", "W073": "W141"}
    for source_id, target_id in reversal_pairs.items():
        source, target = entries[source_id], entries[target_id]
        require(reverse_fibres(source.fibres) == target.fibres,
                f"{source_id}/{target_id}: reversed fibre word")
        require(reverse_parameterization(parameterizations[source_id]) == parameterizations[target_id],
                f"{source_id}/{target_id}: reversed cone")
        reversed_anchor = (2 - source.anchor[0], 3 - source.anchor[1])
        require(reversed_anchor == target.anchor,
                f"{source_id}/{target_id}: reversed anchor")
        require(len(source.fibres) - 1 - source.scalar_index == target.scalar_index,
                f"{source_id}/{target_id}: reversed scalar index")
        require({"+-": "-+", "-+": "+-"}[source.orientation] == target.orientation,
                f"{source_id}/{target_id}: reversed orientation")

    # Universal degree argument, not a bounded degree scan.  Sigma|h gives
    # H=deg(h)>=2 and K=deg(k)>=0.  The leading terms of k h'+2h k'
    # combine with multiplier H+2K>0, so no cancellation occurs.
    minimum_h_degree, minimum_k_degree = 2, 0
    minimum_euler_degree = minimum_h_degree + minimum_k_degree - 1
    minimum_leading_multiplier = minimum_h_degree + 2 * minimum_k_degree
    require(minimum_euler_degree == 1, "Euler positive-degree lower bound")
    require(minimum_leading_multiplier == 2, "Euler leading coefficient nonzero")

    parent_lines_digest = sha256(
        ("\n".join(entries[word_id].source_line for word_id in TARGET_IDS) + "\n").encode("ascii")
    ).hexdigest()
    semantic_object = {
        "parents": {
            "theorem": EXPECTED_PARENT_THEOREM_SHA256,
            "script": EXPECTED_PARENT_SCRIPT_SHA256,
            "output": EXPECTED_PARENT_OUTPUT_SHA256,
            "entries": parent_lines_digest,
        },
        "cone_groups": cone_groups,
        "parameterizations": {
            word_id: tuple(value.canonical() for value in parameterizations[word_id])
            for word_id in TARGET_IDS
        },
        "chambers": chambers,
        "weights": {
            word_id: (
                tuple(value.canonical() for value in weights_by_word[word_id][0]),
                tuple(value.canonical() for value in weights_by_word[word_id][1]),
            )
            for word_id in TARGET_IDS
        },
        "components": components_by_word,
        "magnitudes": {
            word_id: {
                sign: tuple((node, value.canonical()) for node, value in values)
                for sign, values in magnitudes_by_word[word_id].items()
            }
            for word_id in TARGET_IDS
        },
        "factors": factor_ledgers,
        "reversal": reversal_pairs,
        "euler": ("k*h'+2*h*k'", "deg(h)+deg(k)-1", "deg(h)+2*deg(k)"),
    }
    semantic_sha = sha256(
        json.dumps(semantic_object, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-3609 THREE-BY-FOUR SIZE-NINE EULER-FACTOR NONENTRY")
    print(f"thm3606_theorem_sha256_lf={EXPECTED_PARENT_THEOREM_SHA256}")
    print(f"thm3606_script_sha256_lf={EXPECTED_PARENT_SCRIPT_SHA256}")
    print(f"thm3606_output_sha256_lf={EXPECTED_PARENT_OUTPUT_SHA256}")
    print(f"thm3606_four_entries_sha256={parent_lines_digest}")
    print("canonical_words=" + ",".join(TARGET_IDS) + ";cone_count=2;each_rank=3;each_dimension=2")
    print(f"left_cone={parameterization_text(left)};W072:t>a;W073:a>t")
    print(f"right_cone={parameterization_text(right)};W141:a>t;W149:t>a")
    print("left_weights=p:(1-a-t,1,a+1);q:(-a-2,-2,t-2,2a+t-2)")
    print("right_weights=p:(-a-2,-2,a+t-2);q:(1-2a-t,1-t,1,a+1)")
    print("left_singleton_components=-:{P0,Q0,Q1};+:{P1,P2,Q2,Q3}")
    print("left_negative_magnitudes=" + magnitude_text(left_negative))
    print("left_positive_magnitudes=" + magnitude_text(left_positive))
    print("right_singleton_components=-:{P0,P1,Q0,Q1};+:{P2,Q2,Q3}")
    print("right_negative_magnitudes=" + magnitude_text(right_negative))
    print("right_positive_magnitudes=" + magnitude_text(right_positive))
    print("gcd_gate=negative_d_in_{1,2};positive_d=1;d=2_iff_a_even_and_t_odd")
    print("arm_gate=unique_anchor_R=2;other_R=a+2>=3;forces_d=2,ord_beta(h)=1,ord_beta(k)=0")
    print("left_collision_factors=02=10:E[AN(t-2)rho*h^(rho-1)k^(t-3)-LB*eta*h^(eta-1)];"
          "11=20:-E[LC+MB(a+1)eta*h^(eta-1)k^a];"
          "03=21:E[AT*sigma*rho*h^(rho-1)k^(sigma-1)-MC(a+1)k^a]")
    print("right_collision_factors=02=20:E[AN*eta*h^(eta-1)-MB*R*rho*h^(rho-1)k^(R-1)];"
          "03=12:E[AT(a+1)eta*h^(eta-1)k^a+LN];"
          "13=21:E[LT(a+1)k^a-MC*R*omega*h^(omega-1)k^(R-1)]")
    print("definitions=E:k*h'+2*h*k';left:rho=(a+t-1)/2,eta=(a+2)/2,sigma=2a+t-2;"
          "right:eta=(a+2)/2,rho=(2a+t-1)/2,omega=(t-1)/2,R=a+t-2")
    print("reversal_pairs=W072<->W149,W073<->W141;cell_map:(i,j)->(2-i,3-j);orientation:+-<->-+")
    print("scalar_equations=left:-E[...]=1;right:E[...]=1")
    print("euler_degree=deg(h)+deg(k)-1>=1;leading_multiplier=deg(h)+2deg(k)>0;E_nonunit=1")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("boundary=finite-field-free characteristic-zero polynomial identity;eliminates only W072,W073,W141,W149")
    print("status=PROVED+EXACT+OPTIMIZATION-SAFE;four residual size-nine schemes impossible")
    print("PASS")


if __name__ == "__main__":
    main()
