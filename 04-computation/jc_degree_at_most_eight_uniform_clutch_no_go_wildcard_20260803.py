#!/usr/bin/env python3
"""Exact probe for all polynomial B-clutches of degree at most eight.

For each of the two THM-3212 accessory response pairs, this script works in
the exact cubic field K and in the rational-function field K(c), where c is
the nonzero value of B at the simple S-root.  It recursively tunes the eight
remaining Taylor coefficients of B to maximize the S-contact of the
THM-3225 critical resultant.  The first two untunable coefficients have
coprime numerators, so the saturated residual can have S-order at most nine.

This is a critical-point obstruction for the displayed first-coordinate
family.  It is not a polynomial inverse cover and proves no instance of the
planar Jacobian conjecture.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
N = 14
EXPECTED_DIGESTS = (
    "6241cc936e35e6a2999f0ea2e1d251cd4b5f313d0a53331723aa1397e8a47a05",
    "bace7cc5c6c0123742a242cbbddc147595c87f307507b5f67d9dcec6e97951e4",
)
DEPENDENCIES = {
    "04-computation/jc_heptic_affine_B_source_obstruction_thm3225.py":
        "0ada8f35c1523ee802c52ca5e090f990909b829c7ebb767ac1b1f744607ad631",
    "05-knowledge/results/jc_heptic_affine_B_source_obstruction_thm3225.out":
        "7237121880bb8e38841f969dce420595a05980437eaaae3682395b7a6b854652",
    "04-computation/jc_heptic_degree_nine_infinity_wall_thm3237.py":
        "92518f258afeca233e90790fa2f713fcfd375c295271ff85dc4c5c66c0057d81",
    "05-knowledge/results/jc_heptic_degree_nine_infinity_wall_thm3237.out":
        "38599c5d2a9b527098274d0dfccd427b5f091528671df168ca3a7ffb31c3b9cb",
    "04-computation/jc_degree8_tuned_cubic_infinity_wall_thm3257.py":
        "f7c9b6d92204ab0af0271311bfbcfa11fb51014a1e3d875d271ff406898da06d",
    "05-knowledge/results/jc_degree8_tuned_cubic_infinity_wall_thm3257.out":
        "1ade6b505bf92f7cb1a17395a3fcff22ee9ac9692ea1ac449b2bee812831ab1b",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


@dataclass(frozen=True)
class Laurent:
    """A Laurent polynomial in c over one exact cubic field."""

    base: object
    terms: tuple

    @staticmethod
    def make(base, entries) -> "Laurent":
        cleaned = {
            int(exponent): coefficient
            for exponent, coefficient in dict(entries).items()
            if coefficient != base.zero
        }
        return Laurent(base, tuple(sorted(cleaned.items())))

    @staticmethod
    def constant(base, value) -> "Laurent":
        coefficient = value if type(value) is type(base.one) else base.convert(value)
        return Laurent.make(base, {} if coefficient == base.zero else {0: coefficient})

    def as_dict(self):
        return dict(self.terms)

    def __bool__(self):
        return bool(self.terms)

    def _coerce(self, other):
        if isinstance(other, Laurent):
            require(other.base == self.base, "Laurent base mismatch")
            return other
        return Laurent.constant(self.base, other)

    def __eq__(self, other):
        try:
            other = self._coerce(other)
        except Exception:
            return False
        return self.terms == other.terms

    def __add__(self, other):
        other = self._coerce(other)
        entries = self.as_dict()
        for exponent, coefficient in other.terms:
            entries[exponent] = entries.get(exponent, self.base.zero) + coefficient
        return Laurent.make(self.base, entries)

    __radd__ = __add__

    def __neg__(self):
        return Laurent.make(
            self.base, {exponent: -coefficient for exponent, coefficient in self.terms}
        )

    def __sub__(self, other):
        return self + (-self._coerce(other))

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        other = self._coerce(other)
        entries = {}
        for exponent_left, coefficient_left in self.terms:
            for exponent_right, coefficient_right in other.terms:
                exponent = exponent_left + exponent_right
                entries[exponent] = entries.get(exponent, self.base.zero) + (
                    coefficient_left * coefficient_right
                )
        return Laurent.make(self.base, entries)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        require(exponent >= 0, "negative Laurent exponent")
        answer = Laurent.constant(self.base, 1)
        power = self
        while exponent:
            if exponent & 1:
                answer = answer * power
            power = power * power
            exponent //= 2
        return answer

    def __truediv__(self, other):
        other = self._coerce(other)
        require(len(other.terms) == 1, "Laurent division requires a monomial")
        exponent, coefficient = other.terms[0]
        inverse = coefficient**-1
        return Laurent.make(
            self.base,
            {
                source_exponent - exponent: source_coefficient * inverse
                for source_exponent, source_coefficient in self.terms
            },
        )

    def rational_degrees(self):
        require(self.terms, "zero Laurent degree")
        minimum = self.terms[0][0]
        maximum = self.terms[-1][0]
        denominator_degree = max(0, -minimum)
        return maximum + denominator_degree, denominator_degree

    def shifted_polynomial(self, symbol):
        require(self.terms, "zero Laurent polynomial")
        shift = max(0, -self.terms[0][0])
        ring = self.base.poly_ring(symbol)
        C = ring.gens[0]
        answer = ring.zero
        for exponent, coefficient in self.terms:
            answer += coefficient * C ** (exponent + shift)
        return answer, shift


class LaurentDomain:
    def __init__(self, base):
        self.base = base
        self.zero = Laurent.constant(base, 0)
        self.one = Laurent.constant(base, 1)

    def lift(self, value):
        return Laurent.constant(self.base, value)


def convolution(left, right, field):
    answer = [field.zero] * N
    for i, a in enumerate(left[:N]):
        if not a:
            continue
        for j, b in enumerate(right[: N - i]):
            if b:
                answer[i + j] += a * b
    return answer


def series_add(*summands, field):
    return [
        sum(
            (summand[i] if i < len(summand) else field.zero for summand in summands),
            field.zero,
        )
        for i in range(N)
    ]


def series_scale(series, scalar, field):
    return [scalar * entry for entry in series[:N]] + [field.zero] * max(
        0, N - len(series)
    )


def series_power(series, exponent, field):
    answer = [field.one] + [field.zero] * (N - 1)
    base = series[:N] + [field.zero] * max(0, N - len(series))
    while exponent:
        if exponent & 1:
            answer = convolution(answer, base, field)
        base = convolution(base, base, field)
        exponent //= 2
    return answer


def series_derivative(series, field):
    return [(i + 1) * series[i + 1] for i in range(len(series) - 1)] + [
        field.zero
    ]


def series_product(*factors, field):
    answer = [field.one] + [field.zero] * (N - 1)
    for factor in factors:
        answer = convolution(answer, factor, field)
    return answer


def critical_resultant_series(V, A, B, field):
    J = series_add(
        series_scale(
            convolution(V, series_derivative(B, field), field), 2, field
        ),
        series_scale(
            convolution(B, series_derivative(V, field), field), -1, field
        ),
        field=field,
    )

    def term(coefficient, *factors):
        return series_scale(
            series_product(*factors, field=field), coefficient, field
        )

    return series_add(
        term(-1, series_power(A, 3, field), series_power(J, 3, field)),
        term(
            12,
            series_power(A, 2, field),
            series_power(J, 2, field),
            series_power(V, 2, field),
        ),
        term(
            -4,
            A,
            series_power(B, 3, field),
            series_power(J, 2, field),
            V,
        ),
        term(
            4,
            A,
            series_power(B, 2, field),
            J,
            series_power(V, 2, field),
        ),
        term(24, A, B, J, series_power(V, 3, field)),
        term(-48, A, J, series_power(V, 4, field)),
        term(-16, A, series_power(V, 4, field)),
        term(
            -8,
            series_power(B, 4, field),
            J,
            series_power(V, 2, field),
        ),
        term(8, series_power(B, 3, field), J, series_power(V, 3, field)),
        term(32, series_power(B, 2, field), series_power(V, 4, field)),
        term(-96, B, series_power(V, 5, field)),
        term(64, series_power(V, 6, field)),
        field=field,
    )


def shifted_coefficients(poly, root, field, symbol):
    shifted_ring = field.poly_ring(symbol)
    t = shifted_ring.gens[0]
    shifted = shifted_ring.zero
    for (exponent,), coefficient in poly.to_dict().items():
        shifted += coefficient * (t + root) ** exponent
    coefficient_dict = shifted.to_dict()
    return [coefficient_dict.get((i,), field.zero) for i in range(N)]


def canonical_poly_text(poly) -> str:
    rows = []
    for (degree,), coefficient in sorted(poly.to_dict().items()):
        rows.append(f"{degree}:{coefficient}")
    return "|".join(rows)


def build_case(name: str):
    u, x, c, t = sp.symbols("u x c t")
    if name == "4111":
        accessory = sp.Poly(100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ)
        exponent_a, exponent_b = 4, 1
    else:
        accessory = sp.Poly(75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ)
        exponent_a, exponent_b = 3, 2

    cubic_field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = cubic_field.ext
    x_ring = cubic_field.poly_ring(x)
    X = x_ring.gens[0]

    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
        extras = (9, 0)
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)

    gamma = -7 * A0
    q2 = X**2 - alpha * X + accessory_v
    D = X**exponent_a * (X - 1) ** exponent_b * q2
    T = X * (X - 1) * q2
    E = (
        exponent_a * (X - 1) * q2
        + exponent_b * X * q2
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D * T**2 / gamma**2
    A = 2 * S * E * T / gamma

    require(V.degree() == 16 and A.degree() == 8, f"{name} response degrees")
    require(
        2 * V * A.diff(X) - A * V.diff(X) == 2 * V,
        f"{name} response identity",
    )

    boundary = S**3 * T**8 * X ** extras[0] * (X - 1) ** extras[1]
    require(boundary.degree() == 44, f"{name} boundary degree")

    root = -shift
    V_shift = shifted_coefficients(V, root, cubic_field, t)
    A_shift = shifted_coefficients(A, root, cubic_field, t)
    require(
        V_shift[0] == cubic_field.zero and V_shift[1] != cubic_field.zero,
        f"{name} simple S root",
    )
    require(
        A_shift[0] == cubic_field.zero and A_shift[1] == cubic_field.convert(2),
        f"{name} normalized response",
    )

    laurent_field = LaurentDomain(cubic_field)
    C = Laurent.make(cubic_field, {1: cubic_field.one})
    V_series = [laurent_field.lift(value) for value in V_shift]
    A_series = [laurent_field.lift(value) for value in A_shift]
    B_series = [C] + [laurent_field.zero] * (N - 1)

    slopes = []
    tuned_coefficients = []
    for jet in range(1, 9):
        target_degree = jet + 2
        B_series[jet] = laurent_field.zero
        constant = critical_resultant_series(
            V_series, A_series, B_series, laurent_field
        )[target_degree]
        B_series[jet] = laurent_field.one
        unit_value = critical_resultant_series(
            V_series, A_series, B_series, laurent_field
        )[target_degree]
        slope = unit_value - constant
        B_series[jet] = 2 * laurent_field.one
        two_value = critical_resultant_series(
            V_series, A_series, B_series, laurent_field
        )[target_degree]
        expected_slope = 16 * jet * C**4 * V_series[1] ** 3
        require(
            slope == expected_slope,
            f"{name} jet {jet} exact triangular slope",
        )
        require(
            two_value - constant == 2 * slope,
            f"{name} jet {jet} triangular affine dependence",
        )
        require(slope != 0, f"{name} jet {jet} nonzero slope")
        require(
            len(slope.terms) == 1,
            f"{name} jet {jet} slope supported only at c!=0",
        )
        B_series[jet] = -constant / slope
        slopes.append(slope)
        tuned_coefficients.append(B_series[jet])

    K_series = critical_resultant_series(
        V_series, A_series, B_series, laurent_field
    )
    require(all(K_series[degree] == 0 for degree in range(3, 11)),
            f"{name} tuned contact through degree ten")
    F11 = K_series[11]
    F12 = K_series[12]
    require(F11 != 0 and F12 != 0, f"{name} untunable pair nonzero")
    F11_poly, F11_denominator_degree = F11.shifted_polynomial(c)
    F12_poly, F12_denominator_degree = F12.shifted_polynomial(c)
    gcd_poly = F11_poly.gcd(F12_poly)
    expected_jet_degrees = (
        (1, 0),
        (3, 2),
        (4, 3),
        (6, 5),
        (7, 6),
        (9, 8),
        (10, 9),
        (12, 11),
    )
    require(
        tuple(value.rational_degrees() for value in tuned_coefficients)
        == expected_jet_degrees,
        f"{name} tuned jet degree profile",
    )
    require(
        tuple(value.rational_degrees() for value in slopes) == ((4, 0),) * 8,
        f"{name} slope degree profile",
    )
    require(
        (F11_poly.degree(), F11_denominator_degree) == (13, 8)
        and (F12_poly.degree(), F12_denominator_degree) == (15, 10),
        f"{name} untunable degree profile",
    )
    require(gcd_poly.degree() == 0, f"{name} coprime untunable pair")

    digest_text = canonical_poly_text(F11_poly.monic()) + "\n" + canonical_poly_text(
        F12_poly.monic()
    )
    digest = sha256(digest_text.encode("ascii")).hexdigest()

    print(
        f"case={name} degrees=(V,A,boundary)=(16,8,44) "
        f"S_response=(ordV=1,A1=2)"
    )
    print(
        f"case={name} tuned_B_jets=8 "
        f"jet_degrees={tuple(value.rational_degrees() for value in tuned_coefficients)}"
    )
    print(
        f"case={name} slope_degrees="
        f"{tuple(value.rational_degrees() for value in slopes)}"
    )
    print(
        f"case={name} untunable="
        f"(F11_num_degree={F11_poly.degree()},F11_den_degree={F11_denominator_degree},"
        f"F12_num_degree={F12_poly.degree()},F12_den_degree={F12_denominator_degree},"
        f"gcd_degree={gcd_poly.degree()}) digest={digest}"
    )
    print(
        f"case={name} consequence=ord_S(K)<=12;ord_S(H)<=9;"
        "off_boundary_resultant_multiplicity>=43"
    )

    return digest


def universal_degree_ledger() -> None:
    ledger = []
    for degree_B in range(9):
        degree_J = degree_B + 15 if degree_B < 8 else 22
        bounds = (
            24 + 3 * degree_J,
            48 + 2 * degree_J,
            24 + 3 * degree_B + 2 * degree_J,
            40 + 2 * degree_B + degree_J,
            56 + degree_B + degree_J,
            72 + degree_J,
            72,
            32 + 4 * degree_B + degree_J,
            48 + 3 * degree_B + degree_J,
            64 + 2 * degree_B,
            80 + degree_B,
            96,
        )
        ledger.append(max(bounds))
    require(tuple(ledger) == (96,) * 9, "degree-at-most-eight ledger")
    print(f"universal_degree_ledger_d=0..8={tuple(ledger)};unique_top_term=64*V^6")


def universal_boundary_checks() -> None:
    t, v, b, e = sp.symbols("t v b e", nonzero=True)

    def critical_polynomial(V, A, B):
        J = 2 * V * sp.diff(B, t) - B * sp.diff(V, t)
        return sp.expand(
            -A**3 * J**3
            + 12 * A**2 * J**2 * V**2
            - 4 * A * B**3 * J**2 * V
            + 4 * A * B**2 * J * V**2
            + 24 * A * B * J * V**3
            - 48 * A * J * V**4
            - 16 * A * V**4
            - 8 * B**4 * J * V**2
            + 8 * B**3 * J * V**3
            + 32 * B**2 * V**4
            - 96 * B * V**5
            + 64 * V**6
        )

    simple_K = sp.Poly(critical_polynomial(v * t, 2 * t, b + e * t), t)
    require(
        all(simple_K.nth(index) == 0 for index in range(3)),
        "simple S boundary has order at least three",
    )

    leading_rows = []
    for multiplicity in (3, 4, 5, 6):
        A_lead = sp.Rational(2, 2 - multiplicity) * t
        K_lead = sp.Poly(
            critical_polynomial(v * t**multiplicity, A_lead, b + e * t), t
        )
        coefficient = sp.factor(K_lead.nth(3 * multiplicity - 1))
        expected = (
            sp.Rational(16 * multiplicity * (multiplicity - 1), multiplicity - 2)
            * b**5
            * v**3
        )
        require(coefficient == expected, f"T boundary multiplicity {multiplicity}")
        leading_rows.append((multiplicity, 3 * multiplicity - 1))
    print(
        "universal_boundary="
        f"S_order>=3;T_rows={tuple(leading_rows)};"
        "T_lead=16*m*(m-1)/(m-2)*b^5*v^3"
    )


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")


def main() -> None:
    print("degree-at-most-eight uniform polynomial clutch no-go wildcard")
    for dependency, expected_hash in DEPENDENCIES.items():
        require(lf_hash(ROOT / dependency) == expected_hash,
                f"dependency drift: {dependency}")
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    universal_degree_ledger()
    universal_boundary_checks()
    digests = tuple(build_case(name) for name in ("4111", "3211"))
    require(digests == EXPECTED_DIGESTS, "exact case digest drift")
    print(f"case_digests={digests}")
    print("scope=critical-resultant-first-coordinate-family-not-JC2")
    source_audit()


if __name__ == "__main__":
    main()
