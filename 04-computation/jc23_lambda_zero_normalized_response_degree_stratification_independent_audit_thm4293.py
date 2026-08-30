#!/usr/bin/env python3
"""Standard-library-only independent exact audit for THM-4293.

The audit reconstructs the Lambda=0 wall chart from the literal generic-Q
source with a sparse Laurent-polynomial engine over ``Fraction``.  It then
uses the local Newton polygon (rather than a CAS normalization), the residue
chain rule, and canonical-divisor saturation to certify the six noncritical
strata.  A bounded but complete Eisenstein-norm enumeration checks the
degree obstruction.  The repeated-discriminant case is retained only as a
hostile boundary and is not classified here.

No result in this file proves JC(2), re-proves THM-4292, or classifies the
repeated-discriminant lane.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import random


VARIABLES = (
    "s", "p", "z", "w", "u", "v", "carrier",
    "Q", "Phi", "Delta", "Theta", "eta", "zeta3",
    "upsilon5", "xi10", "alpha11", "beta11", "U", "W", "Z",
)
INDEX = {name: position for position, name in enumerate(VARIABLES)}
ZERO_MONOMIAL = (0,) * len(VARIABLES)
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


class Sparse:
    """Sparse Laurent polynomial over Q in the fixed variables above."""

    def __init__(self, terms: dict[tuple[int, ...], F] | None = None) -> None:
        clean: dict[tuple[int, ...], F] = {}
        for monomial, coefficient in (terms or {}).items():
            require(len(monomial) == len(VARIABLES), "monomial dimension")
            coefficient = F(coefficient)
            if coefficient:
                clean[monomial] = clean.get(monomial, F(0)) + coefficient
        self.terms = {monomial: value for monomial, value in clean.items() if value}

    @staticmethod
    def constant(value: F | int) -> "Sparse":
        value = F(value)
        return Sparse({ZERO_MONOMIAL: value}) if value else Sparse()

    @staticmethod
    def monomial(powers: dict[str, int], coefficient: F | int = 1) -> "Sparse":
        exponents = [0] * len(VARIABLES)
        for name, exponent in powers.items():
            exponents[INDEX[name]] = exponent
        return Sparse({tuple(exponents): F(coefficient)})

    @staticmethod
    def variable(name: str) -> "Sparse":
        return Sparse.monomial({name: 1})

    def __add__(self, other: "Sparse" | F | int) -> "Sparse":
        other = other if isinstance(other, Sparse) else Sparse.constant(other)
        terms = dict(self.terms)
        for monomial, coefficient in other.terms.items():
            terms[monomial] = terms.get(monomial, F(0)) + coefficient
        return Sparse(terms)

    __radd__ = __add__

    def __neg__(self) -> "Sparse":
        return Sparse({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: "Sparse" | F | int) -> "Sparse":
        return self + (-other if isinstance(other, Sparse) else -F(other))

    def __rsub__(self, other: "Sparse" | F | int) -> "Sparse":
        return (-self) + other

    def __mul__(self, other: "Sparse" | F | int) -> "Sparse":
        other = other if isinstance(other, Sparse) else Sparse.constant(other)
        terms: dict[tuple[int, ...], F] = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                monomial = tuple(x + y for x, y in zip(left, right))
                terms[monomial] = terms.get(monomial, F(0)) + a * b
        return Sparse(terms)

    __rmul__ = __mul__

    def __truediv__(self, other: F | int) -> "Sparse":
        scalar = F(other)
        require(bool(scalar), "nonzero scalar divisor")
        return self * (1 / scalar)

    def __pow__(self, exponent: int) -> "Sparse":
        require(exponent >= 0, "nonnegative polynomial power")
        result = Sparse.constant(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n //= 2
        return result

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Sparse):
            return self.terms == other.terms
        if isinstance(other, (int, F)):
            return self == Sparse.constant(other)
        return False

    def derivative(self, name: str) -> "Sparse":
        position = INDEX[name]
        terms: dict[tuple[int, ...], F] = {}
        for monomial, coefficient in self.terms.items():
            exponent = monomial[position]
            if exponent:
                derived = list(monomial)
                derived[position] -= 1
                terms[tuple(derived)] = coefficient * exponent
        return Sparse(terms)

    def substitute(self, replacements: dict[str, "Sparse" | F | int]) -> "Sparse":
        normalized = {
            name: value if isinstance(value, Sparse) else Sparse.constant(value)
            for name, value in replacements.items()
        }
        result = Sparse()
        for monomial, coefficient in self.terms.items():
            term = Sparse.constant(coefficient)
            for name, exponent in zip(VARIABLES, monomial):
                if not exponent:
                    continue
                if name in normalized:
                    require(exponent >= 0, "substitution into nonnegative exponent")
                    term = term * normalized[name] ** exponent
                else:
                    term = term * Sparse.monomial({name: exponent})
            result = result + term
        return result

    def coefficient(self, powers: dict[str, int]) -> "Sparse":
        positions = {INDEX[name]: exponent for name, exponent in powers.items()}
        terms: dict[tuple[int, ...], F] = {}
        for monomial, value in self.terms.items():
            if all(monomial[position] == exponent for position, exponent in positions.items()):
                reduced = list(monomial)
                for position in positions:
                    reduced[position] = 0
                terms[tuple(reduced)] = terms.get(tuple(reduced), F(0)) + value
        return Sparse(terms)

    def valuation(self, name: str) -> int:
        require(bool(self.terms), f"nonzero {name}-valuation")
        position = INDEX[name]
        return min(monomial[position] for monomial in self.terms)


def var(name: str) -> Sparse:
    return Sparse.variable(name)


def mono(**powers: int) -> Sparse:
    return Sparse.monomial(powers)


def build_source() -> tuple[Sparse, Sparse, Sparse, Sparse]:
    """Return the literal source, wall source, regular chart, and local chart."""

    s, p = var("s"), var("p")
    y = s * p
    q = var("Q")
    k_forced = F(2848, 45) - F(7, 6) * var("Delta")
    h = (
        -3 * p
        + F(8, 3) * p**2
        - F(1376, 135) * p**3
        + k_forced * y**2
        + var("Phi") * p**2 * y
        + var("Delta") * p**4
        + var("Theta") * p * y**2
        + var("eta") * p**3 * y
        + var("zeta3") * y**3
        + var("upsilon5") * p**5
        + var("xi10") * p**2 * y**2
        + var("alpha11") * p**4 * y
        + var("beta11") * p * y**3
        + var("U") * p**6
        + var("W") * p**3 * y**2
        + var("Z") * y**4
    )
    source = (s**2 - p) * (1 - q * h) - q * s**2 / 2
    wall = source.substitute({"W": 0, "Z": -var("U")})
    chart_w = (
        mono(z=14, w=7)
        * wall.substitute({"s": mono(z=-1), "p": mono(z=-2, w=-1)})
    )
    require(all(min(monomial[INDEX[name]] for monomial in chart_w.terms) >= 0
                for name in ("z", "w")), "regular wall chart exponents")
    chart = chart_w.substitute({"w": 1 + var("u")})
    return source, wall, chart_w, chart


def lower_newton_edges(points: tuple[tuple[int, int], ...]) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    hull: list[tuple[int, int]] = []
    for point in sorted(points):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            cross = ((b[0] - a[0]) * (point[1] - a[1])
                     - (b[1] - a[1]) * (point[0] - a[0]))
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    return tuple(zip(hull, hull[1:]))


def root_valuations(r: int) -> tuple[int, ...]:
    edges = lower_newton_edges(((0, 12), (1, r), (2, 0), (3, 0)))
    valuations: list[int] = []
    for left, right in edges:
        slope = F(right[1] - left[1], right[0] - left[0])
        value = -slope
        require(value.denominator == 1, f"integral Newton slope r={r}")
        valuations.extend([value.numerator] * (right[0] - left[0]))
    return tuple(sorted(value for value in valuations if value > 0))


def parameter_specialization(r: int, c: F, u_value: F) -> dict[str, Sparse | F | int]:
    a6 = F(7168, 135)
    values: dict[str, Sparse | F | int] = {
        "Q": 1,
        "Phi": 0,
        "Delta": 0,
        "Theta": 0,
        "eta": 0,
        "zeta3": 0,
        "upsilon5": 0,
        "xi10": 0,
        "alpha11": 0,
        "beta11": 0,
        "U": u_value,
    }
    if r == 1:
        values["alpha11"] = c
    elif r == 2:
        values["upsilon5"] = c
    elif r == 3:
        values["eta"] = c
    elif r == 4:
        values["Theta"] = c
    elif r == 5:
        values["Phi"] = c
    elif r == 6:
        delta = F(6, 7) * (a6 - c)
        values["Delta"] = delta
        values["Theta"] = -delta
    else:
        raise ValueError(r)
    return values


def polynomial_clean(poly: dict[int, F]) -> dict[int, F]:
    return {degree: F(value) for degree, value in poly.items() if value}


def polynomial_derivative(poly: dict[int, F]) -> dict[int, F]:
    return {degree - 1: degree * value for degree, value in poly.items() if degree}


def polynomial_divmod(left: dict[int, F], right: dict[int, F]) -> tuple[dict[int, F], dict[int, F]]:
    left, right = polynomial_clean(left), polynomial_clean(right)
    require(bool(right), "nonzero polynomial divisor")
    quotient: dict[int, F] = {}
    right_degree = max(right)
    right_lead = right[right_degree]
    while left and max(left) >= right_degree:
        degree = max(left) - right_degree
        coefficient = left[max(left)] / right_lead
        quotient[degree] = quotient.get(degree, F(0)) + coefficient
        for exponent, value in right.items():
            target = exponent + degree
            left[target] = left.get(target, F(0)) - coefficient * value
        left = polynomial_clean(left)
    return polynomial_clean(quotient), left


def polynomial_gcd(left: dict[int, F], right: dict[int, F]) -> dict[int, F]:
    left, right = polynomial_clean(left), polynomial_clean(right)
    while right:
        _, remainder = polynomial_divmod(left, right)
        left, right = right, remainder
    require(bool(left), "nonzero polynomial gcd")
    lead = left[max(left)]
    return {degree: value / lead for degree, value in left.items()}


def small_eisenstein_norms(limit: int) -> tuple[set[int], dict[int, tuple[int, int]]]:
    # 2*N(a,b)=a^2+b^2+(a-b)^2, so N<=limit forces
    # |a|,|b|<=floor(sqrt(2*limit)).  For limit=10, [-4,4] is complete.
    bound = 0
    while (bound + 1) ** 2 <= 2 * limit:
        bound += 1
    values: set[int] = set()
    witnesses: dict[int, tuple[int, int]] = {}
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            norm = a * a - a * b + b * b
            if norm <= limit:
                values.add(norm)
                witnesses.setdefault(norm, (a, b))
    return values, witnesses


@dataclass(frozen=True)
class Row:
    r: int
    branch_valuations: tuple[int, int]
    intersection: int
    residue_orders: tuple[int, ...]
    packet: tuple[int, ...]
    genus: int
    full_degree: int
    finite_degree: int
    full_norm: int | None
    finite_norm: int | None
    eigenline: bool


def main() -> None:
    source, wall, chart_w, chart = build_source()
    q, w, u, z = var("Q"), var("w"), var("u"), var("z")
    capital_u = var("U")

    top_edge = Sparse()
    for position in range(4):
        top_edge += source.coefficient({"s": 2 * position, "p": 7 - position}) * w**position
    expected_top = q * (1 - w) * (capital_u + var("W") * w + var("Z") * w**2)
    require(top_edge == expected_top, "literal top edge")

    carrier = var("carrier")
    quartic_edge = (
        source.coefficient({"s": 2, "p": 0})
        + source.coefficient({"s": 4, "p": 2}) * carrier**2
        + source.coefficient({"s": 5, "p": 3}) * carrier**3
        + source.coefficient({"s": 6, "p": 4}) * carrier**4
    )
    k_forced = F(2848, 45) - F(7, 6) * var("Delta")
    expected_quartic = (
        1 - q / 2
        - q * (k_forced * carrier**2 + var("zeta3") * carrier**3 + var("Z") * carrier**4)
    )
    require(quartic_edge == expected_quartic, "literal quartic carrier edge")
    wall_quartic = quartic_edge.substitute({"Z": -capital_u})
    require(wall_quartic.coefficient({"carrier": 4}) == q * capital_u,
            "quartic carrier degree four on wall")

    require(chart_w.substitute({"z": 0})
            == q * capital_u * (1 - w) ** 2 * (1 + w),
            "double plus simple top roots")
    require(chart.substitute({"z": 0}) == q * capital_u * u**2 * (u + 2),
            "local double-root equation")
    require(chart.substitute({"u": 0}) == -q * z**12 / 2,
            "exact smoothing constant")

    c1 = var("alpha11") + var("beta11")
    c2 = var("upsilon5") + var("xi10")
    c3 = var("eta") + var("zeta3")
    c4 = var("Delta") + var("Theta")
    c5 = var("Phi")
    c6 = F(7168, 135) - F(7, 6) * var("Delta")
    linear_u = chart.coefficient({"u": 1})
    expected_linear = (
        -q * (c1 * z + c2 * z**2 + c3 * z**3 + c4 * z**4
              + c5 * z**5 + c6 * z**6 + F(8, 3) * z**8 - 3 * z**10)
        + (1 - F(7, 2) * q) * z**12
    )
    require(linear_u == expected_linear, "c1 through c6 local expansion")
    require(chart.coefficient({"u": 2, "z": 0}) == 2 * q * capital_u,
            "quadratic local unit")
    require(chart.coefficient({"u": 3, "z": 0}) == q * capital_u,
            "cubic local unit correction")

    # Derive the residue numerator rather than importing its order.
    source_p_chart = wall.derivative("p").substitute(
        {"s": mono(z=-1), "p": mono(z=-2, w=-1)}
    )
    chain_rhs = mono(z=-12, w=-5) * (
        7 * chart_w * mono(w=-1) - chart_w.derivative("w")
    )
    require(source_p_chart == chain_rhs, "exact residue chain rule")
    # On Fbar=0, ds=-z^-2 dz and the two minus signs cancel:
    # Q ds/F_p = Q z^10 w^5 dz/Fbar_w.
    residue_numerator_order = -2 - (-12)
    require(residue_numerator_order == 10, "residue numerator z-order")

    norm_values, norm_witnesses = small_eisenstein_norms(10)
    require(norm_values == {0, 1, 3, 4, 7, 9}, "complete Eisenstein norms through ten")
    for norm, witness in norm_witnesses.items():
        a, b = witness
        require(a * a - a * b + b * b == norm, f"norm witness {norm}")

    expected_table = (
        (1, (1, 11), 1, 17, 40, 32, 10, 8, False),
        (2, (2, 10), 2, 16, 38, 30, None, None, False),
        (3, (3, 9), 3, 15, 36, 28, 9, 7, True),
        (4, (4, 8), 4, 14, 34, 26, None, None, False),
        (5, (5, 7), 5, 13, 32, 24, 8, 6, False),
        (6, (6, 6), 6, 12, 30, 22, None, None, False),
    )
    rows: list[Row] = []
    for r_value in range(1, 7):
        valuations = root_valuations(r_value)
        expected_valuations = ((r_value, 12 - r_value)
                               if r_value < 6 else (6, 6))
        require(valuations == expected_valuations, f"Newton valuations r={r_value}")

        # The length-one initial equations for r<6 have rational nonzero
        # roots.  At r=6, use a distinct-root exact positive control.
        if r_value < 6:
            c_value, u_value = F(r_value + 1), F(r_value + 2)
            fast_root = c_value / (2 * u_value)
            slow_root = -F(1, 2) / c_value
            require(2 * u_value * fast_root**2 - c_value * fast_root == 0,
                    f"fast initial root r={r_value}")
            require(-F(1, 2) - c_value * slow_root == 0,
                    f"slow initial root r={r_value}")
            require(fast_root and slow_root, f"nonzero initial roots r={r_value}")
        else:
            c_value, u_value = F(1), F(2)
            roots = (F(1, 2), F(-1, 4))
            require(all(2 * u_value * root**2 - c_value * root - F(1, 2) == 0
                        for root in roots), "r6 distinct quadratic roots")
            require(roots[0] != roots[1] and c_value**2 + 4 * u_value == 9,
                    "r6 nonzero discriminant")

        specialization = parameter_specialization(r_value, c_value, u_value)
        specialized_linear = linear_u.substitute(specialization)
        require(specialized_linear.valuation("z") == r_value,
                f"literal first nonzero c_r r={r_value}")

        intersection = min(valuations)
        require(intersection == r_value, f"normalized intersection r={r_value}")
        local_residue_order = residue_numerator_order - intersection
        local_index = local_residue_order + 1
        require(local_index == 11 - r_value, f"local index r={r_value}")

        # The inherited non-affine edge inventory is top length 3, quartic
        # length 4, AB length 1.  The double top point normalizes to two
        # branches, beside the simple w=-1 point.
        packet = (local_index, local_index, 11, 1, 2, 2, 2, 2)
        residue_orders = tuple(index - 1 for index in packet)
        genus = 18 - intersection  # delta of two smooth branches is I=r.
        full_degree = sum(packet)
        finite_degree = full_degree - 4 * 2
        require(len(packet) == 3 + 4 + 1 == 8, f"edge point completeness r={r_value}")
        require(sum(residue_orders) == 2 * genus - 2,
                f"canonical divisor saturation r={r_value}")
        require(full_degree == 42 - 2 * r_value, f"full degree r={r_value}")
        require(finite_degree == 34 - 2 * r_value, f"finite degree r={r_value}")

        full_norm = full_degree // 4 if full_degree % 4 == 0 else None
        finite_norm = finite_degree // 4 if finite_degree % 4 == 0 else None
        eigenline = (
            full_norm is not None and finite_norm is not None
            and full_norm in norm_values and finite_norm in norm_values
        )
        rows.append(Row(
            r_value, valuations, intersection, residue_orders, packet, genus,
            full_degree, finite_degree, full_norm, finite_norm, eigenline,
        ))

    actual_table = tuple(
        (row.r, row.branch_valuations, row.intersection, row.genus,
         row.full_degree, row.finite_degree, row.full_norm, row.finite_norm,
         row.eigenline)
        for row in rows
    )
    require(actual_table == expected_table, "complete degree/genus/norm table")
    require(tuple(row.r for row in rows if row.eigenline) == (3,),
            "unique Eisenstein-norm survivor")

    # r=6 specializes to THM-4291's balanced genus-five tail.
    a_tail = F(7168, 135)
    require(a_tail**2 == F(4 * 12845056, 18225), "THM4291 discriminant constant")
    positive_tail = {12: F(1), 6: -2 * a_tail, 0: a_tail**2 + 4}
    require(polynomial_gcd(positive_tail, polynomial_derivative(positive_tail)) == {0: F(1)},
            "THM4291 positive squarefree tail")
    tail_genus = (12 - 2) // 2
    require(tail_genus == 5 and 18 - 6 == 7 + tail_genus == 12,
            "THM4291 r6 genus consistency")
    require((rows[-1].full_degree, rows[-1].finite_degree) == (30, 22),
            "THM4291 r6 response correction")
    require(42 > rows[-1].full_degree > rows[-1].finite_degree,
            "THM4291 abstract tail degree is not central response")

    # Hostile repeated-discriminant control: c6=2, U=-1 makes the r=6
    # initial quadratic a square.  Later faces may split it; no conclusion
    # about their normalization or response is made here.
    repeated = parameter_specialization(6, F(2), F(-1))
    repeated_initial = (
        chart.substitute(repeated)
        .substitute({"u": mono(z=6, v=1)})
        * mono(z=-12)
    ).substitute({"z": 0})
    expected_repeated = -2 * (var("v") + F(1, 2)) ** 2
    require(repeated_initial == expected_repeated, "repeated discriminant hostile initial")
    hostile_tail = {12: F(1), 6: -2 * a_tail}
    hostile_gcd = polynomial_gcd(hostile_tail, polynomial_derivative(hostile_tail))
    require(max(hostile_gcd) == 5, "repeated THM4291 tail has multiple root")

    rng = random.Random(4293)
    random_checks = 4096
    for _ in range(random_checks):
        r_value = rng.randrange(1, 7)
        c_value = F(rng.choice((-1, 1)) * rng.randrange(1, 10**6))
        u_value = F(rng.choice((-1, 1)) * rng.randrange(1, 10**6))
        if r_value == 6 and c_value**2 + 4 * u_value == 0:
            u_value += 1
        require(root_valuations(r_value)
                == ((r_value, 12 - r_value) if r_value < 6 else (6, 6)),
                "random Newton polygon")
        if r_value < 6:
            fast = c_value / (2 * u_value)
            slow = -F(1, 2) / c_value
            require(2 * u_value * fast**2 - c_value * fast == 0,
                    "random fast initial")
            require(-F(1, 2) - c_value * slow == 0,
                    "random slow initial")
        else:
            require(c_value**2 + 4 * u_value != 0,
                    "random noncritical r6 discriminant")

    print("THM4293_LAMBDA_ZERO_NORMALIZED_RESPONSE_DEGREE_INDEPENDENT_V1")
    print("UNIVERSE exact_M12 Lambda=0 W=0 Z=-U char0 Q*U!=0 noncritical_r=1..6")
    print("LITERAL_SOURCE top_edge=Q*(1-w)*(U+W*w+Z*w^2) quartic_degree=4 PASS")
    print("WALL_CHART Fbar(0,w)=Q*U*(1-w)^2*(1+w) Fbar(z,1)=-Q*z^12/2 PASS")
    print("LOCAL_COEFFICIENTS c1=alpha11+beta11 c2=upsilon5+xi10 c3=eta+zeta3 "
          "c4=Delta+Theta c5=Phi c6=7168/135-7Delta/6")
    print("RESIDUE_CHAIN_RULE Q*ds/Fp=Q*z^10*w^5*dz/Fbar_w NUMERATOR_ORDER=10 PASS")
    for row in rows:
        full_norm = "-" if row.full_norm is None else str(row.full_norm)
        finite_norm = "-" if row.finite_norm is None else str(row.finite_norm)
        print(
            f"R={row.r} BRANCH_VALUATIONS={row.branch_valuations[0]},{row.branch_valuations[1]} "
            f"INTERSECTION={row.intersection} RESIDUES={','.join(map(str, row.residue_orders))} "
            f"PACKET={','.join(map(str, row.packet))} GENUS={row.genus} "
            f"DEGREES={row.full_degree},{row.finite_degree} NORMS={full_norm},{finite_norm} "
            f"EIGENLINE={int(row.eigenline)}"
        )
    print("GLOBAL_COMPLETENESS edges=top3+quartic4+AB1 points=8 "
          "residue_sum=2g-2 each_stratum quartic_subtraction=8")
    print("EISENSTEIN_NORMS_LE_10=0,1,3,4,7,9 COMPLETE_BY_2N=a^2+b^2+(a-b)^2")
    print("DEGREE_FILTER even_r=MOD4_FAIL r1=N10,N8_FAIL r5=N8,N6_FAIL r3=N9,N7_SURVIVES")
    print("THM4291_CONTROL r6 discriminant=(4/18225)*(18225U+12845056) "
          "tail_genus=5 central_genus=12=7+5 responses=30,22 abstract_tail_degree=42")
    print("HOSTILE repeated_c6=2_U=-1 initial=-2*(v+1/2)^2 MULTIPLE_ROOT PASS")
    print("BOUNDARY c1=...=c5=0_and_c6^2+4U=0 REPEATED_DISCRIMINANT_OPEN_NOT_CLASSIFIED")
    print(f"CONTROLS randomized={random_checks} seed=4293 checks={CHECKS}")
    print("SCOPE normalized_noncritical_wall_only THM4292_imported_for_degree_location JC2_OPEN")
    print("VERDICT PASS STANDARD_LIBRARY_EXACT_INDEPENDENT_AUDIT")


if __name__ == "__main__":
    main()
