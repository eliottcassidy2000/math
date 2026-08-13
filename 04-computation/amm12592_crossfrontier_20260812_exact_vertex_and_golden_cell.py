#!/usr/bin/env python3
"""Exact AMM-12592 cross-frontier audit.

This file complements the numerical sparse-real-LP probe kps-S169 in two
directions.

1. It gives a rational, nonintegral vertex of the *same* halved sparse
   junk-flow polytope at the positive-control epoch R=8.  Forty-two active
   integer constraints have determinant 1648, while the vertex has
   denominator 103.  Thus this polytope is not integral (and its displayed
   constraint system is not totally unimodular), even at an epoch known to
   possess integer witnesses.

2. It computes the exact finite-profile cell containing
      gamma_* = log_5(phi^2)
   at R=512.  The whole degree word d_n=floor(gamma_* n), 512<=n<1024,
   is unchanged when gamma_* is replaced by 653/1092.  Hence the golden
   constant enters a finite epoch only through its mechanical degree word;
   its special closed form is an asymptotic/all-scale datum.

Only Python's standard library is used.  Every decision is over integers or
fractions; there is no LP solver or floating-point arithmetic.

Reproduction:
  python 04-computation/amm12592_crossfrontier_20260812_exact_vertex_and_golden_cell.py
  python -O 04-computation/amm12592_crossfrontier_20260812_exact_vertex_and_golden_cell.py
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import comb, gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fib_pair(n: int) -> tuple[int, int]:
    """Return (F_n,F_{n+1}) by exact fast doubling."""
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def compare_five_power_to_phi2(d: int, m: int) -> int:
    """Sign of 5^d-phi^(2m), exactly.

    phi^(2m)=(L_{2m}+sqrt(5) F_{2m})/2.  If
    A=2*5^d-L_{2m} is positive, comparison with sqrt(5)F_{2m}
    can be squared without changing the sign.
    """
    require(d >= 0 and m >= 0, "negative exponent")
    f, f_next = fib_pair(2 * m)
    lucas = 2 * f_next - f
    a = 2 * 5**d - lucas
    if a < 0:
        return -1
    if a == 0:
        return -1 if f else 0
    delta = a * a - 5 * f * f
    return (delta > 0) - (delta < 0)


def floor_gamma_star(m: int) -> int:
    """floor(m log_5(phi^2)), using exact integer comparisons."""
    d = 3 * m // 5
    while compare_five_power_to_phi2(d + 1, m) <= 0:
        d += 1
    while compare_five_power_to_phi2(d, m) > 0:
        d -= 1
    return d


def two_g_coeffs(r: int) -> list[int]:
    values = [2]
    b = 1
    for j in range(1, r):
        b = b * (r - j) // j
        values.append((-b if j & 1 else b) - 1)
    return values


def t4_row_load(r: int, d: int) -> list[int]:
    return [
        (-1) ** (d - t) * comb(r - 2 - t, d - t)
        - comb(d + 1, t + 1)
        + 2 * comb(d, t)
        for t in range(d + 1)
    ]


def sparse_half_polytope(r: int):
    """Build the exact y=j/2 version of the kps-S169 sparse LP."""
    profile = [floor_gamma_star(r + i) for i in range(r)]
    degrees = profile[:-1]
    offsets = [0]
    for d in degrees:
        offsets.append(offsets[-1] + d)

    def variable(i: int, t: int) -> int:
        return offsets[i] + t

    variable_meta: list[tuple[int, int] | None] = [None] * offsets[-1]
    for i, d in enumerate(degrees):
        for t in range(d):
            variable_meta[variable(i, t)] = (i, t)

    inequalities: list[tuple[tuple[int, ...], int, tuple]] = []
    bounds: list[tuple[int | None, int | None]] = [
        (None, None) for _ in variable_meta
    ]

    d0 = degrees[0]
    initial = t4_row_load(r, d0)
    require(all(value % 2 == 0 for value in initial), "initial parity")
    for t in range(d0):
        lower_u = -comb(d0 - 1, t)
        upper_u = comb(d0 - 1, t - 1) if t else 0
        load = initial[t] // 2
        bounds[variable(0, t)] = (load - upper_u, load - lower_u)
    require(0 <= initial[d0] // 2 <= 1, "initial top-cell gate")

    feed = two_g_coeffs(r)
    require(all(value % 2 == 0 for value in feed), "dyadic feed parity")

    for i in range(1, len(degrees)):
        d = degrees[i]
        previous_degree = degrees[i - 1]
        kernel = (1, 1) if d == previous_degree else (1, 2, 1)
        require(d - previous_degree in (0, 1), "non-Pascal degree jump")

        feed0 = 0
        feed1 = 0
        if d + i <= r - 1:
            feed0 += feed[d + i] // 2
        if d > previous_degree and d - 1 + i <= r - 1:
            z = feed[d - 1 + i] // 2
            feed0 += z
            feed1 += z

        for t in range(d + 1):
            row = [0] * len(variable_meta)
            for a, coefficient in enumerate(kernel):
                if 0 <= t - a < previous_degree:
                    row[variable(i - 1, t - a)] += coefficient
            if t < d:
                row[variable(i, t)] -= 1

            f = feed0 if t == 0 else (feed1 if t == 1 else 0)
            lower_u = -comb(d - 1, t) if t <= d - 1 else 0
            upper_u = comb(d - 1, t - 1) if t else 0
            inequalities.append(
                (tuple(row), upper_u - f, ("row-upper", i, t))
            )
            inequalities.append(
                (tuple(-z for z in row), f - lower_u, ("row-lower", i, t))
            )

    previous_degree = degrees[-1]
    last_degree = profile[-1]
    kernel = (1, 1) if last_degree == previous_degree else (1, 2, 1)
    require(last_degree - previous_degree in (0, 1), "terminal degree jump")
    last_row = len(degrees) - 1
    for t in range(last_degree + 1):
        row = [0] * len(variable_meta)
        for a, coefficient in enumerate(kernel):
            if 0 <= t - a < previous_degree:
                row[variable(last_row, t - a)] += coefficient
        inequalities.append(
            (tuple(row), comb(last_degree, t), ("terminal-upper", t))
        )
        inequalities.append(
            (tuple(-z for z in row), 0, ("terminal-lower", t))
        )

    require(all(meta is not None for meta in variable_meta), "variable map")
    return profile, variable_meta, inequalities, bounds


def dot(row: tuple[int, ...] | list[int], vector: list[Fraction]) -> Fraction:
    return sum((Fraction(a) * x for a, x in zip(row, vector)), Fraction(0))


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    a = [row[:] for row in matrix]
    n = len(a)
    require(all(len(row) == n for row in a), "nonsquare determinant")
    sign = 1
    previous_pivot = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            require(swap is not None, "singular active basis")
            a[k], a[swap] = a[swap], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous_pivot == 0, "Bareiss division")
                a[i][j] = numerator // previous_pivot
            a[i][k] = 0
        previous_pivot = pivot
    return sign * a[-1][-1]


def certify_fractional_vertex() -> dict[str, object]:
    profile, variable_meta, inequalities, bounds = sparse_half_polytope(8)
    require(profile == [4, 5, 5, 6, 7, 7, 8, 8], "R=8 profile")
    require(len(variable_meta) == 42 and len(inequalities) == 106, "R=8 size")

    raw = [
        "7", "-7", "1", "1",
        "0", "-490/103", "0", "0", "0",
        "-1", "-169/103", "128/103", "-2", "1",
        "-1", "140/103", "-828/103", "911/103", "-3", "1",
        "-1", "-169/103", "894/103", "-2150/103", "751/103", "7", "0",
        "0", "-375/103", "2270/103", "-2801/103", "146/103", "-73/103", "1",
        "1", "346/103", "799/103", "4969/103", "419/103", "-419/103", "8", "0",
    ]
    vertex = [Fraction(value) for value in raw]
    require(len(vertex) == len(variable_meta), "vertex length")

    active: dict[tuple, tuple[tuple[int, ...], int]] = {}
    for row, rhs, name in inequalities:
        value = dot(row, vertex)
        require(value <= rhs, f"violated inequality {name}")
        if value == rhs:
            active[name] = (row, rhs)
    for j, ((lower, upper), value) in enumerate(zip(bounds, vertex)):
        if lower is not None:
            require(value >= lower, f"violated lower bound {variable_meta[j]}")
            if value == lower:
                row = [0] * len(vertex)
                row[j] = 1
                active[("bound-lower",) + variable_meta[j]] = (tuple(row), lower)
        if upper is not None:
            require(value <= upper, f"violated upper bound {variable_meta[j]}")
            if value == upper:
                row = [0] * len(vertex)
                row[j] = 1
                active[("bound-upper",) + variable_meta[j]] = (tuple(row), upper)

    basis_names = [
        ("row-lower", 1, 0), ("row-lower", 1, 2),
        ("row-lower", 1, 3), ("row-upper", 1, 5),
        ("row-upper", 2, 0), ("row-lower", 2, 2),
        ("row-lower", 2, 4), ("row-lower", 2, 5),
        ("row-upper", 3, 0), ("row-lower", 3, 1),
        ("row-upper", 3, 2), ("row-lower", 3, 3),
        ("row-lower", 3, 5), ("row-upper", 3, 6),
        ("row-upper", 4, 0), ("row-upper", 4, 1),
        ("row-lower", 4, 2), ("row-upper", 4, 3),
        ("row-lower", 4, 6), ("row-upper", 4, 7),
        ("row-lower", 5, 0), ("row-upper", 5, 1),
        ("row-lower", 5, 2), ("row-upper", 5, 3),
        ("row-lower", 5, 4), ("row-upper", 5, 5),
        ("row-upper", 5, 6), ("row-lower", 5, 7),
        ("row-lower", 6, 0), ("row-lower", 6, 1),
        ("row-upper", 6, 2), ("row-lower", 6, 3),
        ("row-lower", 6, 4), ("row-lower", 6, 5),
        ("row-lower", 6, 6), ("row-upper", 6, 8),
        ("terminal-upper", 0), ("terminal-upper", 3),
        ("terminal-lower", 5), ("terminal-upper", 7),
        ("terminal-lower", 8), ("bound-lower", 0, 1),
    ]
    require(len(active) == 45, "unexpected active count")
    require(len(basis_names) == len(vertex), "basis size")
    require(all(name in active for name in basis_names), "inactive basis row")

    basis = [list(active[name][0]) for name in basis_names]
    basis_rhs = [active[name][1] for name in basis_names]
    require(
        all(dot(row, vertex) == rhs for row, rhs in zip(basis, basis_rhs)),
        "basis does not contain vertex",
    )
    determinant = bareiss_determinant(basis)
    require(determinant == 1648, "active determinant changed")

    denominators = sorted({value.denominator for value in vertex})
    fractional_count = sum(value.denominator != 1 for value in vertex)
    require(denominators == [1, 103] and fractional_count == 20, "denominator")
    require(1648 == 16 * 103, "determinant factorization")

    digest_text = ";".join(
        f"{i},{t}:{value.numerator}/{value.denominator}"
        for (i, t), value in zip(variable_meta, vertex)
    )
    return {
        "profile": profile,
        "variables": len(vertex),
        "inequalities": len(inequalities),
        "active": len(active),
        "basis_determinant": determinant,
        "fractional_count": fractional_count,
        "denominator": 103,
        "vertex_sha256": sha256(digest_text.encode("ascii")).hexdigest(),
    }


def certify_golden_profile_cell(r: int = 512) -> dict[str, object]:
    profile = [floor_gamma_star(n) for n in range(r, 2 * r)]
    lower, lower_n = max((Fraction(d, n), n) for n, d in zip(range(r, 2 * r), profile))
    upper, upper_n = min(
        (Fraction(d + 1, n), n) for n, d in zip(range(r, 2 * r), profile)
    )

    require(compare_five_power_to_phi2(lower.numerator, lower.denominator) < 0,
            "lower endpoint is not below gamma*")
    require(compare_five_power_to_phi2(upper.numerator, upper.denominator) > 0,
            "upper endpoint is not above gamma*")
    farey_cross = upper.numerator * lower.denominator - lower.numerator * upper.denominator
    require(farey_cross == 1, "profile endpoints are not Farey neighbours")
    require(upper - lower == Fraction(1, lower.denominator * upper.denominator),
            "Farey width")

    mediant = Fraction(
        lower.numerator + upper.numerator,
        lower.denominator + upper.denominator,
    )
    require(lower < mediant < upper, "mediant not in cell")
    rational_profile = [mediant.numerator * n // mediant.denominator for n in range(r, 2 * r)]
    require(rational_profile == profile, "rational replacement changed profile")

    # The determinant-one/Farey lemma says every reduced rational strictly
    # between a/b and c/d has denominator at least b+d.  Exhaustion below
    # b+d is included as an independent finite hostile check.
    for denominator in range(1, mediant.denominator):
        numerator = lower.numerator * denominator // lower.denominator + 1
        require(not (lower < Fraction(numerator, denominator) < upper),
                "smaller-denominator rational found")

    increments = "".join(str(profile[i + 1] - profile[i]) for i in range(len(profile) - 1))
    require(set(increments) <= {"0", "1"}, "degree word is not mechanical")
    return {
        "r": r,
        "first_degree": profile[0],
        "last_degree": profile[-1],
        "increment_ones": increments.count("1"),
        "lower": lower,
        "lower_binder": lower_n,
        "upper": upper,
        "upper_binder": upper_n,
        "width": upper - lower,
        "farey_cross": farey_cross,
        "mediant": mediant,
        "profile_sha256": sha256(
            ",".join(map(str, profile)).encode("ascii")
        ).hexdigest(),
        "increment_sha256": sha256(increments.encode("ascii")).hexdigest(),
    }


def main() -> None:
    vertex = certify_fractional_vertex()
    cell = certify_golden_profile_cell()
    print("AMM12592 CROSS-FRONTIER EXACT AUDIT 2026-08-12")
    print("status=FINITE-EXACT")
    print(
        "R8_fractional_vertex="
        f"profile={vertex['profile']};variables={vertex['variables']};"
        f"inequalities={vertex['inequalities']};active={vertex['active']};"
        f"basis_determinant={vertex['basis_determinant']}=16*103;"
        f"fractional_coordinates={vertex['fractional_count']};"
        f"common_nontrivial_denominator={vertex['denominator']}"
    )
    print(f"R8_vertex_sha256={vertex['vertex_sha256']}")
    print(
        "R512_golden_cell="
        f"[{cell['lower']},{cell['upper']});width={cell['width']};"
        f"binders=({cell['lower_binder']},{cell['upper_binder']});"
        f"farey_cross={cell['farey_cross']};"
        f"smallest_denominator_interior_rational={cell['mediant']}"
    )
    print(
        "R512_profile="
        f"d512={cell['first_degree']};d1023={cell['last_degree']};"
        f"increment_ones={cell['increment_ones']};"
        f"profile_sha256={cell['profile_sha256']};"
        f"increment_sha256={cell['increment_sha256']}"
    )
    print(
        "conclusion_R8=the_sparse_real_polytope_is_not_integral_and_the_"
        "displayed_constraint_system_is_not_totally_unimodular"
    )
    print(
        "conclusion_golden=at_fixed_R512_gamma_star_has_no_role_beyond_"
        "selecting_the_finite_degree_word;653/1092_gives_the_identical_polytope"
    )
    print(
        "nonconsequence=fractional_vertex_does_not_exclude_integer_points;"
        "finite_rational_replacement_does_not_replace_the_all_scale_golden_rate;"
        "no_R512_integer_witness_or_AMM_bound_is_claimed"
    )


if __name__ == "__main__":
    main()
