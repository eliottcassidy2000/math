#!/usr/bin/env python3
"""Primary exact certificate for proposed THM-4368.

The script independently evaluates THM-4364's displayed binomial consumer
and monomial packet law.  It checks the prefix/difference equivalence of all
simplex-order streams, the fixed-row Pascal Gram basis and consecutive unit
minors, the exact annihilator clock, packet boundary generating series,
packet feasibility, triangular natural addresses, reciprocal reflection,
and the named quotient hostiles.  It imports no repository computation.

All verification uses ``require`` rather than ``assert``, so the audit stays
active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import comb, gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + label)


def ceil_div(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError("positive divisor required")
    return -((-a) // b)


def triangular(n: int) -> int:
    if n < 0:
        raise ValueError("nonnegative triangular index required")
    return n * (n + 1) // 2


def integer_triangular(n: int) -> int:
    return n * (n + 1) // 2


def bounds(ell: int) -> tuple[int, int]:
    return ceil_div(ell, 2), ceil_div(ell, 3)


def prefix(values: list[int]) -> list[int]:
    total = 0
    result = []
    for value in values:
        total += value
        result.append(total)
    return result


def difference(values: list[int]) -> list[int]:
    previous = 0
    result = []
    for value in values:
        result.append(value - previous)
        previous = value
    return result


def iterate(function, values: list[int], count: int) -> list[int]:
    result = values[:]
    for _ in range(count):
        result = function(result)
    return result


def direct_stream_consumer(raw: list[int], q: int) -> list[int]:
    """Rows of L_q from a signed diagonal stream beginning at row s."""
    result = []
    for i in range(len(raw)):
        result.append(sum(comb(i - j + q, q) * raw[j]
                          for j in range(i + 1)))
    return result


def matrix_product(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left or not right:
        return []
    width = len(right[0])
    inner = len(right)
    require(all(len(row) == inner for row in left), "matrix product left shape")
    require(all(len(row) == width for row in right), "matrix product right shape")
    return [[sum(left[i][j] * right[j][k] for j in range(inner))
             for k in range(width)]
            for i in range(len(left))]


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*matrix)]


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "determinant square shape")
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for k in range(size - 1):
        if work[k][k] == 0:
            swap = next((i for i in range(k + 1, size) if work[i][k]), None)
            require(swap is not None, "determinant nonzero pivot")
            work[k], work[swap] = work[swap], work[k]
            sign = -sign
        pivot = work[k][k]
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                numerator = work[i][j] * pivot - work[i][k] * work[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                work[i][j] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def full_stream_checks() -> tuple[int, int]:
    cases = 0
    coordinates = 0
    for ell in range(2, 33):
        s, _ = bounds(ell)
        for length in range(1, 17):
            raw = [((ell + 3) * (j + 1) ** 2 - (2 * s + 1) * j - 17) % 31 - 15
                   for j in range(length)]
            # Include exact zero coordinates and alternating signs.
            if length >= 3:
                raw[length // 2] = 0
            if length >= 5:
                raw[-2] *= -1
            for q in range(9):
                direct = direct_stream_consumer(raw, q)
                summed = iterate(prefix, raw, q + 1)
                require(direct == summed, "full stream prefix formula")
                require(iterate(difference, direct, q + 1) == raw,
                        "full stream finite-difference inverse")
                if q == 0:
                    require(difference(direct) == raw, "Delta L0 equals diagonal")
                else:
                    require(difference(direct) == direct_stream_consumer(raw, q - 1),
                            "Delta Lq equals prior simplex order")
                for i in range(length):
                    require(direct[i] == summed[i], "full stream coordinate")
                    coordinates += 1
                cases += 1
    return cases, coordinates


def fixed_row_pascal_checks() -> tuple[int, int, int]:
    gram_cases = 0
    gram_entries = 0
    for r in range(16):
        pascal = [[comb(q, j) if j <= q else 0 for j in range(r + 1)]
                  for q in range(r + 1)]
        gram = [[comb(q + k, q) for k in range(r + 1)]
                for q in range(r + 1)]
        product = matrix_product(pascal, transpose(pascal))
        require(product == gram, "Pas0 fixed-row Pascal Gram factorization")
        require(bareiss_determinant(pascal) == 1, "Pascal determinant one")
        require(bareiss_determinant(gram) == 1, "Gram determinant one")
        for q in range(r + 1):
            for k in range(r + 1):
                require(product[q][k] == comb(q + k, q), "Gram entry")
                gram_entries += 1
        gram_cases += 1

    unit_minors = 0
    for q_start in range(31):
        for size in range(1, 13):
            minor = [[comb(q_start + i + k, k) for k in range(size)]
                     for i in range(size)]
            require(bareiss_determinant(minor) == 1,
                    "consecutive simplex-order unit minor")
            unit_minors += 1
    return gram_cases, gram_entries, unit_minors


def annihilator_orders(ell: int, m: int, depth: int) -> list[int]:
    _, rho = bounds(ell)
    return [q for q in range(rho + 4)
            if q < rho and m + q >= ell and depth <= m + q - ell]


def annihilator_clock_checks() -> tuple[int, int, int]:
    parameter_cases = 0
    clock_rows = 0
    triangular_positions = 0
    for ell in range(2, 121):
        s, rho = bounds(ell)
        for depth in range(16):
            # General fixed-row formula over a window containing onset and saturation.
            for m in range(s, ell + depth + 5):
                q0 = max(0, ell + depth - m)
                expected = list(range(q0, rho)) if q0 < rho else []
                actual = annihilator_orders(ell, m, depth)
                require(actual == expected, "exact THM-4364 annihilator order interval")
                rank = max(0, rho - q0)
                require(len(actual) == rank, "exact annihilator rank")
                require(rank <= m - s + 1, "annihilator rank fits row jet")
                parameter_cases += 1

            m0 = ell + depth - rho + 1
            require(annihilator_orders(ell, m0 - 1, depth) == [],
                    "annihilator clock empty predecessor")
            local_positions = 0
            for j in range(rho):
                m = m0 + j
                expected = list(range(rho - 1 - j, rho))
                require(annihilator_orders(ell, m, depth) == expected,
                        "annihilator clock unsaturated row")
                clock_rows += 1
                local_positions += len(expected)
            require(local_positions == triangular(rho),
                    "annihilator triangular block count")
            triangular_positions += local_positions
            require(annihilator_orders(ell, ell + depth, depth) == list(range(rho)),
                    "annihilator clock saturation")

    named = {
        "row8_A": (8, 8, 2, [2]),
        "row9_A": (10, 9, 2, [3]),
        "row10_C": (10, 10, 3, [3]),
        "row10_A_opposite": (11, 10, 2, [3]),
    }
    for ell, m, depth, expected in named.values():
        require(annihilator_orders(ell, m, depth) == expected,
                "named unique first-entry hierarchy consumer")
    require(annihilator_orders(11, 10, 2) == [3], "THM-4366 A uniqueness")
    require(annihilator_orders(10, 10, 3) == [3], "THM-4366 C uniqueness")
    weights_a = [(-1) ** (n - 6) * comb(13 - n, 3) for n in range(6, 11)]
    weights_c = [(-1) ** (n - 5) * comb(13 - n, 3) for n in range(5, 11)]
    require(weights_a == [35, -20, 10, -4, 1], "THM-4366 A stencil")
    require(weights_c == [56, -35, 20, -10, 4, -1], "THM-4366 C stencil")
    return parameter_cases, clock_rows, triangular_positions


def packet_type_valid(ell: int, degree: int, n0: int) -> bool:
    _, rho = bounds(ell)
    return degree >= rho and n0 >= degree and n0 + degree >= ell


def construct_monomial(ell: int, degree: int, n0: int) -> tuple[int, int, int, int]:
    if not packet_type_valid(ell, degree, n0):
        raise ValueError("invalid packet type")
    e = max(0, ell - 2 * degree)
    c = degree - e
    a = 2 * degree + e - ell
    b = n0 - degree - e
    require(min(a, b, c, e) >= 0, "constructed exponents nonnegative")
    require(a == 2 * c + 3 * e - ell, "constructed packet on diagonal")
    require(c + e == degree, "constructed packet degree")
    require(b + c + 2 * e == n0, "constructed packet start")
    return a, b, c, e


def realizing_monomials(ell: int, degree: int, n0: int) -> list[tuple[int, int, int, int]]:
    lower = max(0, ell - 2 * degree)
    upper = min(degree, n0 - degree)
    if lower > upper:
        return []
    return [(2 * degree + e - ell, n0 - degree - e, degree - e, e)
            for e in range(lower, upper + 1)]


def expanded_diagonal_packet(
        monomial: tuple[int, int, int, int]) -> list[tuple[int, int, int]]:
    a, b, c, e = monomial
    degree = c + e
    n0 = b + c + 2 * e
    r0 = a + 2 * b + e
    return [(n0 + k, r0 + 2 * k, comb(degree, k))
            for k in range(degree + 1)]


def q0_trace_from_expanded_packet(
        ell: int, packet: list[tuple[int, int, int]]) -> list[int]:
    s, _ = bounds(ell)
    terminal = packet[-1][0]
    result = []
    for m in range(s, terminal):
        result.append(sum(((-1) ** (n - s)) * coefficient
                          for n, r, coefficient in packet
                          if n <= m and r == 2 * n - ell))
    return result


def monomial_consumer(ell: int, m: int, q: int,
                      monomial: tuple[int, int, int, int]) -> int:
    a, b, c, e = monomial
    s, _ = bounds(ell)
    degree = c + e
    n0 = b + c + 2 * e
    r0 = a + 2 * b + e
    value = 0
    for k in range(degree + 1):
        n = n0 + k
        r = r0 + 2 * k
        if s <= n <= m and r == 2 * n - ell:
            value += ((-1) ** (n - s)
                      * comb(m + q - n, q)
                      * comb(degree, k))
    return value


def one_minus_z_coefficient(exponent: int, k: int) -> int:
    """Coefficient of z^k in (1-z)^exponent, including negative exponent."""
    if k < 0:
        return 0
    if exponent >= 0:
        return (-1) ** k * comb(exponent, k) if k <= exponent else 0
    return comb(-exponent + k - 1, k)


def packet_stream_formula(ell: int, m: int, q: int,
                          degree: int, n0: int) -> int:
    s, _ = bounds(ell)
    start_order = n0 - s
    return ((-1) ** start_order
            * one_minus_z_coefficient(degree - q - 1, m - n0))


def q0_packet_trace(ell: int, degree: int, n0: int) -> list[int]:
    s, _ = bounds(ell)
    return [packet_stream_formula(ell, m, 0, degree, n0)
            for m in range(s, n0 + degree)]


def address(u: int, v: int) -> int:
    if u < 1 or v < 1:
        raise ValueError("positive address coordinates required")
    return triangular(u + v - 2) + u


def decode_address(value: int) -> tuple[int, int]:
    if value < 1:
        raise ValueError("positive address required")
    total = 1
    while triangular(total) < value:
        total += 1
    u = value - triangular(total - 1)
    v = total + 1 - u
    return u, v


def reciprocal_type(ell: int, degree: int, n0: int) -> tuple[int, int]:
    s, _ = bounds(ell)
    u = n0 - s + 1
    v = degree
    return u, s + v - 1


def packet_and_address_checks() -> tuple[int, int, int, int, int, int, int, tuple[int, int, int]]:
    feasibility_cases = 0
    valid_types = 0
    packet_coordinates = 0
    reciprocal_valid = 0
    reciprocal_invalid = 0
    total_realizers = 0
    max_multiplicity = 0
    max_multiplicity_witness = (0, 0, 0)

    for ell in range(2, 23):
        s, rho = bounds(ell)
        trace_owners: dict[tuple[int, ...], tuple[int, int]] = {}
        for degree in range(1, ell + 7):
            for n0 in range(s, 2 * ell + 11):
                brute_es = []
                for e in range(degree + 1):
                    c = degree - e
                    a = 2 * degree + e - ell
                    b = n0 - degree - e
                    if min(a, b, c, e) >= 0:
                        brute_es.append(e)
                lower = max(0, ell - 2 * degree)
                upper = min(degree, n0 - degree)
                predicted_es = list(range(lower, upper + 1)) if lower <= upper else []
                require(brute_es == predicted_es, "exact exponent-fibre interval")
                brute = bool(brute_es)
                predicted = packet_type_valid(ell, degree, n0)
                require(brute == predicted, "packet feasibility iff")
                feasibility_cases += 1
                if not predicted:
                    continue

                valid_types += 1
                monomial = construct_monomial(ell, degree, n0)
                realizers = realizing_monomials(ell, degree, n0)
                multiplicity = upper - lower + 1
                require(len(realizers) == multiplicity == len(brute_es),
                        "exact exponent-fibre multiplicity")
                require([entry[3] for entry in realizers] == predicted_es,
                        "exponent fibre parametrized by e")
                total_realizers += multiplicity
                if multiplicity > max_multiplicity:
                    max_multiplicity = multiplicity
                    max_multiplicity_witness = (ell, degree, n0)
                n1 = n0 + degree
                q_values = sorted({0, 1, 2, 3, 4, 5,
                                   max(0, degree - 2), degree - 1,
                                   degree, degree + 1})
                for q in q_values:
                    previous = 0
                    for m in range(s, n1 + 6):
                        direct = monomial_consumer(ell, m, q, monomial)
                        formula = packet_stream_formula(ell, m, q, degree, n0)
                        require(direct == formula, "packet Fq coefficient formula")
                        if q == 0:
                            expected_delta = 0
                            k = m - n0
                            if 0 <= k <= degree:
                                expected_delta = ((-1) ** (m - s) * comb(degree, k))
                        else:
                            expected_delta = packet_stream_formula(
                                ell, m, q - 1, degree, n0)
                        require(direct - previous == expected_delta,
                                "packet Delta simplex recurrence")
                        previous = direct
                        packet_coordinates += 1

                # The sharp Fq boundary cases.
                sign = (-1) ** (n0 - s)
                if degree >= 2:
                    require(monomial_consumer(ell, n0, degree - 2, monomial) == sign,
                            "q=N-2 first coefficient")
                    require(monomial_consumer(ell, n0 + 1, degree - 2, monomial) == -sign,
                            "q=N-2 second coefficient")
                    require(monomial_consumer(ell, n0 + 2, degree - 2, monomial) == 0,
                            "q=N-2 compact end")
                require(monomial_consumer(ell, n0, degree - 1, monomial) == sign,
                        "q=N-1 impulse")
                require(monomial_consumer(ell, n0 + 1, degree - 1, monomial) == 0,
                        "q=N-1 impulse end")
                for extra in range(5):
                    require(monomial_consumer(ell, n0 + extra, degree, monomial) == sign,
                            "q=N constant tail")
                    require(monomial_consumer(ell, n0 + extra, degree + 1, monomial)
                            == sign * (extra + 1), "q=N+1 affine tail")

                # q=0 support and boundary orders recover (degree,n0).
                trace = q0_packet_trace(ell, degree, n0)
                support = [s + i for i, value in enumerate(trace) if value]
                require(support == list(range(n0, n1)), "q0 trace exact support interval")
                require(support[0] == n0, "q0 trace zero-boundary valuation")
                require(len(support) == degree, "q0 trace one-boundary valuation")
                require(n1 - ell == monomial[0] + monomial[1],
                        "packet endpoint reconstructs depth")
                trace_key = tuple(trace)
                require(trace_key not in trace_owners,
                        "q0 trace injective on normalized packet types")
                trace_owners[trace_key] = (degree, n0)

                canonical_packet = expanded_diagonal_packet(monomial)
                canonical_trace = q0_trace_from_expanded_packet(ell, canonical_packet)
                require(canonical_trace == trace, "expanded packet gives q0 trace")
                for realizer in realizers:
                    require(expanded_diagonal_packet(realizer) == canonical_packet,
                            "all exponent realizers give same full diagonal packet")
                    require(q0_trace_from_expanded_packet(
                        ell, expanded_diagonal_packet(realizer)) == trace,
                        "all exponent realizers give same full q0 trace")

                # Triangular address and reciprocal boundary reflection.
                u = n0 - s + 1
                v = degree
                value = address(u, v)
                require(decode_address(value) == (u, v), "triangular address decoder")
                total = u + v - 1
                require(triangular(total - 1) < value <= triangular(total),
                        "triangular address block")
                require(value + address(v, u) == total * total + 1,
                        "reciprocal square-center identity")
                degree_star, n0_star = reciprocal_type(ell, degree, n0)
                require(n0_star + degree_star == n1, "reciprocal terminal row")
                require(reciprocal_type(ell, degree_star, n0_star) == (degree, n0),
                        "reciprocal type involution")
                predicted_star = u >= rho and u - v <= s - 1
                actual_star = packet_type_valid(ell, degree_star, n0_star)
                require(actual_star == predicted_star, "reciprocal source-validity iff")
                if actual_star:
                    reciprocal_valid += 1
                    construct_monomial(ell, degree_star, n0_star)
                    reflected = q0_packet_trace(ell, degree_star, n0_star)
                    expected = [0] * len(trace)
                    start_order = u - 1
                    for i in range(start_order + 1):
                        expected[v - 1 + i] += ((-1) ** start_order
                                                * (-1) ** i
                                                * comb(start_order, i))
                    reflection_sign = -1 if (v - u) % 2 else 1
                    require(reflected == [reflection_sign * x for x in expected],
                            "trace z to one-minus-z reflection")
                else:
                    reciprocal_invalid += 1

    # Global address bijection, fixed-point parity, and honest ties.
    seen: set[int] = set()
    for u in range(1, 151):
        for v in range(1, 151):
            value = address(u, v)
            require(value not in seen, "triangular address injective rectangle")
            seen.add(value)
            require(decode_address(value) == (u, v), "triangular address global inverse")
            total = u + v - 1
            require(value + address(v, u) == total * total + 1,
                    "global reciprocal address reflection")
    for value in range(1, triangular(300) + 1):
        u, v = decode_address(value)
        require(u >= 1 and v >= 1, "triangular decoder positive pair")
        require(address(u, v) == value, "triangular address global surjection")
    for total in range(1, 301):
        fixed = [(u, total + 1 - u) for u in range(1, total + 1)
                 if u == total + 1 - u]
        require(len(fixed) == (1 if total % 2 else 0),
                "reciprocal fixed-point parity")
        require(integer_triangular(-total) == triangular(total - 1),
                "negative triangular reflection")

    return (feasibility_cases, valid_types, packet_coordinates,
            reciprocal_valid, reciprocal_invalid, total_realizers,
            max_multiplicity, max_multiplicity_witness)


def named_hostile_checks() -> list[str]:
    ell = 10
    s, rho = bounds(ell)

    named_specs = {
        "orientation_45": (4, 5),
        "orientation_46": (4, 6),
        "reciprocal_54": (5, 4),
        "fixed_44": (4, 4),
        "same_ray_810": (8, 10),
    }
    records: dict[str, dict[str, object]] = {}
    for label, (u, v) in named_specs.items():
        degree = v
        n0 = s + u - 1
        require(packet_type_valid(ell, degree, n0), "named hostile source-valid")
        monomial = construct_monomial(ell, degree, n0)
        records[label] = {
            "u": u,
            "v": v,
            "degree": degree,
            "n0": n0,
            "monomial": monomial,
            "address": address(u, v),
            "signature": ((-1) ** (u - 1), u - 1, v - 1),
        }

    left = records["orientation_45"]
    right = records["orientation_46"]
    require(left["u"] < left["v"] and right["u"] < right["v"],
            "same order orientation hostile")
    require(left["signature"] != right["signature"],
            "order orientation does not determine trace")

    scaled = records["same_ray_810"]
    require(left["u"] * scaled["v"] == left["v"] * scaled["u"],
            "same primitive rational ray")
    require((left["u"] // gcd(left["u"], left["v"]),
             left["v"] // gcd(left["u"], left["v"]))
            == (scaled["u"] // gcd(scaled["u"], scaled["v"]),
                scaled["v"] // gcd(scaled["u"], scaled["v"])),
            "same reduced ray pair")
    require(left["signature"] != scaled["signature"],
            "primitive ray does not determine trace scale")

    reciprocal = records["reciprocal_54"]
    require(reciprocal_type(ell, left["degree"], left["n0"])
            == (reciprocal["degree"], reciprocal["n0"]),
            "named reciprocal packet")
    require(left["address"] == 32 and reciprocal["address"] == 33,
            "named adjacent reciprocal addresses")
    require(left["address"] + reciprocal["address"] == 8 * 8 + 1,
            "named square-center address sum")

    fixed = records["fixed_44"]
    require(reciprocal_type(ell, fixed["degree"], fixed["n0"])
            == (fixed["degree"], fixed["n0"]), "named honest fixed tie")
    require(fixed["address"] == 25, "named fixed address")

    # Address reciprocity is not automatically a source symmetry.
    invalid_monomial = (2, 1, 6, 0)  # x^2 u p^6
    invalid_degree = invalid_monomial[2] + invalid_monomial[3]
    invalid_n0 = invalid_monomial[1] + invalid_monomial[2] + 2 * invalid_monomial[3]
    require(invalid_monomial[0] == 2 * invalid_monomial[2]
            + 3 * invalid_monomial[3] - ell, "source-closure hostile diagonal")
    require(packet_type_valid(ell, invalid_degree, invalid_n0),
            "source-closure hostile valid source")
    invalid_star = reciprocal_type(ell, invalid_degree, invalid_n0)
    require(invalid_star == (3, 10), "source-closure hostile reflected type")
    require(invalid_star[0] < rho and not packet_type_valid(ell, *invalid_star),
            "reciprocal reflection leaves source cone")

    # The diagonal-incidence bit and the erased exponent fibre are both real.
    diagonal_a = (0, 3, 5, 0)
    diagonal_b = (1, 2, 4, 1)
    named_fibre = [(0, 3, 5, 0), (1, 2, 4, 1),
                   (2, 1, 3, 2), (3, 0, 2, 3)]
    require(realizing_monomials(10, 5, 8) == named_fibre,
            "named exponent fibre e=0 through 3")
    require((diagonal_a[2] + diagonal_a[3],
             diagonal_a[1] + diagonal_a[2] + 2 * diagonal_a[3])
            == (diagonal_b[2] + diagonal_b[3],
                diagonal_b[1] + diagonal_b[2] + 2 * diagonal_b[3])
            == (5, 8), "distinct exponents share one diagonal packet type")
    for q in range(6):
        for m in range(s, 17):
            require(monomial_consumer(ell, m, q, diagonal_a)
                    == monomial_consumer(ell, m, q, diagonal_b),
                    "same packet type gives same hierarchy stream")
    off_diagonal = (1, 3, 5, 0)
    require(off_diagonal[0] != 2 * off_diagonal[2] + 3 * off_diagonal[3] - ell,
            "off-diagonal hostile misses intercept")
    for q in range(6):
        for m in range(s, 17):
            require(monomial_consumer(ell, m, q, off_diagonal) == 0,
                    "off-diagonal packet invisible to hierarchy")

    # THM-4367 arithmetic contrast: one reduced metric ray, distinct binders.
    A = 3371
    M = 14 * A
    S = 1303
    a, b, kappa = 34, 75, 29
    require(a + S * b == A * kappa, "THM-4367 primitive-ray equation")
    require(gcd(a, b) == 1 and gcd(kappa, 14 * b) == 1,
            "THM-4367 reduced metric pair")
    lrc_rows = []
    for scale in (155, 169):
        require((scale * kappa - 1) % 14 == 0,
                "THM-4367 scale residue")
        P = b * scale
        rho_value = A - a * scale
        tooth = (scale * kappa - 1) // 14
        require(P >= 11019 and P % 2 == 1 and abs(rho_value) < A,
                "THM-4367 strict active tail cell")
        require(S * P == M * tooth + rho_value,
                "THM-4367 centered residue reconstruction")
        exit_value = Fraction(S, M) + Fraction(a, M * b)
        require(exit_value == Fraction(kappa, 14 * b),
                "THM-4367 metric exit ray factorization")
        lrc_rows.append((P, rho_value, tooth, exit_value))
    require(lrc_rows[0][3] == lrc_rows[1][3], "LRC same metric ray consumer")
    require(lrc_rows[0][:3] != lrc_rows[1][:3], "LRC physical binder needs scale")

    ledger = []
    for label in named_specs:
        record = records[label]
        ledger.append(
            f"{label}:pair=({record['u']},{record['v']}),"
            f"type=({record['degree']},{record['n0']}),"
            f"mono={record['monomial']},addr={record['address']},"
            f"trace_sig={record['signature']}"
        )
    ledger.append(
        "source_reflection_hostile:pair=(3,6),reflected_type=(3,10),valid=false"
    )
    ledger.append(
        "packet_fibre:ell=10,N=5,n0=8,e=0..3,mu=4;"
        "off_diagonal=(1,3,5,0)->zero"
    )
    ledger.append(
        "lrc_ray_contrast:"
        + ";".join(f"P={P},rho={rho_value},tooth={tooth},exit={exit_value}"
                   for P, rho_value, tooth, exit_value in lrc_rows)
    )
    return ledger


def main() -> None:
    full_cases, full_coordinates = full_stream_checks()
    gram_cases, gram_entries, unit_minors = fixed_row_pascal_checks()
    clock_cases, clock_rows, clock_positions = annihilator_clock_checks()
    (feasibility_cases, valid_types, packet_coordinates,
     reciprocal_valid, reciprocal_invalid, total_realizers,
     max_multiplicity, max_multiplicity_witness) = packet_and_address_checks()
    ledger = named_hostile_checks()
    semantic_digest = sha256(("\n".join(ledger) + "\n").encode("ascii")).hexdigest()

    print("THM-4368 diagonal boundary-address/simplex-stream primary: PASS")
    print(f"checks={CHECKS}")
    print(f"full_stream_cases={full_cases} coordinates={full_coordinates}")
    print(f"fixed_row_gram_cases={gram_cases} entries={gram_entries}")
    print(f"consecutive_unit_minors={unit_minors}")
    print(f"annihilator_parameter_cases={clock_cases}")
    print(f"annihilator_clock_rows={clock_rows} positions={clock_positions}")
    print("thm4366_unique_consumers=A:(ell=11,m=10,d=2,q=3);"
          "C:(ell=10,m=10,d=3,q=3)")
    print(f"packet_feasibility_cases={feasibility_cases} valid_types={valid_types}")
    print("exponent_fibre=e:max(0,ell-2N)..min(N,n0-N);mu=upper-lower+1")
    print(f"packet_exponent_realizers={total_realizers} "
          f"bounded_max_multiplicity={max_multiplicity} "
          f"bounded_first_max_witness={max_multiplicity_witness}")
    print(f"packet_Fq_coordinates={packet_coordinates}")
    print(f"reciprocal_source_valid={reciprocal_valid} invalid={reciprocal_invalid}")
    for line in ledger:
        print(line)
    print(f"semantic_digest={semantic_digest}")


if __name__ == "__main__":
    main()
