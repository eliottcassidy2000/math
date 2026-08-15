#!/usr/bin/env python3
"""Exact referee for THM-3423's odd-interval ratio and dyadic clique law."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for divisor in range(3, isqrt(n) + 1, 2):
        if n % divisor == 0:
            return False
    return True


def first_prime(start: int, residue: int | None = None, modulus: int = 1) -> int:
    candidate = max(3, start)
    if candidate % 2 == 0:
        candidate += 1
    while True:
        if (residue is None or candidate % modulus == residue) and is_prime(candidate):
            return candidate
        candidate += 2


def largest_odd_interval_radius(p: int, h: int) -> int:
    radius = (p - 1) // h
    if radius % 2 == 0:
        radius -= 1
    require(radius >= 1 and radius % 2 == 1, ("radius", p, h, radius))
    require(h * radius < p < h * (radius + 2), ("maximal radius", p, h, radius))
    return radius


def rational_connection(h: int) -> set[Fraction]:
    out: set[Fraction] = set()
    for a in range(1, h):
        for b in range(1, h - a + 1):
            if a + b <= h and gcd(a, b) == 1 and (a - b) % 2:
                out.add(Fraction(a, b))
                out.add(Fraction(-a, b))
    return out


def mod_fraction(value: Fraction, p: int) -> int:
    return (value.numerator % p) * pow(value.denominator % p, -1, p) % p


def odd_interval(radius: int, p: int) -> set[int]:
    return {value % p for value in range(-radius, radius + 1, 2)}


def ratio_set(values: set[int], p: int) -> set[int]:
    return {left * pow(right, -1, p) % p for left in values for right in values}


def v2_integer(value: int) -> int:
    require(value != 0, ("zero v2", value))
    value = abs(value)
    exponent = 0
    while value % 2 == 0:
        exponent += 1
        value //= 2
    return exponent


def v2_fraction(value: Fraction) -> int:
    return v2_integer(value.numerator) - v2_integer(value.denominator)


def clique_tariff(h: int) -> int:
    return (h - 1).bit_length()


def threshold_one(h: int, c: int) -> int:
    return (h + 1) * (h * (h + 1) + c)


def threshold_two(h: int) -> int:
    return 2 * (h - 1) ** 3 + 1


def find_control_prime(h: int, start: int = 3) -> tuple[int, int, int]:
    candidate = first_prime(max(start, 2 * h + 1))
    while True:
        radius = largest_odd_interval_radius(candidate, h)
        c = candidate - h * radius
        if candidate >= threshold_one(h, c) and candidate >= threshold_two(h):
            return candidate, radius, c
        candidate = first_prime(candidate + 2)


def check_rational_graphs() -> tuple[int, int, int]:
    connection_cells = 0
    color_cells = 0
    clique_cells = 0
    for h in range(2, 41):
        connection = rational_connection(h)
        colors = clique_tariff(h)
        for value in connection:
            valuation = v2_fraction(value)
            require(valuation != 0, ("zero connection valuation", h, value))
            require(abs(valuation) <= colors - 1, ("valuation range", h, value, valuation))
            require(valuation % colors != 0, ("color collision", h, value, colors))
            connection_cells += 1
            color_cells += 1

        sharp = [Fraction(2**power, 1) for power in range(colors)]
        for i, left in enumerate(sharp):
            for right in sharp[i + 1 :]:
                require(left / right in connection, ("missing sharp edge", h, left, right))
                clique_cells += 1
    return connection_cells, color_cells, clique_cells


def check_rounding_cells(controls: list[tuple[int, int, int, int]]) -> int:
    checks = 0
    for h, p, radius, c in controls:
        require(p >= threshold_one(h, c), ("threshold one", h, p, c))
        require(Fraction(p, h + 1) + (h + 1) <= radius, ("radius inequality", h, p))
        upper_a = p // (h + 1)
        for b in range(1, h + 1):
            for a in range(1, upper_a + 1):
                if gcd(a, b) != 1 or (a - b) % 2 == 0 or a + b < h + 1:
                    continue
                bound = Fraction(p, a + b) + max(a, b)
                require(bound <= Fraction(p, h + 1) + (h + 1), ("convex bound", h, p, a, b))
                require(bound <= radius, ("odd lift radius", h, p, a, b, bound, radius))
                checks += 1
    return checks


def check_finite_fields() -> tuple[list[tuple[int, int, int, int, int, int]], int, int]:
    records: list[tuple[int, int, int, int, int, int]] = []
    ratio_cells = 0
    graph_cells = 0
    for h in range(2, 17):
        p, radius, c = find_control_prime(h)
        connection = rational_connection(h)
        reduced = {mod_fraction(value, p) for value in connection}
        require(len(reduced) == len(connection), ("connection collision", h, p))
        values = odd_interval(radius, p)
        ratios = ratio_set(values, p)
        complement = set(range(1, p)) - ratios
        require(complement == reduced, ("ratio complement", h, p, len(complement), len(reduced)))
        ratio_cells += len(values) ** 2

        mod_to_rational = {mod_fraction(value, p): value for value in connection}
        vertices = [1] + sorted(reduced)
        rational_vertex = {1: Fraction(1)} | mod_to_rational
        colors = clique_tariff(h)
        for i, left in enumerate(vertices):
            for right in vertices[i + 1 :]:
                modular_edge = left * pow(right, -1, p) % p in reduced
                rational_edge = rational_vertex[left] / rational_vertex[right] in connection
                require(modular_edge == rational_edge, ("graph descent", h, p, left, right))
                if modular_edge:
                    left_color = v2_fraction(rational_vertex[left]) % colors
                    right_color = v2_fraction(rational_vertex[right]) % colors
                    require(left_color != right_color, ("modular color", h, p, left, right))
                graph_cells += 1

        records.append((h, p, radius, c, len(connection), colors))
    return records, ratio_cells, graph_cells


def check_h7_residue_controls() -> list[tuple[int, int, int, int]]:
    controls: list[tuple[int, int, int, int]] = []
    for residue in (1, 9, 11, 13):
        p = first_prime(512, residue, 14)
        radius = largest_odd_interval_radius(p, 7)
        c = p - 7 * radius
        require(c == {1: 8, 9: 2, 11: 4, 13: 6}[residue], ("h7 c", residue, p, c))
        require(p >= threshold_one(7, c) and p >= threshold_two(7), ("h7 threshold", p))
        connection = rational_connection(7)
        complement = set(range(1, p)) - ratio_set(odd_interval(radius, p), p)
        require(complement == {mod_fraction(value, p) for value in connection}, ("h7 complement", p))
        controls.append((residue, p, radius, len(complement)))
    require(len(rational_connection(7)) == 24, "h7 connection count")
    require(clique_tariff(7) == 3, "h7 clique tariff")
    return controls


def main() -> None:
    rational_cells, color_cells, sharp_edges = check_rational_graphs()
    records, ratio_cells, graph_cells = check_finite_fields()
    h7_controls = check_h7_residue_controls()
    rounding_controls = [(h, p, radius, c) for h, p, radius, c, _, _ in records]
    for h in range(17, 25):
        p, radius, c = find_control_prime(h)
        rounding_controls.append((h, p, radius, c))
    rounding_controls += [(7, p, radius, p - 7 * radius) for _, p, radius, _ in h7_controls]
    rounding_cells = check_rounding_cells(rounding_controls)

    semantic_payload = (tuple(records), tuple(h7_controls), rational_cells, sharp_edges)
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    print(f"rational_connection_cells={rational_cells}")
    print(f"dyadic_color_cells={color_cells}")
    print(f"sharp_clique_edges={sharp_edges}")
    print(f"rounding_cells={rounding_cells}")
    print(f"finite_ratio_pair_cells={ratio_cells}")
    print(f"finite_normalized_graph_cells={graph_cells}")
    print("finite_controls=" + ";".join(f"h{h}:p{p}:L{radius}:c{c}:D{size}:w{omega}" for h, p, radius, c, size, omega in records))
    print("h7_residue_controls=" + ";".join(f"r{r}:p{p}:L{radius}:D{size}" for r, p, radius, size in h7_controls))
    print(f"semantic_sha256={semantic}")
    print("odd_interval_ratio_complement=PASS")


if __name__ == "__main__":
    main()
