#!/usr/bin/env python3
"""Exact referee for THM-3426's rough-composite odd-interval law."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def least_prime_factor(n: int) -> int:
    require(n >= 2, ("least prime factor domain", n))
    if n % 2 == 0:
        return 2
    for divisor in range(3, isqrt(n) + 1, 2):
        if n % divisor == 0:
            return divisor
    return n


def largest_odd_interval_radius(n: int, h: int) -> int:
    radius = (n - 1) // h
    if radius % 2 == 0:
        radius -= 1
    require(radius >= 1 and radius % 2 == 1, ("radius", n, h, radius))
    require(h * radius < n < h * (radius + 2), ("maximal radius", n, h, radius))
    return radius


def rational_connection(h: int) -> set[Fraction]:
    out: set[Fraction] = set()
    for a in range(1, h):
        for b in range(1, h - a + 1):
            if a + b <= h and gcd(a, b) == 1 and (a - b) % 2:
                out.add(Fraction(a, b))
                out.add(Fraction(-a, b))
    return out


def mod_fraction(value: Fraction, modulus: int) -> int:
    require(gcd(value.denominator, modulus) == 1, ("nonunit denominator", value, modulus))
    return value.numerator % modulus * pow(value.denominator, -1, modulus) % modulus


def odd_interval(radius: int, modulus: int) -> set[int]:
    return {value % modulus for value in range(-radius, radius + 1, 2)}


def units(modulus: int) -> list[int]:
    return [value for value in range(1, modulus) if gcd(value, modulus) == 1]


def collision_set(values: set[int], modulus: int) -> tuple[set[int], int]:
    collisions: set[int] = set()
    cells = 0
    for multiplier in units(modulus):
        for value in values:
            cells += 1
            if multiplier * value % modulus in values:
                collisions.add(multiplier)
                break
    return collisions, cells


def v2_integer(value: int) -> int:
    require(value != 0, ("zero v2", value))
    value = abs(value)
    exponent = 0
    while value % 2 == 0:
        value //= 2
        exponent += 1
    return exponent


def v2_fraction(value: Fraction) -> int:
    return v2_integer(value.numerator) - v2_integer(value.denominator)


def clique_tariff(h: int) -> int:
    return (h - 1).bit_length()


def threshold_one(h: int, c: int) -> int:
    return (h + 1) * (h * (h + 1) + c)


def threshold_two(h: int) -> int:
    return 2 * (h - 1) ** 3 + 1


def find_rough_controls(h: int, count: int) -> list[tuple[int, int, int]]:
    controls: list[tuple[int, int, int]] = []
    candidate = max(9, threshold_two(h))
    if candidate % 2 == 0:
        candidate += 1
    while len(controls) < count:
        factor = least_prime_factor(candidate)
        if factor > h and factor < candidate:
            radius = largest_odd_interval_radius(candidate, h)
            c = candidate - h * radius
            if 0 < c < 2 * h and candidate >= threshold_one(h, c):
                controls.append((candidate, radius, c))
        candidate += 2
    return controls


def check_rational_graphs() -> tuple[int, int, int]:
    connection_cells = 0
    color_cells = 0
    sharp_edges = 0
    for h in range(2, 41):
        connection = rational_connection(h)
        colors = clique_tariff(h)
        for value in connection:
            valuation = v2_fraction(value)
            require(valuation != 0, ("zero connection valuation", h, value))
            require(abs(valuation) <= colors - 1, ("valuation range", h, value))
            require(valuation % colors != 0, ("color collision", h, value))
            connection_cells += 1
            color_cells += 1
        sharp = [Fraction(2**power, 1) for power in range(colors)]
        for i, left in enumerate(sharp):
            for right in sharp[i + 1 :]:
                require(left / right in connection, ("missing sharp edge", h, left, right))
                sharp_edges += 1
    return connection_cells, color_cells, sharp_edges


def check_rounding_cells(controls: list[tuple[int, int, int, int]]) -> int:
    cells = 0
    for h, modulus, radius, c in controls:
        require(modulus >= threshold_one(h, c), ("threshold one", h, modulus, c))
        require(Fraction(modulus, h + 1) + (h + 1) <= radius, ("radius bound", h, modulus))
        for b in range(1, h + 1):
            for a in range(1, modulus // (h + 1) + 1):
                if gcd(a, b) != 1 or (a - b) % 2 == 0 or a + b < h + 1:
                    continue
                bound = Fraction(modulus, a + b) + max(a, b)
                require(bound <= Fraction(modulus, h + 1) + (h + 1),
                        ("convex bound", h, modulus, a, b))
                require(bound <= radius, ("odd lift radius", h, modulus, a, b))
                cells += 1
    return cells


def check_rough_composites() -> tuple[list[tuple[int, int, int, int, int, int, int]], int, int]:
    records: list[tuple[int, int, int, int, int, int, int]] = []
    action_cells = 0
    graph_cells = 0
    for h in range(2, 17):
        for modulus, radius, c in find_rough_controls(h, 2):
            factor = least_prime_factor(modulus)
            require(factor > h and factor < modulus, ("rough composite", h, modulus, factor))
            require(modulus >= threshold_one(h, c), ("first threshold", h, modulus))
            require(modulus >= threshold_two(h), ("second threshold", h, modulus))
            connection = rational_connection(h)
            reduced = {mod_fraction(value, modulus) for value in connection}
            require(len(reduced) == len(connection), ("connection collision", h, modulus))
            values = odd_interval(radius, modulus)
            collisions, cells = collision_set(values, modulus)
            action_cells += cells
            complement = set(units(modulus)) - collisions
            require(complement == reduced,
                    ("collision complement", h, modulus, sorted(complement ^ reduced)))

            nonunit_entries = sum(gcd(value, modulus) != 1 for value in values)
            require(nonunit_entries > 0, ("missing nonunit interval control", h, modulus))
            mod_to_rational = {mod_fraction(value, modulus): value for value in connection}
            vertices = [1] + sorted(reduced)
            rational_vertex = {1: Fraction(1)} | mod_to_rational
            colors = clique_tariff(h)
            for i, left in enumerate(vertices):
                for right in vertices[i + 1 :]:
                    modular_edge = left * pow(right, -1, modulus) % modulus in reduced
                    rational_edge = rational_vertex[left] / rational_vertex[right] in connection
                    require(modular_edge == rational_edge,
                            ("graph descent", h, modulus, left, right))
                    if modular_edge:
                        require(v2_fraction(rational_vertex[left]) % colors
                                != v2_fraction(rational_vertex[right]) % colors,
                                ("modular color", h, modulus, left, right))
                    graph_cells += 1

            sharp = [pow(2, power, modulus) for power in range(colors)]
            for i, left in enumerate(sharp):
                for right in sharp[i + 1 :]:
                    require(left * pow(right, -1, modulus) % modulus in reduced,
                            ("sharp modular edge", h, modulus, left, right))
            records.append((h, modulus, factor, radius, c, len(connection), nonunit_entries))
    return records, action_cells, graph_cells


def check_small_prime_hostile() -> tuple[int, int, int, tuple[int, ...]]:
    h = 9
    modulus = 969
    radius = 107
    c = modulus - h * radius
    require(c == 6 and modulus >= threshold_one(h, c), ("hostile threshold", c))
    require(least_prime_factor(modulus) == 3, ("hostile factor", modulus))
    values = odd_interval(radius, modulus)
    collisions, _ = collision_set(values, modulus)
    complement = set(units(modulus)) - collisions
    reduced_unit_fractions = {
        mod_fraction(value, modulus)
        for value in rational_connection(h)
        if gcd(value.numerator * value.denominator, modulus) == 1
    }
    extras = tuple(sorted(complement - reduced_unit_fractions))
    require(extras == (161, 325, 644, 808), ("roughness hostile", extras))
    require(not (reduced_unit_fractions - complement), "hostile reverse difference")
    return h, modulus, radius, extras


def main() -> None:
    rational_cells, color_cells, sharp_edges = check_rational_graphs()
    records, action_cells, graph_cells = check_rough_composites()
    rounding_cells = check_rounding_cells(
        [(h, modulus, radius, c) for h, modulus, _, radius, c, _, _ in records]
    )
    hostile = check_small_prime_hostile()
    semantic_payload = (tuple(records), hostile, rational_cells, sharp_edges)
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    print(f"rational_connection_cells={rational_cells}")
    print(f"dyadic_color_cells={color_cells}")
    print(f"sharp_clique_edges={sharp_edges}")
    print(f"rounding_cells={rounding_cells}")
    print(f"rough_unit_action_cells={action_cells}")
    print(f"rough_normalized_graph_cells={graph_cells}")
    print("rough_composite_controls=" + ";".join(
        f"h{h}:n{modulus}:spf{factor}:L{radius}:c{c}:D{size}:nonunit{nonunits}"
        for h, modulus, factor, radius, c, size, nonunits in records
    ))
    print(f"small_prime_hostile=h{hostile[0]}:n{hostile[1]}:L{hostile[2]}:extras{hostile[3]}")
    print(f"semantic_sha256={semantic}")
    print("rough_composite_odd_interval_collision=PASS")


if __name__ == "__main__":
    main()
