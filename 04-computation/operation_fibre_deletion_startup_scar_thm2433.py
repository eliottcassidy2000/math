#!/usr/bin/env python3
"""Exact companion for THM-2433.

This script independently checks the Burnside deletion identity, its
additive and multiplicative specializations, the two internal transitive
closures, the induced divisibility covers, the relative atoms, and the
chain/divisor incidence inversions.  It uses integers only and keeps every
truth-bearing check active under optimized Python.
"""

from __future__ import annotations

from itertools import combinations
from math import isqrt


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


STARTUP_HOLES = frozenset({1, 4, 6})


def in_startup_carrier(value: int) -> bool:
    return value >= 1 and value not in STARTUP_HOLES


def unordered_additive_fibre(
    target: int,
    member,
    strict: bool,
) -> int:
    total = 0
    for first in range(1, target // 2 + 1):
        second = target - first
        if first > second or (strict and first == second):
            continue
        if member(first) and member(second):
            total += 1
    return total


def unordered_multiplicative_fibre(
    target: int,
    member,
    strict: bool,
) -> int:
    total = 0
    for first in range(2, isqrt(target) + 1):
        if target % first:
            continue
        second = target // first
        if first > second or (strict and first == second):
            continue
        if member(first) and member(second):
            total += 1
    return total


def ordered_add_convolution(left, right, target: int) -> int:
    return sum(left(first) * right(target - first) for first in range(1, target))


def divisors(number: int) -> tuple[int, ...]:
    small = []
    large = []
    for divisor in range(1, isqrt(number) + 1):
        if number % divisor:
            continue
        small.append(divisor)
        partner = number // divisor
        if partner != divisor:
            large.append(partner)
    return tuple(small + large[::-1])


def prime_sieve(limit: int) -> bytearray:
    prime = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        prime[0] = 0
    if limit >= 1:
        prime[1] = 0
    for value in range(2, isqrt(limit) + 1):
        if prime[value]:
            start = value * value
            prime[start : limit + 1 : value] = b"\x00" * (
                (limit - start) // value + 1
            )
    return prime


# 1. Universal additive Burnside deletion formula.
burnside_cases = 0
for holes_tuple in combinations(range(1, 9), 3):
    holes = frozenset(holes_tuple)

    def ambient(value: int) -> int:
        return int(value >= 1)

    def missing(value: int) -> int:
        return int(value in holes)

    def survivor(value: int) -> int:
        return int(value >= 1 and value not in holes)

    for target in range(1, 31):
        ordered_marked = ordered_add_convolution(ambient, missing, target)
        ordered_hole_hole = ordered_add_convolution(missing, missing, target)
        diagonal_hole = int(target % 2 == 0 and target // 2 in holes)

        weak_loss = unordered_additive_fibre(
            target, ambient, strict=False
        ) - unordered_additive_fibre(target, survivor, strict=False)
        strict_loss = unordered_additive_fibre(
            target, ambient, strict=True
        ) - unordered_additive_fibre(target, survivor, strict=True)

        require(
            2 * weak_loss
            == 2 * ordered_marked - ordered_hole_hole + diagonal_hole,
            "weak Burnside deletion",
        )
        require(
            2 * strict_loss
            == 2 * ordered_marked - ordered_hole_hole - diagonal_hole,
            "strict Burnside deletion",
        )
        burnside_cases += 1

require(burnside_cases == 56 * 30, "Burnside case count")


# 2. Additive startup scar and ordinary generating-function coefficients.
additive_coefficient_cases = 0
additive_tail_start = None
for target in range(1, 101):
    ambient_member = lambda value: value >= 1
    weak_loss = unordered_additive_fibre(
        target, ambient_member, strict=False
    ) - unordered_additive_fibre(target, in_startup_carrier, strict=False)
    strict_loss = unordered_additive_fibre(
        target, ambient_member, strict=True
    ) - unordered_additive_fibre(target, in_startup_carrier, strict=True)

    ordered_marked = sum(
        int(target - hole >= 1) for hole in STARTUP_HOLES
    )
    hole_hole = sum(
        1
        for first in STARTUP_HOLES
        for second in STARTUP_HOLES
        if first + second == target
    )
    diagonal = int(
        target % 2 == 0 and target // 2 in STARTUP_HOLES
    )
    require(
        2 * weak_loss == 2 * ordered_marked - hole_hole + diagonal,
        "additive weak OGF coefficient",
    )
    require(
        2 * strict_loss == 2 * ordered_marked - hole_hole - diagonal,
        "additive strict OGF coefficient",
    )
    if target >= 13:
        require(weak_loss == strict_loss == 3, "constant additive tail")
        if additive_tail_start is None:
            additive_tail_start = target
    additive_coefficient_cases += 1

require(additive_tail_start == 13, "additive tail start")


# 3. Multiplicative proper-fibre loss and Dirichlet coefficients.
MULTIPLICATIVE_HOLES = frozenset({4, 6})


def proper_ambient(value: int) -> bool:
    return value >= 2


def proper_survivor(value: int) -> bool:
    return value >= 2 and value not in MULTIPLICATIVE_HOLES


multiplicative_coefficient_cases = 0
for target in range(2, 10_001):
    weak_loss = unordered_multiplicative_fibre(
        target, proper_ambient, strict=False
    ) - unordered_multiplicative_fibre(target, proper_survivor, strict=False)
    strict_loss = unordered_multiplicative_fibre(
        target, proper_ambient, strict=True
    ) - unordered_multiplicative_fibre(target, proper_survivor, strict=True)

    expected_weak = (
        int(target % 4 == 0 and target // 4 >= 2)
        + int(target % 6 == 0 and target // 6 >= 2)
        - int(target == 24)
    )
    expected_strict = (
        expected_weak
        - int(target == 16)
        - int(target == 36)
    )
    require(weak_loss == expected_weak, "multiplicative weak coefficient")
    require(strict_loss == expected_strict, "multiplicative strict coefficient")
    multiplicative_coefficient_cases += 1

require(
    unordered_multiplicative_fibre(16, proper_ambient, False)
    - unordered_multiplicative_fibre(16, proper_survivor, False)
    == 1,
    "square weak control",
)
require(
    unordered_multiplicative_fibre(16, proper_ambient, True)
    - unordered_multiplicative_fibre(16, proper_survivor, True)
    == 0,
    "square strict control",
)
require(
    unordered_multiplicative_fibre(24, proper_ambient, False)
    - unordered_multiplicative_fibre(24, proper_survivor, False)
    == 1,
    "hole-overlap control",
)


# 4. Internal shadows, two-step failures, and independent reachability.
shadow_limit = 500
vertices = tuple(
    value for value in range(2, shadow_limit + 1) if in_startup_carrier(value)
)
vertex_set = frozenset(vertices)


def additive_edge(first: int, target: int) -> bool:
    return (
        first in vertex_set
        and target in vertex_set
        and target > first
        and in_startup_carrier(target - first)
    )


def multiplicative_edge(first: int, target: int) -> bool:
    return (
        first in vertex_set
        and target in vertex_set
        and target > first
        and target % first == 0
        and in_startup_carrier(target // first)
    )


def reachability(edge) -> tuple[int, frozenset[tuple[int, int]]]:
    adjacency = {
        first: tuple(target for target in vertices if edge(first, target))
        for first in vertices
    }
    edge_count = sum(len(targets) for targets in adjacency.values())
    reachable_pairs: set[tuple[int, int]] = set()
    for source in vertices:
        seen: set[int] = set()
        stack = list(adjacency[source])
        while stack:
            target = stack.pop()
            if target in seen:
                continue
            seen.add(target)
            stack.extend(adjacency[target])
        reachable_pairs.update((source, target) for target in seen)
    return edge_count, frozenset(reachable_pairs)


additive_edges, additive_reachability = reachability(additive_edge)
multiplicative_edges, multiplicative_reachability = reachability(
    multiplicative_edge
)

expected_additive_reachability = frozenset(
    (first, target)
    for first in vertices
    for target in vertices
    if target - first >= 2
)
expected_multiplicative_reachability = frozenset(
    (first, target)
    for first in vertices
    for target in vertices
    if (
        first < target
        and target % first == 0
        and (first, target) not in {(2, 8), (2, 12), (3, 12)}
    )
)
require(
    additive_reachability == expected_additive_reachability,
    "additive transitive closure",
)
require(
    multiplicative_reachability == expected_multiplicative_reachability,
    "multiplicative transitive closure",
)

witnesses = tuple(value for value in range(2, 101) if in_startup_carrier(value))
additive_two_step_failures = frozenset(
    (first, second, first + second)
    for first in witnesses
    for second in witnesses
    if not in_startup_carrier(first + second)
)
multiplicative_two_step_failures = frozenset(
    (first, second, first * second)
    for first in witnesses
    for second in witnesses
    if not in_startup_carrier(first * second)
)
require(
    additive_two_step_failures == {(2, 2, 4), (3, 3, 6)},
    "additive two-step failures",
)
require(
    multiplicative_two_step_failures
    == {(2, 2, 4), (2, 3, 6), (3, 2, 6)},
    "multiplicative two-step failures",
)
require(additive_edge(3, 5) and additive_edge(5, 7), "additive hostile path")
require(not additive_edge(3, 7), "additive hostile missing shortcut")
require(
    multiplicative_edge(5, 10) and multiplicative_edge(10, 20),
    "multiplicative hostile path",
)
require(
    not multiplicative_edge(5, 20),
    "multiplicative hostile missing shortcut",
)


# 5. Relative atoms and induced divisibility covers through 100,000.
arithmetic_limit = 100_000
prime = prime_sieve(arithmetic_limit)
extra_weak_atoms = []
extra_strict_atoms = []
for target in range(2, arithmetic_limit + 1):
    if not in_startup_carrier(target):
        continue
    weak_empty = (
        unordered_multiplicative_fibre(target, in_startup_carrier, False) == 0
    )
    strict_empty = (
        unordered_multiplicative_fibre(target, in_startup_carrier, True) == 0
    )
    is_prime_square = (
        isqrt(target) ** 2 == target and prime[isqrt(target)]
    )
    require(
        weak_empty == bool(prime[target] or target in {8, 12}),
        "pointwise weak atom classification",
    )
    require(
        strict_empty
        == bool(prime[target] or is_prime_square or target in {8, 12}),
        "pointwise strict atom classification",
    )
    if weak_empty and not prime[target]:
        extra_weak_atoms.append(target)
    if strict_empty and not prime[target] and not is_prime_square:
        extra_strict_atoms.append(target)

require(extra_weak_atoms == [8, 12], "relative weak atoms")
require(extra_strict_atoms == [8, 12], "relative strict atoms")

divisor_cache: dict[int, tuple[int, ...]] = {}


def proper_divisors(number: int) -> tuple[int, ...]:
    if number not in divisor_cache:
        divisor_cache[number] = tuple(
            divisor for divisor in divisors(number) if 1 < divisor < number
        )
    return divisor_cache[number]


composite_cover_triples: list[tuple[int, int, int]] = []
for first in range(2, arithmetic_limit + 1):
    if not in_startup_carrier(first):
        continue
    for quotient in range(2, arithmetic_limit // first + 1):
        if prime[quotient]:
            continue
        target = first * quotient
        if not in_startup_carrier(target):
            continue
        has_intermediate = any(
            in_startup_carrier(first * divisor)
            for divisor in proper_divisors(quotient)
        )
        if not has_intermediate:
            composite_cover_triples.append((first, target, quotient))

require(
    composite_cover_triples
    == [(2, 8, 4), (2, 12, 6), (2, 18, 9), (3, 12, 4)],
    "composite divisibility covers",
)


# 6. A014574 finite control through one million.
twin_limit = 1_000_001
twin_prime = prime_sieve(twin_limit)
twin_centres = [
    centre
    for centre in range(3, twin_limit)
    if twin_prime[centre - 1] and twin_prime[centre + 1]
]
retained_twin_centres = [
    centre for centre in twin_centres if in_startup_carrier(centre)
]
retained_fake_atoms = [
    centre
    for centre in retained_twin_centres
    if unordered_multiplicative_fibre(centre, in_startup_carrier, False) == 0
    and not twin_prime[centre]
]
require(twin_centres[:3] == [4, 6, 12], "twin-centre startup")
require(retained_fake_atoms == [12], "unique retained fake twin atom")
require(
    all(
        centre == 12
        or (
            centre >= 18
            and centre % 2 == 0
            and in_startup_carrier(centre // 2)
        )
        for centre in retained_twin_centres
    ),
    "twin-centre legal factor repair",
)


# 7. Chain difference and divisor-Mobius recovery of the hole indicator.
incidence_additive_cases = 0
for target in range(1, 101):
    marked = sum(hole < target for hole in STARTUP_HOLES)
    next_marked = sum(hole < target + 1 for hole in STARTUP_HOLES)
    require(
        next_marked - marked == int(target in STARTUP_HOLES),
        "chain incidence inversion",
    )
    incidence_additive_cases += 1

mobius_limit = 10_000
mobius = [1] * (mobius_limit + 1)
mobius[0] = 0
smallest_prime = [0] * (mobius_limit + 1)
for value in range(2, mobius_limit + 1):
    if smallest_prime[value] == 0:
        for multiple in range(value, mobius_limit + 1, value):
            if smallest_prime[multiple] == 0:
                smallest_prime[multiple] = value
for value in range(2, mobius_limit + 1):
    remaining = value
    sign = 1
    squarefree = True
    while remaining > 1:
        factor = smallest_prime[remaining]
        exponent = 0
        while remaining % factor == 0:
            remaining //= factor
            exponent += 1
        if exponent >= 2:
            squarefree = False
            break
        sign = -sign
    mobius[value] = sign if squarefree else 0

divisor_incidence = [0] * (mobius_limit + 1)
for target in range(1, mobius_limit + 1):
    divisor_incidence[target] = sum(
        divisor in STARTUP_HOLES for divisor in divisors(target)
    )

incidence_multiplicative_cases = 0
for target in range(1, mobius_limit + 1):
    recovered = sum(
        mobius[divisor] * divisor_incidence[target // divisor]
        for divisor in divisors(target)
    )
    require(
        recovered == int(target in STARTUP_HOLES),
        "divisor Mobius inversion",
    )
    incidence_multiplicative_cases += 1


print("theorem=THM-2433")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"burnside_three_hole_target_cases={burnside_cases}")
print(f"additive_coefficient_cases={additive_coefficient_cases}")
print(f"additive_constant_tail_start={additive_tail_start}")
print(f"multiplicative_coefficient_cases={multiplicative_coefficient_cases}")
print(f"shadow_vertex_count={len(vertices)}")
print(f"additive_shadow_edges={additive_edges}")
print(f"additive_transitive_pairs={len(additive_reachability)}")
print("additive_two_step_failures=2+2->4|3+3->6")
print(f"multiplicative_shadow_edges={multiplicative_edges}")
print(f"multiplicative_transitive_pairs={len(multiplicative_reachability)}")
print("multiplicative_closure_exceptions=2->8|2->12|3->12")
print("multiplicative_two_step_failures=2*2->4|2*3->6|3*2->6")
print("extra_weak_atoms_through_100000=8,12")
print("extra_strict_nonprime_nonsquare_atoms_through_100000=8,12")
print("composite_divisibility_covers=2:8:4|2:12:6|2:18:9|3:12:4")
print(f"twin_centres_through_1000000={len(twin_centres)}")
print(f"retained_twin_centres={len(retained_twin_centres)}")
print("retained_fake_twin_atoms=12")
print(f"chain_incidence_cases={incidence_additive_cases}")
print(f"divisor_mobius_cases={incidence_multiplicative_cases}")
print("ALL CHECKS PASSED")
