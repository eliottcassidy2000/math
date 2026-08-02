#!/usr/bin/env python3
"""Exact companion for THM-3055 (grid-witness divisor upset).

The proof is elementary; this companion exhausts every nonempty speed-residue
support modulo q for 2 <= q <= 12, checks the adjoining-speed update law, and
rebuilds the six corrected THM-3043 witness rows.  It deliberately computes
the geometric grid predicate and the divisor-upset predicate independently.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def divisors(q: int) -> tuple[int, ...]:
    return tuple(d for d in range(1, q + 1) if q % d == 0)


def phi(d: int) -> int:
    return sum(1 for a in range(1, d + 1) if gcd(a, d) == 1)


def safe_grid_numerators(q: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        a
        for a in range(q)
        if all((a * v) % q != 0 for v in speeds)
    )


def unhit_divisors(q: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        d
        for d in divisors(q)
        if all(v % d != 0 for v in speeds)
    )


def predicted_grid_numerators(q: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    unhit = set(unhit_divisors(q, speeds))
    return tuple(a for a in range(q) if q // gcd(a, q) in unhit)


def residue_supports(q: int):
    # The representatives 1,...,q encode all nonempty value supports modulo q.
    residues = tuple(range(1, q + 1))
    for size in range(1, q + 1):
        yield from combinations(residues, size)


set_rows = 0
grid_cells = 0
upset_cells = 0
pure_top_rows = 0
for q in range(2, 13):
    q_divisors = divisors(q)
    for speeds in residue_supports(q):
        set_rows += 1
        actual = safe_grid_numerators(q, speeds)
        predicted = predicted_grid_numerators(q, speeds)
        require(actual == predicted, f"grid/upset mismatch q={q}, speeds={speeds}")
        grid_cells += q

        unhit = unhit_divisors(q, speeds)
        require(
            len(actual) == sum(phi(d) for d in unhit),
            f"totient count mismatch q={q}, speeds={speeds}",
        )

        # If d is unhit and d|e|q, then e is unhit: the denominator carrier is
        # an upward-closed upset in the divisor poset.
        for d in unhit:
            for e in q_divisors:
                if e % d == 0:
                    require(e in unhit, f"upset failure q={q}, d={d}, e={e}")
                upset_cells += 1

        pure_top = unhit == (q,)
        boundary = (
            all(v % q != 0 for v in speeds)
            and all(
                any(v % d == 0 for v in speeds)
                for d in q_divisors
                if 1 < d < q
            )
        )
        require(pure_top == boundary, f"pure-top iff failure q={q}, speeds={speeds}")
        pure_top_rows += int(pure_top)


extension_cells = 0
for q in range(2, 11):
    q_divisors = set(divisors(q))
    for speeds in residue_supports(q):
        before = set(unhit_divisors(q, speeds))
        for v in range(1, q + 1):
            after = set(unhit_divisors(q, speeds + (v,)))
            removed = {d for d in q_divisors if gcd(v, q) % d == 0}
            require(
                after == before - removed,
                f"adjoining law failure q={q}, speeds={speeds}, v={v}",
            )
            extension_cells += 1


examples = (
    (4, (1, 2, 3), (Fraction(1, 4), Fraction(3, 4))),
    (5, (1, 2, 3, 4), tuple(Fraction(a, 5) for a in range(1, 5))),
    (5, (1, 3, 4, 7), tuple(Fraction(a, 5) for a in range(1, 5))),
    (7, (1, 2, 3, 4, 5, 6), tuple(Fraction(a, 7) for a in range(1, 7))),
    (6, (1, 3, 4, 5, 9), (Fraction(1, 6), Fraction(5, 6))),
    (8, (1, 2, 3, 4, 5, 6, 7), tuple(Fraction(a, 8) for a in (1, 3, 5, 7))),
)
for q, speeds, expected in examples:
    actual = tuple(Fraction(a, q) for a in safe_grid_numerators(q, speeds))
    require(actual == expected, f"THM-3043 witness mismatch q={q}, speeds={speeds}")
    require(unhit_divisors(q, speeds) == (q,), f"THM-3043 top-upset mismatch q={q}")


# Hostile to the tempting claim that every safe q-grid numerator must be a
# unit: with q=8 and speed 1, all denominator layers 2,4,8 survive.
hostile_q = 8
hostile_speeds = (1,)
hostile_unhit = unhit_divisors(hostile_q, hostile_speeds)
hostile_grid = tuple(Fraction(a, hostile_q) for a in safe_grid_numerators(hostile_q, hostile_speeds))
require(hostile_unhit == (2, 4, 8), "composite hostile divisor layers changed")
require(
    hostile_grid == tuple(Fraction(a, 8) for a in range(1, 8)),
    "composite hostile grid changed",
)


print("theorem=grid-safe iff reduced denominator divides no speed")
print(
    f"exhaustive_universe=q=2..12;nonempty residue supports={set_rows};"
    f"grid cells={grid_cells};upset cells={upset_cells}"
)
print(f"pure_top_rows={pure_top_rows};pure_top_iff=q unhit and every proper divisor >1 hit")
print(f"adjoining_speed_cells={extension_cells};update=remove divisors of gcd(v,q)")
print("thm3043_rows=6/6;unhit_upsets={q};reduced denominators exactly q")
print("hostile=q8,S={1};unhit={2,4,8};nonunit grid witnesses survive")
print("conclusion=grid witness carrier is an upward-closed unhit-divisor upset")
