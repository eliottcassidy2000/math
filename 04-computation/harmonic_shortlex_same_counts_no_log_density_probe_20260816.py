#!/usr/bin/env python3
"""Exact controls for equal shortlex counts with different/no log density.

Binary shortlex level n is [2^n,2^(n+1)-1].  Its first and last halves have
the same cardinality and Kraft mass, but limiting harmonic masses log(3/2)
and log(4/3).  Alternating enormously dominant blocks of first/last-half
levels therefore gives an arbitrary language with a_n=2^(n-1) at every
positive level and no logarithmic density.

The asymptotic proof is elementary and recorded in the companion reflection;
this file freezes the exact finite interval, harmonic, Kraft, K4, and block-
dominance controls used there.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from itertools import combinations


EXPECTED_LEDGER_SHA256 = "e135ce689fe04866568fd98961a45f3ebc9fabca7c03c9bf92c50be762f2cb45"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def interval_mass(first: int, last: int) -> Fraction:
    require(1 <= first <= last, "invalid harmonic interval")
    return sum((Fraction(1, value) for value in range(first, last + 1)), Fraction())


rows = []
for level in range(1, 13):
    first = 2**level
    middle = 3 * 2 ** (level - 1)
    last = 2 ** (level + 1) - 1
    half_count = 2 ** (level - 1)

    early = interval_mass(first, middle - 1)
    late = interval_mass(middle, last)
    full = interval_mass(first, last)

    require(middle - first == half_count, "early half count changed")
    require(last - middle + 1 == half_count, "late half count changed")
    require(early + late == full, "two halves ceased to partition the level")
    require(early > late, "shortlex order ceased to favour the early half")

    # THM-3499 equation (21), specialized to q=2 and a_n=2^(n-1).
    endpoint_lower = Fraction(half_count, 2 ** (level + 1) - 1)
    endpoint_upper = Fraction(half_count, 2**level)
    require(endpoint_lower <= late <= early <= endpoint_upper, "Kraft endpoint squeeze failed")
    require(Fraction(half_count, 2**level) == Fraction(1, 2), "level Kraft mass changed")

    rows.append(
        (
            level,
            first,
            middle,
            last,
            half_count,
            early.numerator,
            early.denominator,
            late.numerator,
            late.denominator,
        )
    )


# Binary words of length two are four K4 vertices.  Their six unordered
# comparisons are the six K4 edges; they are not six tournament vertices.
words = ("00", "01", "10", "11")
edges = tuple(combinations(words, 2))
require(len(words) == 4 and len(edges) == 6, "K4 vertex/edge census changed")
require(len(tuple(combinations(edges, 2))) == 15, "a tournament on six vertices would not have 15 comparisons")

# The three nonzero F_2-linear forms give the three 2+2 cuts/perfect
# matchings.  The first-bit cut is exactly early versus late shortlex halves.
bit_vectors = ((0, 0), (0, 1), (1, 0), (1, 1))
linear_forms = ((1, 0), (0, 1), (1, 1))
matching_partitions = []
for form in linear_forms:
    fibres = tuple(
        tuple(index for index, vector in enumerate(bit_vectors) if sum(form[j] * vector[j] for j in range(2)) % 2 == value)
        for value in (0, 1)
    )
    require(tuple(sorted(map(len, fibres))) == (2, 2), "Walsh form lost its 2+2 split")
    matching_partitions.append(fibres)
require(len(set(matching_partitions)) == 3, "the three Walsh matchings collapsed")
require(matching_partitions[0] == ((0, 1), (2, 3)), "first-bit cut no longer matches shortlex halves")


# Stage k has 2^(2^k) levels and alternates between early and late halves.
# Each new stage asymptotically swamps the complete preceding history.
stage_rows = []
past = 0
for stage in range(1, 9):
    length = 2 ** (2**stage)
    ratio = Fraction(past, length)
    if stage >= 2:
        require(ratio < Fraction(1, 2 ** (2 ** (stage - 1) - 1)), "stage dominance weakened")
    stage_rows.append((stage, "early" if stage % 2 else "late", length, ratio.numerator, ratio.denominator))
    past += length


ledger = repr((tuple(rows), tuple(edges), tuple(matching_partitions), tuple(stage_rows)))
ledger_hash = hashlib.sha256(ledger.encode("ascii")).hexdigest()
if EXPECTED_LEDGER_SHA256 != "TO_BE_PINNED":
    require(ledger_hash == EXPECTED_LEDGER_SHA256, "semantic ledger changed")

last_row = rows[-1]
last_early = Fraction(last_row[5], last_row[6])
last_late = Fraction(last_row[7], last_row[8])
print("== same counts, different/no shortlex logarithmic density ==")
print("binary level n=[2^n,2^(n+1)-1]; each selected half has a_n=2^(n-1), Kraft mass=1/2")
print("early limiting level mass=log(3/2); late limiting level mass=log(4/3)")
print(
    "level-12 exact masses: early~"
    f"{float(last_early):.12f} ({len(str(last_early.numerator))}/{len(str(last_early.denominator))} digits), "
    f"late~{float(last_late):.12f} ({len(str(last_late.numerator))}/{len(str(last_late.denominator))} digits); early>late: PASS"
)
print("alternating stage lengths 2^(2^k) make prior/current length ratio tend to zero: PASS")
print("therefore the staged language has identical level counts but two distinct normalized harmonic subsequential limits")
print("depth-2 carrier: 4 binary words=K4 vertices; 6 comparisons=K4 edges; 3 Walsh forms=3 perfect matchings")
print("an order-6 tournament would have 15 comparisons, so the six-edge carrier is not T6")
print(f"semantic ledger sha256={ledger_hash}")
print("scope: arbitrary binary shortlex language boundary; no automaton, LRC current, JC flux, or ancestry recovery")
print("all exact controls passed")
