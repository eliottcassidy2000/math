#!/usr/bin/env python3
"""Exact arithmetic audit for the complete-packet monoid and branch transplant.

The geometric input is THM-3528: the fixed cleared norm sends every complete
packet polynomial to another polynomial.  This companion checks packet-grade
additivity, linearity of the renewal matrix, quotient grades after an old-L
return, and the diagonal ancestry shift forced by norm multiplicativity.  It
also contrasts the faithful formal return series with lossy harmonic sums.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json


EXPECTED_SEMANTIC_SHA256 = (
    "f004cd7643933e81a2fbce73a3df6d72c9e4943702cca7ee17b32d6690be3874"
)
MATRIX = ((7, -2), (3, -2))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: tuple[int, int], right: tuple[int, int],
             scale: int = 1) -> tuple[int, int]:
    return left[0] - scale * right[0], left[1] - scale * right[1]


def step(row: tuple[int, int]) -> tuple[int, int]:
    e, m = row
    return 7 * e - 2 * m, 3 * e - 2 * m


def iterate(row: tuple[int, int], count: int) -> tuple[int, int]:
    for _ in range(count):
        row = step(row)
    return row


def packet_signature(row: tuple[int, int]) -> tuple[tuple[int, ...], ...]:
    """Exponent vectors in the five complete packet faces."""
    e, m = row
    require(m % 3 == 0, ("packet divisibility", row))
    return (
        (e, m),
        (3 * e - 2 * m, e),
        (3 * e - 2 * m, e - m, 2 * m // 3, m // 3),
        (2 * e - 4 * m // 3, 2 * e - 2 * m // 3),
        (e, e - 2 * m // 3),
    )


def signature_add(left: tuple[tuple[int, ...], ...],
                  right: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(a + b for a, b in zip(lface, rface))
                 for lface, rface in zip(left, right))


def admissible(row: tuple[int, int]) -> bool:
    e, m = row
    return e >= 0 and m >= 0 and m % 3 == 0 and all(
        exponent >= 0 for face in packet_signature(row) for exponent in face
    )


def harmonic_subset_census() -> tuple[object, ...]:
    buckets: dict[Fraction, list[int]] = defaultdict(list)
    for mask in range(1 << 13):
        total = sum((Fraction(1, denominator)
                     for denominator in range(1, 14)
                     if (mask >> (denominator - 1)) & 1), Fraction(0))
        buckets[total].append(mask)
    collision_values = sum(len(masks) > 1 for masks in buckets.values())
    maximum_multiplicity = max(map(len, buckets.values()))
    left = 1 << (2 - 1)
    right = (1 << (3 - 1)) | (1 << (6 - 1))
    require(Fraction(1, 2) == Fraction(1, 3) + Fraction(1, 6),
            "harmonic hostile identity")
    half_fibre = tuple(buckets[Fraction(1, 2)])
    require(left in half_fibre and right in half_fibre and len(half_fibre) == 3,
            ("canonical harmonic collision", half_fibre))
    return (len(buckets), collision_values, maximum_multiplicity,
            left, right, half_fibre)


def main() -> None:
    rows = [(1, 0)]
    for _ in range(15):
        rows.append(step(rows[-1]))
    rows = tuple(rows)

    # Complete packet faces multiply, so all exponent vectors add.  The norm
    # transform grade M is additive on the same packet monoid.
    pair_checks = 0
    for left in rows:
        for right in rows:
            summed = add(left, right)
            require(packet_signature(summed)
                    == signature_add(packet_signature(left), packet_signature(right)),
                    ("packet product", left, right))
            require(step(summed) == add(step(left), step(right)),
                    ("renewal additivity", left, right))
            pair_checks += 1

    # If P with packet A(e,m) contains L^s, division by L^s removes s from e
    # and none from m.  The exact sharp packet bound is s<=e-m; at equality
    # the beta-face z exponent reaches zero.
    quotient_record = []
    for index, row in enumerate(rows[1:11], start=1):
        e, m = row
        boundary = e - m
        require(boundary >= 1, ("positive quotient boundary", index, row))
        samples = tuple(sorted({0, 1, 2, boundary}))
        sample_record = []
        for multiplicity in samples:
            if multiplicity > boundary:
                continue
            quotient = subtract(row, rows[0], multiplicity)
            require(admissible(quotient),
                    ("L quotient packet", index, multiplicity, quotient))
            require(packet_signature(row)
                    == signature_add(packet_signature((multiplicity, 0)),
                                     packet_signature(quotient)),
                    ("L packet factorization", index, multiplicity))
            sample_record.append((multiplicity, quotient))
        require(not admissible((m - 1, m)) if m else True,
                ("beyond-boundary hostile", index, row))
        quotient_record.append((index, row, boundary, tuple(sample_record)))

    # Formal branch transplant: an L^s factor first present in P_n becomes a
    # P_k^s factor in P_(n+k).  The residual packet is transported by M^k.
    transplant_record = []
    for event_index, multiplicity in ((8, 1), (8, 2), (10, 3)):
        seed_residual = subtract(rows[event_index], rows[0], multiplicity)
        require(admissible(seed_residual), ("event residual", event_index, multiplicity))
        descendants = []
        for offset in range(0, 6):
            target_grade = rows[event_index + offset]
            ancestor_grade = rows[offset]
            residual_grade = subtract(target_grade, ancestor_grade, multiplicity)
            require(residual_grade == iterate(seed_residual, offset),
                    ("transplant linearity", event_index, multiplicity, offset))
            require(admissible(residual_grade),
                    ("transplant residual packet", event_index, multiplicity, offset))
            descendants.append((event_index + offset, offset, multiplicity,
                                target_grade, ancestor_grade, residual_grade))
        transplant_record.append((event_index, multiplicity, tuple(descendants)))

    # A finitely supported defect sequence is retained exactly by its formal
    # coefficient word.  Applying the ancestry shift moves every coefficient
    # one place, i.e. multiplies the formal series by t.
    defect_word = (0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1)
    shifted_word = (0,) + defect_word
    require(shifted_word[1:] == defect_word and shifted_word[0] == 0,
            "formal return-series shift")

    harmonic_record = harmonic_subset_census()
    require(harmonic_record[:3] == (3712, 2944, 3), harmonic_record)

    record = (
        MATRIX,
        rows,
        pair_checks,
        tuple(quotient_record),
        tuple(transplant_record),
        defect_word,
        shifted_word,
        harmonic_record,
        ("operator_identity", "T(PQ)=T(P)T(Q); T(cP)=c^3 T(P)"),
        ("scope", "formal packet/factor transport; later defects and primality OPEN"),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== fixed Keller packet monoid and branch-transplant audit ==")
    print(f"orbit_rows_0_to_15={rows}")
    print(f"packet_product_and_M_additivity_checks={pair_checks}: PASS")
    print("cleared_norm_monoid_identity=T(PQ)=T(P)T(Q);scalar_law=T(cP)=c^3T(P)")
    print("L_factor_quotient=A(e,m)/A(s,0)=A(e-s,m);sharp_packet_bound=s<=e-m: PASS")
    print(f"quotient_boundary_record={tuple(quotient_record)}")
    print(f"synthetic_branch_transplant_controls={tuple(transplant_record)}")
    print(f"formal_return_word={defect_word};one_step_shift={shifted_word}")
    print(f"harmonic_subset_census_1_to_13={harmonic_record[:3]};collision=1/2=1/3+1/6")
    print("encoding_boundary=formal_series_retains indexed multiplicities;harmonic scalar does not")
    print(f"semantic_sha256={semantic}")
    print("scope=exact arithmetic corollary of complete-packet closure;synthetic defects are hostiles,not observed Keller returns")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
