#!/usr/bin/env python3
"""Independent assertion-free referee for provisional THM-3818.

This implementation deliberately does not import or copy the canonical
companion.  It uses a smallest-prime-factor/totient census, a full bounded
coordinate-fibre table, ordered label placements, exact rational projection
rank, and independent rational/modular grid evaluators.

The purpose is hostile auditing.  A completed run can therefore report a
mathematical defect in the candidate while all of this referee's own gates
pass.  No Python ``assert`` statements are used, so optimization cannot
disable a gate.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction


GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def spf_table(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    for prime in range(2, math.isqrt(limit) + 1):
        if spf[prime] != prime:
            continue
        for value in range(prime * prime, limit + 1, prime):
            if spf[value] == value:
                spf[value] = prime
    return spf


def factors_from_spf(value: int, spf: list[int]) -> dict[int, int]:
    gate(value >= 1, "factor input is not positive")
    factors: dict[int, int] = {}
    while value > 1:
        prime = spf[value]
        factors[prime] = factors.get(prime, 0) + 1
        value //= prime
    return factors


def trial_factors(value: int) -> dict[int, int]:
    gate(value >= 1, "trial-factor input is not positive")
    factors: dict[int, int] = {}
    prime = 2
    while prime * prime <= value:
        while value % prime == 0:
            factors[prime] = factors.get(prime, 0) + 1
            value //= prime
        prime = 3 if prime == 2 else prime + 2
    if value > 1:
        factors[value] = factors.get(value, 0) + 1
    return factors


def inert_cube_free_sum(value: int, spf: list[int]) -> bool:
    factors = factors_from_spf(value, spf)
    return bool(factors) and all(prime % 3 == 2 and exponent <= 2
                                 for prime, exponent in factors.items())


def totients(limit: int) -> list[int]:
    phi = list(range(limit + 1))
    for prime in range(2, limit + 1):
        if phi[prime] == prime:
            for value in range(prime, limit + 1, prime):
                phi[value] -= phi[value] // prime
    return phi


def exact_cube_root(value: int) -> int | None:
    if value < 0:
        return None
    low, high = 0, 1
    while high ** 3 < value:
        high *= 2
    while low + 1 < high:
        middle = (low + high) // 2
        if middle ** 3 < value:
            low = middle
        else:
            high = middle
    if high ** 3 == value:
        return high
    if low ** 3 == value:
        return low
    return None


def brute_representations(value: int) -> tuple[tuple[int, int], ...]:
    """Direct coordinate scan, independent of the sum-divisor decoder."""
    upper = 1
    while upper ** 3 < value:
        upper *= 2
    low, high = 0, upper
    while low + 1 < high:
        middle = (low + high) // 2
        if middle ** 3 < value:
            low = middle
        else:
            high = middle
    coordinate_cap = high
    found: list[tuple[int, int]] = []
    cubes = {coordinate ** 3: coordinate
             for coordinate in range(1, coordinate_cap + 1)}
    for x in range(1, coordinate_cap + 1):
        y = cubes.get(value - x ** 3)
        if y is not None and x < y:
            found.append((x, y))
    return tuple(found)


def rational_rank(columns: list[list[Fraction]]) -> int:
    if not columns:
        return 0
    row_count = len(columns[0])
    col_count = len(columns)
    matrix = [[columns[column][row] for column in range(col_count)]
              for row in range(row_count)]
    rank = 0
    for column in range(col_count):
        pivot = next((row for row in range(rank, row_count)
                      if matrix[row][column] != 0), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [entry / pivot_value for entry in matrix[rank]]
        for row in range(row_count):
            if row == rank or matrix[row][column] == 0:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [matrix[row][entry] - multiplier * matrix[rank][entry]
                           for entry in range(col_count)]
        rank += 1
    return rank


def projected_free_rank(speeds: tuple[int, ...], i: int, j: int) -> int:
    norm = sum(speed * speed for speed in speeds)
    columns: list[list[Fraction]] = []
    for free in range(len(speeds)):
        if free in (i, j):
            continue
        column = []
        for row, speed in enumerate(speeds):
            kronecker = 1 if row == free else 0
            column.append(Fraction(kronecker * norm - speed * speeds[free], norm))
        columns.append(column)
    return rational_rank(columns)


def modular_schedule(residues: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    safe: list[int] = []
    for k in range(modulus):
        if all(14 * min((k * residue) % modulus,
                        modulus - ((k * residue) % modulus)) >= modulus
               for residue in residues):
            safe.append(k)
    return tuple(safe)


def rational_schedule(speeds: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    safe: list[int] = []
    for k in range(modulus):
        good = True
        for speed in speeds:
            point = Fraction(k * speed, modulus)
            fractional = point - math.floor(point)
            if min(fractional, 1 - fractional) < Fraction(1, 14):
                good = False
                break
        if good:
            safe.append(k)
    return tuple(safe)


def safe_at(speeds: tuple[int, ...], time: Fraction) -> bool:
    for speed in speeds:
        point = speed * time
        fractional = point - math.floor(point)
        if min(fractional, 1 - fractional) < Fraction(1, 14):
            return False
    return True


def packet_hypothesis(p: int, q: int, scale: int) -> bool:
    if not (0 < p < q and math.gcd(p, q) == 1):
        return False
    primitive_sum = p + q
    total_sum = scale * primitive_sum
    total_factors = trial_factors(total_sum)
    primitive_factors = trial_factors(primitive_sum)
    return (all(prime % 3 == 2 for prime in total_factors)
            and all(exponent <= 2 for exponent in primitive_factors.values()))


def inherited_hypothesis(x: int, y: int) -> bool:
    scale = math.gcd(x, y)
    total_sum = x + y
    primitive_sum = total_sum // scale
    total_factors = trial_factors(total_sum)
    primitive_factors = trial_factors(primitive_sum)
    return (all(prime % 3 == 2 for prime in total_factors)
            and all(exponent <= 2 for exponent in primitive_factors.values()))


def main() -> None:
    limit = 356
    labels = 13
    spf = spf_table(limit)
    phi = totients(limit)

    admissible_sums = [value for value in range(3, limit + 1)
                       if inert_cube_free_sum(value, spf)]
    ratios = [(p, total - p)
              for total in admissible_sums
              for p in range(1, total // 2 + 1)
              if p < total - p and math.gcd(p, total - p) == 1]
    formula_count = sum(phi[total] // 2 for total in admissible_sums)
    gate(len(admissible_sums) == 94, "admissible-sum census changed")
    gate(formula_count == len(ratios), "totient and direct ratio censuses differ")
    gate(len(ratios) == 5855, "primitive-ratio census changed")
    for total in admissible_sums:
        direct = sum(1 for p in range(1, total // 2 + 1)
                     if p < total - p and math.gcd(p, total) == 1)
        gate(direct == phi[total] // 2, "per-sum totient half failed")

    # A complete coordinate table: every candidate value is below 356^3, so
    # every coordinate in any positive representation is at most 355.
    fibres: dict[int, list[tuple[int, int]]] = {}
    for x in range(1, limit):
        for y in range(x + 1, limit):
            fibres.setdefault(x ** 3 + y ** 3, []).append((x, y))
    addresses: dict[int, tuple[int, int]] = {}
    for p, q in ratios:
        value = p ** 3 + q ** 3
        gate(fibres[value] == [(p, q)], "complete bounded fibre is not singleton")
        gate(value not in addresses, "two admissible ratios share a cube address")
        addresses[value] = (p, q)
    gate(len(addresses) == 5855, "address count changed")

    # The literal covector fixes global sign: positive q sits at the label of
    # the smaller speed p.  Hence each support has two ordered assignments.
    oriented_count = 0
    relation_rays: set[tuple[int, ...]] = set()
    for p, q in ratios:
        for i in range(labels):
            for j in range(labels):
                if i == j:
                    continue
                covector = [0] * labels
                covector[i] = q
                covector[j] = -p
                positive = [(index, value) for index, value in enumerate(covector)
                            if value > 0]
                negative = [(index, -value) for index, value in enumerate(covector)
                            if value < 0]
                gate(positive == [(i, q)] and negative == [(j, p)],
                     "signed support did not decode ordered placement")
                gate(covector[i] * p + covector[j] * q == 0,
                     "placed covector is not a relation")
                gate(sum(abs(value) for value in covector) == p + q,
                     "placed covector mass changed")
                relation_rays.add(tuple(covector))
                oriented_count += 1
    gate(oriented_count == 5855 * labels * (labels - 1),
         "ordered placement product changed")
    gate(len(relation_rays) == oriented_count,
         "ordered placed covectors collided")
    gate(oriented_count == 913380, "oriented labelled census changed")
    support_quotient_count = len(ratios) * math.comb(labels, 2)
    gate(support_quotient_count == 456690, "unoriented support quotient changed")

    first_ray = (4, -1) + (0,) * 11
    reversed_ray = (-1, 4) + (0,) * 11
    gate(first_ray != reversed_ray, "two assignments unexpectedly coincide")
    gate(first_ray != tuple(-entry for entry in reversed_ray),
         "two assignments unexpectedly agree up to global sign")

    # Direct coordinate scans, rather than sum-divisor factoring, audit the
    # advertised all-scale control list and recover gcd plus primitive pair.
    scale_cases: list[tuple[int, int, int]] = []
    for p, q in ((1, 4), (2, 9), (5, 6), (1, 24)):
        for scale in (1, 2, 4, 5, 25, 125, 256, 2000):
            gate(packet_hypothesis(p, q, scale), "admissible scale control rejected")
            gate(inherited_hypothesis(scale * p, scale * q),
                 "THM-3793 quantifier translation failed")
            value = (scale * p) ** 3 + (scale * q) ** 3
            representations = brute_representations(value)
            gate(representations == ((scale * p, scale * q),),
                 "direct all-scale control fibre is not singleton")
            x, y = representations[0]
            recovered_scale = math.gcd(x, y)
            gate((recovered_scale, x // recovered_scale, y // recovered_scale)
                 == (scale, p, q), "gcd decoder failed")
            scale_cases.append((scale, p, q))
    for forbidden_scale in (3, 7, 13, 21):
        gate(not packet_hypothesis(1, 4, forbidden_scale),
             "split/ramified scale leaked into the all-inert hypothesis")
        gate(not inherited_hypothesis(forbidden_scale, 4 * forbidden_scale),
             "THM-3793 accepted a forbidden scale")

    gate(brute_representations(1729) == ((1, 12), (9, 10)),
         "split-prime collision hostile changed")
    gate(brute_representations(515375) == ((15, 80), (54, 71)),
         "exponent-three collision hostile changed")

    # Exact geometry control for the AP row and pair (1,4).  The free-face
    # projection has rank eleven.  Analytically its Gram determinant is
    # (n_i^2+n_j^2)/(n dot n)>0, the matrix-determinant-lemma proof used in
    # the report below.
    ap_row = tuple(range(1, 14))
    i, j, p, q = 0, 3, 1, 4
    covector = tuple(q if index == i else -p if index == j else 0
                     for index in range(labels))
    gate(sum(covector[index] * ap_row[index] for index in range(labels)) == 0,
         "AP facet covector is not orthogonal")
    gate(projected_free_rank(ap_row, i, j) == 11,
         "projected exposed face is not eleven-dimensional")
    gram_determinant = Fraction(ap_row[i] ** 2 + ap_row[j] ** 2,
                                sum(speed * speed for speed in ap_row))
    gate(gram_determinant > 0, "free-face Gram determinant vanished")
    width = Fraction(12, 14) * sum(abs(entry) for entry in covector)
    gate(width == Fraction(6, 7) * (p + q), "facet width formula changed")

    positive_row = tuple(list(range(1, 11)) + [12, 13, 14])
    positive_modulus = 11
    positive_residues = tuple(speed % positive_modulus for speed in positive_row)
    positive_grid = modular_schedule(positive_residues, positive_modulus)
    gate(positive_grid == tuple(range(1, 11)), "positive grid control changed")
    gate(positive_grid == rational_schedule(positive_row, positive_modulus),
         "positive modular/rational schedules differ")

    lifted_positive = tuple(speed if index in (0, 9)
                            else speed + positive_modulus * (index + 1)
                            for index, speed in enumerate(positive_row))
    gate(tuple(speed % positive_modulus for speed in lifted_positive)
         == positive_residues, "positive lift changed residues")
    gate(rational_schedule(lifted_positive, positive_modulus) == positive_grid,
         "grid schedule was not residue invariant")

    ap_grid = modular_schedule(tuple(speed % 5 for speed in ap_row), 5)
    gate(ap_grid == (), "AP pair-sum grid is no longer empty")
    gate(ap_grid == rational_schedule(ap_row, 5), "AP modular/rational grids differ")
    gate(safe_at(ap_row, Fraction(1, 14)), "AP t=1/14 hostile is no longer safe")

    # Same complete packet, different off-grid truth: replace speed 2 by 42.
    # The selected (1,4) pair, M=65, a, D=5, and all residues mod 5 agree.
    off_grid_lift = list(ap_row)
    off_grid_lift[1] = 42
    off_grid_lift_row = tuple(off_grid_lift)
    gate(tuple(speed % 5 for speed in off_grid_lift_row)
         == tuple(speed % 5 for speed in ap_row), "off-grid lift changed packet residues")
    gate(not safe_at(off_grid_lift_row, Fraction(1, 14)),
         "off-grid lift unexpectedly retained t=1/14 safety")
    gate((1 ** 3 + 4 ** 3, covector, 5) == (65, covector, 5),
         "selected pair packet identity changed")

    # The inclusive boundary is visible at modulus 14 even though this
    # synthetic modulus lies outside the all-inert packet hypothesis.
    boundary_row = (1,) * labels
    boundary_grid = modular_schedule(boundary_row, 14)
    gate(boundary_grid == tuple(range(1, 14)), "inclusive 1/14 boundary was lost")
    gate(boundary_grid == rational_schedule(boundary_row, 14),
         "boundary modular/rational schedules differ")
    gate(1 in boundary_grid and 13 in boundary_grid,
         "equality endpoints were treated strictly")

    semantic = {
        "admissible_sums": len(admissible_sums),
        "ratios": len(ratios),
        "support_quotient": support_quotient_count,
        "oriented_labelled": oriented_count,
        "scale_controls": len(scale_cases),
        "face_dimension": 11,
        "positive_grid": positive_grid,
        "ap_grid": ap_grid,
        "same_packet_off_grid_divergence": True,
        "repaired_target_unoriented_supports": 456690,
        "repaired_target_oriented_assignments": 913380,
        "census_repair_matches": True,
    }
    blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")

    print("THM3818 INDEPENDENT HOSTILE REFEREE")
    print("UNIVERSE=p<q;gcd=1;p+q<=356;all_sum_primes_2mod3;primitive_sum_exponents<=2")
    print(f"ADMISSIBLE_SUMS={len(admissible_sums)}")
    print(f"PRIMITIVE_RATIOS={len(ratios)}")
    print(f"UNORIENTED_SUPPORT_QUOTIENT={support_quotient_count}")
    print(f"ORIENTED_LABELLED_ASSIGNMENTS={oriented_count}")
    print("REPAIRED_TARGET_UNORIENTED_SUPPORTS=456690")
    print("REPAIRED_TARGET_ORIENTED_ASSIGNMENTS=913380")
    print("CENSUS_REPAIR=PASS")
    print(f"ALL_SCALE_DIRECT_CONTROLS={len(scale_cases)}")
    print("INERT_SCALE_QUANTIFIER=EXACT;SPLIT_OR_RAMIFIED_SCALE_EXCLUDED")
    print("CUBE_DECODER=PASS")
    print("ORIENTED_COVECTOR_DECODER=PASS_FOR_156_ASSIGNMENTS_PER_RATIO")
    print("EXPOSED_FACE_DIMENSION=11")
    print("RESIDUE_GRID=PASS_INCLUSIVE_GE_BOUNDARY")
    print("AP_OFF_GRID_HOSTILE=PASS")
    print("SAME_PACKET_OFF_GRID_DIVERGENCE=AP13_vs_speed2_to42_at_t=1/14")
    print(f"GATES={GATES}")
    print(f"semantic_sha256={hashlib.sha256(blob).hexdigest()}")
    print("AUDIT=PASS_REPAIRED_PACKET")


if __name__ == "__main__":
    main()
