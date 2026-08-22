#!/usr/bin/env python3
"""Finite exact low-weight probe for the prime-even half-twist tail.

This unnumbered companion does four bounded jobs:

1. check exact strict-endpoint formulas for the height-eight A/A and B/E
   tangencies;
2. decide the first uncovered boundary prime p=601 with the independently
   frozen profile DFS;
3. scan every coefficient against normalized A/B/E roots for every prime
   607 <= p <= 997; and
4. build the height-at-most-seven candidate zero graph on one prime in each
   invertible residue class modulo 420.

The fourth calculation is conditional telemetry for the still-open uniform
height-at-most-seven lemma.  Nothing here extrapolates the finite scans.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.py"
EXPECTED_BASE_SHA256 = "ed54b01f9bf155643b6407c6af8ee6f15c7a099f8ac3b60399dd49b496fb1d12"
EXPECTED_SEMANTIC_SHA256 = "17411c3f578f369ca6ab424db0db6730b547db5deb598e87fa341da49dd6473b"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def primes_between(lower: int, upper: int) -> tuple[int, ...]:
    return tuple(value for value in range(max(2, lower), upper + 1) if is_prime(value))


def coefficient_type(coefficient: int) -> str:
    if coefficient % 2:
        return "A"
    if coefficient % 4:
        return "B"
    return "E"


def danger(prime: int, coefficient: int, odd_sheet: int) -> bool:
    period = 4 * prime
    phase = coefficient * odd_sheet % period
    return 14 * min(phase, period - phase) < period


def danger_mask(prime: int, coefficient: int) -> int:
    modulus = 2 * prime
    return sum(
        1 << sheet
        for sheet in range(modulus)
        if danger(prime, coefficient, 2 * sheet + 1)
    )


def literal_weight(prime: int, first: int, second: int) -> int:
    return (danger_mask(prime, first) & danger_mask(prime, second)).bit_count()


def count_residue(lower: int, upper: int, residue: int, modulus: int) -> int:
    if lower > upper:
        return 0
    residue %= modulus
    return (
        (upper - residue) // modulus
        - ((lower - 1) - residue) // modulus
    )


def aa_tangent_count(prime: int) -> int:
    """One-sign count for 3x-5y=+/-2p in the two A intervals."""
    lower = 4 * prime // 21 + 1
    upper = (2 * prime - 1) // 7
    residues = tuple(
        residue
        for residue in range(10)
        if residue % 2 and (3 * residue - 2 * prime) % 5 == 0
    )
    require(len(residues) == 1, ("A/A residue", prime, residues))
    return count_residue(lower, upper, residues[0], 10)


def be_35_tangent_count(prime: int) -> int:
    """One-sign count for 3x-10y=+/-p (raw B-to-E ratio 3/5)."""
    lower = 2 * prime // 21 + 1
    upper = (prime - 1) // 7
    residue = 7 * prime % 10
    require((3 * residue - prime) % 10 == 0, ("B/E 3/5 residue", prime))
    return count_residue(lower, upper, residue, 10)


def be_53_tangent_count(prime: int) -> int:
    """One-sign count for 5x-6y=+/-p (raw B-to-E ratio 5/3)."""
    lower = 4 * prime // 35 + 1
    upper = (prime - 1) // 7
    residues = tuple(
        residue
        for residue in range(6)
        if residue % 2 and (5 * residue - prime) % 6 == 0
    )
    require(len(residues) == 1, ("B/E 5/3 residue", prime, residues))
    return count_residue(lower, upper, residues[0], 6)


def coefficient_solutions(prime: int, root: int, numerator: int,
                          denominator: int) -> tuple[int, ...]:
    residue = numerator * root * pow(denominator, -1, prime) % prime
    return tuple(
        coefficient
        for coefficient in (residue, residue + prime)
        if 0 < coefficient < 2 * prime and coefficient != prime
    )


H8_FRACTIONS = ((3, 5), (-3, 5), (5, 3), (-5, 3))


def h8_formula_control(prime: int):
    cache = {}

    def weight(first, second):
        key = tuple(sorted((first, second)))
        if key not in cache:
            cache[key] = literal_weight(prime, first, second)
        return cache[key]

    aa_edges = {}
    for numerator, denominator in H8_FRACTIONS:
        for candidate in coefficient_solutions(prime, 1, numerator, denominator):
            if (coefficient_type(candidate) == "A"
                    and (denominator * candidate - numerator) % (4 * prime) == 2 * prime):
                aa_edges[candidate] = weight(1, candidate)
    require(aa_edges, ("missing A/A H8 tangent", prime))
    require(set(aa_edges.values()) == {2 * aa_tangent_count(prime)},
            ("A/A formula", prime, aa_edges, aa_tangent_count(prime)))

    be_edges = {}
    for numerator, denominator in H8_FRACTIONS:
        for candidate in coefficient_solutions(prime, 2, numerator, denominator):
            if coefficient_type(candidate) != "E":
                continue
            expected = (
                4 * be_35_tangent_count(prime)
                if abs(numerator) == 3
                else 4 * be_53_tangent_count(prime)
            )
            be_edges[candidate] = (weight(2, candidate), expected)
    require(be_edges and all(actual == expected for actual, expected in be_edges.values()),
            ("B/E formula", prime, be_edges))

    eb_edges = {}
    for numerator, denominator in H8_FRACTIONS:
        for candidate in coefficient_solutions(prime, 4, numerator, denominator):
            if coefficient_type(candidate) != "B":
                continue
            # Reversing the ordered edge reciprocates the raw coefficient ratio.
            expected = (
                4 * be_53_tangent_count(prime)
                if abs(numerator) == 3
                else 4 * be_35_tangent_count(prime)
            )
            eb_edges[candidate] = (weight(4, candidate), expected)
    require(eb_edges and all(actual == expected for actual, expected in eb_edges.values()),
            ("E/B formula", prime, eb_edges))

    return (
        prime,
        2 * aa_tangent_count(prime),
        4 * be_35_tangent_count(prime),
        4 * be_53_tangent_count(prime),
        len(aa_edges),
        len(be_edges),
        len(eb_edges),
    )


def rational_residues(prime: int, height: int) -> frozenset[int]:
    values = set()
    for numerator in range(1, height):
        for denominator in range(1, height - numerator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            value = numerator * pow(denominator, -1, prime) % prime
            values.add(value)
            values.add(-value % prime)
    return frozenset(values)


def scan_normalized_roots(prime: int):
    modulus = 2 * prime
    period = 4 * prime
    roots = (1, 2, 4)
    root_sheets = {
        root: tuple(
            2 * sheet + 1
            for sheet in range(modulus)
            if danger(prime, root, 2 * sheet + 1)
        )
        for root in roots
    }
    height_seven = rational_residues(prime, 7)
    zero_degrees = []
    minima = []
    low_edges = []
    height_misses = []
    for root in roots:
        zero_count = 0
        positive_minimum = modulus + 1
        inverse = pow(root, -1, prime)
        for candidate in range(1, modulus):
            if candidate in {prime, root}:
                continue
            weight = 0
            for odd_sheet in root_sheets[root]:
                phase = candidate * odd_sheet % period
                weight += 14 * min(phase, period - phase) < period
            if weight == 0:
                zero_count += 1
                ratio = candidate * inverse % prime
                if ratio not in height_seven:
                    height_misses.append((root, candidate, ratio))
            else:
                positive_minimum = min(positive_minimum, weight)
                if weight <= 8:
                    low_edges.append((root, candidate, weight))
        zero_degrees.append(zero_count)
        minima.append(positive_minimum)
    return tuple(zero_degrees), tuple(minima), tuple(low_edges), tuple(height_misses)


def height_seven_fractions() -> tuple[tuple[int, int], ...]:
    values = []
    for numerator in range(1, 7):
        for denominator in range(1, 8 - numerator):
            if gcd(numerator, denominator) == 1:
                values.extend(((numerator, denominator), (-numerator, denominator)))
    return tuple(values)


HEIGHT_SEVEN_FRACTIONS = height_seven_fractions()


def candidate_coefficients(prime: int, root: int) -> tuple[int, ...]:
    coefficients = {root}
    for numerator, denominator in HEIGHT_SEVEN_FRACTIONS:
        coefficients.update(coefficient_solutions(
            prime, root, numerator, denominator
        ))
    return tuple(sorted(coefficients))


def maximum_cliques_with_root(prime: int, root: int):
    coefficients = candidate_coefficients(prime, root)
    masks = {coefficient: danger_mask(prime, coefficient) for coefficient in coefficients}
    neighbors = tuple(
        coefficient
        for coefficient in coefficients
        if coefficient != root and not (masks[root] & masks[coefficient])
    )
    best_size = 0
    best = set()

    def extend(chosen, candidates):
        nonlocal best_size, best
        if len(chosen) + len(candidates) < best_size:
            return
        if len(chosen) > best_size:
            best_size = len(chosen)
            best = {tuple(sorted(chosen))}
        elif len(chosen) == best_size:
            best.add(tuple(sorted(chosen)))

        remaining = candidates
        while remaining:
            if len(chosen) + len(remaining) < best_size:
                return
            vertex = remaining[0]
            tail = remaining[1:]
            next_candidates = tuple(
                other for other in tail if not (masks[vertex] & masks[other])
            )
            extend(chosen + (vertex,), next_candidates)
            remaining = tail

    extend((root,), neighbors)
    best = {clique for clique in best if len(clique) == best_size}
    return len(neighbors), best_size, best


def first_prime_in_class(residue: int, modulus: int, lower: int) -> int:
    candidate = lower + (residue - lower) % modulus
    while not is_prime(candidate):
        candidate += modulus
    return candidate


def mod_420_candidate_atlas():
    residues = tuple(residue for residue in range(420) if gcd(residue, 420) == 1)
    reports = []
    for residue in residues:
        prime = first_prime_in_class(residue, 420, 607)
        root_reports = []
        maximum = 0
        maximum_cliques = set()
        for root in (1, 2, 4):
            degree, clique_size, cliques = maximum_cliques_with_root(prime, root)
            root_reports.append((root, degree, clique_size, len(cliques)))
            if clique_size > maximum:
                maximum = clique_size
                maximum_cliques = set(cliques)
            elif clique_size == maximum:
                maximum_cliques.update(cliques)

        pairable = tuple(
            all(2 * prime - vertex in clique for vertex in clique)
            for clique in maximum_cliques
        )
        type_profiles = {
            tuple(sorted(coefficient_type(vertex) for vertex in clique))
            for clique in maximum_cliques
        }
        require(maximum == 6, ("candidate clique", residue, prime, maximum))
        require(pairable and all(pairable), ("unpaired candidate clique", residue, prime))
        require(type_profiles == {("A", "A", "A", "A", "B", "E")},
                ("candidate clique types", residue, prime, type_profiles))
        reports.append((residue, prime, tuple(root_reports), len(maximum_cliques)))
    return tuple(reports)


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")
    require(lf_sha256(BASE_PATH) == EXPECTED_BASE_SHA256, "base companion changed")

    # The exact count functions have affine recurrences on the common modulus
    # 210.  These finite residue checks certify the stated global thresholds.
    for prime in range(1, 421, 2):
        require(aa_tangent_count(prime + 210) == aa_tangent_count(prime) + 2,
                ("A/A recurrence", prime))
        require(be_35_tangent_count(prime + 210) == be_35_tangent_count(prime) + 1,
                ("B/E 3/5 recurrence", prime))
        require(be_53_tangent_count(prime + 210) == be_53_tangent_count(prime) + 1,
                ("B/E 5/3 recurrence", prime))
    require(max(p for p in range(1, 1000, 2) if aa_tangent_count(p) <= 4) == 521,
            "A/A threshold")
    require(max(p for p in range(1, 1000, 2) if be_35_tangent_count(p) <= 2) == 601,
            "B/E 3/5 threshold")
    require(max(p for p in range(1, 1000, 2) if be_53_tangent_count(p) <= 2) == 609,
            "B/E 5/3 odd threshold")
    require(max(p for p in primes_between(3, 1000) if be_53_tangent_count(p) <= 2) == 593,
            "B/E 5/3 prime threshold")
    require(all(
        min(4 * be_35_tangent_count(p), 4 * be_53_tangent_count(p)) >= 12
        for p in primes_between(607, 2000)
    ), "prime H8 threshold")

    formula_controls = tuple(
        h8_formula_control(prime)
        for prime in (401, 521, 547, 571, 593, 601, 607, 1009)
    )

    base = runpy.run_path(str(BASE_PATH))
    p601 = base["decide_prime"](601)
    require(p601[0] is None, ("p=601 cover", p601[0]))
    p601_tally = p601[3]
    p601_report = (p601_tally.profiles, p601_tally.nodes, p601_tally.branches)
    require(p601_report == (15, 65, 50), ("p=601 tally", p601_report))

    scan_primes = primes_between(607, 997)
    scan_rows = []
    for prime in scan_primes:
        row = scan_normalized_roots(prime)
        require(row[0] == (19, 31, 23), ("zero degrees", prime, row[0]))
        require(not row[2], ("low positive edge", prime, row[2]))
        require(not row[3], ("height-seven zero miss", prime, row[3]))
        scan_rows.append((prime, row[0], row[1]))
    scan_rows = tuple(scan_rows)
    minimum_histogram = tuple(sorted(Counter(row[2] for row in scan_rows).items()))
    require(min(value for row in scan_rows for value in row[2]) == 12,
            "scan minimum")

    mod_420 = mod_420_candidate_atlas()
    require(len(mod_420) == 96, ("mod-420 classes", len(mod_420)))
    require(max(row[1] for row in mod_420) == 3319, "mod-420 prime boundary")
    require(all(
        tuple(root_row[1] for root_row in row[2]) == (19, 31, 23)
        for row in mod_420
    ), "mod-420 candidate degrees")
    require(all(row[3] == 4 for row in mod_420), "mod-420 equality count")

    semantic_payload = (
        formula_controls,
        p601_report,
        scan_rows,
        minimum_histogram,
        mod_420,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("prime-even half-twist low-weight tail probe")
    print("status=FINITE-EXACT controls; uniform height<=7 lemma OPEN")
    print("H8 formulas:")
    print("  w_AA=2*#{odd x: 4p/21<x<2p/7, 3x=2p (mod 5)}")
    print("  w_BE_3_5=4*#{x: 2p/21<x<p/7, 3x=p (mod 10)}")
    print("  w_BE_5_3=4*#{odd x: 4p/35<x<p/7, 5x=p (mod 6)}")
    print("H8 thresholds: last_AA_w<=8_odd=521 last_BE_3_5_w<=8_odd=601 "
          "last_BE_5_3_w<=8_odd=609 last_BE_5_3_w<=8_prime=593")
    print(f"formula_controls={formula_controls}")
    print(f"p601_exact_negative profiles_nodes_branches={p601_report}")
    print(f"full_root_scan primes={len(scan_primes)} range=({scan_primes[0]},{scan_primes[-1]}) "
          f"global_min_positive=12 zero_degrees=(19,31,23) low_edges=0 height_misses=0")
    print(f"full_root_scan_minimum_histogram={minimum_histogram}")
    print("mod420_height7_candidate_atlas: classes=96 first_prime_range=(607,3319) "
          "degrees=(19,31,23) clique=6 equality_cliques_per_class=4 "
          "all_equalities=three_complementary_pairs types=A4BE")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
