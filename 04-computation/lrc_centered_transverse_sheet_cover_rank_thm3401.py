#!/usr/bin/env python3
"""Exact companion for THM-3401's centered transverse sheet-cover rank.

The verifier is deliberately self-contained and standard-library only.  It
computes every danger set from the strict integer inequality, classifies all
residue-speed owners for 15 <= q <= 28, checks the explicit unit/prime-kernel
construction, and independently solves each finite set-cover problem by
exhausting combinations of the distinct literal coverage masks.  Duplicate
coverage masks may be discarded because using the same set twice never
enlarges a union.

The q=14 strict-endpoint failure, q=28 endpoint exclusion, and q=29 cycle
boundary are checked from the same literal danger-set routine.  Runtime gates
survive ``python -O``.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


EXPECTED_RANKS = (6, 5, 8, 5, 9, 6, 8, 7, 11, 6, 11, 8, 10, 8)
EXPECTED_SEMANTIC_DIGEST = "41c31d22b1a1b6fa3637d129fad3a311565d21bf45ad9f322bd550d595c157cd"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def prime_divisors(n):
    answer = []
    divisor = 2
    remaining = n
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            answer.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        answer.append(remaining)
    return tuple(answer)


def phi(n):
    answer = n
    for prime in prime_divisors(n):
        answer = answer // prime * (prime - 1)
    return answer


def rank_formula(q):
    return phi(q) // 2 + sum(prime < q for prime in prime_divisors(q))


def cyclic_distance_numerator(residue, q):
    residue %= q
    return min(residue, q - residue)


def danger_set(q, speed):
    """Return D_(q,speed)(0) using 14*distance < q, with no fractions."""
    require(q >= 2 and speed > 0 and speed % q != 0, (q, speed))
    return frozenset(
        sheet
        for sheet in range(q)
        if 14 * cyclic_distance_numerator(speed * sheet, q) < q
    )


def quotient_block_set(q, speed):
    """Return pi_m^(-1)(a^(-1) R_m) in the general gcd quotient model."""
    g = gcd(speed, q)
    m = q // g
    if m == 1:
        allowed = frozenset((0,))
    else:
        a = speed // g % m
        inverse = pow(a, -1, m)
        short_block = frozenset(
            residue
            for residue in range(m)
            if 14 * cyclic_distance_numerator(residue, m) < m
        )
        allowed = frozenset(inverse * residue % m for residue in short_block)
    return frozenset(sheet for sheet in range(q) if sheet % m in allowed)


def predicted_set(q, speed):
    g = gcd(speed, q)
    m = q // g
    if g == 1:
        inverse = pow(speed, -1, q)
        return frozenset((0, inverse, (-inverse) % q))
    return frozenset(range(0, q, m))


def sign_representative(residue, q):
    residue %= q
    return min(residue, (-residue) % q)


def canonical_construction(q):
    unit_speeds = {
        sign_representative(pow(sheet, -1, q), q)
        for sheet in range(1, q // 2 + 1)
        if gcd(sheet, q) == 1
    }
    kernel_speeds = {q // prime for prime in prime_divisors(q) if prime < q}
    require(unit_speeds.isdisjoint(kernel_speeds), (q, unit_speeds, kernel_speeds))
    return tuple(sorted(unit_speeds | kernel_speeds))


def mask_of(sheets):
    answer = 0
    for sheet in sheets:
        answer |= 1 << sheet
    return answer


def literal_coverage_types(q):
    """Distinct literal masks from the complete transverse residue universe."""
    by_mask = {}
    for speed in range(1, q):
        mask = mask_of(danger_set(q, speed))
        if mask not in by_mask or speed < by_mask[mask]:
            by_mask[mask] = speed
    return tuple(sorted((speed, mask) for mask, speed in by_mask.items()))


def minimum_set_cover_search(q):
    """Exact combination search, independent of the classification formula."""
    coverage_types = literal_coverage_types(q)
    full = (1 << q) - 1
    tested = 0
    for rank in range(1, len(coverage_types) + 1):
        for chosen in combinations(coverage_types, rank):
            tested += 1
            union = 0
            for _, mask in chosen:
                union |= mask
            if union == full:
                return rank, tuple(speed for speed, _ in chosen), len(coverage_types), tested
    return None, (), len(coverage_types), tested


def unit_sign_pairs(q):
    return tuple(
        (sheet, (-sheet) % q)
        for sheet in range(1, q // 2 + 1)
        if gcd(sheet, q) == 1
    )


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    rows = []
    quotient_records = []
    quotient_block_checks = 0
    classification_checks = 0
    construction_checks = 0
    total_search_subsets = 0

    for q in range(14, 30):
        for speed in range(1, q):
            literal = danger_set(q, speed)
            quotient = quotient_block_set(q, speed)
            require(literal == quotient, ("quotient block", q, speed, literal, quotient))
            quotient_records.append((q, speed, tuple(sorted(quotient))))
            quotient_block_checks += 1
    require(quotient_block_checks == 328, quotient_block_checks)

    for q, expected_rank in zip(range(15, 29), EXPECTED_RANKS):
        require(rank_formula(q) == expected_rank, (q, rank_formula(q), expected_rank))

        for speed in range(1, q):
            literal = danger_set(q, speed)
            predicted = predicted_set(q, speed)
            require(literal == predicted, ("classification", q, speed, literal, predicted))
            g = gcd(speed, q)
            m = q // g
            if g == 1:
                require(literal - {0} <= {sheet for pair in unit_sign_pairs(q) for sheet in pair}, (q, speed))
            else:
                require(2 <= m <= 14, ("nonunit quotient", q, speed, m))
                require(literal == frozenset(range(0, q, m)), (q, speed, m, literal))
            classification_checks += 1

        construction = canonical_construction(q)
        construction_union = frozenset().union(*(danger_set(q, speed) for speed in construction))
        require(construction_union == frozenset(range(q)), ("construction", q, construction))
        require(len(construction) == expected_rank, ("construction rank", q, construction))
        require(
            all(
                frozenset().union(
                    *(danger_set(q, other) for other in construction if other != speed)
                )
                != frozenset(range(q))
                for speed in construction
            ),
            ("construction irredundancy", q, construction),
        )
        construction_checks += len(construction) + 1

        search_rank, search_witness, type_count, tested = minimum_set_cover_search(q)
        require(search_rank == expected_rank, ("set-cover rank", q, search_rank, expected_rank))
        require(
            frozenset().union(*(danger_set(q, speed) for speed in search_witness))
            == frozenset(range(q)),
            ("set-cover witness", q, search_witness),
        )
        total_search_subsets += tested

        rows.append(
            (
                q,
                phi(q),
                prime_divisors(q),
                expected_rank,
                type_count,
                tested,
                construction,
                search_witness,
            )
        )

    q14_units = tuple(speed for speed in range(1, 14) if gcd(speed, 14) == 1)
    require(all(danger_set(14, speed) == frozenset((0,)) for speed in q14_units), q14_units)
    q14_union = frozenset().union(*(danger_set(14, speed) for speed in range(1, 14)))
    q14_nonunits = frozenset(sheet for sheet in range(14) if gcd(sheet, 14) != 1)
    q14_missing = tuple(sorted(set(range(14)) - q14_union))
    require(q14_union == q14_nonunits, (q14_union, q14_nonunits))
    require(q14_missing == q14_units, (q14_missing, q14_units))
    q14_search = minimum_set_cover_search(14)
    require(q14_search[0] is None, q14_search)

    q28_endpoint_checks = 0
    for speed in range(1, 28):
        if gcd(speed, 28) != 1:
            continue
        inverse = pow(speed, -1, 28)
        for signed_two in (2, -2):
            sheet = signed_two * inverse % 28
            require(14 * cyclic_distance_numerator(speed * sheet, 28) == 28, (speed, sheet))
            require(sheet not in danger_set(28, speed), (speed, sheet))
            q28_endpoint_checks += 1
    require(q28_endpoint_checks == 24, q28_endpoint_checks)

    for speed in range(1, 29):
        inverse = pow(speed, -1, 29)
        expected = frozenset(
            (0, inverse, (-inverse) % 29, 2 * inverse % 29, -2 * inverse % 29)
        )
        require(danger_set(29, speed) == expected, ("q29 block", speed))

    q29_edges = set()
    for speed in range(1, 29):
        inverse_class = sign_representative(pow(speed, -1, 29), 29)
        doubled_class = sign_representative(2 * inverse_class, 29)
        q29_edges.add(tuple(sorted((inverse_class, doubled_class))))
    q29_edges = tuple(sorted(q29_edges))
    require(len(q29_edges) == 14, q29_edges)

    adjacency = {vertex: set() for vertex in range(1, 15)}
    for left, right in q29_edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    require(all(len(neighbours) == 2 for neighbours in adjacency.values()), adjacency)
    reached = {1}
    frontier = [1]
    while frontier:
        vertex = frontier.pop()
        for neighbour in adjacency[vertex]:
            if neighbour not in reached:
                reached.add(neighbour)
                frontier.append(neighbour)
    require(reached == set(range(1, 15)), reached)

    q29_orbit = []
    vertex = 1
    while vertex not in q29_orbit:
        q29_orbit.append(vertex)
        vertex = sign_representative(2 * vertex, 29)
    require(vertex == 1 and len(q29_orbit) == 14, (q29_orbit, vertex))
    require(pow(2, 14, 29) == 28, pow(2, 14, 29))
    require(all(pow(2, exponent, 29) not in (1, 28) for exponent in range(1, 14)), "q29 order")

    q29_search = minimum_set_cover_search(29)
    require(q29_search[0] == 7, q29_search)
    require(rank_formula(29) == 14, rank_formula(29))

    edge_cover_rank = None
    edge_cover_witness = ()
    edge_cover_tests = 0
    for rank in range(1, len(q29_edges) + 1):
        for chosen in combinations(q29_edges, rank):
            edge_cover_tests += 1
            vertices = {vertex for edge in chosen for vertex in edge}
            if vertices == set(range(1, 15)):
                edge_cover_rank = rank
                edge_cover_witness = chosen
                break
        if edge_cover_rank is not None:
            break
    require(edge_cover_rank == 7, (edge_cover_rank, edge_cover_witness))

    semantic = ExactDigest()
    semantic.add(("quotient blocks", tuple(quotient_records)))
    semantic.add(("rows", tuple(rows)))
    semantic.add(("q14", q14_units, tuple(sorted(q14_union)), q14_missing, q14_search))
    semantic.add(("q28", q28_endpoint_checks))
    semantic.add(
        (
            "q29",
            q29_edges,
            tuple(q29_orbit),
            q29_search,
            edge_cover_rank,
            edge_cover_witness,
            edge_cover_tests,
        )
    )
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("THM-3401 CENTERED TRANSVERSE SHEET-COVER RANK EXACT COMPANION")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=PROVED analytic classification plus VERIFIED-EXACT independent literal set-cover search")
    print("universe=all q=15..28 and all transverse speed residues 1..q-1 at source centre zero")
    print(f"rank_table={tuple((row[0], row[3]) for row in rows)}")
    print(
        f"quotient_block_checks_q14_to_q29={quotient_block_checks};"
        f"classification_checks={classification_checks};construction_checks={construction_checks}"
    )
    print(f"independent_search_subsets_tested={total_search_subsets}")
    for row in rows:
        q, phi_q, primes, rank, type_count, tested, construction, witness = row
        print(
            f"q={q};phi={phi_q};prime_divisors={primes};rank={rank};"
            f"literal_types={type_count};search_subsets={tested};"
            f"construction={construction};search_witness={witness}"
        )
    print(
        f"q14_strict_hostile=unit_owners_cover_only_zero;all_transverse_union="
        f"{tuple(sorted(q14_union))};missing_units={q14_missing};cover_rank=NONE"
    )
    print(f"q28_endpoint_hostile=+-2 strict_equalities_excluded;checks={q28_endpoint_checks}")
    print(
        f"q29_boundary=unit_blocks_have_two_sign_pairs;times2_orbit={tuple(q29_orbit)};"
        f"edge_count={len(q29_edges)};cycle_min_edge_cover={edge_cover_rank};"
        f"literal_cover_rank={q29_search[0]};old_formula={rank_formula(29)}"
    )
    print(f"q29_edge_cover_witness={edge_cover_witness};literal_speed_witness={q29_search[1]}")
    print("q15_canonical_construction=(1,2,3,4,5,7);rank=6;zero_cochain_triphase_edge_explained")
    print("scope=centered_zero_cochain_transverse_cover_only;nonzero_centres_and_LRC14_remain_open")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
