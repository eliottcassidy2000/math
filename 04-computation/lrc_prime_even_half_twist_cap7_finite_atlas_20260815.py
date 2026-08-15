#!/usr/bin/env python3
"""Finite-exact prime-even literal half-twist cap-seven atlas.

For every target-free odd prime p <= 599, this companion constructs the
literal masks on Q=2p directly and decides whether at most seven distinct
transverse masks cover Z/QZ.  The search is organized by the three dyadic
coefficient types and uses an exact total-overlap budget.  It is a bounded
calculation, not an all-prime classification.

The final section records finite weighted-disjointness telemetry on one large
prime in every nonzero residue class modulo 14.  That telemetry motivates a
bounded-rational tail argument, but no such tail argument is assumed here.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
from math import gcd, isqrt
from pathlib import Path


LIMIT = 599
LOWER_BASE_PRIMES = frozenset({5, 11, 13, 23, 29})
POSITIVE_CONTROLS = {
    7: (1, 3, 4, 5, 9, 11, 13),
    19: (1, 9, 17, 20, 21, 29, 37),
}
TAIL_CONTROLS = (547, 521, 593, 541, 599, 587)
EXPECTED_SEMANTIC_SHA256 = "19f846db72803c683f608149aac8c7da2015eaca6e1e3524d9954eaeef0fa826"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def primes_through(limit: int) -> tuple[int, ...]:
    values = []
    for value in range(2, limit + 1):
        if all(value % divisor for divisor in range(2, isqrt(value) + 1)):
            values.append(value)
    return tuple(values)


def coefficient_type(coefficient: int) -> str:
    if coefficient % 2:
        return "A"
    if coefficient % 4:
        return "B"
    return "E"


def danger_mask(modulus: int, coefficient: int) -> int:
    """Return the strict radius-1/14 mask as a Python integer bitset."""
    period = 2 * modulus
    mask = 0
    for sheet in range(modulus):
        phase = coefficient * (2 * sheet + 1) % period
        distance = min(phase, period - phase)
        if 14 * distance < period:
            mask |= 1 << sheet
    return mask


def expected_type_size(prime: int, owner_type: str) -> int:
    if prime == 7:
        return 0 if owner_type == "B" else 2
    quotient, residue = divmod(prime, 14)
    require(residue in {1, 3, 5, 9, 11, 13}, (prime, residue))
    offsets = {
        "A": {1: 0, 3: 0, 5: 2, 9: 2, 11: 4, 13: 4},
        "B": {1: 0, 3: 0, 5: 0, 9: 4, 11: 4, 13: 4},
        "E": {1: 2, 3: 2, 5: 2, 9: 2, 11: 2, 13: 2},
    }
    return 4 * quotient + offsets[owner_type][residue]


def canonical_coefficients(prime: int) -> tuple[int, ...]:
    modulus = 2 * prime
    # Sign gives representatives 1,...,Q.  The transverse endpoint Q and the
    # identically empty coefficient p can both be discarded.
    return tuple(
        coefficient
        for coefficient in range(1, modulus)
        if coefficient != prime
    )


def build_masks(prime: int):
    modulus = 2 * prime
    coefficients = canonical_coefficients(prime)
    masks = {
        coefficient: danger_mask(modulus, coefficient)
        for coefficient in coefficients
    }
    groups = {
        owner_type: tuple(
            coefficient
            for coefficient in coefficients
            if coefficient_type(coefficient) == owner_type
        )
        for owner_type in "ABE"
    }

    for owner_type in "ABE":
        observed = {masks[coefficient].bit_count() for coefficient in groups[owner_type]}
        require(observed == {expected_type_size(prime, owner_type)},
                ("type size", prime, owner_type, observed))
    return coefficients, masks, groups


@dataclass
class SearchTally:
    profiles: int = 0
    nodes: int = 0
    branches: int = 0


def profiles_at_most_seven(groups, masks, modulus: int):
    sizes = {
        owner_type: masks[groups[owner_type][0]].bit_count()
        for owner_type in "ABE"
    }
    for owner_count in range(1, 8):
        for a_count in range(owner_count + 1):
            for b_count in range(owner_count - a_count + 1):
                o_count = owner_count - a_count - b_count
                counts = {"A": o_count, "B": b_count, "E": a_count}
                if any(counts[t] > len(groups[t]) for t in "ABE"):
                    continue
                mass = sum(counts[t] * sizes[t] for t in "ABE")
                if mass >= modulus:
                    yield counts, mass - modulus


def search_profile(modulus: int, masks, groups, counts, omega: int,
                   tally: SearchTally):
    """Search one exact type profile after an odd-unit orbit normalization."""
    roots = {"A": 1, "B": 2, "E": 4}
    available_types = tuple(t for t in "ABE" if counts[t])

    # Every family with this profile is odd-unit equivalent to one containing
    # the canonical root of any selected nonempty type.  Choose the root whose
    # low-weight neighborhood is smallest; this changes speed, not coverage.
    root_type = min(
        available_types,
        key=lambda t: (
            sum(
                1
                for candidate in canonical_coefficients(modulus // 2)
                if candidate != roots[t]
                and (masks[candidate] & masks[roots[t]]).bit_count() <= omega
            ),
            t,
        ),
    )
    root = roots[root_type]
    require(root in groups[root_type], ("missing root", modulus, root_type))

    remaining = dict(counts)
    remaining[root_type] -= 1
    last = {"A": 0, "B": 0, "E": 0}
    last[root_type] = root
    root_mask = masks[root]

    def visit(selected, union_mask: int, overlap_used: int, rem, lower):
        tally.nodes += 1
        slots = sum(rem.values())
        if slots == 0:
            if union_mask.bit_count() == modulus:
                require(overlap_used == omega,
                        ("budget iff", modulus, selected, overlap_used, omega))
                return selected
            return None

        remaining_mass = sum(
            rem[t] * masks[groups[t][0]].bit_count() for t in "ABE"
        )
        if union_mask.bit_count() + remaining_mass < modulus:
            return None

        options_by_type = {}
        budget_left = omega - overlap_used
        for owner_type in "ABE":
            if rem[owner_type] == 0:
                continue
            options = tuple(
                coefficient
                for coefficient in groups[owner_type]
                if coefficient > lower[owner_type]
                and coefficient not in selected
                and (masks[coefficient] & union_mask).bit_count() <= budget_left
            )
            if len(options) < rem[owner_type]:
                return None
            options_by_type[owner_type] = options

        owner_type = min(
            options_by_type,
            key=lambda t: (len(options_by_type[t]), -rem[t], t),
        )
        options = options_by_type[owner_type]
        need = rem[owner_type]
        for index, coefficient in enumerate(options):
            if len(options) - index < need:
                break
            tally.branches += 1
            mask = masks[coefficient]
            increment = (mask & union_mask).bit_count()
            new_rem = dict(rem)
            new_rem[owner_type] -= 1
            new_lower = dict(lower)
            new_lower[owner_type] = coefficient
            answer = visit(
                selected + (coefficient,),
                union_mask | mask,
                overlap_used + increment,
                new_rem,
                new_lower,
            )
            if answer is not None:
                return answer
        return None

    return visit((root,), root_mask, 0, remaining, last)


def decide_prime(prime: int):
    modulus = 2 * prime
    coefficients, masks, groups = build_masks(prime)
    full = (1 << modulus) - 1
    tally = SearchTally()
    found = None
    found_profile = None
    found_omega = None
    for counts, omega in profiles_at_most_seven(groups, masks, modulus):
        tally.profiles += 1
        answer = search_profile(modulus, masks, groups, counts, omega, tally)
        if answer is not None:
            union = 0
            mass = 0
            for coefficient in answer:
                union |= masks[coefficient]
                mass += masks[coefficient].bit_count()
            require(union == full, ("false cover", prime, answer))
            require(mass - union.bit_count() == omega,
                    ("false overlap", prime, answer))
            found = tuple(sorted(answer))
            found_profile = tuple(counts[t] for t in "ABE")
            found_omega = omega
            break
    return found, found_profile, found_omega, tally, masks, groups


def verify_positive_control(prime: int, witness: tuple[int, ...], masks) -> tuple:
    modulus = 2 * prime
    union = 0
    multiplicities = [0] * modulus
    for coefficient in witness:
        mask = masks[coefficient]
        union |= mask
        for sheet in range(modulus):
            multiplicities[sheet] += bool(mask & (1 << sheet))
    require(union == (1 << modulus) - 1, ("control does not cover", prime))
    types = Counter(coefficient_type(r) for r in witness)
    sizes = tuple(masks[r].bit_count() for r in witness)
    orders = tuple(modulus // gcd(modulus, r) for r in witness)
    return (
        prime,
        witness,
        tuple(sorted(types.items())),
        sizes,
        orders,
        tuple(sorted(Counter(multiplicities).items())),
    )


def maximum_disjoint_clique(vertices, masks, allowed_type=None):
    roots = (1,) if allowed_type == "A" else (1, 2, 4)
    best = ()

    def extend(chosen, candidates):
        nonlocal best
        if len(chosen) + len(candidates) <= len(best):
            return
        if len(chosen) > len(best):
            best = chosen
        remaining = candidates
        while remaining:
            if len(chosen) + len(remaining) <= len(best):
                return
            vertex = remaining[0]
            tail = remaining[1:]
            neighbors = tuple(
                other for other in tail if not (masks[vertex] & masks[other])
            )
            extend(chosen + (vertex,), neighbors)
            remaining = tail

    for root in roots:
        if allowed_type is not None and coefficient_type(root) != allowed_type:
            continue
        neighbors = tuple(
            vertex
            for vertex in vertices
            if vertex != root
            and (allowed_type is None or coefficient_type(vertex) == allowed_type)
            and not (masks[root] & masks[vertex])
        )
        extend((root,), neighbors)
    return best


def rational_height(residue: int, prime: int, maximum: int = 7):
    """Find a signed reduced a/b of height a+b <= maximum, if present."""
    possibilities = []
    for numerator in range(1, maximum):
        for denominator in range(1, maximum - numerator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            inverse = pow(denominator, -1, prime)
            value = numerator * inverse % prime
            for sign in (1, -1):
                if sign * value % prime == residue:
                    possibilities.append((numerator + denominator, sign, numerator, denominator))
    if not possibilities:
        return None
    return min(possibilities)


def tail_telemetry(prime: int):
    modulus = 2 * prime
    vertices, masks, groups = build_masks(prime)
    roots = {"A": 1, "B": 2, "E": 4}
    zero_degrees = []
    positive_minima = []
    low_weight_histograms = []
    atlas_heights = []
    for owner_type in "ABE":
        root = roots[owner_type]
        weights = tuple(
            (masks[root] & masks[vertex]).bit_count()
            for vertex in vertices if vertex != root
        )
        zero_neighbors = tuple(
            vertex for vertex in vertices
            if vertex != root and not (masks[root] & masks[vertex])
        )
        zero_degrees.append((owner_type, len(zero_neighbors)))
        positive_minima.append((owner_type, min(weight for weight in weights if weight)))
        weight_histogram = Counter(weights)
        low_weight_histograms.append(
            (owner_type, tuple(
                (weight, weight_histogram[weight])
                for weight in sorted(weight_histogram)
                if weight <= 8
            ))
        )

        # The raw coefficient ratio modulo p is only a finite diagnostic.  Its
        # bounded rational representative does not restore the A-sheet
        # orientation or the dyadic type of either endpoint.
        for vertex in zero_neighbors:
            ratio = vertex * pow(root, -1, prime) % prime
            height = rational_height(ratio, prime)
            require(height is not None, ("height-seven miss", prime, root, vertex))
            atlas_heights.append(height[0])

    clique = maximum_disjoint_clique(vertices, masks)
    odd_clique = maximum_disjoint_clique(vertices, masks, "A")
    require(len(clique) == 6, ("mixed clique", prime, clique))
    require(len(odd_clique) == 4, ("odd clique", prime, odd_clique))
    partner_pairs = tuple(
        sorted(
            (min(vertex, modulus - vertex), max(vertex, modulus - vertex))
            for vertex in clique
            if vertex < modulus - vertex
        )
    )
    require(len(partner_pairs) == 3, ("pair packet", prime, clique, partner_pairs))
    require({vertex for pair in partner_pairs for vertex in pair} == set(clique),
            ("unpaired equality witness", prime, clique))

    # Profile budgets are reported independently of the DFS.  The raw maximum
    # includes profiles forbidden by THM-3435's proved target-free gates.  On
    # Q=2p those gates require exactly seven owners, at least two A owners, and
    # at least one E owner.  The common E fixed fibre also spends 2(a-1) sheets
    # before any other overlap.  The literal DFS above does not use these gates.
    raw_omegas = []
    gated_omegas = []
    for counts, omega in profiles_at_most_seven(groups, masks, modulus):
        if sum(counts.values()) != 7:
            continue
        raw_omegas.append(omega)
        fixed_fibre_cost = 2 * (counts["E"] - 1)
        if (counts["A"] >= 2 and counts["E"] >= 1
                and omega >= fixed_fibre_cost):
            gated_omegas.append(omega)
    gated_maximum = max(gated_omegas) if gated_omegas else None
    require(gated_maximum is None or gated_maximum <= 8,
            ("gated budget", prime, gated_omegas))

    return (
        prime,
        prime % 14,
        tuple((t, masks[groups[t][0]].bit_count()) for t in "ABE"),
        tuple(zero_degrees),
        tuple(positive_minima),
        tuple(low_weight_histograms),
        len(clique),
        clique,
        tuple(coefficient_type(vertex) for vertex in clique),
        partner_pairs,
        len(odd_clique),
        max(raw_omegas),
        gated_maximum,
        len(atlas_heights),
        max(atlas_heights),
    )


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert found")

    target_free_primes = tuple(
        prime for prime in primes_through(LIMIT)
        if prime % 2 and prime not in LOWER_BASE_PRIMES
    )
    require(target_free_primes[0] == 3 and target_free_primes[-1] == LIMIT,
            "prime universe endpoints")

    support = []
    controls = []
    scan_rows = []
    aggregate = SearchTally()
    for prime in target_free_primes:
        found, profile, omega, tally, masks, _groups = decide_prime(prime)
        aggregate.profiles += tally.profiles
        aggregate.nodes += tally.nodes
        aggregate.branches += tally.branches
        scan_rows.append((prime, tally.profiles, tally.nodes, tally.branches, found is not None))
        if found is not None:
            support.append((prime, found, profile, omega))
        if prime in POSITIVE_CONTROLS:
            controls.append(verify_positive_control(prime, POSITIVE_CONTROLS[prime], masks))

    require(tuple(prime for prime, *_ in support) == (7, 19), support)
    negatives = tuple(
        prime for prime in target_free_primes
        if prime not in {row[0] for row in support}
    )
    late_negatives = tuple(prime for prime in negatives if 191 <= prime <= LIMIT)
    require(late_negatives and late_negatives[0] == 191 and late_negatives[-1] == 599,
            ("late negatives", late_negatives[:1], late_negatives[-1:]))

    telemetry = tuple(tail_telemetry(prime) for prime in TAIL_CONTROLS)
    require(tuple(row[1] for row in telemetry) == (1, 3, 5, 9, 11, 13),
            "tail residue controls")
    require(all(row[3] == (("A", 19), ("B", 31), ("E", 23)) for row in telemetry),
            "zero-degree atlas")
    require(tuple(row[11] for row in telemetry) == (12, 8, 4, 10, 6, 2),
            "raw tail overlap budgets")
    require(tuple(row[12] for row in telemetry) == (8, None, 4, 4, 4, 0),
            "gated tail overlap budgets")
    require(all(row[13] == 73 and row[14] <= 7 for row in telemetry),
            "height-seven atlas")

    semantic_payload = (
        target_free_primes,
        tuple(support),
        tuple(controls),
        negatives,
        late_negatives,
        tuple(scan_rows),
        (aggregate.profiles, aggregate.nodes, aggregate.branches),
        telemetry,
    )
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("prime-even literal half-twist cap-seven finite atlas")
    print(f"universe=target-free odd primes p<=599; primes={len(target_free_primes)}")
    print(f"excluded_lower_base_primes={tuple(sorted(LOWER_BASE_PRIMES))}")
    print(f"search_profiles={aggregate.profiles} nodes={aggregate.nodes} branches={aggregate.branches}")
    print(f"finite_support={tuple(prime for prime, *_ in support)}")
    for prime, found, profile, omega in support:
        print(f"search_witness p={prime} Q={2*prime}: {found}; profile_ABE={profile}; omega={omega}")
    for control in controls:
        print(f"positive_control={control}")
    print(f"finite_negatives={negatives}")
    print(f"new_negative_tail_191_599={late_negatives}")
    print("tail_weighted_disjointness_controls:")
    for row in telemetry:
        print(row)
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
