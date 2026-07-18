#!/usr/bin/env python3
"""Exact scale-22 AP-centred Hamming-six owner-deficit certificate.

The certificate traverses all 924 labelled six-supports and every hereditary
leave-one-out-lcm order word over {1,2,11,22}.  A unit-independent scalar gate
is followed by literal owner-local sheet-union reachability.  The calculation
never quotients supports or unit words; multiplication orbits and tournaments
are telemetry only.

Methodological independence for a future referee: this primary solves CRT
representatives algebraically (while auditing them by literal search), stores
reachable unions as immutable Python sets, and hashes every reachable mask.
It uses only the Python standard library.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 22
ORDERS = (1, 2, 11, 22)
ORDER_INDEX = {order: i for i, order in enumerate(ORDERS)}
UNITS = {
    1: (0,),
    2: (1,),
    11: tuple(range(1, 11)),
    22: tuple(u for u in range(1, 22) if gcd(u, 22) == 1),
}
FULL_MASK = (1 << C) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def inverse_mod_13(value: int) -> int:
    return pow(value, -1, P)


def crt_base_algebraic(label: int, order: int, unit: int) -> int:
    """Solve x=D*label (mod 13), x=unit (mod D)."""
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def crt_base_literal(label: int, order: int, unit: int) -> int:
    for value in range(P * order):
        if value % P == order * label % P and value % order == unit % order:
            return value
    raise AssertionError("literal CRT search failed")


def local_mask(label: int, order: int, unit: int, owner: int) -> int:
    base = crt_base_algebraic(label, order, unit)
    inverse = inverse_mod_13(owner)
    result = 0
    for sheet in range(C):
        value = centered(base * (inverse + P * sheet), P * order)
        if -order < value <= order:
            result |= 1 << sheet
    return result


def analytic_cardinality(label: int, order: int, owner: int) -> int:
    """Independent one-period count, repeated C/D times."""
    ratio = label * inverse_mod_13(owner) % P
    target = order * ratio % P
    period_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * period_count


def build_tables():
    masks: dict[tuple[int, int, int, int], int] = {}
    cards: dict[tuple[int, int, int], int] = {}
    mask_hash = sha256()
    for label in range(1, P):
        for order in ORDERS:
            for owner in range(1, P):
                card = analytic_cardinality(label, order, owner)
                cards[label, order, owner] = card
                for unit in UNITS[order]:
                    algebraic = crt_base_algebraic(label, order, unit)
                    literal = crt_base_literal(label, order, unit)
                    require(algebraic == literal, "algebraic/literal CRT mismatch")
                    mask = local_mask(label, order, unit, owner)
                    require(mask.bit_count() == card, "mask/cardinality mismatch")
                    masks[label, order, unit, owner] = mask
                    mask_hash.update(bytes((label, order, unit, owner)))
                    mask_hash.update(mask.to_bytes(4, "little"))
    return masks, cards, mask_hash.hexdigest()


def hereditary(word: tuple[int, ...]) -> bool:
    return all(
        lcm(*(word[j] for j in range(6) if j != omitted)) == C
        for omitted in range(6)
    )


def hereditary_prime_power(word: tuple[int, ...]) -> bool:
    return (
        sum(order % 2 == 0 for order in word) >= 2
        and sum(order % 11 == 0 for order in word) >= 2
    )


def build_order_words():
    result = []
    digest = sha256()
    for word in product(ORDERS, repeat=6):
        h_lcm = hereditary(word)
        require(
            h_lcm == hereditary_prime_power(word),
            "lcm/prime-power hereditary mismatch",
        )
        if h_lcm:
            result.append(word)
            digest.update(bytes(word))
    return result, digest.hexdigest()


def scalar_capacities(
    support: tuple[int, ...], word: tuple[int, ...], cards
) -> tuple[int, ...]:
    return tuple(
        sum(cards[label, order, owner] for label, order in zip(support, word))
        for owner in support
    )


def build_scalar_bank(supports, words, cards):
    bank = []
    contexts_per_support = Counter()
    multiplicities = Counter()
    capacity_vectors = set()
    minimum_slack = Counter()
    maximum_slack = Counter()
    tight_owners = Counter()
    digest = sha256()
    checked = 0
    for support in supports:
        support_count = 0
        for word in words:
            checked += 1
            capacities = scalar_capacities(support, word, cards)
            if min(capacities) < C:
                continue
            row = (support, word, capacities)
            bank.append(row)
            support_count += 1
            multiplicities[tuple(word.count(order) for order in ORDERS)] += 1
            capacity_vectors.add(capacities)
            minimum_slack[min(capacities) - C] += 1
            maximum_slack[max(capacities) - C] += 1
            tight_owners[sum(value == C for value in capacities)] += 1
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes(capacities))
        if support_count:
            contexts_per_support[support_count] += 1
    require(checked == 924 * len(words), "scalar traversal count mismatch")
    return {
        "bank": bank,
        "contexts_per_support": contexts_per_support,
        "multiplicities": multiplicities,
        "capacity_vectors": capacity_vectors,
        "minimum_slack": minimum_slack,
        "maximum_slack": maximum_slack,
        "tight_owners": tight_owners,
        "digest": digest.hexdigest(),
    }


def owner_local_audit(support, word, owner, masks):
    reachable = frozenset({0})
    for label, order in zip(support, word):
        choices = frozenset(
            masks[label, order, unit, owner] for unit in UNITS[order]
        )
        reachable = frozenset(partial | choice for partial in reachable for choice in choices)
    maximum = max(mask.bit_count() for mask in reachable)
    feasible = FULL_MASK in reachable
    require(feasible == (maximum == C), "feasibility/maximum mismatch")
    return feasible, maximum, reachable


def strongly_connected_components(adjacency):
    n = len(adjacency)

    def reach(start: int, reverse: bool) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adjacency[w][v] if reverse else adjacency[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(range(n))
    count = 0
    while remaining:
        vertex = min(remaining)
        component = reach(vertex, False) & reach(vertex, True)
        remaining -= component
        count += 1
    return count


def hamiltonian_path_count(adjacency) -> int:
    n = len(adjacency)
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for subset in range(1 << n):
        for last in range(n):
            if not dp[subset][last]:
                continue
            for nxt in range(n):
                if not (subset >> nxt) & 1 and adjacency[last][nxt]:
                    dp[subset | (1 << nxt)][nxt] += dp[subset][last]
    return sum(dp[-1])


def tournament_audit(capacities, local_rows):
    n = 6
    adjacency = [[False] * n for _ in range(n)]
    ties = 0
    flips = 0
    for left in range(n):
        for right in range(left + 1, n):
            left_key = (local_rows[left][0], local_rows[left][1], capacities[left])
            right_key = (local_rows[right][0], local_rows[right][1], capacities[right])
            if left_key == right_key:
                ties += 1
                winner, loser = left, right
            elif left_key > right_key:
                winner, loser = left, right
            else:
                winner, loser = right, left
                flips += 1
            adjacency[winner][loser] = True
    scores = tuple(sorted(sum(row) for row in adjacency))
    triangles = sum(
        adjacency[a][b] and adjacency[b][c] and adjacency[c][a]
        or adjacency[a][c] and adjacency[c][b] and adjacency[b][a]
        for a, b, c in combinations(range(n), 3)
    )
    return (
        ties,
        flips,
        scores,
        triangles,
        strongly_connected_components(adjacency),
        hamiltonian_path_count(adjacency),
    )


def multiplication_orbit_histogram(support_bank):
    remaining = set(support_bank)
    histogram = Counter()
    while remaining:
        support = min(remaining)
        orbit = {
            tuple(sorted((multiplier * label) % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= support_bank, "support bank is not multiplication-invariant")
        remaining -= orbit
        histogram[len(orbit)] += 1
    return histogram


def format_counter(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    masks, cards, mask_hash = build_tables()
    words, order_hash = build_order_words()
    supports = tuple(combinations(range(1, P), 6))
    state_words_per_support = sum(
        prod(len(UNITS[order]) for order in word) for word in words
    )

    scalar = build_scalar_bank(supports, words, cards)
    bank = scalar["bank"]
    support_bank = {support for support, _, _ in bank}

    feasible_contexts = Counter()
    maximum_union = Counter()
    minimum_owner_maximum = Counter()
    owner_vectors = set()
    reachable_counts = Counter()
    reachable_hash = sha256()
    tie_histogram = Counter()
    flip_histogram = Counter()
    all_tournaments_ok = True
    feasible_rows = 0

    for support, word, capacities in bank:
        local_rows = []
        for owner in support:
            feasible, maximum, reachable = owner_local_audit(
                support, word, owner, masks
            )
            local_rows.append((feasible, maximum))
            feasible_rows += feasible
            maximum_union[maximum] += 1
            reachable_counts[len(reachable)] += 1
            reachable_hash.update(bytes(support))
            reachable_hash.update(bytes(word))
            reachable_hash.update(bytes((owner,)))
            for mask in sorted(reachable):
                reachable_hash.update(mask.to_bytes(4, "little"))
        feasible_count = sum(feasible for feasible, _ in local_rows)
        feasible_contexts[feasible_count] += 1
        owner_vector = tuple(maximum for _, maximum in local_rows)
        owner_vectors.add(owner_vector)
        minimum_owner_maximum[min(owner_vector)] += 1
        tournament = tournament_audit(capacities, local_rows)
        ties, flips, scores, triangles, sccs, paths = tournament
        tie_histogram[ties] += 1
        flip_histogram[flips] += 1
        all_tournaments_ok &= (
            scores == (0, 1, 2, 3, 4, 5)
            and triangles == 0
            and sccs == 6
            and paths == 1
        )

    expected_multiplicities = Counter(
        {
            (0, 2, 0, 4): 36,
            (0, 2, 1, 3): 144,
            (0, 2, 2, 2): 216,
            (0, 2, 3, 1): 144,
            (0, 2, 4, 0): 36,
            (0, 3, 1, 2): 288,
            (0, 3, 2, 1): 96,
            (0, 3, 3, 0): 24,
        }
    )
    require(len(supports) == 924, "support count")
    require(len(words) == 3249, "hereditary word count")
    require(state_words_per_support == 100_975_500, "literal state count")
    require(len(bank) == 984 and len(support_bank) == 180, "scalar bank")
    require(scalar["multiplicities"] == expected_multiplicities, "multiplicities")
    require(
        scalar["contexts_per_support"] == Counter({2: 96, 3: 24, 6: 24, 16: 36}),
        "contexts/support",
    )
    require(feasible_contexts == Counter({0: 792, 1: 192}), "owner deficit")
    require(
        maximum_union == Counter({16: 864, 17: 1584, 18: 2784, 19: 480, 22: 192}),
        "maximum-union histogram",
    )
    require(all_tournaments_ok, "tournament fingerprint")

    orbit_histogram = multiplication_orbit_histogram(support_bank)

    print("scale-twenty-two AP-centred Hamming-six owner-deficit certificate")
    print("divisor grammar 1,2,11,22; literal states 22")
    print(
        f"supports {len(supports)}; hereditary order words {len(words)}; "
        f"labelled order contexts {len(supports) * len(words)}"
    )
    print(
        f"state words/support {state_words_per_support}; raw labelled states "
        f"{len(supports) * state_words_per_support}"
    )
    print(f"mask SHA256 {mask_hash}; order SHA256 {order_hash}")
    print(
        f"scalar contexts {len(bank)} on {len(support_bank)} supports; "
        f"multiplicity patterns {len(scalar['multiplicities'])}; "
        f"scalar-bank SHA256 {scalar['digest']}"
    )
    print("scalar multiplicities n1,n2,n11,n22", end=" ")
    print(" ".join(f"{','.join(map(str, key))}:{value}" for key, value in sorted(scalar["multiplicities"].items())))
    print("contexts-per-support histogram", format_counter(scalar["contexts_per_support"]))
    print("multiplication orbit-size histogram", format_counter(orbit_histogram), "(telemetry; no quotient)")
    print(f"capacity vectors {len(scalar['capacity_vectors'])}")
    print("minimum scalar-slack histogram", format_counter(scalar["minimum_slack"]))
    print("maximum scalar-slack histogram", format_counter(scalar["maximum_slack"]))
    print("tight-owner/context histogram", format_counter(scalar["tight_owners"]))
    print(f"owner-local rows {6 * len(bank)}; feasible rows {feasible_rows}")
    print("feasible-owner/context histogram", format_counter(feasible_contexts))
    print("maximum reachable sheet-union histogram", format_counter(maximum_union))
    print("minimum owner maximum/context histogram", format_counter(minimum_owner_maximum))
    print(f"distinct owner max-union vectors {len(owner_vectors)}")
    print(f"reachable-union-bank SHA256 {reachable_hash.hexdigest()}")
    print("reachable-count histogram", format_counter(reachable_counts))
    print("owner-local all-six contexts", feasible_contexts.get(6, 0))
    print("tournament pair observable exact ordered (feasible,max-union,capacity) owner summaries; lexicographic switch and coordinate tie Hamiltonian path")
    print("tournament fingerprints all 984 transitive: scores 0,1,2,3,4,5; cycles 0; SCCs 6; Hamiltonian paths 1")
    print("tournament tie-edge histogram", format_counter(tie_histogram))
    print("tournament edge-flip histogram", format_counter(flip_histogram))
    print("challenged vertices owner obligations preserve the terminal deficit through their feasibility/max-union vector and the scalar gate through capacities; the tournament loses thresholds and magnitudes, while provider, divisor, residue, isolated sheet, and wall-event vertices lose shared-unit incidence")
    print("frontier verdict scalar-empty no; owner-local all-six empty yes; next legal common scale 24 (23 is prime-excluded)")
    print("local D22 mask table at owner one (units 1,3,5,7,9,13,15,17,19,21; ratios 1..12 in hex)")
    for unit in UNITS[22]:
        row = " ".join(f"{masks[label, 22, unit, 1]:x}" for label in range(1, P))
        print(f"  e={unit}: {row}")


if __name__ == "__main__":
    main()
