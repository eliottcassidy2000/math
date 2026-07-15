#!/usr/bin/env python3
"""Exact joint-cut certificate for the primitive AP-centred H6 scale cap.

This script uses only integer arithmetic.  It enumerates the complete
order-only relaxation supplied by THM-860's leave-one-out lcm law and all 63
relative cuts, then tests the unit-independent scalar owner-capacity condition
at the largest arithmetic scales.

Tournament Analysis is deliberately included as a loss audit.  Its vertices
are the six order slots.  A raw order comparison and a complement-conditioned
pair comparison are completed by the slot-order tie gauge.  Both completions
are transitive on the extremal rows, while the decisive object is the full
63-cut prime-power hypergraph; the pair completion forgets that object.
"""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations
from math import gcd, lcm


N = 6
ALL = (1 << N) - 1
PRIME_RANGES = ((2, 5), (3, 3), (5, 1), (7, 1))


def rho_num(m):
    """Numerator of rho(m)=ceil(2m/13)/m."""
    return (2 * m + 12) // 13


def rho(m):
    return Fraction(rho_num(m), m)


def lcm_many(values):
    answer = 1
    for value in values:
        answer = lcm(answer, value)
    return answer


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def relative_orders(orders, mask):
    complement_lcm = lcm_many(
        orders[i] for i in range(N) if not ((mask >> i) & 1)
    )
    return tuple(
        orders[i] // gcd(orders[i], complement_lcm)
        for i in range(N)
        if (mask >> i) & 1
    )


def cut_holds(orders, mask):
    ms = relative_orders(orders, mask)
    denominator = lcm_many(ms)
    numerator = sum(rho_num(m) * (denominator // m) for m in ms)
    return numerator >= denominator


def all_cuts_hold(orders):
    return all(cut_holds(orders, mask) for mask in range(1, ALL + 1))


def unique_permutations(word, cache={}):
    if word not in cache:
        cache[word] = tuple(sorted(set(permutations(word))))
    return cache[word]


def admissible_prime_profiles(prime, exponent_range):
    profiles = []
    for exponents in combinations_with_replacement(range(exponent_range + 1), N):
        # gcd-normalization and leave-one-out lcm, prime by prime.
        if exponents[0] != 0 or exponents.count(exponents[-1]) < 2:
            continue
        orders = tuple(prime**exponent for exponent in exponents)
        if all_cuts_hold(orders):
            profiles.append(exponents)
    return tuple(profiles)


def align_prime(words, prime, profiles):
    aligned = set()
    for word in words:
        for profile in profiles:
            for exponents in unique_permutations(profile):
                aligned.add(
                    tuple(sorted(word[i] * prime**exponents[i] for i in range(N)))
                )
    return aligned


def enumerate_normalized_words():
    profiles = {
        prime: admissible_prime_profiles(prime, exponent_range)
        for prime, exponent_range in PRIME_RANGES
    }

    # Simultaneous colour permutation lets us fix the dyadic profile in sorted
    # order.  Every later prime profile is aligned in every distinct way.
    words = set()
    alignments_23 = 0
    for profile_2 in profiles[2]:
        for profile_3 in profiles[3]:
            for exponents_3 in unique_permutations(profile_3):
                alignments_23 += 1
                words.add(
                    tuple(
                        sorted(
                            2 ** profile_2[i] * 3 ** exponents_3[i]
                            for i in range(N)
                        )
                    )
                )
    raw_23 = len(words)
    words = {word for word in words if all_cuts_hold(word)}
    kept_23 = len(words)

    stage_counts = [("2x3", raw_23, kept_23)]
    for prime in (5, 7):
        words = align_prime(words, prime, profiles[prime])
        raw = len(words)
        words = {word for word in words if all_cuts_hold(word)}
        stage_counts.append((f"through {prime}", raw, len(words)))

    for word in words:
        assert gcd_many(word) == 1
        common_scale = lcm_many(word)
        assert all(
            lcm_many(word[j] for j in range(N) if j != i) == common_scale
            for i in range(N)
        )

    return profiles, alignments_23, stage_counts, tuple(sorted(words))


def enumerate_multiplier_rows(normalized_words):
    rows = []
    for word in normalized_words:
        # THM-860's exact all-six pivot is min(D_i)<=72.
        for common_gcd in range(1, 72 // min(word) + 1):
            if common_gcd % 13 == 0:
                continue
            orders = tuple(common_gcd * value for value in word)
            # Every proper cut is unchanged by a common multiplier and was
            # already checked on ``word``.  Only the all-six cut changes.
            if cut_holds(orders, ALL):
                rows.append((lcm_many(orders), orders, word, common_gcd))
    return tuple(sorted(rows))


def interval_residue_count(order, provider, owner):
    target = order * provider * pow(owner, -1, 13) % 13
    return sum(1 for z in range(-order + 1, order + 1) if z % 13 == target)


def scalar_owner_scan(orders):
    """Necessary, unit-independent scalar common-sheet capacity scan."""
    common_scale = lcm_many(orders)
    capacity = {}
    for i, order in enumerate(orders):
        for provider in range(1, 13):
            for owner in range(1, 13):
                capacity[i, provider, owner] = (
                    common_scale
                    // order
                    * interval_residue_count(order, provider, owner)
                )

    assignments = 0
    survivors = 0
    best_floor = -1
    best_labels = None
    for tail in permutations(range(2, 13), 5):
        labels = (1,) + tail
        assignments += 1
        owner_floor = min(
            sum(capacity[i, labels[i], owner] for i in range(N))
            for owner in labels
        )
        if owner_floor > best_floor:
            best_floor = owner_floor
            best_labels = labels
        if owner_floor >= common_scale:
            survivors += 1
    assert assignments == 55_440
    return assignments, survivors, best_floor, best_labels


def tournament_edges(orders, conditioned):
    edges = set()
    ties = 0
    for i, j in combinations(range(N), 2):
        if conditioned:
            other_lcm = lcm_many(
                orders[k] for k in range(N) if k not in (i, j)
            )
            left = rho(orders[i] // gcd(orders[i], other_lcm))
            right = rho(orders[j] // gcd(orders[j], other_lcm))
        else:
            left = orders[i]
            right = orders[j]
        if left == right:
            ties += 1
            winner = i  # fixed slot-order tie Hamiltonian path
        elif left > right:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        edges.add((winner, loser))
    return frozenset(edges), ties


def tournament_fingerprint(edges):
    scores = [sum((i, j) in edges for j in range(N) if j != i) for i in range(N)]
    score_histogram = {score: scores.count(score) for score in sorted(set(scores))}
    triangles = 0
    for i, j, k in combinations(range(N), 3):
        if (
            ((i, j) in edges and (j, k) in edges and (k, i) in edges)
            or ((j, i) in edges and (k, j) in edges and (i, k) in edges)
        ):
            triangles += 1

    reach = [[i == j or (i, j) in edges for j in range(N)] for i in range(N)]
    for k in range(N):
        for i in range(N):
            for j in range(N):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    seen = set()
    scc_sizes = []
    for i in range(N):
        if i not in seen:
            component = {j for j in range(N) if reach[i][j] and reach[j][i]}
            seen.update(component)
            scc_sizes.append(len(component))

    hamiltonian_paths = sum(
        all((path[i], path[i + 1]) in edges for i in range(N - 1))
        for path in permutations(range(N))
    )
    return score_histogram, triangles, tuple(sorted(scc_sizes)), hamiltonian_paths


def fraction_string(value):
    return str(value.numerator) if value.denominator == 1 else str(value)


def main():
    profiles, alignments_23, stage_counts, normalized_words = enumerate_normalized_words()
    rows = enumerate_multiplier_rows(normalized_words)
    arithmetic_cap = max(scale for scale, *_ in rows)
    extremals = tuple(orders for scale, orders, *_ in rows if scale == arithmetic_cap)
    distinct_scales = tuple(sorted({scale for scale, *_ in rows}, reverse=True))

    expected_extremals = (
        (20, 35, 56, 80, 1120, 1120),
        (20, 35, 80, 224, 280, 1120),
        (20, 35, 80, 280, 1120, 1120),
        (20, 35, 112, 160, 280, 1120),
    )
    assert tuple(sorted(extremals)) == expected_extremals
    assert distinct_scales[:2] == (1120, 1008)

    print("THM-860 JOINT RAMIFICATION CAP -- EXACT INTEGER CERTIFICATE")
    print("method: gcd-normalized prime profiles + all 63 relative cuts")
    print("rho(m)=ceil(2m/13)/m; every comparison uses a common integer denominator")
    print("monotone prune: rho(qk)<=rho(q), so a failed partial-prime cut cannot recover")
    print()
    print("prime-profile counts:")
    for prime, _ in PRIME_RANGES:
        print(f"  p={prime}: {len(profiles[prime])}")
    print(f"2x3 profile alignments represented: {alignments_23}")
    for stage, raw, kept in stage_counts:
        print(f"{stage}: raw distinct={raw}, all-cut survivors={kept}")
    print(f"final normalized words: {len(normalized_words)}")
    print(f"feasible (normalized word, common-gcd) rows: {len(rows)}")
    print(f"largest two arithmetic scales: {distinct_scales[0]}, {distinct_scales[1]}")
    print(f"arithmetic relaxation cap: {arithmetic_cap}")
    print()
    print("arithmetic extremals:")
    for orders in extremals:
        all_sum = sum((rho(order) for order in orders), Fraction())
        leave_one_out = all(
            lcm_many(orders[j] for j in range(N) if j != i) == arithmetic_cap
            for i in range(N)
        )
        print(
            f"  {orders}: cuts=63/63, leave-one-out={leave_one_out}, "
            f"all-six-sum={fraction_string(all_sum)}"
        )

    print()
    print("unit-independent scalar owner-capacity scan at c=1120:")
    scalar_results = []
    for orders in extremals:
        result = scalar_owner_scan(orders)
        scalar_results.append(result)
        assignments, survivors, best_floor, best_labels = result
        print(
            f"  {orders}: assignments={assignments}, covers={survivors}, "
            f"best-floor={best_floor}/1120, deficit={1120-best_floor}, "
            f"first-best-labels={best_labels}"
        )
    assert all(survivors == 0 for _, survivors, _, _ in scalar_results)
    print("common-sheet consequence: c=1120 impossible; hence c<=1008")

    print()
    print("Tournament Analysis loss audit on the four arithmetic extremals:")
    print("  vertices=six sorted order slots")
    print("  raw gauge=larger D wins; conditioned gauge=larger relative rho wins")
    print("  tie gauge=lower slot index, giving the tie Hamiltonian path 0>1>...>5")
    for orders in extremals:
        raw_edges, raw_ties = tournament_edges(orders, conditioned=False)
        conditioned_edges, conditioned_ties = tournament_edges(orders, conditioned=True)
        raw_fp = tournament_fingerprint(raw_edges)
        conditioned_fp = tournament_fingerprint(conditioned_edges)
        flips = sum(
            ((i, j) in raw_edges) != ((i, j) in conditioned_edges)
            for i, j in combinations(range(N), 2)
        )
        print(f"  {orders}:")
        print(
            f"    raw ties={raw_ties}, scores={raw_fp[0]}, 3-cycles={raw_fp[1]}, "
            f"SCC-sizes={raw_fp[2]}, HP={raw_fp[3]}"
        )
        print(
            f"    conditioned ties={conditioned_ties}, scores={conditioned_fp[0]}, "
            f"3-cycles={conditioned_fp[1]}, SCC-sizes={conditioned_fp[2]}, "
            f"HP={conditioned_fp[3]}, edge-flips={flips}"
        )
    print("  verdict: every conditioned pair is tied, while the 63-cut hypergraph decides")
    print("  preserved by pair completion: a scalar total order and tie path")
    print("  destroyed by pair completion: prime-power upper sets and multi-colour cut debt")
    print("  challenged vertices: cut obligations/prime-power fibres, not only runners or arcs")


if __name__ == "__main__":
    main()
