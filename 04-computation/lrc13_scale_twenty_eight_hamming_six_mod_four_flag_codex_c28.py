#!/usr/bin/env python3
"""Exact scale-28 AP-centred Hamming-six mod-four flag certificate.

This certificate reconstructs the hereditary divisor grammar, ratio-cardinality
table, complete scalar survivor bank, and all owner-local union obligations for
the primitive proper common-scale c=28 face.  Its proof-facing quotient keeps
the four sheet fibres modulo 4 and only their saturated capacities (seven per
fibre).  On the scalar survivor bank this quotient is equivalent to literal
28-bit cover feasibility.  A parallel mod-seven flag is included to show that
the choice of quotient is substantive: it has false positives.

Tournament Analysis is telemetry only.  The vertices are owner obligations,
not runners: pairwise lexicographic comparison of exact or flag observables is
completed along the coordinate-order tie Hamiltonian path.  This deliberately
lossy total-order switch is audited by scores, triangles, SCCs, edge flips, and
Hamiltonian-path counts; it is not used in the obstruction.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

import numpy as np


P = 13
C = 28
LABELS = tuple(range(1, P))
ORDERS = (1, 2, 4, 7, 14, 28)
UNITS = {
    order: (
        (0,)
        if order == 1
        else tuple(unit for unit in range(1, order) if gcd(unit, order) == 1)
    )
    for order in ORDERS
}
FULL = (1 << C) - 1


def require(condition: bool, message: str) -> None:
    """Optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    """Representative in (-modulus/2, modulus/2]."""
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def literal_crt_base(label: int, order: int, unit: int) -> int:
    """Find x = order*label (mod 13), x = unit (mod order) literally."""
    candidates = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(len(candidates) == 1, "CRT representative is not unique")
    return candidates[0]


def algebraic_crt_base(label: int, order: int, unit: int) -> int:
    if order == 1:
        return label
    step = ((unit - order * label) * pow(P, -1, order)) % order
    return (order * label + P * step) % (P * order)


def literal_mask(label: int, order: int, unit: int, owner: int) -> int:
    """Strict-left/closed-right local sheet mask."""
    base = literal_crt_base(label, order, unit)
    owner_inverse = pow(owner, -1, P)
    mask = 0
    for sheet in range(C):
        value = centered(base * (owner_inverse + P * sheet), P * order)
        if -order < value <= order:
            mask |= 1 << sheet
    return mask


def rotate(mask: int, amount: int) -> int:
    amount %= C
    if amount == 0:
        return mask
    return ((mask << amount) | (mask >> (C - amount))) & FULL


def analytic_cardinality(ratio: int, order: int) -> int:
    target = order * ratio % P
    one_period = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * one_period


def projection_signature(mask: int, modulus: int) -> tuple[int, ...]:
    return tuple(
        sum((mask >> sheet) & 1 for sheet in range(residue, C, modulus))
        for residue in range(modulus)
    )


def build_tables():
    """Build normalized masks and audit every literal owner gauge."""
    normalized = {}
    cards = np.zeros((len(ORDERS), P), dtype=np.int16)
    base_digest = sha256()
    mask_digest = sha256()
    signature_digest = sha256()

    for order_index, order in enumerate(ORDERS):
        for ratio in LABELS:
            options = frozenset(
                literal_mask(ratio, order, unit, 1) for unit in UNITS[order]
            )
            expected = analytic_cardinality(ratio, order)
            require(
                {mask.bit_count() for mask in options} == {expected},
                "normalized mask cardinality depends on exact-order residue",
            )
            normalized[ratio, order] = options
            cards[order_index, ratio] = expected
            signature_digest.update(bytes((ratio, order)))
            for mask in sorted(options):
                signature_digest.update(bytes(projection_signature(mask, 4)))
                signature_digest.update(bytes(projection_signature(mask, 7)))

    for label in LABELS:
        for order in ORDERS:
            for unit in UNITS[order]:
                literal = literal_crt_base(label, order, unit)
                algebraic = algebraic_crt_base(label, order, unit)
                require(literal == algebraic, "literal/algebraic CRT mismatch")
                base_digest.update(bytes((label, order, unit)))
                base_digest.update(literal.to_bytes(2, "little"))
            for owner in LABELS:
                ratio = label * pow(owner, -1, P) % P
                shift = next(
                    candidate
                    for candidate in range(C)
                    if (pow(owner, -1, P) + P * candidate) % C == 1
                )
                actual_options = []
                for unit in UNITS[order]:
                    actual = literal_mask(label, order, unit, owner)
                    expected = rotate(literal_mask(ratio, order, unit, 1), shift)
                    require(actual == expected, "cyclic owner gauge mismatch")
                    require(
                        actual.bit_count() == cards[ORDERS.index(order), ratio],
                        "literal mask/cardinality mismatch",
                    )
                    actual_options.append(actual)
                    mask_digest.update(bytes((label, order, unit, owner)))
                    mask_digest.update(actual.to_bytes(4, "little"))
                require(
                    len(set(actual_options)) == len(normalized[ratio, order]),
                    "owner gauge changed the option-bank cardinality",
                )

    expected_rows = {
        1: (28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        2: (14, 0, 0, 0, 0, 14, 14, 0, 0, 0, 0, 0),
        4: (7, 0, 7, 7, 0, 7, 7, 0, 7, 7, 0, 0),
        7: (8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
        14: (6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4),
        28: (5, 4, 4, 4, 4, 5, 5, 4, 4, 4, 4, 4),
    }
    actual_rows = {
        order: tuple(int(cards[ORDERS.index(order), ratio]) for ratio in LABELS)
        for order in ORDERS
    }
    require(actual_rows == expected_rows, "ratio-cardinality table changed")
    return {
        "normalized": normalized,
        "cards": cards,
        "rows": actual_rows,
        "base_digest": base_digest.hexdigest(),
        "mask_digest": mask_digest.hexdigest(),
        "signature_digest": signature_digest.hexdigest(),
    }


def hereditary(word: tuple[int, ...]) -> bool:
    return all(
        lcm(*(word[index] for index in range(6) if index != omitted)) == C
        for omitted in range(6)
    )


def build_grammar():
    words = []
    weighted_states = 0
    digest = sha256()
    for index_word in product(range(len(ORDERS)), repeat=6):
        word = tuple(ORDERS[index] for index in index_word)
        via_lcm = hereditary(word)
        via_valuations = (
            sum(order % 4 == 0 for order in word) >= 2
            and sum(order % 7 == 0 for order in word) >= 2
        )
        require(via_lcm == via_valuations, "hereditary valuation grammar mismatch")
        if not via_lcm:
            continue
        fibre = prod(len(UNITS[order]) for order in word)
        words.append(index_word)
        weighted_states += fibre
        digest.update(bytes(index_word))
        digest.update(fibre.to_bytes(8, "little"))
    require(len(words) == 26_961, "hereditary word census changed")
    require(weighted_states == 429_048_576, "weighted grammar census changed")
    return np.asarray(words, dtype=np.int8), weighted_states, digest.hexdigest()


def scalar_bank(words: np.ndarray, cards: np.ndarray):
    rows = []
    feasible_owner_histogram = np.zeros(7, dtype=np.int64)
    support_histogram = Counter()
    profile_histogram = Counter()
    minimum_slack = Counter()
    maximum_slack = Counter()
    digest = sha256()

    for support in combinations(LABELS, 6):
        ratio_cards = np.empty((6, len(ORDERS), 6), dtype=np.int16)
        for provider, label in enumerate(support):
            for order_index in range(len(ORDERS)):
                for owner_index, owner in enumerate(support):
                    ratio = label * pow(owner, -1, P) % P
                    ratio_cards[provider, order_index, owner_index] = cards[
                        order_index, ratio
                    ]

        capacities = np.zeros((len(words), 6), dtype=np.int16)
        for provider in range(6):
            capacities += ratio_cards[provider, words[:, provider], :]
        feasible_counts = (capacities >= C).sum(axis=1)
        feasible_owner_histogram += np.bincount(feasible_counts, minlength=7)
        survivor_indices = np.flatnonzero(feasible_counts == 6)
        support_histogram[len(survivor_indices)] += 1

        for row_index in survivor_indices:
            word = tuple(ORDERS[int(index)] for index in words[row_index])
            capacity_word = tuple(int(value) for value in capacities[row_index])
            rows.append((support, word, capacity_word))
            profile_histogram[
                tuple(word.count(order) for order in ORDERS)
            ] += 1
            minimum_slack[min(capacity_word) - C] += 1
            maximum_slack[max(capacity_word) - C] += 1
            digest.update(bytes(support))
            digest.update(bytes(word))
            digest.update(bytes(capacity_word))

    expected_feasible = (120_024, 5_260_824, 10_675_332, 6_969_052,
                         1_724_910, 158_652, 3_170)
    require(
        tuple(int(value) for value in feasible_owner_histogram) == expected_feasible,
        "scalar feasible-owner histogram changed",
    )
    require(len(rows) == 3_170, "scalar survivor census changed")
    require(sum(count for survivors, count in support_histogram.items() if survivors) == 206,
            "scalar support census changed")
    return {
        "rows": rows,
        "feasible_owner_histogram": tuple(int(x) for x in feasible_owner_histogram),
        "support_histogram": support_histogram,
        "profile_histogram": profile_histogram,
        "minimum_slack": minimum_slack,
        "maximum_slack": maximum_slack,
        "digest": digest.hexdigest(),
    }


def normalized_owner_key(support, word, owner):
    inverse = pow(owner, -1, P)
    return tuple(
        sorted((label * inverse % P, order) for label, order in zip(support, word))
    )


def exact_owner_bank(key, masks):
    reachable = frozenset((0,))
    for ratio, order in key:
        reachable = frozenset(
            partial | choice
            for partial in reachable
            for choice in masks[ratio, order]
        )
    maximum = max(mask.bit_count() for mask in reachable)
    digest = sha256()
    for mask in sorted(reachable):
        digest.update(mask.to_bytes(4, "little"))
    return maximum == C, maximum, len(reachable), digest.digest()


def capacity_flag_bank(key, masks, modulus):
    threshold = C // modulus
    target = (threshold,) * modulus
    states = frozenset(((0,) * modulus,))
    for ratio, order in key:
        options = frozenset(
            projection_signature(mask, modulus) for mask in masks[ratio, order]
        )
        states = frozenset(
            tuple(min(threshold, left + right) for left, right in zip(state, option))
            for state in states
            for option in options
        )
    maximum = max(sum(state) for state in states)
    digest = sha256()
    for state in sorted(states):
        digest.update(bytes(state))
    return target in states, maximum, len(states), digest.digest()


def tournament(keys):
    out = [0] * 6
    ties = 0
    natural_flips = 0
    for left in range(6):
        for right in range(left + 1, 6):
            if keys[left] == keys[right]:
                winner = left
                ties += 1
            else:
                winner = left if keys[left] > keys[right] else right
            loser = left + right - winner
            out[winner] |= 1 << loser
            natural_flips += int(winner == right)

    scores = sorted(mask.bit_count() for mask in out)
    require(scores == list(range(6)), "owner tournament is not transitive")
    triangles = 0
    for a, b, c in combinations(range(6), 3):
        triangles += int(
            ((out[a] >> b) & 1 and (out[b] >> c) & 1 and (out[c] >> a) & 1)
            or ((out[a] >> c) & 1 and (out[c] >> b) & 1 and (out[b] >> a) & 1)
        )
    require(triangles == 0, "owner tournament acquired a directed triangle")

    paths = [[0] * 6 for _ in range(64)]
    for vertex in range(6):
        paths[1 << vertex][vertex] = 1
    for subset in range(1, 64):
        for last in range(6):
            if not (subset >> last) & 1:
                continue
            previous_subset = subset ^ (1 << last)
            for previous in range(6):
                if (previous_subset >> previous) & 1 and (out[previous] >> last) & 1:
                    paths[subset][last] += paths[previous_subset][previous]
    require(sum(paths[63]) == 1, "owner tournament Hamiltonian-path count changed")
    # The transitive score word also certifies six singleton SCCs.
    return tuple(out), ties, natural_flips


def multiplication_orbits(rows):
    row_set = {(support, word, capacities) for support, word, capacities in rows}
    seen = set()
    histogram = Counter()
    for row in sorted(row_set):
        if row in seen:
            continue
        support, word, capacities = row
        keyed = {label: (order, capacity) for label, order, capacity in zip(support, word, capacities)}
        orbit = set()
        for multiplier in LABELS:
            transported = sorted(
                (multiplier * label % P, order, capacity)
                for label, (order, capacity) in keyed.items()
            )
            orbit.add(
                (
                    tuple(item[0] for item in transported),
                    tuple(item[1] for item in transported),
                    tuple(item[2] for item in transported),
                )
            )
        require(orbit <= row_set, "scalar bank is not multiplication-invariant")
        seen |= orbit
        histogram[len(orbit)] += 1
    require(seen == row_set, "multiplication orbit partition is incomplete")
    return histogram


def audit_owners(scalar, masks):
    exact_memo = {}
    mod4_memo = {}
    mod7_memo = {}
    feasible_histogram = Counter()
    maximum_histogram = Counter()
    flag4_maximum_histogram = Counter()
    maximum_by_order = Counter()
    profile_feasible = Counter()
    confusion4 = Counter()
    confusion7 = Counter()
    exact_ties = Counter()
    flag_ties = Counter()
    exact_natural_flips = Counter()
    flag_natural_flips = Counter()
    gauge_edge_flips = Counter()
    exact_incidence_total = 0
    exact_digest = sha256()
    flag4_digest = sha256()
    flag7_digest = sha256()

    for support, word, capacities in scalar["rows"]:
        exact_summaries = []
        flag4_summaries = []
        feasible_owners = 0
        profile = tuple(word.count(order) for order in ORDERS)
        for owner_index, owner in enumerate(support):
            key = normalized_owner_key(support, word, owner)
            if key not in exact_memo:
                exact_memo[key] = exact_owner_bank(key, masks)
                mod4_memo[key] = capacity_flag_bank(key, masks, 4)
                mod7_memo[key] = capacity_flag_bank(key, masks, 7)
            exact = exact_memo[key]
            flag4 = mod4_memo[key]
            flag7 = mod7_memo[key]
            require(exact[0] == flag4[0], "mod-four flag changed exact feasibility")
            confusion4[exact[0], flag4[0]] += 1
            confusion7[exact[0], flag7[0]] += 1
            feasible_owners += int(exact[0])
            maximum_histogram[exact[1]] += 1
            flag4_maximum_histogram[flag4[1]] += 1
            maximum_by_order[word[owner_index], exact[1]] += 1
            exact_incidence_total += exact[2]
            exact_digest.update(bytes((owner_index,)))
            exact_digest.update(bytes((int(exact[0]), exact[1])))
            exact_digest.update(exact[2].to_bytes(4, "little"))
            exact_digest.update(exact[3])
            flag4_digest.update(bytes((owner_index, int(flag4[0]), flag4[1])))
            flag4_digest.update(flag4[2].to_bytes(4, "little"))
            flag4_digest.update(flag4[3])
            flag7_digest.update(bytes((owner_index, int(flag7[0]), flag7[1])))
            flag7_digest.update(flag7[2].to_bytes(4, "little"))
            flag7_digest.update(flag7[3])
            exact_summaries.append(exact)
            flag4_summaries.append(flag4)

        require(feasible_owners <= 2, "a scalar row has three feasible owners")
        feasible_histogram[feasible_owners] += 1
        profile_feasible[profile, feasible_owners] += 1

        exact_keys = tuple(
            (int(summary[0]), summary[1], capacities[index], summary[2])
            for index, summary in enumerate(exact_summaries)
        )
        flag_keys = tuple(
            (int(summary[0]), summary[1], capacities[index], summary[2])
            for index, summary in enumerate(flag4_summaries)
        )
        exact_out, exact_tie, exact_flip = tournament(exact_keys)
        flag_out, flag_tie, flag_flip = tournament(flag_keys)
        exact_ties[exact_tie] += 1
        flag_ties[flag_tie] += 1
        exact_natural_flips[exact_flip] += 1
        flag_natural_flips[flag_flip] += 1
        gauge_edge_flips[
            sum((exact_out[i] ^ flag_out[i]).bit_count() for i in range(6)) // 2
        ] += 1

    require(feasible_histogram == Counter({0: 2_018, 1: 912, 2: 240}),
            "owner feasible-row histogram changed")
    require(confusion4 == Counter({(False, False): 17_628, (True, True): 1_392}),
            "mod-four exact/flag contingency changed")
    require(confusion7 == Counter({(False, False): 11_556,
                                  (False, True): 6_072,
                                  (True, True): 1_392}),
            "mod-seven exact/flag contingency changed")
    require(all(order == 2 for (order, maximum), count in maximum_by_order.items()
                if maximum == C and count),
            "a non-order-two owner became feasible")
    key_orders = Counter()
    feasible_key_orders = Counter()
    for key, exact in exact_memo.items():
        owner_orders = tuple(order for ratio, order in key if ratio == 1)
        require(len(owner_orders) == 1, "normalized key lost its unique owner")
        key_orders[owner_orders[0]] += 1
        if exact[0]:
            feasible_key_orders[owner_orders[0]] += 1
    return {
        "exact_keys": len(exact_memo),
        "key_orders": key_orders,
        "feasible_key_orders": feasible_key_orders,
        "exact_bank_sizes": Counter(summary[2] for summary in exact_memo.values()),
        "mod4_state_sizes": Counter(summary[2] for summary in mod4_memo.values()),
        "mod7_state_sizes": Counter(summary[2] for summary in mod7_memo.values()),
        "feasible_histogram": feasible_histogram,
        "maximum_histogram": maximum_histogram,
        "flag4_maximum_histogram": flag4_maximum_histogram,
        "maximum_by_order": maximum_by_order,
        "profile_feasible": profile_feasible,
        "confusion4": confusion4,
        "confusion7": confusion7,
        "exact_incidence_total": exact_incidence_total,
        "exact_ties": exact_ties,
        "flag_ties": flag_ties,
        "exact_natural_flips": exact_natural_flips,
        "flag_natural_flips": flag_natural_flips,
        "gauge_edge_flips": gauge_edge_flips,
        "exact_digest": exact_digest.hexdigest(),
        "flag4_digest": flag4_digest.hexdigest(),
        "flag7_digest": flag7_digest.hexdigest(),
    }


def format_counter(counter):
    return " ".join(f"{key}:{value}" for key, value in sorted(counter.items()))


def format_profiles(counter):
    return " ".join(
        f"{profile}:{count}" for profile, count in sorted(counter.items())
    )


def main():
    tables = build_tables()
    grammar, weighted_states, grammar_digest = build_grammar()
    scalar = scalar_bank(grammar, tables["cards"])
    orbits = multiplication_orbits(scalar["rows"])
    audit = audit_owners(scalar, tables["normalized"])

    print("scale-twenty-eight AP-centred Hamming-six mod-four flag certificate")
    print("scope: primitive proper common-scale H6 owner-local gate only")
    print(f"orders={ORDERS} exact-order-residues={tuple(len(UNITS[d]) for d in ORDERS)}")
    for order in ORDERS:
        print(f"D{order} ratio-cardinalities={tables['rows'][order]}")
    print(
        "hereditary-grammar=at-least-two-v4-maximal-and-at-least-two-v7-maximal "
        f"words={len(grammar)} weighted-unit-words={weighted_states}"
    )
    print(
        f"literal-state-contexts={924 * weighted_states} "
        f"support-order-rows={924 * len(grammar)}"
    )
    print(f"scalar-feasible-owner-counts={scalar['feasible_owner_histogram']}")
    print(
        f"scalar-survivors={len(scalar['rows'])} supports=206 "
        f"support-histogram={format_counter(scalar['support_histogram'])}"
    )
    print(f"scalar-multiplicities={format_profiles(scalar['profile_histogram'])}")
    print(
        f"scalar-slack min={format_counter(scalar['minimum_slack'])} "
        f"max={format_counter(scalar['maximum_slack'])}"
    )
    print(f"multiplication-orbits={format_counter(orbits)}")
    print(
        f"normalized-owner-keys={audit['exact_keys']} "
        f"exact-feasible-owners-per-row={format_counter(audit['feasible_histogram'])}"
    )
    print(
        f"normalized-owner-key-orders={format_counter(audit['key_orders'])} "
        f"feasible-key-orders={format_counter(audit['feasible_key_orders'])}"
    )
    print(
        f"normalized-bank-size-range={min(audit['exact_bank_sizes'])}..{max(audit['exact_bank_sizes'])} "
        f"mod4-state-counts={format_counter(audit['mod4_state_sizes'])} "
        f"mod7-state-range={min(audit['mod7_state_sizes'])}..{max(audit['mod7_state_sizes'])}"
    )
    print(f"profile-feasible-owners={format_profiles(audit['profile_feasible'])}")
    print(f"exact-maximum-union={format_counter(audit['maximum_histogram'])}")
    print(f"mod4-flag-maximum={format_counter(audit['flag4_maximum_histogram'])}")
    print(f"exact-maximum-by-owner-order={format_counter(audit['maximum_by_order'])}")
    print(f"mod4-exact-contingency={format_counter(audit['confusion4'])}")
    print(f"mod7-exact-contingency={format_counter(audit['confusion7'])}")
    print(
        "owner-obstruction=mod4 saturated-capacity flag is exact on all 19020 "
        "owner rows; only D2 owners are feasible; every scalar row has at most "
        "two feasible owners; all-six-owner-row=0"
    )
    print(
        f"exact-reachable-mask-incidences={audit['exact_incidence_total']} "
        f"flag-equivalence-checks={6 * len(scalar['rows'])}"
    )
    print(
        "owner-tournament pair-observable exact=(feasible,max-union,capacity,bank-size) "
        "flag=(feasible,mod4-score,capacity,state-count)"
    )
    print(
        "owner-tournament switch=lexicographic; tie-Hamiltonian-path=coordinate-order; "
        "scores=(0,1,2,3,4,5) directed-triangles=0 SCCs=(1,1,1,1,1,1) "
        "Hamiltonian-paths=1 rows=3170"
    )
    print(
        f"owner-tournament exact-ties={format_counter(audit['exact_ties'])} "
        f"flag-ties={format_counter(audit['flag_ties'])} "
        f"gauge-edge-flips={format_counter(audit['gauge_edge_flips'])}"
    )
    print(
        f"owner-tournament natural-order-flips exact={format_counter(audit['exact_natural_flips'])} "
        f"flag={format_counter(audit['flag_natural_flips'])}"
    )
    print(
        "alternate-carrier=four Z/4 sheet fibres of size seven with saturated "
        "provider-capacity flags; D1,D2,D4 masks are whole-fibre hyperedges"
    )
    print(
        "flag-preserves=owner-cover-feasibility on the complete scalar survivor bank; "
        "destroys=within-fibre points, overlap multiplicity, exact units, and witnesses"
    )
    print(
        "challenged-assumptions=vertices need not be runners or owners; residues, "
        "provider masks, and proof obligations were compared; mod7 has 6072 false positives"
    )
    print(
        f"SHA256 crt-bases={tables['base_digest']} literal-masks={tables['mask_digest']} "
        f"projection-signatures={tables['signature_digest']}"
    )
    print(f"SHA256 grammar={grammar_digest} scalar={scalar['digest']}")
    print(
        f"SHA256 exact-owner-banks={audit['exact_digest']} "
        f"mod4-flags={audit['flag4_digest']} mod7-flags={audit['flag7_digest']}"
    )
    print("verdict=scale-28 primitive proper AP-centred common-sheet H6 face empty")


if __name__ == "__main__":
    main()
