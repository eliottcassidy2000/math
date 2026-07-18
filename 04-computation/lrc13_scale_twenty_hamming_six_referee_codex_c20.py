#!/usr/bin/env python3
"""Independent exact referee for the common-scale-twenty Hamming-six face.

The effective-order/unit grammar and every CRT representative are constructed
algebraically.  Scalar owner capacity filters the complete support/order bank;
the surviving owner projections are then rebuilt as immutable sets of exact
twenty-sheet masks.  SHA-256 certificates use only explicitly ordered byte
serializations and therefore do not depend on Python hash iteration order.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 20
LABELS = tuple(range(1, P))
DIVISORS = (1, 2, 4, 5, 10, 20)
FULL = (1 << C) - 1
STATES = tuple(
    (divisor, unit)
    for divisor in DIVISORS
    for unit in range(divisor)
    if (divisor == 1 and unit == 0)
    or (unit != 0 and gcd(unit, divisor) == 1)
)
STATE_INDICES = tuple(
    tuple(
        state_index
        for state_index, (order, _unit) in enumerate(STATES)
        if order == divisor
    )
    for divisor in DIVISORS
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label, state_index):
    """Solve x = D*label (mod 13), x = unit (mod D) directly."""
    order, unit = STATES[state_index]
    coefficient = (label - unit * pow(order, -1, P)) % P
    answer = unit + order * coefficient
    require(answer % P == order * label % P, "mod-thirteen CRT failure")
    require(answer % order == unit % order, "effective-order CRT failure")
    require(0 <= answer < P * order, "CRT representative outside range")
    return answer


def build_masks():
    masks = [[[0] * P for _state in STATES] for _label in range(P)]
    digest = sha256()
    for label in LABELS:
        for state_index, (order, _unit) in enumerate(STATES):
            base = algebraic_crt_base(label, state_index)
            for owner in LABELS:
                owner_inverse = pow(owner, -1, P)
                mask = 0
                for sheet in range(C):
                    value = centered(
                        base * (owner_inverse + P * sheet), P * order
                    )
                    if -order < value <= order:
                        mask |= 1 << sheet
                masks[label][state_index][owner] = mask
                digest.update(mask.to_bytes(3, "little"))
    return masks, digest.hexdigest()


def hereditary(word):
    return all(
        lcm(*(DIVISORS[word[index]] for index in range(6) if index != omitted))
        == C
        for omitted in range(6)
    )


def grammar_census():
    words = tuple(
        word for word in product(range(len(DIVISORS)), repeat=6)
        if hereditary(word)
    )
    order_digest = sha256()
    for word in words:
        order_digest.update(bytes(word))

    weighted = sum(
        prod(len(STATE_INDICES[order_index]) for order_index in word)
        for word in words
    )
    literal_count = 0
    literal_digest = sha256()
    buffer = bytearray()
    for word in words:
        fibres = tuple(STATE_INDICES[order_index] for order_index in word)
        for states in product(*fibres):
            literal_count += 1
            buffer.extend(states)
            if len(buffer) >= 1 << 20:
                literal_digest.update(buffer)
                buffer.clear()
    literal_digest.update(buffer)
    return (
        words,
        weighted,
        literal_count,
        order_digest.hexdigest(),
        literal_digest.hexdigest(),
    )


def owner_reachable(masks, labels, word, owner):
    reachable = frozenset((0,))
    for provider_index, label in enumerate(labels):
        options = frozenset(
            masks[label][state_index][owner]
            for state_index in STATE_INDICES[word[provider_index]]
        )
        reachable = frozenset(
            partial | option
            for partial in reachable
            for option in options
        )
    return reachable


def strongly_connected_component_sizes(out):
    reach = list(out)
    for vertex in range(6):
        reach[vertex] |= 1 << vertex
    for middle in range(6):
        middle_bit = 1 << middle
        for source in range(6):
            if reach[source] & middle_bit:
                reach[source] |= reach[middle]
    unused = set(range(6))
    sizes = []
    while unused:
        root = min(unused)
        component = {
            vertex
            for vertex in unused
            if (reach[root] >> vertex) & 1 and (reach[vertex] >> root) & 1
        }
        require(component, "empty tournament SCC")
        sizes.append(len(component))
        unused -= component
    return tuple(sorted(sizes))


def tournament_fingerprint(summaries):
    """Lexicographic owner-summary switch with the coordinate tie path."""
    out = [0] * 6
    scores = [0] * 6
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        winner = left
        if summaries[left] == summaries[right]:
            ties += 1
        elif summaries[right] > summaries[left]:
            winner = right
            flips += 1
        loser = left + right - winner
        out[winner] |= 1 << loser
        scores[winner] += 1

    triangles = 0
    for first, second, third in combinations(range(6), 3):
        forward = (
            (out[first] >> second)
            & (out[second] >> third)
            & (out[third] >> first)
            & 1
        )
        reverse = (
            (out[first] >> third)
            & (out[third] >> second)
            & (out[second] >> first)
            & 1
        )
        triangles += bool(forward or reverse)

    paths = [[0] * 6 for _mask in range(1 << 6)]
    for last in range(6):
        paths[1 << last][last] = 1
    for mask in range(1, 1 << 6):
        for last in range(6):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            for previous in range(6):
                if (
                    (previous_mask >> previous) & 1
                    and (out[previous] >> last) & 1
                ):
                    paths[mask][last] += paths[previous_mask][previous]
    return (
        ties,
        flips,
        tuple(scores),
        triangles,
        strongly_connected_component_sizes(out),
        sum(paths[-1]),
    )


def multiply_support(labels, multiplier):
    return tuple(sorted(multiplier * label % P for label in labels))


def histogram_text(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def tuple_histogram_digest(counter):
    digest = sha256()
    for key, count in sorted(counter.items()):
        digest.update(bytes(key))
        digest.update(count.to_bytes(8, "little"))
    return digest.hexdigest()


def main():
    require(len(STATES) == 20, "literal state alphabet mismatch")
    for order_index, divisor in enumerate(DIVISORS):
        actual = {STATES[state][1] for state in STATE_INDICES[order_index]}
        expected = (
            {0}
            if divisor == 1
            else {unit for unit in range(1, divisor) if gcd(unit, divisor) == 1}
        )
        require(actual == expected, "unit grammar mismatch")
        require(
            all(STATES[state][0] == divisor for state in STATE_INDICES[order_index]),
            "state/order index mismatch",
        )

    masks, mask_digest = build_masks()
    cardinality = [[[0] * P for _order in DIVISORS] for _label in range(P)]
    for label in LABELS:
        for order_index in range(len(DIVISORS)):
            for owner in LABELS:
                sizes = {
                    masks[label][state][owner].bit_count()
                    for state in STATE_INDICES[order_index]
                }
                require(len(sizes) == 1, "mask cardinality depends on unit")
                size = next(iter(sizes))
                cardinality[label][order_index][owner] = size
                ratio = label * pow(owner, -1, P) % P
                require(
                    {size} == {
                        masks[ratio][state][1].bit_count()
                        for state in STATE_INDICES[order_index]
                    },
                    "provider/owner ratio reduction failure",
                )

    (
        order_words,
        weighted_states,
        literal_states,
        order_digest,
        state_digest,
    ) = grammar_census()
    require(weighted_states == literal_states, "weighted/literal grammar mismatch")

    supports = tuple(combinations(LABELS, 6))
    require(len(supports) == 924, "support census mismatch")
    scalar_bank = []
    scalar_supports = set()
    contexts_per_support = Counter()
    scalar_patterns = Counter()
    scalar_digest = sha256()

    for labels in supports:
        # Each seven-bit field is safe: six providers contribute at most 120.
        packed = [[0] * len(DIVISORS) for _provider in range(6)]
        for provider_index, label in enumerate(labels):
            for order_index in range(len(DIVISORS)):
                packed[provider_index][order_index] = sum(
                    cardinality[label][order_index][owner]
                    << (7 * owner_index)
                    for owner_index, owner in enumerate(labels)
                )

        for word in order_words:
            capacity = sum(
                packed[provider_index][word[provider_index]]
                for provider_index in range(6)
            )
            if any(
                ((capacity >> (7 * owner_index)) & 127) < C
                for owner_index in range(6)
            ):
                continue
            scalar_bank.append((labels, word))
            scalar_supports.add(labels)
            contexts_per_support[labels] += 1
            pattern = tuple(word.count(index) for index in range(len(DIVISORS)))
            scalar_patterns[pattern] += 1
            scalar_digest.update(bytes(labels + word))

    all_support_context_hist = Counter(
        contexts_per_support.get(labels, 0) for labels in supports
    )
    remaining = set(scalar_supports)
    orbit_size_hist = Counter()
    orbit_context_hist = Counter()
    while remaining:
        representative = min(remaining)
        orbit = {
            multiply_support(representative, multiplier) for multiplier in LABELS
        }
        require(orbit <= scalar_supports, "multiplication leaves scalar support bank")
        require(orbit <= remaining, "multiplication support orbits overlap")
        multiplicities = {contexts_per_support[support] for support in orbit}
        require(
            len(multiplicities) == 1,
            "context multiplicity changes within multiplication orbit",
        )
        orbit_size_hist[len(orbit)] += 1
        orbit_context_hist[(len(orbit), multiplicities.pop())] += 1
        remaining -= orbit

    feasible_rows = 0
    feasible_owner_hist = Counter()
    maximum_union_hist = Counter()
    minimum_owner_union_hist = Counter()
    reachable_count_hist = Counter()
    maximum_mask_count_hist = Counter()
    owner_vectors = Counter()
    owner_digest = sha256()
    tournament_ties = Counter()
    tournament_flips = Counter()
    tournament_scores = Counter()
    tournament_score_vectors = Counter()
    tournament_triangles = Counter()
    tournament_sccs = Counter()
    tournament_paths = Counter()

    for labels, word in scalar_bank:
        summaries = []
        maxima = []
        feasible_mask = 0
        for owner_index, owner in enumerate(labels):
            reachable = owner_reachable(masks, labels, word, owner)
            maximum = max(mask.bit_count() for mask in reachable)
            maximum_mask_count = sum(
                mask.bit_count() == maximum for mask in reachable
            )
            feasible = FULL in reachable
            require(feasible == (maximum == C), "owner threshold mismatch")
            feasible_mask |= int(feasible) << owner_index
            feasible_rows += feasible
            maximum_union_hist[maximum] += 1
            reachable_count_hist[len(reachable)] += 1
            maximum_mask_count_hist[maximum_mask_count] += 1
            maxima.append(maximum)
            summaries.append(
                (int(feasible), maximum, len(reachable), maximum_mask_count)
            )

            reachable_digest = sha256()
            for mask in sorted(reachable):
                reachable_digest.update(mask.to_bytes(3, "little"))
            owner_digest.update(bytes(labels + word + (owner_index, owner)))
            owner_digest.update(bytes((int(feasible), maximum)))
            owner_digest.update(len(reachable).to_bytes(4, "little"))
            owner_digest.update(maximum_mask_count.to_bytes(4, "little"))
            owner_digest.update(reachable_digest.digest())

        feasible_owner_hist[feasible_mask.bit_count()] += 1
        minimum_owner_union_hist[min(maxima)] += 1
        owner_vectors[tuple(maxima)] += 1
        owner_digest.update(bytes((feasible_mask,) + tuple(maxima)))

        ties, flips, scores, triangles, sccs, paths = tournament_fingerprint(
            tuple(summaries)
        )
        tournament_ties[ties] += 1
        tournament_flips[flips] += 1
        tournament_scores.update(scores)
        tournament_score_vectors[scores] += 1
        tournament_triangles[triangles] += 1
        tournament_sccs[sccs] += 1
        tournament_paths[paths] += 1

    scalar_digest_hex = scalar_digest.hexdigest()
    owner_digest_hex = owner_digest.hexdigest()
    pattern_digest = tuple_histogram_digest(scalar_patterns)
    score_vector_digest = tuple_histogram_digest(tournament_score_vectors)

    require(len(order_words) == 26_961, "hereditary order-word census mismatch")
    require(
        weighted_states == 56_908_800 and literal_states == weighted_states,
        "literal state-word census mismatch",
    )
    require(
        mask_digest
        == "db45c23efceded1efcfe72556ab64337ecc96c016ad5427caf3a28cad02b7892",
        "algebraic CRT mask-table digest mismatch",
    )
    require(
        order_digest
        == "788eab466cbe0b5f2b0baee84bcb259357341af8c280426bce39b45b7679620e",
        "order-grammar digest mismatch",
    )
    require(
        state_digest
        == "e5fcb44c7c1a0100505ce672ff590fe671774850ab666573b46cb090c149f6f8",
        "literal-state-grammar digest mismatch",
    )
    require(
        len(scalar_bank) == 12_584
        and len(scalar_supports) == 830
        and len(scalar_patterns) == 65,
        "scalar bank census mismatch",
    )
    require(
        pattern_digest
        == "362f3f85efa803f064dba65056f81b39ac50e94de39ceeac352b8afd77b1d6e6",
        "scalar multiplicity digest mismatch",
    )
    require(
        all_support_context_hist
        == Counter(
            {
                0: 94,
                1: 12,
                2: 96,
                3: 24,
                4: 96,
                6: 48,
                7: 96,
                13: 216,
                15: 6,
                17: 24,
                18: 24,
                19: 24,
                25: 24,
                29: 24,
                30: 6,
                38: 48,
                44: 24,
                61: 24,
                65: 12,
                85: 2,
            }
        ),
        "all-support scalar-context histogram mismatch",
    )
    require(
        orbit_size_hist == Counter({2: 1, 6: 2, 12: 68}),
        "support multiplication-orbit census mismatch",
    )
    require(
        orbit_context_hist
        == Counter(
            {
                (2, 85): 1,
                (6, 15): 1,
                (6, 30): 1,
                (12, 1): 1,
                (12, 2): 8,
                (12, 3): 2,
                (12, 4): 8,
                (12, 6): 4,
                (12, 7): 8,
                (12, 13): 18,
                (12, 17): 2,
                (12, 18): 2,
                (12, 19): 2,
                (12, 25): 2,
                (12, 29): 2,
                (12, 38): 4,
                (12, 44): 2,
                (12, 61): 2,
                (12, 65): 1,
            }
        ),
        "support-orbit context-multiplicity histogram mismatch",
    )
    require(
        scalar_digest_hex
        == "68aab474423d01352a20bf5a1c769634e9bbfd9bb030845d0c04972f31648167",
        "scalar-bank digest mismatch",
    )
    require(feasible_rows == 3_720, "feasible owner-row census mismatch")
    require(
        feasible_owner_hist == Counter({0: 9_632, 1: 2_184, 2: 768}),
        "feasible-owner histogram mismatch",
    )
    require(
        maximum_union_hist
        == Counter(
            {
                12: 1_320,
                13: 2_904,
                14: 9_720,
                15: 15_924,
                16: 22_800,
                17: 12_876,
                18: 5_700,
                19: 540,
                20: 3_720,
            }
        ),
        "maximum-union histogram mismatch",
    )
    require(
        minimum_owner_union_hist
        == Counter(
            {
                12: 1_194,
                13: 1_422,
                14: 4_282,
                15: 3_748,
                16: 1_668,
                17: 222,
                18: 48,
            }
        ),
        "minimum owner-maximum histogram mismatch",
    )
    require(len(owner_vectors) == 4_061, "owner maximum-vector census mismatch")
    require(max(feasible_owner_hist) == 2, "owner-local all-six context survived")
    require(
        owner_digest_hex
        == "16af5ce1ce6743d7ca5dcd5cdd349776df3cfe3b00a23b36370e394952888dc4",
        "owner reachable-bank digest mismatch",
    )
    require(
        tournament_ties
        == Counter({0: 2_112, 1: 7_272, 2: 2_796, 3: 300, 4: 24, 6: 32, 7: 48}),
        "tournament tie-edge histogram mismatch",
    )
    require(
        tournament_flips
        == Counter(
            {
                0: 34,
                1: 91,
                2: 237,
                3: 472,
                4: 955,
                5: 1_430,
                6: 2_151,
                7: 2_294,
                8: 1_899,
                9: 1_418,
                10: 805,
                11: 461,
                12: 217,
                13: 88,
                14: 32,
            }
        ),
        "tournament edge-flip histogram mismatch",
    )
    require(
        tournament_scores == Counter({score: 12_584 for score in range(6)}),
        "tournament aggregate score histogram mismatch",
    )
    require(
        len(tournament_score_vectors) == 710
        and score_vector_digest
        == "73eaa91302f457d399b83966da894a03c1457fd406bc3a4d5207eb0bf50eb6f3",
        "tournament score-vector certificate mismatch",
    )
    require(
        tournament_triangles == Counter({0: 12_584}),
        "tournament directed-cycle census mismatch",
    )
    require(
        tournament_sccs == Counter({(1, 1, 1, 1, 1, 1): 12_584}),
        "tournament SCC fingerprint mismatch",
    )
    require(
        tournament_paths == Counter({1: 12_584}),
        "tournament Hamiltonian-path census mismatch",
    )

    print("scale-twenty independent algebraic-CRT Python referee")
    print("divisor grammar 1,2,4,5,10,20; literal states 20; supports 924")
    print(
        f"hereditary order words {len(order_words)}; "
        f"state words/support {literal_states}; raw {len(supports) * literal_states}"
    )
    print(f"mask-table SHA256 {mask_digest}")
    print(f"order-grammar SHA256 {order_digest}")
    print(f"literal-state-grammar SHA256 {state_digest}")
    print(
        f"scalar contexts {len(scalar_bank)} on {len(scalar_supports)} supports; "
        f"multiplicity patterns {len(scalar_patterns)}"
    )
    print(f"scalar-multiplicity SHA256 {pattern_digest}")
    print("all-support contexts histogram " + histogram_text(all_support_context_hist))
    print("multiplication orbit-size histogram " + histogram_text(orbit_size_hist))
    print(
        "multiplication orbit (size,contexts/support) histogram "
        + histogram_text(orbit_context_hist)
    )
    print(f"scalar-bank SHA256 {scalar_digest_hex}")
    print(f"owner rows {len(scalar_bank) * 6}; feasible {feasible_rows}")
    print("feasible-owner/context histogram " + histogram_text(feasible_owner_hist))
    print("maximum-union histogram " + histogram_text(maximum_union_hist))
    print(
        "minimum-of-owner-maxima histogram "
        + histogram_text(minimum_owner_union_hist)
    )
    print(
        f"distinct owner maximum vectors {len(owner_vectors)}; "
        f"owner-reachable-bank SHA256 {owner_digest_hex}"
    )
    print(
        "exact reachable masks are hashed in sorted order for every owner row; "
        f"owner-local all-six contexts {feasible_owner_hist.get(6, 0)}; "
        "hence global literal unit fibres 0"
    )
    print(
        "tournament vertices owner obligations; pair observable exact "
        "(feasible,maximum,reachable-count,maximum-mask-count) summaries; "
        "lexicographic switch, coordinate tie Hamiltonian path"
    )
    print("tournament tie-edge histogram " + histogram_text(tournament_ties))
    print("tournament edge-flip histogram " + histogram_text(tournament_flips))
    print("tournament aggregate score histogram " + histogram_text(tournament_scores))
    print(
        f"tournament score vectors {len(tournament_score_vectors)}; "
        f"score-vector SHA256 {score_vector_digest}"
    )
    print(
        "tournament directed-cycle (triangle) histogram "
        + histogram_text(tournament_triangles)
    )
    print("tournament SCC-size histogram " + histogram_text(tournament_sccs))
    print(
        "tournament Hamiltonian-path-count histogram "
        + histogram_text(tournament_paths)
    )
    print(
        "preserved audit: the immutable owner banks preserve each FULL-mask "
        "predicate and maximum deficit; retained summaries preserve their "
        "feasibility, maximum, and reachability-count ordering"
    )
    print(
        "lost audit: the oriented tournament discards absolute values, exact "
        "masks, sheet identities, unit witnesses, and simultaneous witness "
        "compatibility, so it cannot prove the all-owner LRC predicate"
    )
    print(
        "challenged vertex assumption: owners/proof obligations were selected; "
        "providers, gaps, fixed sections and boundaries, wall events, residues, "
        "cover arcs, Fourier modes, and matroid circuits lose the labelled "
        "owner-local terminal predicate"
    )


if __name__ == "__main__":
    main()
