#!/usr/bin/env python3
"""Independent exact referee for common-scale-twenty-one Hamming six.

The divisor/unit grammar and CRT representatives are built algebraically with
the Python standard library.  The complete hereditary grammar is enumerated,
scalar owner capacity is tested on all 924 labelled supports, and every
surviving owner row is replayed as an immutable set of exact 21-sheet union
masks.  All certificates use explicitly ordered SHA-256 serializations.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 21
LABELS = tuple(range(1, P))
DIVISORS = (1, 3, 7, 21)
FULL = (1 << C) - 1
MASK_BYTES = (C + 7) // 8
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
    """Optimization-stable replacement for a proof-bearing assertion."""
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    """Return the representative in (-modulus/2, modulus/2]."""
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label, state_index):
    """Solve x=D*label (mod 13), x=unit (mod D) without search."""
    order, unit = STATES[state_index]
    coefficient = (label - unit * pow(order, -1, P)) % P
    answer = unit + order * coefficient
    require(answer % P == order * label % P, "mod-thirteen CRT failure")
    require(answer % order == unit % order, "effective-order CRT failure")
    require(0 <= answer < P * order, "CRT representative outside range")
    return answer


def build_masks():
    """Construct and hash every provider/state/owner 21-sheet mask."""
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
                digest.update(mask.to_bytes(MASK_BYTES, "little"))
    return masks, digest.hexdigest()


def hereditary(word):
    """Every leave-one-out lcm must still equal the common scale."""
    return all(
        lcm(
            *(
                DIVISORS[word[index]]
                for index in range(6)
                if index != omitted
            )
        )
        == C
        for omitted in range(6)
    )


def grammar_census():
    """Enumerate and hash order words and all literal unit refinements."""
    words = tuple(
        word
        for word in product(range(len(DIVISORS)), repeat=6)
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
    """Exact owner-local union DP, immutable after every provider."""
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
    """Return sorted SCC sizes for a six-vertex tournament."""
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
    """Orient lexicographic owner summaries; coordinate order breaks ties."""
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
        values = key if isinstance(key, tuple) else (key,)
        digest.update(bytes(values))
        digest.update(count.to_bytes(8, "little"))
    return digest.hexdigest()


def tuple_histogram_text(counter):
    return " ".join(
        f"{','.join(str(value) for value in key)}:{counter[key]}"
        for key in sorted(counter)
    )


def main():
    require(len(STATES) == C, "literal state alphabet mismatch")
    for order_index, divisor in enumerate(DIVISORS):
        actual = {STATES[state][1] for state in STATE_INDICES[order_index]}
        expected = (
            {0}
            if divisor == 1
            else {
                unit
                for unit in range(1, divisor)
                if gcd(unit, divisor) == 1
            }
        )
        require(actual == expected, "unit grammar mismatch")
        require(
            all(
                STATES[state][0] == divisor
                for state in STATE_INDICES[order_index]
            ),
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
                    {size}
                    == {
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
        # Six seven-bit fields are safe: six providers contribute at most 126.
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
            pattern = tuple(
                word.count(order_index)
                for order_index in range(len(DIVISORS))
            )
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
            multiply_support(representative, multiplier)
            for multiplier in LABELS
        }
        require(orbit <= scalar_supports, "multiplication leaves scalar bank")
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
    scalar_capacity_hist = Counter()
    reachable_count_hist = Counter()
    maximum_mask_count_hist = Counter()
    owner_vectors = Counter()
    reachable_bank_digest = sha256()
    owner_profile_digest = sha256()
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
            capacity = sum(
                cardinality[label][word[provider_index]][owner]
                for provider_index, label in enumerate(labels)
            )
            require(capacity >= C, "scalar survivor has deficient owner")
            reachable = owner_reachable(masks, labels, word, owner)
            maximum = max(mask.bit_count() for mask in reachable)
            maximum_mask_count = sum(
                mask.bit_count() == maximum for mask in reachable
            )
            feasible = FULL in reachable
            require(feasible == (maximum == C), "owner threshold mismatch")

            row_digest = sha256()
            row_digest.update(len(reachable).to_bytes(4, "little"))
            for mask in sorted(reachable):
                row_digest.update(mask.to_bytes(MASK_BYTES, "little"))
            row_key = bytes(labels + word + (owner_index, owner))
            reachable_bank_digest.update(row_key)
            reachable_bank_digest.update(row_digest.digest())

            feasible_mask |= int(feasible) << owner_index
            feasible_rows += int(feasible)
            maximum_union_hist[maximum] += 1
            scalar_capacity_hist[capacity] += 1
            reachable_count_hist[len(reachable)] += 1
            maximum_mask_count_hist[maximum_mask_count] += 1
            maxima.append(maximum)
            summary = (
                int(feasible),
                maximum,
                capacity,
            )
            summaries.append(summary)
            owner_profile_digest.update(row_key)
            owner_profile_digest.update(bytes(summary))
            owner_profile_digest.update(len(reachable).to_bytes(4, "little"))
            owner_profile_digest.update(
                maximum_mask_count.to_bytes(4, "little")
            )
            owner_profile_digest.update(row_digest.digest())

        feasible_owner_hist[feasible_mask.bit_count()] += 1
        minimum_owner_union_hist[min(maxima)] += 1
        owner_vectors[tuple(maxima)] += 1
        owner_profile_digest.update(bytes((feasible_mask,) + tuple(maxima)))

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
    multiplicity_digest = tuple_histogram_digest(scalar_patterns)
    reachable_bank_digest_hex = reachable_bank_digest.hexdigest()
    owner_profile_digest_hex = owner_profile_digest.hexdigest()
    score_vector_digest = tuple_histogram_digest(tournament_score_vectors)

    # These expectations were frozen from the algebraic-CRT Python replay
    # before the pre-existing scale-twenty-one C++ source or output was read.
    require(len(order_words) == 3_249, "hereditary order-word census mismatch")
    require(
        weighted_states == 84_210_192 and literal_states == weighted_states,
        "literal state-word census mismatch",
    )
    require(
        len(supports) * literal_states == 77_810_217_408,
        "raw labelled state-context census mismatch",
    )
    require(
        mask_digest
        == "6677373c8360f1d89c751e9fe1a22ed4a7fe062c883efcadd2e782297902520e",
        "algebraic CRT mask-table digest mismatch",
    )
    require(
        order_digest
        == "f7c0254d8ac9108d318f4a9a21d0d2e5b244be91087b22c162bb563956e9b474",
        "order-grammar digest mismatch",
    )
    require(
        state_digest
        == "d67b07c19a008fa8ae6f292f65e071fb93f00412adb7d74db20d1480e93f1024",
        "literal-state-grammar digest mismatch",
    )
    require(
        len(scalar_bank) == 350
        and len(scalar_supports) == 80
        and len(scalar_patterns) == 9,
        "scalar bank census mismatch",
    )
    require(
        multiplicity_digest
        == "1416f868fd8b1e43ccff24303670993610479194a29ddd360c7d853b767b32f1",
        "scalar multiplicity digest mismatch",
    )
    require(
        all_support_context_hist
        == Counter({0: 844, 1: 14, 2: 48, 8: 12, 24: 6}),
        "all-support scalar-context histogram mismatch",
    )
    require(
        orbit_size_hist == Counter({2: 1, 6: 1, 12: 6}),
        "support multiplication-orbit census mismatch",
    )
    require(
        orbit_context_hist
        == Counter(
            {
                (2, 1): 1,
                (6, 24): 1,
                (12, 1): 1,
                (12, 2): 4,
                (12, 8): 1,
            }
        ),
        "support-orbit context-multiplicity histogram mismatch",
    )
    require(
        scalar_digest_hex
        == "6af0546c482c0d7a0981b5450c4fef51f3f2293566cc6d1df69abc9dfd86d2b5",
        "scalar-bank digest mismatch",
    )
    require(feasible_rows == 384, "feasible owner-row census mismatch")
    require(
        scalar_capacity_hist
        == Counter({21: 84, 22: 720, 23: 264, 24: 144,
                    26: 504, 27: 288, 30: 96}),
        "scalar owner-capacity histogram mismatch",
    )
    require(
        feasible_owner_hist == Counter({0: 182, 1: 96, 4: 72}),
        "feasible-owner histogram mismatch",
    )
    require(
        maximum_union_hist
        == Counter({17: 1_056, 18: 648, 20: 12, 21: 384}),
        "maximum-union histogram mismatch",
    )
    require(
        minimum_owner_union_hist == Counter({17: 348, 20: 2}),
        "minimum owner-maximum histogram mismatch",
    )
    require(
        reachable_count_hist
        == Counter(
            {
                43: 288,
                63: 36,
                82: 96,
                83: 96,
                112: 1_080,
                123: 36,
                124: 96,
                153: 24,
                165: 48,
                177: 12,
                297: 24,
                509: 60,
                515: 24,
                521: 24,
                857: 24,
                863: 12,
                1_167: 12,
                1_215: 24,
                1_253: 24,
                1_289: 48,
                81_986: 12,
            }
        ),
        "reachable-mask-count histogram mismatch",
    )
    require(
        maximum_mask_count_hist
        == Counter({1: 384, 15: 144, 18: 12, 20: 288,
                    30: 1_080, 40: 96, 70: 96}),
        "maximum-mask-count histogram mismatch",
    )
    require(len(owner_vectors) == 34, "owner maximum-vector census mismatch")
    require(max(feasible_owner_hist) == 4, "all-owner local context survived")
    require(
        reachable_bank_digest_hex
        == "27bb8bf24629b0bafb836d231a9753011c582711c5ff5f0c8a54315e5d0919c2",
        "sorted reachable-mask bank digest mismatch",
    )
    require(
        owner_profile_digest_hex
        == "33574db3cb36927439d77b19f0fadfedfc0ab8b8638a5d66368039f8464cbb82",
        "owner-profile digest mismatch",
    )
    require(
        tournament_ties == Counter({2: 72, 4: 24, 6: 36, 7: 216, 15: 2}),
        "tournament tie-edge histogram mismatch",
    )
    require(
        tournament_flips
        == Counter(
            {0: 19, 1: 20, 2: 27, 3: 29, 4: 69,
             5: 56, 6: 47, 7: 45, 8: 30, 9: 8}
        ),
        "tournament edge-flip histogram mismatch",
    )
    require(
        tournament_scores == Counter({score: 350 for score in range(6)}),
        "tournament aggregate score histogram mismatch",
    )
    require(
        len(tournament_score_vectors) == 79
        and score_vector_digest
        == "a294eded1ff8ae2cdf009b283757a7253b58b4969db7f64e94e2f487bdeae7e9",
        "tournament score-vector certificate mismatch",
    )
    require(
        tournament_triangles == Counter({0: 350}),
        "tournament directed-cycle census mismatch",
    )
    require(
        tournament_sccs == Counter({(1, 1, 1, 1, 1, 1): 350}),
        "tournament SCC fingerprint mismatch",
    )
    require(
        tournament_paths == Counter({1: 350}),
        "tournament Hamiltonian-path census mismatch",
    )

    print("scale-twenty-one independent algebraic-CRT Python referee")
    print("divisor grammar 1,3,7,21; literal states 21; supports 924")
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
    print(
        "scalar multiplicity (n1,n3,n7,n21) histogram "
        + tuple_histogram_text(scalar_patterns)
    )
    print(f"scalar-multiplicity SHA256 {multiplicity_digest}")
    print("all-support contexts histogram " + histogram_text(all_support_context_hist))
    print("multiplication orbit-size histogram " + histogram_text(orbit_size_hist))
    print(
        "multiplication orbit (size,contexts/support) histogram "
        + histogram_text(orbit_context_hist)
    )
    print(f"scalar-bank SHA256 {scalar_digest_hex}")
    print(f"owner rows {len(scalar_bank) * 6}; feasible {feasible_rows}")
    print("scalar owner-capacity histogram " + histogram_text(scalar_capacity_hist))
    print("feasible-owner/context histogram " + histogram_text(feasible_owner_hist))
    print("maximum-union histogram " + histogram_text(maximum_union_hist))
    print(
        "minimum-of-owner-maxima histogram "
        + histogram_text(minimum_owner_union_hist)
    )
    print("reachable-mask-count histogram " + histogram_text(reachable_count_hist))
    print(
        "maximum-mask-count histogram "
        + histogram_text(maximum_mask_count_hist)
    )
    print(
        f"distinct owner maximum vectors {len(owner_vectors)}; "
        f"sorted-reachable-banks SHA256 {reachable_bank_digest_hex}; "
        f"owner-profile SHA256 {owner_profile_digest_hex}"
    )
    print(
        "every owner row hashes its complete mask set in increasing-mask order; "
        f"owner-local all-six contexts {feasible_owner_hist.get(6, 0)}; "
        "hence global literal unit fibres 0"
    )
    print(
        "tournament vertices owner obligations; pair observable exact "
        "(feasible,maximum,scalar-capacity) summaries; lexicographic switch, "
        "coordinate tie Hamiltonian path"
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
        "preserved audit: the owner-obligation vector preserves each labelled "
        "FULL-mask predicate, maximum union deficit, and scalar capacity; the "
        "summary tournament preserves only their lexicographic ordering"
    )
    print(
        "lost audit: tournament orientation discards absolute values, exact "
        "masks, sheet identities, unit witnesses, and simultaneous witness "
        "compatibility, so it cannot prove the all-owner LRC predicate"
    )
    print(
        "challenged vertex assumption: providers, gaps, fixed sections and "
        "boundaries, wall events, residues, cover arcs, Fourier modes, matroid "
        "circuits, and proof obligations were considered; owner obligations "
        "alone retain the labelled terminal predicate, while all listed "
        "quotients destroy some owner/unit or sheet incidence"
    )


if __name__ == "__main__":
    main()
