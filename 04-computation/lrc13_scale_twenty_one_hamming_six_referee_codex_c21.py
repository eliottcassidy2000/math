#!/usr/bin/env python3
"""Independent exact referee for the common-scale-twenty-one H6 face.

The C++ primary finds CRT representatives by literal search and stores unions
in an epoch-marked packed table.  This referee instead solves every CRT system
algebraically and carries each owner-local bank as an immutable Python set.
Every sorted reachable-mask bank is SHA-256 hashed.  Explicit byte encodings
make normal and ``python -O`` output directly comparable.

The terminal carrier is the labelled six-tuple of owner-local mask banks.  A
tournament on owner obligations is useful telemetry but is not faithful: its
lexicographic switch forgets absolute sheet thresholds and the masks that
realize its pairwise comparisons.
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

FNV_OFFSET = 14_695_981_039_346_656_037
FNV_PRIME = 1_099_511_628_211
FNV_MASK = (1 << 64) - 1


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fnv_bytes(digest, data):
    for value in data:
        digest ^= value
        digest = (digest * FNV_PRIME) & FNV_MASK
    return digest


def centered(value, modulus):
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
    masks = [[[0] * P for _state in STATES] for _label in range(P)]
    digest = sha256()
    primary_fnv = FNV_OFFSET
    for label in LABELS:
        for state_index, (order, _unit) in enumerate(STATES):
            base = algebraic_crt_base(label, state_index)
            for owner in LABELS:
                mask = 0
                owner_inverse = pow(owner, -1, P)
                for sheet in range(C):
                    value = centered(
                        base * (owner_inverse + P * sheet), P * order
                    )
                    if -order < value <= order:
                        mask |= 1 << sheet
                masks[label][state_index][owner] = mask
                encoded = mask.to_bytes(4, "little")
                primary_fnv = fnv_bytes(primary_fnv, encoded)
                digest.update(encoded[:3])
    return masks, digest.hexdigest(), primary_fnv


def analytic_cardinality(label, order_index, owner):
    """Count one D-period arithmetically, then repeat it 21/D times."""
    order = DIVISORS[order_index]
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * period_count


def hereditary(word):
    return all(
        lcm(*(DIVISORS[word[index]] for index in range(6) if index != omitted))
        == C
        for omitted in range(6)
    )


def hereditary_prime_provider_audit(word):
    three_providers = sum(DIVISORS[index] % 3 == 0 for index in word)
    seven_providers = sum(DIVISORS[index] % 7 == 0 for index in word)
    return three_providers >= 2 and seven_providers >= 2


def grammar_census():
    words = []
    order_digest = sha256()
    literal_descriptor_digest = sha256()
    primary_order_fnv = FNV_OFFSET
    weighted = 0
    fibre_histogram = Counter()

    for word in product(range(len(DIVISORS)), repeat=6):
        by_lcm = hereditary(word)
        require(
            by_lcm == hereditary_prime_provider_audit(word),
            "lcm and prime-provider hereditary audits disagree",
        )
        if not by_lcm:
            continue
        words.append(word)
        order_digest.update(bytes(word))
        primary_order_fnv = fnv_bytes(primary_order_fnv, bytes(word) + b"\xff")
        fibre = prod(len(STATE_INDICES[index]) for index in word)
        weighted += fibre
        fibre_histogram[fibre] += 1

        # A compact, independently parsed description of the complete literal
        # fibre.  The exact word count is the product of these option lengths.
        literal_descriptor_digest.update(bytes(word))
        for index in word:
            options = STATE_INDICES[index]
            literal_descriptor_digest.update(bytes((len(options),)))
            literal_descriptor_digest.update(bytes(options))
        literal_descriptor_digest.update(fibre.to_bytes(8, "little"))

    return (
        tuple(words),
        weighted,
        fibre_histogram,
        order_digest.hexdigest(),
        literal_descriptor_digest.hexdigest(),
        primary_order_fnv,
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
    """Primary-compatible owner-summary switch and coordinate tie path."""
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


def integer_histogram_digest(counter):
    digest = sha256()
    for key, count in sorted(counter.items()):
        digest.update(key.to_bytes(4, "little"))
        digest.update(count.to_bytes(8, "little"))
    return digest.hexdigest()


def main():
    require(len(STATES) == C, "literal state alphabet mismatch")
    require(
        tuple(len(indices) for indices in STATE_INDICES) == (1, 2, 6, 12),
        "unit-count grammar mismatch",
    )
    for order_index, divisor in enumerate(DIVISORS):
        actual = {STATES[state][1] for state in STATE_INDICES[order_index]}
        expected = (
            {0}
            if divisor == 1
            else {unit for unit in range(1, divisor) if gcd(unit, divisor) == 1}
        )
        require(actual == expected, "unit grammar mismatch")

    masks, mask_digest, primary_mask_fnv = build_masks()
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
                require(
                    size == analytic_cardinality(label, order_index, owner),
                    "literal and period-count cardinalities disagree",
                )
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
        literal_states,
        fibre_histogram,
        order_digest,
        literal_descriptor_digest,
        primary_order_fnv,
    ) = grammar_census()

    supports = tuple(combinations(LABELS, 6))
    require(len(supports) == 924, "support census mismatch")
    scalar_bank = []
    scalar_supports = set()
    contexts_per_support = Counter()
    scalar_patterns = Counter()
    capacity_vectors = Counter()
    minimum_slack_histogram = Counter()
    maximum_slack_histogram = Counter()
    tight_owner_histogram = Counter()
    scalar_digest = sha256()
    capacity_digest = sha256()
    primary_scalar_fnv = FNV_OFFSET
    primary_capacity_fnv = FNV_OFFSET

    for labels in supports:
        # Six providers contribute at most 126, so seven-bit packed fields do
        # not carry into the adjacent owner field.
        packed = [[0] * len(DIVISORS) for _provider in range(6)]
        for provider_index, label in enumerate(labels):
            for order_index in range(len(DIVISORS)):
                packed[provider_index][order_index] = sum(
                    cardinality[label][order_index][owner]
                    << (7 * owner_index)
                    for owner_index, owner in enumerate(labels)
                )

        for word in order_words:
            packed_capacity = sum(
                packed[provider_index][word[provider_index]]
                for provider_index in range(6)
            )
            capacities = tuple(
                (packed_capacity >> (7 * owner_index)) & 127
                for owner_index in range(6)
            )
            direct = tuple(
                sum(
                    cardinality[label][word[provider_index]][owner]
                    for provider_index, label in enumerate(labels)
                )
                for owner in labels
            )
            require(capacities == direct, "packed/direct capacities disagree")
            if min(capacities) < C:
                continue

            scalar_bank.append((labels, word, capacities))
            scalar_supports.add(labels)
            contexts_per_support[labels] += 1
            pattern = tuple(word.count(index) for index in range(len(DIVISORS)))
            scalar_patterns[pattern] += 1
            capacity_vectors[capacities] += 1
            slacks = tuple(value - C for value in capacities)
            minimum_slack_histogram[min(slacks)] += 1
            maximum_slack_histogram[max(slacks)] += 1
            tight_owner_histogram[slacks.count(0)] += 1

            scalar_row = bytes(labels + word)
            scalar_digest.update(scalar_row)
            capacity_digest.update(bytes(labels + word + capacities))
            primary_scalar_fnv = fnv_bytes(primary_scalar_fnv, scalar_row)
            primary_capacity_fnv = fnv_bytes(
                primary_capacity_fnv, bytes(capacities)
            )

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
        require(orbit <= scalar_supports, "multiplication leaves scalar bank")
        require(orbit <= remaining, "multiplication support orbits overlap")
        multiplicities = {contexts_per_support[support] for support in orbit}
        require(len(multiplicities) == 1, "orbit changes context multiplicity")
        orbit_size_hist[len(orbit)] += 1
        orbit_context_hist[(len(orbit), multiplicities.pop())] += 1
        remaining -= orbit

    primary_pattern_fnv = FNV_OFFSET
    for pattern, count in sorted(scalar_patterns.items()):
        primary_pattern_fnv = fnv_bytes(primary_pattern_fnv, bytes(pattern))
        primary_pattern_fnv = fnv_bytes(
            primary_pattern_fnv, count.to_bytes(8, "little")
        )
    primary_capacity_vector_fnv = FNV_OFFSET
    for vector, count in sorted(capacity_vectors.items()):
        primary_capacity_vector_fnv = fnv_bytes(
            primary_capacity_vector_fnv, bytes(vector)
        )
        primary_capacity_vector_fnv = fnv_bytes(
            primary_capacity_vector_fnv, count.to_bytes(8, "little")
        )

    feasible_rows = 0
    feasible_owner_hist = Counter()
    maximum_union_hist = Counter()
    minimum_owner_union_hist = Counter()
    reachable_count_hist = Counter()
    maximum_mask_count_hist = Counter()
    owner_vectors = Counter()
    owner_digest = sha256()
    primary_owner_fnv = FNV_OFFSET
    reachable_mask_total = 0
    all_order_rows = 0
    tournament_ties = Counter()
    tournament_flips = Counter()
    tournament_scores = Counter()
    tournament_score_vectors = Counter()
    tournament_triangles = Counter()
    tournament_sccs = Counter()
    tournament_paths = Counter()

    for labels, word, capacities in scalar_bank:
        summaries = []
        maxima = []
        reachable_counts = []
        feasible_mask = 0
        all_order = word == (len(DIVISORS) - 1,) * 6
        all_order_rows += all_order

        for owner_index, owner in enumerate(labels):
            reachable = owner_reachable(masks, labels, word, owner)
            ordered_reachable = sorted(reachable)
            maximum = max(mask.bit_count() for mask in ordered_reachable)
            maximum_mask_count = sum(
                mask.bit_count() == maximum for mask in ordered_reachable
            )
            feasible = FULL in reachable
            require(feasible == (maximum == C), "owner threshold mismatch")
            feasible_mask |= int(feasible) << owner_index
            feasible_rows += feasible
            maximum_union_hist[maximum] += 1
            reachable_count_hist[len(reachable)] += 1
            maximum_mask_count_hist[maximum_mask_count] += 1
            reachable_mask_total += len(reachable)
            maxima.append(maximum)
            reachable_counts.append(len(reachable))
            summaries.append((int(feasible), maximum, capacities[owner_index]))

            row_digest = sha256()
            for mask in ordered_reachable:
                row_digest.update(mask.to_bytes(3, "little"))
            owner_digest.update(bytes(labels + word + (owner_index, owner)))
            owner_digest.update(bytes((int(feasible), maximum)))
            owner_digest.update(len(reachable).to_bytes(4, "little"))
            owner_digest.update(maximum_mask_count.to_bytes(4, "little"))
            owner_digest.update(row_digest.digest())

        require(
            not all_order
            or (capacities == (C,) * 6 and tuple(maxima) == (C - 1,) * 6),
            "all-order overlap-only obstruction mismatch",
        )
        feasible_owner_hist[feasible_mask.bit_count()] += 1
        minimum_owner_union_hist[min(maxima)] += 1
        owner_vectors[tuple(maxima)] += 1

        primary_owner_fnv = fnv_bytes(primary_owner_fnv, bytes(labels + word))
        primary_owner_fnv = fnv_bytes(primary_owner_fnv, bytes((feasible_mask,)))
        for maximum, count in zip(maxima, reachable_counts):
            primary_owner_fnv = fnv_bytes(primary_owner_fnv, bytes((maximum,)))
            primary_owner_fnv = fnv_bytes(
                primary_owner_fnv, count.to_bytes(4, "little")
            )

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
    capacity_digest_hex = capacity_digest.hexdigest()
    owner_digest_hex = owner_digest.hexdigest()
    pattern_digest = tuple_histogram_digest(scalar_patterns)
    capacity_vector_digest = tuple_histogram_digest(capacity_vectors)
    score_vector_digest = tuple_histogram_digest(tournament_score_vectors)
    reachable_count_digest = integer_histogram_digest(reachable_count_hist)
    maximum_mask_count_digest = integer_histogram_digest(maximum_mask_count_hist)

    require(len(order_words) == 3_249, "hereditary word census mismatch")
    require(literal_states == 84_210_192, "literal state census mismatch")
    require(924 * literal_states == 77_810_217_408, "raw census mismatch")
    require(
        mask_digest
        == "6677373c8360f1d89c751e9fe1a22ed4a7fe062c883efcadd2e782297902520e",
        "algebraic CRT mask-table SHA mismatch",
    )
    require(
        order_digest
        == "f7c0254d8ac9108d318f4a9a21d0d2e5b244be91087b22c162bb563956e9b474",
        "order-grammar SHA mismatch",
    )
    require(
        literal_descriptor_digest
        == "be9236c7a255afeaf7ec6fe0c6716603d1be59818a8b4f2e71d2b030a6fb12b0",
        "literal-fibre descriptor SHA mismatch",
    )
    require(primary_mask_fnv == 0x88BD211E0C707E3D, "primary mask FNV mismatch")
    require(primary_order_fnv == 0xAE4E6CDBC98C0522, "primary order FNV mismatch")
    require(
        len(scalar_bank) == 350
        and len(scalar_supports) == 80
        and len(scalar_patterns) == 9,
        "scalar bank census mismatch",
    )
    require(primary_pattern_fnv == 0x0C84A8F8C824E11D, "pattern FNV mismatch")
    require(primary_scalar_fnv == 0xE835F3B670CFC659, "scalar-bank FNV mismatch")
    require(
        pattern_digest
        == "1416f868fd8b1e43ccff24303670993610479194a29ddd360c7d853b767b32f1",
        "multiplicity SHA mismatch",
    )
    require(
        scalar_digest_hex
        == "6af0546c482c0d7a0981b5450c4fef51f3f2293566cc6d1df69abc9dfd86d2b5",
        "scalar-bank SHA mismatch",
    )
    require(
        capacity_digest_hex
        == "6598c794af639163d894100f2da7e3b7388429529746b831cc7ff8c0e30ecf3d",
        "capacity-row SHA mismatch",
    )
    require(
        all_support_context_hist == Counter({0: 844, 1: 14, 2: 48, 8: 12, 24: 6}),
        "all-support context histogram mismatch",
    )
    require(orbit_size_hist == Counter({2: 1, 6: 1, 12: 6}), "orbit census mismatch")
    require(
        orbit_context_hist
        == Counter({(2, 1): 1, (6, 24): 1, (12, 1): 1, (12, 2): 4, (12, 8): 1}),
        "orbit/context census mismatch",
    )
    require(len(capacity_vectors) == 131, "capacity-vector census mismatch")
    require(
        capacity_vector_digest
        == "933431170549185541b83f15f5ccab42bd0fe545a52b26c21e8b4bf4c2cfe3a2",
        "capacity-vector SHA mismatch",
    )
    require(
        primary_capacity_vector_fnv == 0x6B9D23AFB54BC2AD,
        "capacity-vector FNV mismatch",
    )
    require(primary_capacity_fnv == 0x2039AF5E9A8EAB41, "row-capacity FNV mismatch")
    require(minimum_slack_histogram == Counter({0: 56, 1: 180, 2: 114}), "minimum slack mismatch")
    require(maximum_slack_histogram == Counter({0: 2, 5: 180, 6: 72, 9: 96}), "maximum slack mismatch")
    require(tight_owner_histogram == Counter({0: 294, 1: 36, 2: 18, 6: 2}), "tight-owner mismatch")
    require(feasible_rows == 384, "feasible owner-row census mismatch")
    require(feasible_owner_hist == Counter({0: 182, 1: 96, 4: 72}), "feasible-owner mismatch")
    require(maximum_union_hist == Counter({17: 1_056, 18: 648, 20: 12, 21: 384}), "maximum-union mismatch")
    require(minimum_owner_union_hist == Counter({17: 348, 20: 2}), "minimum owner-maximum mismatch")
    require(len(owner_vectors) == 34, "owner-vector census mismatch")
    require(primary_owner_fnv == 0x9E5738B21D4A7F95, "owner-profile FNV mismatch")
    require(
        owner_digest_hex
        == "c3e7ce9c46dffd5846bb663f9f490ac9bf3fae76885a7085d5c3172a2907a8aa",
        "exact reachable-bank SHA mismatch",
    )
    require(reachable_mask_total == 1_393_896, "reachable-mask total mismatch")
    require(
        reachable_count_digest
        == "e3f4d0c8ad9de804b14173da3e8ba81687820c3f844ec14e9b58dbc9daef4dc0",
        "reachable-count histogram SHA mismatch",
    )
    require(
        maximum_mask_count_digest
        == "bd0a1ee3a6265722536ff61691840f01729d4586876373fe815485ed800e6858",
        "maximum-mask-count histogram SHA mismatch",
    )
    require(max(feasible_owner_hist) == 4, "all-six owner-local context survived")
    require(all_order_rows == 2, "all-order scalar-row census mismatch")
    require(tournament_ties == Counter({2: 72, 4: 24, 6: 36, 7: 216, 15: 2}), "tournament tie mismatch")
    require(tournament_flips == Counter({0: 19, 1: 20, 2: 27, 3: 29, 4: 69, 5: 56, 6: 47, 7: 45, 8: 30, 9: 8}), "tournament flip mismatch")
    require(tournament_scores == Counter({score: 350 for score in range(6)}), "tournament aggregate scores mismatch")
    require(
        len(tournament_score_vectors) == 79
        and score_vector_digest
        == "a294eded1ff8ae2cdf009b283757a7253b58b4969db7f64e94e2f487bdeae7e9",
        "tournament score-vector certificate mismatch",
    )
    require(tournament_triangles == Counter({0: 350}), "tournament cycle mismatch")
    require(tournament_sccs == Counter({(1, 1, 1, 1, 1, 1): 350}), "tournament SCC mismatch")
    require(tournament_paths == Counter({1: 350}), "tournament path mismatch")

    print("scale-twenty-one independent algebraic-CRT Python referee")
    print("divisor grammar 1,3,7,21; literal states 21; supports 924")
    print(
        f"hereditary order words {len(order_words)}; "
        f"state words/support {literal_states}; raw {len(supports) * literal_states}"
    )
    print("literal-fibre-size histogram " + histogram_text(fibre_histogram))
    print(f"mask-table SHA256 {mask_digest}")
    print(f"order-grammar SHA256 {order_digest}")
    print(f"literal-fibre-descriptor SHA256 {literal_descriptor_digest}")
    print(
        f"primary-compatible FNV64 mask {primary_mask_fnv:016x} "
        f"order {primary_order_fnv:016x}"
    )
    print(
        f"scalar contexts {len(scalar_bank)} on {len(scalar_supports)} supports; "
        f"multiplicity patterns {len(scalar_patterns)}"
    )
    print("scalar multiplicities " + histogram_text(scalar_patterns))
    print(f"scalar-multiplicity SHA256 {pattern_digest}")
    print("all-support contexts histogram " + histogram_text(all_support_context_hist))
    print("multiplication orbit-size histogram " + histogram_text(orbit_size_hist))
    print("multiplication orbit (size,contexts/support) histogram " + histogram_text(orbit_context_hist))
    print(f"scalar-bank SHA256 {scalar_digest_hex}")
    print(f"capacity-row SHA256 {capacity_digest_hex}")
    print(
        f"capacity vectors {len(capacity_vectors)}; "
        f"capacity-vector SHA256 {capacity_vector_digest}"
    )
    print("minimum scalar-slack histogram " + histogram_text(minimum_slack_histogram))
    print("maximum scalar-slack histogram " + histogram_text(maximum_slack_histogram))
    print("tight-owner/context histogram " + histogram_text(tight_owner_histogram))
    print(f"owner rows {len(scalar_bank) * 6}; feasible {feasible_rows}")
    print("feasible-owner/context histogram " + histogram_text(feasible_owner_hist))
    print("maximum-union histogram " + histogram_text(maximum_union_hist))
    print("minimum-of-owner-maxima histogram " + histogram_text(minimum_owner_union_hist))
    print(
        f"distinct owner maximum vectors {len(owner_vectors)}; "
        f"owner-reachable-bank SHA256 {owner_digest_hex}"
    )
    print(
        f"reachable masks total {reachable_mask_total}; "
        f"count-histogram SHA256 {reachable_count_digest}; "
        f"maximum-mask-count-histogram SHA256 {maximum_mask_count_digest}"
    )
    print(
        f"all-order scalar rows {all_order_rows}; each has capacities "
        "21,21,21,21,21,21 and owner maxima 20,20,20,20,20,20"
    )
    print(
        "exact reachable masks are hashed in sorted order for every owner row; "
        f"owner-local all-six contexts {feasible_owner_hist.get(6, 0)}; "
        "hence global literal unit fibres 0"
    )
    print(
        "tournament vertices owner obligations; pair observable exact "
        "(feasible,maximum,capacity) summaries; lexicographic switch, "
        "coordinate tie Hamiltonian path"
    )
    print("tournament tie-edge histogram " + histogram_text(tournament_ties))
    print("tournament edge-flip histogram " + histogram_text(tournament_flips))
    print("tournament aggregate score histogram " + histogram_text(tournament_scores))
    print(
        f"tournament score vectors {len(tournament_score_vectors)}; "
        f"score-vector SHA256 {score_vector_digest}"
    )
    print("tournament directed-cycle histogram " + histogram_text(tournament_triangles))
    print("tournament SCC-size histogram " + histogram_text(tournament_sccs))
    print("tournament Hamiltonian-path histogram " + histogram_text(tournament_paths))
    print(
        "preserved audit: labelled owner-local mask banks preserve the FULL-mask "
        "predicate, sheet identities, maxima, and every owner deficit"
    )
    print(
        "lost audit: the tournament keeps only pairwise ordering of three-number "
        "summaries; it discards thresholds, masks, unit witnesses, and simultaneous "
        "cross-owner compatibility, so it cannot itself prove emptiness"
    )
    print(
        "challenged vertex assumption: owner proof obligations are faithful after "
        "the scalar gate; runners/providers, gaps, fixed sections and boundaries, "
        "wall events, residues, cover arcs, Fourier modes, and matroid circuits all "
        "lose labelled shared-unit incidence"
    )


if __name__ == "__main__":
    main()
