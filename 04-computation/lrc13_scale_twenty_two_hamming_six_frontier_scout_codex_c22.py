#!/usr/bin/env python3
"""Exact scale-22 AP-centred common-scale Hamming-six frontier certificate.

This program constructs the effective-order/unit alphabet and every local CRT
mask from first principles.  It enumerates the complete hereditary divisor
grammar on every labelled six-subset of F_13^*, applies the exact scalar owner
capacity gate, and then performs immutable-set union reachability separately
at all six owner obligations.  All serialized banks have an explicit order;
the resulting SHA-256 certificates are independent of Python hash iteration
order and of ``python`` versus ``python -O`` execution.

The owner-local test is a necessary projection of the global literal-unit
problem: a global word must cover every owner, whereas an empty owner bank
already rules it out.  No converse or cross-owner gluing claim is made here.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 22
LABELS = tuple(range(1, P))
DIVISORS = (1, 2, 11, 22)
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
    """Optimization-stable replacement for assertions."""
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    """The unique representative in (-modulus/2, modulus/2]."""
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label, state_index):
    """Solve x = D*label (mod 13), x = unit (mod D) algebraically."""
    order, unit = STATES[state_index]
    coefficient = (label - unit * pow(order, -1, P)) % P
    answer = unit + order * coefficient
    require(answer % P == order * label % P, "mod-thirteen CRT failure")
    require(answer % order == unit % order, "effective-order CRT failure")
    require(0 <= answer < P * order, "CRT representative outside range")
    # Brute uniqueness is a deliberately differently shaped audit of the
    # algebraic formula.  It is cheap because the largest modulus is 286.
    brute = tuple(
        value
        for value in range(P * order)
        if value % P == order * label % P and value % order == unit % order
    )
    require(brute == (answer,), "algebraic/brute CRT uniqueness mismatch")
    return answer


def build_masks():
    """Build label/state/owner masks and their canonical byte certificate."""
    masks = [[[0] * P for _state in STATES] for _label in range(P)]
    digest = sha256()
    base_digest = sha256()
    for label in LABELS:
        for state_index, (order, _unit) in enumerate(STATES):
            base = algebraic_crt_base(label, state_index)
            base_digest.update(bytes((label, state_index)))
            base_digest.update(base.to_bytes(2, "little"))
            for owner in LABELS:
                owner_inverse = pow(owner, -1, P)
                mask = 0
                for sheet in range(C):
                    value = centered(
                        base * (owner_inverse + P * sheet), P * order
                    )
                    if -order < value <= order:
                        mask |= 1 << sheet
                require(mask & ~FULL == 0, "mask outside twenty-two sheets")
                masks[label][state_index][owner] = mask
                digest.update(mask.to_bytes(MASK_BYTES, "little"))
    return masks, digest.hexdigest(), base_digest.hexdigest()


def analytic_cardinality(label, order, owner):
    """Count one order-period arithmetically and repeat it C/order times."""
    ratio = label * pow(owner, -1, P) % P
    target = order * ratio % P
    period_count = sum(
        value % P == target for value in range(-order + 1, order + 1)
    )
    return (C // order) * period_count


def hereditary(word):
    """Every leave-one-out lcm remains the full common scale."""
    return all(
        lcm(*(DIVISORS[word[index]] for index in range(6) if index != omitted))
        == C
        for omitted in range(6)
    )


def hereditary_prime_provider_audit(word):
    """For C=2*11, each prime must occur at least twice."""
    two_providers = sum(DIVISORS[index] % 2 == 0 for index in word)
    eleven_providers = sum(DIVISORS[index] % 11 == 0 for index in word)
    return two_providers >= 2 and eleven_providers >= 2


def grammar_census():
    """Enumerate the divisor bank and count every literal unit fibre exactly."""
    words = []
    order_digest = sha256()
    weighted_digest = sha256()
    weighted_states = 0
    for word in product(range(len(DIVISORS)), repeat=6):
        by_lcm = hereditary(word)
        by_primes = hereditary_prime_provider_audit(word)
        require(by_lcm == by_primes, "lcm/prime-provider grammar mismatch")
        if not by_lcm:
            continue
        words.append(word)
        fibre = prod(len(STATE_INDICES[index]) for index in word)
        weighted_states += fibre
        order_digest.update(bytes(word))
        weighted_digest.update(bytes(word))
        weighted_digest.update(fibre.to_bytes(8, "little"))
    return (
        tuple(words),
        weighted_states,
        order_digest.hexdigest(),
        weighted_digest.hexdigest(),
    )


def owner_reachable(masks, labels, word, owner, provider_order):
    """Exact projected union bank for one owner in a declared provider order."""
    reachable = frozenset((0,))
    layer_sizes = []
    for provider_index in provider_order:
        label = labels[provider_index]
        options = frozenset(
            masks[label][state_index][owner]
            for state_index in STATE_INDICES[word[provider_index]]
        )
        require(options, "empty state-option fibre")
        reachable = frozenset(
            partial | option
            for partial in reachable
            for option in options
        )
        require(all(mask & ~FULL == 0 for mask in reachable), "bad union mask")
        layer_sizes.append(len(reachable))
    return reachable, tuple(layer_sizes)


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
    """Orient exact owner summaries, breaking equality by coordinate order."""
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
    require(len(STATES) == C, "literal state alphabet mismatch")
    require(tuple(order for order, _unit in STATES) == tuple(sorted(
        order for order, _unit in STATES
    )), "state alphabet is not canonically ordered")
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

    masks, mask_digest, crt_base_digest = build_masks()
    cardinality = [[[0] * P for _order in DIVISORS] for _label in range(P)]
    for label in LABELS:
        for order_index, divisor in enumerate(DIVISORS):
            for owner in LABELS:
                sizes = {
                    masks[label][state][owner].bit_count()
                    for state in STATE_INDICES[order_index]
                }
                require(len(sizes) == 1, "mask cardinality depends on unit")
                size = next(iter(sizes))
                require(
                    size == analytic_cardinality(label, divisor, owner),
                    "literal/analytic mask cardinality mismatch",
                )
                cardinality[label][order_index][owner] = size
    require(
        max(
            cardinality[label][order_index][owner]
            for label in LABELS
            for order_index in range(len(DIVISORS))
            for owner in LABELS
        )
        * 6
        < 256,
        "eight-bit scalar packing is unsafe",
    )

    (
        order_words,
        literal_states,
        order_digest,
        weighted_grammar_digest,
    ) = grammar_census()

    supports = tuple(combinations(LABELS, 6))
    require(len(supports) == 924, "support census mismatch")
    scalar_bank = []
    scalar_supports = set()
    contexts_per_support = Counter()
    scalar_patterns = Counter()
    scalar_digest = sha256()
    capacity_digest = sha256()

    for labels in supports:
        # Eight-bit fields cannot carry: six providers contribute at most 132.
        packed = [[0] * len(DIVISORS) for _provider in range(6)]
        for provider_index, label in enumerate(labels):
            for order_index in range(len(DIVISORS)):
                packed[provider_index][order_index] = sum(
                    cardinality[label][order_index][owner]
                    << (8 * owner_index)
                    for owner_index, owner in enumerate(labels)
                )

        for word in order_words:
            packed_capacity = sum(
                packed[provider_index][word[provider_index]]
                for provider_index in range(6)
            )
            capacities = tuple(
                (packed_capacity >> (8 * owner_index)) & 255
                for owner_index in range(6)
            )
            require(max(capacities) <= 6 * C, "packed scalar field overflow")
            if min(capacities) < C:
                continue
            direct = tuple(
                sum(
                    cardinality[labels[provider]][word[provider]][owner]
                    for provider in range(6)
                )
                for owner in labels
            )
            require(capacities == direct, "packed/direct scalar mismatch")
            scalar_bank.append((labels, word, capacities))
            scalar_supports.add(labels)
            contexts_per_support[labels] += 1
            pattern = tuple(word.count(index) for index in range(len(DIVISORS)))
            scalar_patterns[pattern] += 1
            scalar_digest.update(bytes(labels + word))
            capacity_digest.update(bytes(labels + word + capacities))

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
        require(
            len(multiplicities) == 1,
            "context multiplicity changes in multiplication orbit",
        )
        orbit_size_hist[len(orbit)] += 1
        orbit_context_hist[(len(orbit), multiplicities.pop())] += 1
        remaining -= orbit

    feasible_rows = 0
    feasible_owner_hist = Counter()
    feasible_owner_mask_hist = Counter()
    maximum_union_hist = Counter()
    minimum_owner_union_hist = Counter()
    reachable_count_hist = Counter()
    maximum_mask_count_hist = Counter()
    owner_vectors = Counter()
    owner_digest = sha256()
    layer_digest = sha256()
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
        feasible_mask = 0
        for owner_index, owner in enumerate(labels):
            reachable, forward_layers = owner_reachable(
                masks, labels, word, owner, range(6)
            )
            reverse_reachable, reverse_layers = owner_reachable(
                masks, labels, word, owner, range(5, -1, -1)
            )
            require(
                reachable == reverse_reachable,
                "forward/reverse provider reachability mismatch",
            )
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
            summary = (
                int(feasible),
                maximum,
                capacities[owner_index],
                len(reachable),
                maximum_mask_count,
            )
            summaries.append(summary)

            reachable_digest = sha256()
            for mask in sorted(reachable):
                reachable_digest.update(mask.to_bytes(MASK_BYTES, "little"))
            owner_digest.update(bytes(labels + word + (owner_index, owner)))
            owner_digest.update(bytes((int(feasible), maximum)))
            owner_digest.update(capacities[owner_index].to_bytes(2, "little"))
            owner_digest.update(len(reachable).to_bytes(4, "little"))
            owner_digest.update(maximum_mask_count.to_bytes(4, "little"))
            owner_digest.update(reachable_digest.digest())
            layer_digest.update(bytes(labels + word + (owner_index,)))
            for layer_size in forward_layers:
                layer_digest.update(layer_size.to_bytes(4, "little"))
            for layer_size in reverse_layers:
                layer_digest.update(layer_size.to_bytes(4, "little"))

        feasible_count = feasible_mask.bit_count()
        feasible_owner_hist[feasible_count] += 1
        feasible_owner_mask_hist[feasible_mask] += 1
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
    capacity_digest_hex = capacity_digest.hexdigest()
    owner_digest_hex = owner_digest.hexdigest()
    layer_digest_hex = layer_digest.hexdigest()
    pattern_digest = tuple_histogram_digest(scalar_patterns)
    score_vector_digest = tuple_histogram_digest(tournament_score_vectors)

    # Frozen primary certificate.  These constants come from the initial
    # clean first-principles replay and make every subsequent run fail closed.
    require(len(order_words) == 3_249, "hereditary order-word census mismatch")
    require(literal_states == 100_975_500, "literal state census mismatch")
    require(len(supports) * literal_states == 93_301_362_000, "raw census mismatch")
    require(
        crt_base_digest
        == "fe217f797c08702c8f607d3a936321fb7ffb6c6e73770a0097cf71c597297793",
        "CRT-base digest mismatch",
    )
    require(
        mask_digest
        == "54587d940a12b70601943dbe7505d4797363395a2579ea6f3e09583db5a01282",
        "mask-table digest mismatch",
    )
    require(
        order_digest
        == "f7c0254d8ac9108d318f4a9a21d0d2e5b244be91087b22c162bb563956e9b474",
        "order-grammar digest mismatch",
    )
    require(
        weighted_grammar_digest
        == "b4ec74d190864c2a050409126bacfd79fb0fac97ce2952628b17341b1718c4dd",
        "weighted literal-grammar digest mismatch",
    )
    require(
        len(scalar_bank) == 984
        and len(scalar_supports) == 180
        and len(scalar_patterns) == 8,
        "scalar bank census mismatch",
    )
    require(
        scalar_patterns
        == Counter(
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
        ),
        "scalar multiplicity histogram mismatch",
    )
    require(
        pattern_digest
        == "5e1220b62de69e3d493ae5ae6731ffde1dce8ed7728e442abaa42629ef36a80d",
        "scalar multiplicity digest mismatch",
    )
    require(
        all_support_context_hist == Counter({0: 744, 2: 96, 3: 24, 6: 24, 16: 36}),
        "all-support scalar-context histogram mismatch",
    )
    require(
        orbit_size_hist == Counter({12: 15}),
        "support multiplication-orbit census mismatch",
    )
    require(
        orbit_context_hist
        == Counter({(12, 2): 8, (12, 3): 2, (12, 6): 2, (12, 16): 3}),
        "support-orbit context histogram mismatch",
    )
    require(
        scalar_digest_hex
        == "29067f69b228b9956239b27a43af9bc72e8c141acfb47587f536dc557cebb1de",
        "scalar-bank digest mismatch",
    )
    require(
        capacity_digest_hex
        == "5f10732eda4cd0dcf9fe2eb0166e4191673774d40fae03ec4993d143caa3528f",
        "capacity-bank digest mismatch",
    )
    require(feasible_rows == 192, "feasible owner-row census mismatch")
    require(
        feasible_owner_hist == Counter({0: 792, 1: 192}),
        "feasible-owner histogram mismatch",
    )
    require(
        feasible_owner_mask_hist
        == Counter({0: 792, 1: 38, 2: 34, 4: 24, 8: 24, 16: 34, 32: 38}),
        "feasible-owner-mask histogram mismatch",
    )
    require(
        maximum_union_hist
        == Counter({16: 864, 17: 1_584, 18: 2_784, 19: 480, 22: 192}),
        "maximum-union histogram mismatch",
    )
    require(
        minimum_owner_union_hist == Counter({16: 408, 17: 180, 18: 396}),
        "minimum owner-maximum histogram mismatch",
    )
    require(
        reachable_count_hist
        == Counter(
            {
                1: 192,
                255: 24,
                265: 624,
                315: 96,
                317: 408,
                332: 24,
                360: 24,
                395: 72,
                410: 288,
                415: 192,
                432: 192,
                445: 48,
                470: 288,
                570: 48,
                575: 48,
                627: 24,
                660: 240,
                665: 48,
                675: 888,
                677: 1_656,
                705: 24,
                715: 72,
                720: 96,
                730: 48,
                740: 48,
                750: 96,
                760: 48,
                770: 48,
            }
        ),
        "reachable-count histogram mismatch",
    )
    require(
        maximum_mask_count_hist
        == Counter(
            {
                1: 192,
                20: 240,
                25: 240,
                85: 72,
                90: 288,
                95: 48,
                100: 96,
                120: 144,
                150: 96,
                152: 408,
                185: 624,
                200: 1_128,
                205: 240,
                210: 1_848,
                220: 48,
                252: 192,
            }
        ),
        "maximum-mask-count histogram mismatch",
    )
    require(len(owner_vectors) == 127, "owner maximum-vector census mismatch")
    require(
        owner_digest_hex
        == "b881af1af73fe6dde434d92ca0598bac79828881ad32feff87d711fa448a0eef",
        "owner reachable-bank digest mismatch",
    )
    require(
        layer_digest_hex
        == "b08be64149acffb006205465bbcb5825b06d4e2ffb1ece561506f0e1746b8baa",
        "forward/reverse layer-bank digest mismatch",
    )
    require(max(feasible_owner_hist) == 1, "two owner-local projections survived")
    require(
        tournament_ties == Counter({0: 288, 1: 240, 2: 216, 3: 192, 4: 24, 7: 24}),
        "tournament tie-edge histogram mismatch",
    )
    require(
        tournament_flips
        == Counter(
            {
                0: 5,
                1: 3,
                2: 30,
                3: 60,
                4: 95,
                5: 122,
                6: 142,
                7: 157,
                8: 135,
                9: 100,
                10: 66,
                11: 38,
                12: 22,
                13: 8,
                14: 1,
            }
        ),
        "tournament edge-flip histogram mismatch",
    )
    require(
        tournament_scores == Counter({score: 984 for score in range(6)}),
        "tournament aggregate score histogram mismatch",
    )
    require(
        len(tournament_score_vectors) == 411
        and score_vector_digest
        == "fdf170efdaacdab8cc6fb54ae2ca4dad4c03bfa3213f60d52c798905fb45aa0b",
        "tournament score-vector certificate mismatch",
    )
    require(
        tournament_triangles == Counter({0: 984}),
        "tournament directed-cycle census mismatch",
    )
    require(
        tournament_sccs == Counter({(1, 1, 1, 1, 1, 1): 984}),
        "tournament SCC fingerprint mismatch",
    )
    require(
        tournament_paths == Counter({1: 984}),
        "tournament Hamiltonian-path census mismatch",
    )

    print("scale-twenty-two algebraic-CRT exact frontier certificate")
    print("divisor grammar 1,2,11,22; literal states 22; supports 924")
    print(
        f"hereditary order words {len(order_words)}; "
        f"state words/support {literal_states}; raw {len(supports) * literal_states}"
    )
    print(f"CRT-base SHA256 {crt_base_digest}")
    print(f"mask-table SHA256 {mask_digest}")
    print(f"order-grammar SHA256 {order_digest}")
    print(f"weighted-literal-grammar SHA256 {weighted_grammar_digest}")
    print(
        "execution-mode audit fail-closed require checks; frozen normal and "
        "python -O replays byte-identical"
    )
    print(
        f"scalar contexts {len(scalar_bank)} on {len(scalar_supports)} supports; "
        f"multiplicity patterns {len(scalar_patterns)}"
    )
    print(f"scalar-multiplicity SHA256 {pattern_digest}")
    print("scalar multiplicities n1,n2,n11,n22 " + histogram_text(scalar_patterns))
    print("all-support contexts histogram " + histogram_text(all_support_context_hist))
    print("multiplication orbit-size histogram " + histogram_text(orbit_size_hist))
    print(
        "multiplication orbit (size,contexts/support) histogram "
        + histogram_text(orbit_context_hist)
    )
    print(f"scalar-bank SHA256 {scalar_digest_hex}")
    print(f"capacity-bank SHA256 {capacity_digest_hex}")
    print(f"owner rows {len(scalar_bank) * 6}; feasible {feasible_rows}")
    print("feasible-owner/context histogram " + histogram_text(feasible_owner_hist))
    print("feasible-owner-mask histogram " + histogram_text(feasible_owner_mask_hist))
    print("maximum-union histogram " + histogram_text(maximum_union_hist))
    print(
        "minimum-of-owner-maxima histogram "
        + histogram_text(minimum_owner_union_hist)
    )
    print("reachable-count histogram " + histogram_text(reachable_count_hist))
    print("maximum-mask-count histogram " + histogram_text(maximum_mask_count_hist))
    print(
        f"distinct owner maximum vectors {len(owner_vectors)}; "
        f"owner-reachable-bank SHA256 {owner_digest_hex}"
    )
    print(f"forward/reverse layer-bank SHA256 {layer_digest_hex}")
    print(
        "exact reachable masks are hashed in sorted order for every owner row; "
        f"owner-local all-six contexts {feasible_owner_hist.get(6, 0)}; "
        "hence global literal unit fibres 0"
    )
    print(
        "tournament vertices owner obligations; pair observable exact ordered "
        "(feasible,maximum,capacity,reachable-count,maximum-mask-count); "
        "lexicographic switch with the coordinate tie Hamiltonian path"
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
        "preserved audit: exact owner banks preserve the FULL-mask predicate, "
        "maximum deficit, sheet masks, and local unit choices; the tournament "
        "summary retains feasibility, deficit, scalar capacity, and bank-size "
        "ordering"
    )
    print(
        "lost audit: projection forgets simultaneous cross-owner unit gluing; "
        "the tournament further discards exact masks, sheet identities, absolute "
        "threshold meaning, and witness incidence, so it is diagnostic only"
    )
    print(
        "challenged vertex assumption: owners/proof obligations were selected; "
        "providers, gaps, fixed circle sections, section boundaries, wall events, "
        "residues, cover arcs, Fourier modes, and matroid circuits do not retain "
        "the labelled owner-local terminal predicate without extra incidence data"
    )


if __name__ == "__main__":
    main()
