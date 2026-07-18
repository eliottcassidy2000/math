#!/usr/bin/env python3
"""Independent exact referee for the common-scale-eighteen Hamming-six face.

The primary C++ certificate finds CRT bases by literal search and hashes with
FNV-1a.  This referee solves the CRT algebraically, enumerates the hereditary
grammar with Python products, hashes with SHA-256, and computes each owner's
reachable sheet unions with immutable Python sets.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

P = 13
C = 18
FULL = (1 << C) - 1
DIVISORS = (1, 2, 3, 6, 9, 18)
STATES = (
    (1, 0),
    (2, 1),
    (3, 1), (3, 2),
    (6, 1), (6, 5),
    (9, 1), (9, 2), (9, 4), (9, 5), (9, 7), (9, 8),
    (18, 1), (18, 5), (18, 7), (18, 11), (18, 13), (18, 17),
)
STATE_INDICES = (
    (0,),
    (1,),
    (2, 3),
    (4, 5),
    (6, 7, 8, 9, 10, 11),
    (12, 13, 14, 15, 16, 17),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label, state):
    """Solve x = D*label (mod 13), x = e (mod D) algebraically."""
    order, unit = STATES[state]
    coefficient = (label - unit * pow(order, -1, P)) % P
    base = unit + order * coefficient
    require(base % P == order * label % P, "mod-thirteen CRT failure")
    require(base % order == unit % order, "effective-order CRT failure")
    require(0 <= base < P * order, "CRT base outside canonical interval")
    return base


def build_masks():
    masks = [[[0 for _owner in range(P)] for _state in STATES]
             for _label in range(P)]
    digest = sha256()
    for label in range(1, P):
        for state, (order, _unit) in enumerate(STATES):
            base = algebraic_crt_base(label, state)
            for owner in range(1, P):
                inverse = pow(owner, -1, P)
                mask = 0
                for sheet in range(C):
                    value = centered(
                        base * (inverse + P * sheet), P * order)
                    if -order < value <= order:
                        mask |= 1 << sheet
                masks[label][state][owner] = mask
                digest.update(mask.to_bytes(3, "little"))
    return masks, digest.hexdigest()


def hereditary(word):
    return all(
        lcm(*(DIVISORS[word[i]] for i in range(6) if i != omitted)) == C
        for omitted in range(6)
    )


def grammar_census():
    words = tuple(word for word in product(range(6), repeat=6)
                  if hereditary(word))
    order_digest = sha256()
    for word in words:
        order_digest.update(bytes(word))

    weighted = sum(prod(len(STATE_INDICES[index]) for index in word)
                   for word in words)
    literal = 0
    literal_digest = sha256()
    buffer = bytearray()
    for word in words:
        fibres = tuple(STATE_INDICES[index] for index in word)
        for states in product(*fibres):
            literal += 1
            buffer.extend(states)
            if len(buffer) >= 1 << 20:
                literal_digest.update(buffer)
                buffer.clear()
    literal_digest.update(buffer)
    return (words, weighted, literal, order_digest.hexdigest(),
            literal_digest.hexdigest())


def owner_local(masks, labels, word, owner):
    reachable = frozenset((0,))
    for provider in range(6):
        options = frozenset(
            masks[labels[provider]][state][owner]
            for state in STATE_INDICES[word[provider]]
        )
        reachable = frozenset(
            partial | option
            for partial in reachable for option in options
        )
    maximum = max(mask.bit_count() for mask in reachable)
    feasible = FULL in reachable
    require(feasible == (maximum == C),
            "owner feasibility and maximum union disagree")
    return feasible, maximum, len(reachable)


def tournament_fingerprint(locals_):
    """Orient owner summaries; exact ties follow the coordinate path."""
    out = [0] * 6
    scores = [0] * 6
    ties = 0
    flips = 0
    for left, right in combinations(range(6), 2):
        left_key = locals_[left][:2]
        right_key = locals_[right][:2]
        winner = left
        if left_key == right_key:
            ties += 1
        elif right_key > left_key:
            winner = right
            flips += 1
        loser = left + right - winner
        out[winner] |= 1 << loser
        scores[winner] += 1

    triangles = 0
    for i, j, k in combinations(range(6), 3):
        forward = (out[i] >> j) & (out[j] >> k) & (out[k] >> i) & 1
        reverse = (out[i] >> k) & (out[k] >> j) & (out[j] >> i) & 1
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
                if ((previous_mask >> previous) & 1
                        and (out[previous] >> last) & 1):
                    paths[mask][last] += paths[previous_mask][previous]
    require(sorted(scores) == list(range(6)), "tournament score failure")
    require(triangles == 0, "owner tournament has a directed triangle")
    require(sum(paths[-1]) == 1, "owner tournament path-count failure")
    return ties, flips


def histogram_text(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def digest_counter(counter):
    digest = sha256()
    for key, count in sorted(counter.items()):
        values = key if isinstance(key, tuple) else (key,)
        digest.update(bytes(values))
        digest.update(count.to_bytes(8, "little"))
    return digest.hexdigest()


def main():
    for order_index, divisor in enumerate(DIVISORS):
        actual = {STATES[state][1] for state in STATE_INDICES[order_index]}
        expected = ({0} if divisor == 1 else
                    {unit for unit in range(1, divisor)
                     if gcd(unit, divisor) == 1})
        require(actual == expected, "unit grammar mismatch")
        require(all(STATES[state][0] == divisor
                    for state in STATE_INDICES[order_index]),
                "state/order indexing mismatch")

    masks, mask_digest = build_masks()
    require(mask_digest ==
            "67efd99def67e9eec2ff5182ca9ef12a33d12751e725b1a1bf039107d25db8cd",
            "algebraic CRT mask-table digest mismatch")
    cardinality = [[[0] * P for _order in DIVISORS]
                   for _label in range(P)]
    for label in range(1, P):
        for order_index in range(6):
            for owner in range(1, P):
                cardinalities = {
                    masks[label][state][owner].bit_count()
                    for state in STATE_INDICES[order_index]
                }
                require(len(cardinalities) == 1,
                        "unit-dependent scalar cardinality")
                cardinality[label][order_index][owner] = cardinalities.pop()

    (order_words, weighted_states, literal_states, order_digest,
     state_digest) = grammar_census()
    require(len(order_words) == 26_961, "hereditary word census mismatch")
    require(weighted_states == 29_751_948,
            "weighted state census mismatch")
    require(literal_states == weighted_states,
            "literal state census mismatch")
    require(order_digest ==
            "447038aeee836c5f42c98ae8ae90ce0f4d86be83969dc70402b7691c6e1c7f61",
            "order-grammar digest mismatch")
    require(state_digest ==
            "8266e17dfec2534891ec4336ed4601f05310fcf7e2cd16bcdf8017bbe0e52718",
            "literal state-grammar digest mismatch")

    scalar_contexts = 0
    feasible_rows = 0
    scalar_supports = set()
    contexts_per_support = Counter()
    scalar_patterns = Counter()
    feasible_owner_hist = Counter()
    maximum_union_hist = Counter()
    minimum_owner_hist = Counter()
    owner_vectors = Counter()
    tournament_ties = Counter()
    tournament_flips = Counter()
    scalar_digest = sha256()
    owner_digest = sha256()

    supports = tuple(combinations(range(1, P), 6))
    require(len(supports) == 924, "support census mismatch")
    for labels in supports:
        # Six seven-bit fields: every field is at most 6*18 = 108.
        packed = [[0] * 6 for _provider in range(6)]
        for provider, label in enumerate(labels):
            for order_index in range(6):
                packed[provider][order_index] = sum(
                    cardinality[label][order_index][owner] << (7 * owner_index)
                    for owner_index, owner in enumerate(labels)
                )

        for word in order_words:
            capacity = sum(packed[provider][word[provider]]
                           for provider in range(6))
            if any(((capacity >> (7 * owner)) & 127) < C
                   for owner in range(6)):
                continue

            scalar_contexts += 1
            scalar_supports.add(labels)
            contexts_per_support[labels] += 1
            scalar_patterns[tuple(word.count(index) for index in range(6))] += 1
            scalar_digest.update(bytes(labels + word))

            locals_ = tuple(owner_local(masks, labels, word, owner)
                            for owner in labels)
            maxima = tuple(entry[1] for entry in locals_)
            feasible_mask = sum(
                bool(entry[0]) << owner_index
                for owner_index, entry in enumerate(locals_)
            )
            feasible = feasible_mask.bit_count()
            feasible_rows += feasible
            feasible_owner_hist[feasible] += 1
            maximum_union_hist.update(maxima)
            minimum_owner_hist[min(maxima)] += 1
            owner_vectors[maxima] += 1
            owner_digest.update(bytes(labels + word + (feasible_mask,) + maxima))
            for entry in locals_:
                owner_digest.update(entry[2].to_bytes(4, "little"))
            ties, flips = tournament_fingerprint(locals_)
            tournament_ties[ties] += 1
            tournament_flips[flips] += 1

    contexts_hist = Counter(contexts_per_support.values())
    remaining = set(scalar_supports)
    orbit_hist = Counter()
    while remaining:
        support = min(remaining)
        orbit = {
            tuple(sorted((multiplier * label) % P for label in support))
            for multiplier in range(1, P)
        }
        require(orbit <= scalar_supports,
                "multiplication leaves scalar support bank")
        require(orbit <= remaining, "multiplication orbits overlap")
        remaining -= orbit
        orbit_hist[len(orbit)] += 1

    require(924 * len(order_words) == 24_911_964,
            "labelled order-context census mismatch")
    require(924 * weighted_states == 27_490_799_952,
            "raw labelled state-context census mismatch")
    require(scalar_contexts == 13_098 and len(scalar_supports) == 684
            and len(scalar_patterns) == 63,
            "scalar classification mismatch")
    require(feasible_rows == 5_268, "feasible owner-row census mismatch")
    require(feasible_owner_hist == Counter({
        0: 8_922, 1: 3_108, 2: 1_056, 4: 12,
    }), "feasible-owner histogram mismatch")
    require(maximum_union_hist == Counter({
        10: 48, 11: 96, 12: 2_472, 13: 6_732, 14: 21_696,
        15: 19_584, 16: 16_068, 17: 6_624, 18: 5_268,
    }), "maximum-union histogram mismatch")
    require(minimum_owner_hist == Counter({
        10: 48, 11: 72, 12: 1_800, 13: 4_104, 14: 4_734,
        15: 2_088, 16: 252,
    }), "minimum-owner histogram mismatch")
    require(contexts_hist == Counter({
        1: 96, 2: 24, 3: 24, 5: 96, 6: 48, 10: 120,
        12: 96, 18: 24, 19: 24, 23: 24, 28: 24, 34: 6,
        65: 24, 85: 12, 91: 6, 102: 24, 156: 12,
    }), "contexts/support histogram mismatch")
    require(orbit_hist == Counter({6: 2, 12: 56}),
            "multiplication orbit census mismatch")
    require(len(owner_vectors) == 4_575, "owner-vector census mismatch")
    require(max(feasible_owner_hist) == 4,
            "feasible-owner maximum mismatch")
    require(tournament_ties == Counter({
        1: 1_080, 2: 2_748, 3: 2_634, 4: 3_048,
        6: 1_860, 7: 1_320, 10: 360, 15: 48,
    }), "tournament tie census mismatch")
    require(tournament_flips == Counter({
        0: 314, 1: 363, 2: 772, 3: 1_206, 4: 1_911,
        5: 1_951, 6: 2_016, 7: 1_795, 8: 1_237, 9: 767,
        10: 415, 11: 234, 12: 93, 13: 20, 14: 4,
    }), "tournament flip census mismatch")

    scalar_digest_hex = scalar_digest.hexdigest()
    owner_digest_hex = owner_digest.hexdigest()
    multiplicity_digest = digest_counter(scalar_patterns)
    require(multiplicity_digest ==
            "e5d1c2ce7fbb048cd424458d3cb0bb55e24b5152d8268d46dd9c01eeddc8fa7e",
            "scalar multiplicity digest mismatch")
    require(scalar_digest_hex ==
            "f0ad58ab958e47ca6cd4d81eac525c1f7540b430686e22c59ac005832863d97f",
            "scalar-bank digest mismatch")
    require(owner_digest_hex ==
            "4edd78570e677687594feb93298600bef46402c4f1ac4a3fe7a286bcfa2f5c52",
            "owner-profile digest mismatch")

    print("scale-eighteen independent algebraic-CRT Python referee")
    print("divisor grammar 1,2,3,6,9,18; supports 924")
    print(f"hereditary order words {len(order_words)}; "
          f"state words/support {literal_states}; raw {924 * literal_states}")
    print(f"mask-table SHA256 {mask_digest}")
    print(f"order-grammar SHA256 {order_digest}")
    print(f"literal-state-grammar SHA256 {state_digest}")
    print(f"scalar contexts {scalar_contexts} on {len(scalar_supports)} supports; "
          f"multiplicity patterns {len(scalar_patterns)}")
    print(f"scalar-multiplicity SHA256 {multiplicity_digest}")
    print("contexts-per-support histogram " + histogram_text(contexts_hist))
    print("multiplication orbit-size histogram " + histogram_text(orbit_hist))
    print(f"scalar-bank SHA256 {scalar_digest_hex}")
    print(f"owner rows {scalar_contexts * 6}; feasible {feasible_rows}")
    print("feasible-owner/context histogram "
          + histogram_text(feasible_owner_hist))
    print("maximum-union histogram " + histogram_text(maximum_union_hist))
    print("minimum-owner histogram " + histogram_text(minimum_owner_hist))
    print(f"distinct owner vectors {len(owner_vectors)}; "
          f"owner-profile SHA256 {owner_digest_hex}")
    print("owner-local all-six contexts 0; global literal unit fibres 0")
    print("tournament pair observable is the exact ordered owner-summary pair; "
          "lexicographic switch and coordinate tie path yield transitive telemetry")
    print("tournament tie-edge histogram " + histogram_text(tournament_ties))
    print("tournament edge-flip histogram " + histogram_text(tournament_flips))
    print("challenged vertices owner obligations preserve the terminal deficit; "
          "providers, residues, divisor words, and isolated owner-sheet sets lose "
          "the exact six-owner feasibility/max-union vector")


if __name__ == "__main__":
    main()
