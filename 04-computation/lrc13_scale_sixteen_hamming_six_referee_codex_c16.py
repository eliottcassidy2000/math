#!/usr/bin/env python3
"""Independent exact referee for the common-scale-sixteen Hamming-six face."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod

P = 13
C = 16
FULL = (1 << C) - 1
DIVISORS = (1, 2, 4, 8, 16)
STATES = (
    (1, 0),
    (2, 1),
    (4, 1), (4, 3),
    (8, 1), (8, 3), (8, 5), (8, 7),
    (16, 1), (16, 3), (16, 5), (16, 7),
    (16, 9), (16, 11), (16, 13), (16, 15),
)
STATE_INDICES = (
    (0,),
    (1,),
    (2, 3),
    (4, 5, 6, 7),
    (8, 9, 10, 11, 12, 13, 14, 15),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value, modulus):
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label, state):
    """Solve x=D*label (mod 13), x=e (mod D) without CRT search."""
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
                digest.update(mask.to_bytes(2, "little"))
    return masks, digest.hexdigest()


def hereditary(word):
    for omitted in range(6):
        residual = 1
        for coordinate, order_index in enumerate(word):
            if coordinate != omitted:
                residual = lcm(residual, DIVISORS[order_index])
        if residual != C:
            return False
    return True


def grammar_census():
    words = [word for word in product(range(5), repeat=6)
             if hereditary(word)]
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
    reachable = {0}
    for provider in range(6):
        options = {
            masks[labels[provider]][state][owner]
            for state in STATE_INDICES[word[provider]]
        }
        reachable = {partial | option
                     for partial in reachable for option in options}
    maximum = max(mask.bit_count() for mask in reachable)
    feasible = FULL in reachable
    require(feasible == (maximum == C),
            "owner feasibility and maximum union disagree")
    return feasible, maximum, len(reachable)


def tournament_fingerprint(locals_):
    """Orient exact owner-summary pairs; ties follow 0->1->...->5."""
    out = [0] * 6
    scores = [0] * 6
    ties = 0
    flips = 0
    for left in range(6):
        for right in range(left + 1, 6):
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
    hamiltonian_paths = sum(paths[-1])
    require(sorted(scores) == list(range(6)), "tournament score failure")
    require(triangles == 0, "owner tournament has a directed triangle")
    require(hamiltonian_paths == 1, "owner tournament path count failure")
    return ties, flips


def histogram_text(counter):
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def multiplicity_text(counter):
    return " ".join(
        f"{','.join(map(str, key))}:{counter[key]}" for key in sorted(counter)
    )


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
            "328249e6cdd6f5fbb53d746ea3c476c4f1a7c38023cfcb1533769a8b40ab2b1e",
            "algebraic CRT mask-table digest mismatch")
    for label in range(1, P):
        for owner in range(1, P):
            for order_index in range(5):
                cardinalities = {
                    masks[label][state][owner].bit_count()
                    for state in STATE_INDICES[order_index]
                }
                require(len(cardinalities) == 1,
                        "unit-dependent scalar cardinality")
    require(max(masks[label][state][owner].bit_count()
                for label in range(1, P)
                for state in range(len(STATES))
                for owner in range(1, P)) <= C,
            "mask cardinality exceeds sheet count")

    (order_words, weighted_states, literal_states, order_digest,
     state_digest) = grammar_census()
    require(len(order_words) == 5_385, "hereditary word census mismatch")
    require(weighted_states == 14_942_208,
            "weighted state census mismatch")
    require(literal_states == weighted_states,
            "literal state census mismatch")
    require(order_digest ==
            "13dab484f3c13ea8b6934dc9e6c39dfaaffa681c9214092d9a70037be845433d",
            "order-grammar digest mismatch")
    require(state_digest ==
            "3170258dc4291dd9b45140b3cff808d828fe488c98e23fa2c4b09a94cbd232f7",
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

    supports = list(combinations(range(1, P), 6))
    require(len(supports) == 924, "support census mismatch")
    representative = tuple(indices[0] for indices in STATE_INDICES)
    for labels in supports:
        # Six seven-bit fields allow ordinary integer addition without carries:
        # each of six providers contributes at most sixteen sheets, total 96.
        packed = [[0] * 5 for _provider in range(6)]
        for provider, label in enumerate(labels):
            for order_index in range(5):
                value = 0
                state = representative[order_index]
                for owner_index, owner in enumerate(labels):
                    value |= (masks[label][state][owner].bit_count()
                              << (7 * owner_index))
                packed[provider][order_index] = value

        for word in order_words:
            capacity = sum(packed[provider][word[provider]]
                           for provider in range(6))
            if any(((capacity >> (7 * owner)) & 127) < C
                   for owner in range(6)):
                continue

            scalar_contexts += 1
            scalar_supports.add(labels)
            contexts_per_support[labels] += 1
            pattern = tuple(word.count(index) for index in range(5))
            scalar_patterns[pattern] += 1
            scalar_digest.update(bytes(labels + word))

            locals_ = tuple(owner_local(masks, labels, word, owner)
                            for owner in labels)
            maxima = tuple(entry[1] for entry in locals_)
            feasible = sum(entry[0] for entry in locals_)
            feasible_rows += feasible
            feasible_owner_hist[feasible] += 1
            maximum_union_hist.update(maxima)
            minimum_owner_hist[min(maxima)] += 1
            owner_vectors[maxima] += 1
            owner_digest.update(bytes(labels + word + maxima))
            ties, flips = tournament_fingerprint(locals_)
            tournament_ties[ties] += 1
            tournament_flips[flips] += 1

    expected_patterns = Counter({
        (0, 0, 0, 3, 3): 28, (0, 0, 0, 4, 2): 54,
        (0, 0, 1, 0, 5): 48, (0, 0, 1, 1, 4): 216,
        (0, 0, 1, 2, 3): 360, (0, 0, 1, 3, 2): 444,
        (0, 0, 2, 2, 2): 132, (0, 0, 3, 0, 3): 16,
        (0, 0, 3, 1, 2): 192, (0, 0, 4, 0, 2): 30,
        (0, 1, 1, 2, 2): 48, (0, 1, 2, 0, 3): 120,
        (0, 1, 2, 1, 2): 120, (0, 1, 3, 0, 2): 96,
        (0, 2, 0, 0, 4): 36, (0, 2, 0, 1, 3): 144,
        (0, 2, 0, 2, 2): 216, (0, 2, 1, 0, 3): 48,
        (0, 2, 1, 1, 2): 48, (0, 3, 1, 0, 2): 144,
    })
    require(scalar_contexts == 2_540 and len(scalar_supports) == 404,
            "scalar census mismatch")
    require(scalar_patterns == expected_patterns,
            "scalar pattern census mismatch")
    require(feasible_rows == 636, "feasible row census mismatch")
    require(feasible_owner_hist == Counter({0: 2_006, 1: 432, 2: 102}),
            "feasible-owner histogram mismatch")
    require(maximum_union_hist == Counter({
        9: 144, 10: 468, 11: 876, 12: 2_316,
        13: 4_068, 14: 4_740, 15: 1_992, 16: 636,
    }), "maximum-union histogram mismatch")
    require(minimum_owner_hist == Counter({
        9: 120, 10: 462, 11: 336, 12: 720, 13: 886, 14: 16,
    }), "minimum-owner histogram mismatch")
    require(max(feasible_owner_hist) == 2,
            "some context passes more than two owner gates")
    require(len(owner_vectors) == 1_210,
            "owner-vector census mismatch")

    contexts_hist = Counter(contexts_per_support.values())
    require(contexts_hist == Counter({
        1: 102, 2: 96, 3: 24, 4: 24, 7: 24, 8: 24,
        11: 36, 13: 24, 16: 24, 17: 24, 109: 2,
    }), "contexts/support histogram mismatch")

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
    require(orbit_hist == Counter({2: 1, 6: 1, 12: 33}),
            "multiplication orbit census mismatch")

    require(tournament_ties == Counter({
        1: 120, 2: 648, 3: 390, 4: 744, 6: 344, 7: 270, 10: 24,
    }), "tournament tie census mismatch")
    require(tournament_flips == Counter({
        0: 47, 1: 43, 2: 115, 3: 209, 4: 387, 5: 384,
        6: 510, 7: 360, 8: 256, 9: 151, 10: 48, 11: 25,
        12: 4, 14: 1,
    }), "tournament flip census mismatch")

    scalar_digest_hex = scalar_digest.hexdigest()
    owner_digest_hex = owner_digest.hexdigest()
    require(scalar_digest_hex ==
            "14552a8ca9b63b8efd4961518fa60a5be7694c3295314ab0db7616fa7b9b91e3",
            "scalar-bank digest mismatch")
    require(owner_digest_hex ==
            "89b80ed80e555d84025c6c3279ea8b23f7d58b2c88b55af61e89fa80f40a1e73",
            "owner-profile digest mismatch")

    print("scale-sixteen independent algebraic-CRT Python referee")
    print("divisor grammar 1,2,4,8,16; supports 924")
    print(f"hereditary order words {len(order_words)}; "
          f"state words/support {literal_states}; raw {924 * literal_states}")
    print(f"mask-table SHA256 {mask_digest}")
    print(f"order-grammar SHA256 {order_digest}")
    print(f"literal-state-grammar SHA256 {state_digest}")
    print(f"scalar contexts {scalar_contexts} on {len(scalar_supports)} supports")
    print("scalar multiplicities n1,n2,n4,n8,n16 "
          + multiplicity_text(scalar_patterns))
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
    print("owner-local all-six contexts 0; global unit fibres 0")
    print("tournament pair observable is the exact ordered owner-summary pair; "
          "lexicographic switch and coordinate tie path yield transitive telemetry")
    print("tournament tie-edge histogram " + histogram_text(tournament_ties))
    print("tournament edge-flip histogram " + histogram_text(tournament_flips))
    print("challenged vertices owner obligations preserve the terminal obstruction; "
          "providers, residues, divisor words, and isolated owner-sheet sets lose "
          "the exact six-owner feasibility/max-union vector")


if __name__ == "__main__":
    main()
