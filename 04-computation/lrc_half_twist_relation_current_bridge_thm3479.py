#!/usr/bin/env python3
"""Standalone exact core companion for promoted THM-3479.

The program checks four separately typed layers.

1. One explicit all-91-unit relation is realized by two q=11 owner tuples.
2. U_full has the required half-twist masks, a positive delayed word, and
   all 169 unrestricted THM-2334 target aggregates nonzero.
3. U_clock has a literal common zero centre in the same delayed owner word,
   a positive delayed mass, a two-embedding nonconstant endpoint bank, and
   one explicit THM-2331 atomic address term.
4. Exact q=27 and q=51 congruence lifts preserve their literal covers, with
   the q=51 affine k mod 3 character retained.
5. The natural C13-equivariant attempt to identify a 13-character fibre with
   the proved private-support carrier's 13 edges is obstructed exactly:
   the graph automorphism group has no element of order 13.  This core
   companion does not manufacture a non-equivariant role chart; the promoted
   package records the later explicit 72-chart sidecar separately.

The two endpoint implementations are independent.  The full-bank engine
uses periodic subtraction and a scaled sweep.  The referee engine refines
at every exact boundary, tests exact midpoints, and directly intersects all
169 word preimages.  Nothing here computes a grouped exact-address
coefficient, an all-91-unit target projector, ancestry, a bispectrum, or an
LRC(14) conclusion.
"""

from __future__ import annotations

import ast
from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
OUTPUT = "05-knowledge/results/lrc_half_twist_relation_current_bridge_thm3479.out"

HASH_PINS = (
    (
        "THM-2309",
        ROOT / "01-canon/theorems/THM-2309-owner-aligned-pivot-packets-and-visible-height-separation.md",
        "dd4db397efac078a7c049956491f56c0c1b111daeb17bf9fddc117d39a253bc1",
        "PROVED",
    ),
    (
        "THM-2331",
        ROOT / "01-canon/theorems/THM-2331-two-sided-septimal-address-embedding-in-marked-current.md",
        "1e986d42105f9186ed0809f14798cbf92dcaae1af3f50813a56687ed02dbac35",
        "PROVED",
    ),
    (
        "THM-2334",
        ROOT / "01-canon/theorems/THM-2334-relation-residue-current-and-character-twist-pushforward.md",
        "c80f7bb3d31274a02f046fa6cea3b36cd56c62be611936ca32ee723881cd3899",
        "PROVED",
    ),
    (
        "THM-2349",
        ROOT / "01-canon/theorems/THM-2349-first-depth-one-delayed-shallow-restart.md",
        "fcf6af64b1810d763473bb2d84efb5c76868c533d35bddfaad8b0a39e80ea3d6",
        "PROVED",
    ),
    (
        "THM-3398",
        ROOT / "01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md",
        "01901da2bb382184cfe4466550afe79255598f580f00a761fc32731a52ec9378",
        "PROVED",
    ),
    (
        "THM-3405",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
        "PROVED",
    ),
    (
        "THM-3429",
        ROOT / "01-canon/theorems/THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers.md",
        "58ebf850fc79fc9afed57966b7599e7376a6684fa3bbc5a2aa2e1a8e6e0ca148",
        "PROVED",
    ),
    (
        "THM-3453",
        ROOT / "01-canon/theorems/THM-3453-global-literal-half-twist-cap-seven-support-classification.md",
        "a4444813a2b24aa55613eef5f4cc0538ca0148f16297dff46f80d723f1beb247",
        "PROVED",
    ),
    (
        "THM-3461",
        ROOT / "01-canon/theorems/THM-3461-literal-half-twist-common-centre-lifts-and-q83-rank-nine-boundary.md",
        "371d39375740f97d896d722fbddbe124dc7d91a79d7fe7179e2fcc1bd9ce615f",
        "PROVED",
    ),
)

CARRIER_HASH_PINS = (
    (
        "THM-3473",
        ROOT / "01-canon/theorems/THM-3473-three-times-p-eight-owner-private-sheet-partition-and-irredundancy.md",
        "05236a7338a4c92b443b2798508d407e2e4f82d942111d022064f6ec9fe86ca2",
        "PROVED",
    ),
    (
        "PRIVATE-7x13-SCRIPT",
        ROOT / "04-computation/lrc_private_support_7x13_incidence_h1_probe_20260815.py",
        "6efa87aa9f9b50d57d7a2db3c282ad216057b41c10af08f648e3f3398e457b91",
        "FINITE-EXACT-SIDECAR",
    ),
    (
        "PRIVATE-7x13-OUTPUT",
        ROOT / "05-knowledge/results/lrc_private_support_7x13_incidence_h1_probe_20260815.out",
        "77a3d5786845fdac988e8a4d23867c699368dc6d0a5714224476e24357d7ae20",
        "FINITE-EXACT-SIDECAR",
    ),
)

EXPECTED_SEMANTIC_SHA256 = (
    "1bf53086d3da347dde443175462be716da6e9dac54c96582718d19ec8fddff21"
)
EXPECTED_CARRIER_SEMANTIC_SHA256 = (
    "3f17c2206feec73da48a989ab2150ceb1c7d1bc275c77291df476d882957581a"
)

REL_LABELS = ("c1", "c2", "c3", "H", "q1", "q2", "q3", "q4", "q5")
CUR_LABELS = ("H", "q1", "q2", "q3", "q4", "q5", "c1", "c2", "c3")
REL_TO_CURRENT = (3, 4, 5, 6, 7, 8, 0, 1, 2)
GUARD, OWNER, TARGET_A, TARGET_B = 0, 6, 7, 8
UNITS = (1, 2, 3, 4, 5)
OMITTED_UNIT = 5
R_DILATION = 13**2
X_FREQUENCY = 13
M_DEEP = 1

PRIVATE_OWNER_LABELS = tuple(f"u{index}" for index in range(1, 9))
PRIVATE_PACKETS = ((0, 3, 4, 5), (1, 2, 4, 7), (4, 6))
PRIVATE_EDGES = (
    (0, 3), (0, 4), (0, 5), (1, 2), (1, 4), (1, 7), (2, 4),
    (2, 7), (3, 4), (3, 5), (4, 5), (4, 6), (4, 7),
)

RELATION = (-27, -27, -27, 20110798, -41, -27, -27, -27, 38)
U_FULL_REL = (13, 2197, 742586, 1, 183, 27, 131, 53, 313)
U_CLOCK_REL = (65, 2197, 742586, 5, 661549, 655231, 658533, 661445, 291)
U_Q27_REL = (
    28405, 7599423, 18279868269, 3459, 2016, 2757, 1041, 3693,
    11163142875,
)
U_Q51_REL = (
    70993, 7199569, 30105550319, 5825, 4200, 5214, 7684, 4421,
    18313194875,
)

U_CLOCK_ATOM_U = (-3, -3, 1, 2840374, -3, 3, 2, -3, -48974)
U_CLOCK_ATOM_V = (24, 24, 29, -17270424, 38, 30, 29, 24, -49012)

PATTERN_E = {
    GUARD: "guard_safe",
    1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    OWNER: "in", TARGET_A: "out", TARGET_B: "out",
}
PATTERN_QA = {
    GUARD: "guard_safe",
    1: "out", 2: "out", 3: "out", 4: "out", 5: "out",
    OWNER: "out", TARGET_A: "in", TARGET_B: "out",
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_integers(values: tuple[int, ...]) -> str:
    body = ",".join(str(value) for value in values).encode("ascii")
    return sha256(body).hexdigest()


def to_current(values: tuple[int, ...]) -> tuple[int, ...]:
    require(len(values) == 9, ("coordinate count", values))
    return tuple(values[index] for index in REL_TO_CURRENT)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right))


def valuation(number: int, prime: int) -> int:
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def lcm_tuple(values: tuple[int, ...]) -> int:
    answer = 1
    for value in values:
        answer = answer * value // gcd(answer, value)
    return answer


def product_factorization(factors: tuple[tuple[int, int], ...]) -> int:
    answer = 1
    for prime, exponent in factors:
        answer *= prime**exponent
    return answer


def all_permutations(values: tuple[int, ...]):
    if not values:
        yield ()
        return
    for position, value in enumerate(values):
        remainder = values[:position] + values[position + 1:]
        for suffix in all_permutations(remainder):
            yield (value,) + suffix


def permutation_order(permutation: tuple[int, ...]) -> int:
    seen = set()
    answer = 1
    for start in range(len(permutation)):
        if start in seen:
            continue
        cursor = start
        length = 0
        while cursor not in seen:
            seen.add(cursor)
            cursor = permutation[cursor]
            length += 1
        answer = answer * length // gcd(answer, length)
    return answer


def carrier_map_obstruction_certificate() -> tuple[object, ...]:
    source_hashes = tuple(
        (label, status, lf_sha256(path))
        for label, path, _expected, status in CARRIER_HASH_PINS
    )
    expected_hashes = tuple(
        (label, status, expected)
        for label, _path, expected, status in CARRIER_HASH_PINS
    )
    require(source_hashes == expected_hashes,
            ("carrier source hash drift", source_hashes, expected_hashes))

    derived_edges = set()
    for packet in PRIVATE_PACKETS:
        for left_index in range(len(packet)):
            for right_index in range(left_index + 1, len(packet)):
                derived_edges.add(tuple(sorted(
                    (packet[left_index], packet[right_index])
                )))
    require(tuple(sorted(derived_edges)) == PRIVATE_EDGES,
            ("private carrier edge drift", derived_edges))

    degrees = tuple(
        sum(vertex in edge for edge in PRIVATE_EDGES) for vertex in range(8)
    )
    require(degrees == (3, 3, 3, 3, 7, 3, 1, 3),
            ("private carrier degrees", degrees))

    frozen_edges = frozenset(frozenset(edge) for edge in PRIVATE_EDGES)
    automorphisms = []
    for permutation in all_permutations(tuple(range(8))):
        image = frozenset(
            frozenset((permutation[left], permutation[right]))
            for left, right in PRIVATE_EDGES
        )
        if image == frozen_edges:
            automorphisms.append(tuple(permutation))
    automorphisms = tuple(automorphisms)
    require(len(automorphisms) == 72,
            ("private carrier automorphism count", len(automorphisms)))
    automorphism_orders = tuple(sorted({
        permutation_order(permutation) for permutation in automorphisms
    }))
    require(13 not in automorphism_orders,
            ("unexpected order-13 graph automorphism", automorphism_orders))

    remaining = set(PRIVATE_EDGES)
    edge_orbits = []
    while remaining:
        seed = min(remaining)
        orbit = tuple(sorted({
            tuple(sorted((permutation[seed[0]], permutation[seed[1]])))
            for permutation in automorphisms
        }))
        edge_orbits.append(orbit)
        remaining.difference_update(orbit)
    edge_orbits = tuple(sorted(edge_orbits, key=lambda row: (len(row), row)))
    require(tuple(len(orbit) for orbit in edge_orbits) == (1, 6, 6),
            ("private carrier edge orbits", edge_orbits))
    require(edge_orbits[0] == ((4, 6),),
            ("private carrier bridge orbit", edge_orbits[0]))

    source_line = tuple((0, beta) for beta in range(13))
    translated_line = tuple((alpha, (beta + 1) % 13)
                            for alpha, beta in source_line)
    require(set(translated_line) == set(source_line),
            "source C13 line is not translation-stable")
    source_translation = tuple((index + 1) % 13 for index in range(13))
    require(permutation_order(source_translation) == 13,
            "source translation loses order 13")

    shared_labels = tuple(sorted(set(REL_LABELS) & set(PRIVATE_OWNER_LABELS)))
    require(shared_labels == (), ("unexpected label dictionary", shared_labels))
    equivariant_bijection_exists = any(
        permutation_order(permutation) == 13 for permutation in automorphisms
    )
    require(not equivariant_bijection_exists,
            "C13-equivariant carrier bijection unexpectedly exists")

    return (
        source_hashes,
        REL_LABELS,
        PRIVATE_OWNER_LABELS,
        shared_labels,
        PRIVATE_PACKETS,
        PRIVATE_EDGES,
        degrees,
        len(automorphisms),
        automorphism_orders,
        edge_orbits,
        len(source_line),
        permutation_order(source_translation),
        equivariant_bijection_exists,
        "NO_C13_EQUIVARIANT_13-TWIST-TO-EDGE_BIJECTION",
    )


def small_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def verify_lucas_certificate(
    prime: int,
    factors: tuple[tuple[int, int], ...],
    witnesses: tuple[tuple[int, int], ...],
    label: str,
) -> None:
    require(product_factorization(factors) == prime - 1,
            (label, "incomplete p-1 factorization"))
    witness_map = dict(witnesses)
    require(set(witness_map) == {q for q, _exponent in factors},
            (label, "witness support"))
    for q, _exponent in factors:
        require(small_prime(q), (label, "composite factor", q))
        base = witness_map[q]
        require(pow(base, prime - 1, prime) == 1,
                (label, "Fermat failure", q))
        require(gcd(pow(base, (prime - 1) // q, prime) - 1, prime) == 1,
                (label, "Lucas gcd failure", q))


def verify_embedding(
    prime: int,
    root: int,
    nn: int,
    nn_factors: tuple[tuple[int, int], ...],
    label: str,
) -> None:
    require(pow(root, nn, prime) == 1, (label, "root power"))
    for q, _exponent in nn_factors:
        require(pow(root, nn // q, prime) != 1,
                (label, "root loses prime", q))


FULL_NN_FACTORS = (
    (2, 2), (3, 3), (7, 1), (13, 8), (53, 1), (61, 1), (131, 1),
    (313, 1),
)
FULL_EMBEDDINGS = (
    (
        572252886246508880869,
        279936,
        (
            (2, 2), (3, 3), (7, 2), (13, 8), (53, 1), (61, 1),
            (131, 1), (313, 1),
        ),
        (
            (2, 2), (3, 3), (7, 2), (13, 2), (53, 2), (61, 2),
            (131, 2), (313, 2),
        ),
    ),
    (
        1062755360172087921613,
        34522712143931,
        (
            (2, 2), (3, 3), (7, 1), (13, 9), (53, 1), (61, 1),
            (131, 1), (313, 1),
        ),
        (
            (2, 2), (3, 3), (7, 2), (13, 2), (53, 2), (61, 2),
            (131, 2), (313, 2),
        ),
    ),
)

CLOCK_NN_FACTORS = (
    (2, 2), (3, 1), (5, 1), (7, 3), (13, 8), (17, 1), (23, 1),
    (31, 1), (73, 1), (97, 1), (263, 1), (503, 1), (587, 1),
    (38543, 1),
)
CLOCK_EMBEDDINGS = (
    (
        77625621335503120714074176051966761,
        75047496554032956760519149093721,
        (
            (2, 3), (3, 3), (5, 1), (7, 3), (13, 8), (17, 1),
            (23, 1), (31, 1), (73, 1), (97, 1), (263, 1), (503, 1),
            (587, 1), (38543, 1),
        ),
        (
            (2, 29), (3, 3), (5, 2), (7, 2), (13, 2), (17, 2),
            (23, 2), (31, 2), (73, 2), (97, 2), (263, 2), (503, 2),
            (587, 2), (38543, 2),
        ),
    ),
    (
        112125897484615618809218254297285321,
        36743065895590462951408390936726169,
        (
            (2, 3), (3, 1), (5, 1), (7, 3), (13, 9), (17, 1),
            (23, 1), (31, 1), (73, 1), (97, 1), (263, 1), (503, 1),
            (587, 1), (38543, 1),
        ),
        (
            (2, 29), (3, 3), (5, 2), (7, 2), (13, 2), (17, 2),
            (23, 2), (31, 2), (73, 2), (97, 2), (263, 2), (503, 2),
            (587, 2), (38543, 2),
        ),
    ),
)


def rank_mod_13(vectors: tuple[tuple[int, ...], ...]) -> int:
    matrix = [[entry % 13 for entry in row] for row in vectors]
    rank = 0
    for column in range(9):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], 11, 13)
        matrix[rank] = [(inverse * value) % 13 for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                (matrix[row][col] - multiplier * matrix[rank][col]) % 13
                for col in range(9)
            ]
        rank += 1
    return rank


def owner_packet_rows(word: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    pivot_labels = (OWNER, GUARD, 1, 2, 3, 4)
    rows = []
    for label in pivot_labels:
        row = [0] * 9
        row[OMITTED_UNIT] += word[label]
        row[label] -= word[OMITTED_UNIT]
        if label == 1:
            row[OMITTED_UNIT] += word[TARGET_A]
            row[TARGET_A] -= word[OMITTED_UNIT]
        if label == 2:
            row[OMITTED_UNIT] += word[TARGET_B]
            row[TARGET_B] -= word[OMITTED_UNIT]
        rows.append(tuple(value % 13 for value in row))
    return tuple(rows)


def relation_packet_certificate(
    relation_word: tuple[int, ...],
) -> tuple[object, ...]:
    word = to_current(relation_word)
    relation_current = to_current(RELATION)
    packet = owner_packet_rows(word)
    word_mod = tuple(value % 13 for value in word)
    relation_mod = tuple(value % 13 for value in relation_current)
    target_axis = tuple(
        (relation_mod[index] - (1 if index == TARGET_A else 0)) % 13
        for index in range(9)
    )
    require(dot(RELATION, relation_word) == 0, "exact relation failure")
    require(all(gcd(abs(value), 91) == 1 for value in RELATION),
            "relation is not all-91-unit")
    require(all(dot(row, word_mod) % 13 == 0 for row in packet),
            "packet leaves K_13")
    packet_rank = rank_mod_13(packet)
    target_rank = rank_mod_13(packet + (target_axis,))
    relation_rank = rank_mod_13(packet + (relation_mod,))
    require((packet_rank, target_rank, relation_rank) == (6, 6, 7),
            ("target-axis quotient", packet_rank, target_rank, relation_rank))
    return (
        dot(RELATION, relation_word),
        tuple(value % 91 for value in RELATION),
        tuple(value % 13 for value in RELATION),
        (packet_rank, target_rank, relation_rank),
    )


def in_half_arc(word: int, modulus: int) -> bool:
    residue = word % modulus
    return 14 * min(residue, modulus - residue) < modulus


def half_mask(modulus: int, coefficient: int) -> tuple[int, ...]:
    return tuple(
        sheet
        for sheet in range(modulus)
        if in_half_arc(coefficient * (2 * sheet + 1), 2 * modulus)
    )


def normalized_half_residue(modulus: int, coefficient: int) -> int:
    residue = coefficient % (2 * modulus)
    if residue > modulus:
        residue = 2 * modulus - residue
    return residue


def multiplicity_profile(
    modulus: int,
    coefficients: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    counts = [0] * modulus
    for coefficient in coefficients:
        for sheet in half_mask(modulus, coefficient):
            counts[sheet] += 1
    return tuple(
        (value, counts.count(value))
        for value in sorted(set(counts))
    )


def mask_lift_certificate(
    name: str,
    modulus: int,
    relation_word: tuple[int, ...],
    witness: tuple[int, ...],
    expected_residues: tuple[int, ...],
    expected_mod7: tuple[int, ...] | None = None,
    expected_v13: tuple[int, ...] | None = None,
) -> tuple[object, ...]:
    residues = tuple(value % (2 * modulus) for value in relation_word)
    normalized = tuple(
        normalized_half_residue(modulus, value) for value in relation_word
    )
    require(residues == expected_residues, (name, "residue drift", residues))
    require(dot(RELATION, relation_word) == 0, (name, "exact relation"))
    require(dot(RELATION, residues) % (2 * modulus) == 0,
            (name, "relation congruence"))
    require(gcd(*relation_word) == 1, (name, "not primitive"))
    require(set(witness) <= set(normalized), (name, "missing masks"))
    profile = multiplicity_profile(modulus, witness)
    require(profile[0][0] >= 1, (name, "witness does not cover"))
    if expected_mod7 is not None:
        require(tuple(value % 7 for value in relation_word) == expected_mod7,
                (name, "mod-7 decoration"))
    if expected_v13 is not None:
        require(tuple(valuation(value, 13) for value in relation_word)
                == expected_v13, (name, "13-adic decoration"))
    mask_rows = tuple(half_mask(modulus, value) for value in relation_word)
    return (
        residues,
        normalized,
        mask_rows,
        witness,
        profile,
        tuple(value % 7 for value in relation_word),
        tuple(valuation(value, 13) for value in relation_word),
    )


def q51_affine_certificate() -> tuple[object, ...]:
    residues = tuple(value % 102 for value in U_Q51_REL)
    bars = tuple(value % 34 for value in residues)
    lifts = tuple(((value - bar) // 34) % 3
                  for value, bar in zip(residues, bars))
    require(lifts == (0, 0, 0, 0, 0, 0, 1, 1, 0),
            ("q51 lift characters", lifts))

    def fibre_block(coefficient: int, base: int) -> tuple[int, ...]:
        return tuple(
            t for t in range(17)
            if base + 3 * t in half_mask(51, coefficient)
        )

    rows = tuple(
        (
            base,
            fibre_block(34, base),
            fibre_block(1, base),
            fibre_block(35, base),
        )
        for base in range(3)
    )
    require(rows[0][1] == rows[2][1] == (), "inactive missed fibres")
    require(rows[1][1] == tuple(range(17)), "inactive fixed fibre")
    require(set(rows[0][2]).isdisjoint(rows[0][3]), "q51 y=0 lift collision")
    require(set(rows[2][2]).isdisjoint(rows[2][3]), "q51 y=2 lift collision")
    require(rows[1][2] == rows[1][3], "q51 fixed-fibre equality")
    return residues, bars, lifts, rows


def circle_distance(value: Fraction) -> Fraction:
    floor = value.numerator // value.denominator
    fractional = value - floor
    return min(fractional, 1 - fractional)


def pattern_is_strict(
    distances: tuple[Fraction, ...],
    pattern: dict[int, str],
) -> bool:
    for index, mode in pattern.items():
        if mode == "guard_safe" and not distances[index] > Fraction(1, 7):
            return False
        if mode == "in" and not distances[index] < Fraction(1, 14):
            return False
        if mode == "out" and not distances[index] > Fraction(1, 14):
            return False
    return True


def clock_centre_certificate() -> tuple[object, ...]:
    word = to_current(U_CLOCK_REL)
    centre = Fraction(1, 22)
    base_distances = tuple(circle_distance(value * centre) for value in word)
    future_distances = tuple(
        circle_distance(value * R_DILATION * centre) for value in word
    )
    expected_base = (
        Fraction(5, 22), Fraction(9, 22), Fraction(5, 22),
        Fraction(7, 22), Fraction(7, 22), Fraction(5, 22),
        Fraction(1, 22), Fraction(3, 22), Fraction(1, 11),
    )
    expected_future = (
        Fraction(9, 22), Fraction(3, 22), Fraction(9, 22),
        Fraction(5, 22), Fraction(5, 22), Fraction(9, 22),
        Fraction(7, 22), Fraction(1, 22), Fraction(4, 11),
    )
    require(base_distances == expected_base, "base half-twist distances")
    require(future_distances == expected_future, "future half-twist distances")
    require(pattern_is_strict(base_distances, PATTERN_E), "centre not in E")
    require(pattern_is_strict(future_distances, PATTERN_QA),
            "same clock not in Q_a")

    mode_centres = []
    amplitudes = []
    for speed in word:
        residue = speed % 22
        mode_index = (speed - residue) // 22
        lifted = (Fraction(mode_index, 1) + Fraction(residue, 22)) / speed
        mode_centres.append(lifted)
        amplitudes.append(22 * speed * lifted)
    require(tuple(mode_centres) == (centre,) * 9, "nonzero mode cochain")
    require(tuple(amplitudes) == word, "amplitude ray mismatch")
    wedges = tuple(
        amplitudes[i] * word[j] - word[i] * amplitudes[j]
        for i in range(9) for j in range(i + 1, 9)
    )
    require(set(wedges) == {0}, "projective wedge does not vanish")
    margins = (
        Fraction(1, 14) - base_distances[OWNER],
        base_distances[GUARD] - Fraction(1, 7),
        Fraction(1, 14) - future_distances[TARGET_A],
        future_distances[GUARD] - Fraction(1, 7),
    )
    require(margins == (
        Fraction(2, 77), Fraction(13, 154),
        Fraction(2, 77), Fraction(41, 154),
    ), "strict margin drift")
    return (
        centre,
        base_distances,
        future_distances,
        tuple(mode_centres),
        tuple(amplitudes),
        tuple(sorted(set(wedges))),
        margins,
    )


def atomic_term_certificate() -> tuple[object, ...]:
    e_c3 = tuple(1 if index == 2 else 0 for index in range(9))
    address = tuple(
        U_CLOCK_ATOM_U[index] + e_c3[index] - U_CLOCK_ATOM_V[index]
        for index in range(9)
    )
    left_frequency = dot(U_CLOCK_ATOM_U, U_CLOCK_REL)
    right_frequency = dot(U_CLOCK_ATOM_V, U_CLOCK_REL)
    require(address == RELATION, "atomic address mismatch")
    require(left_frequency == X_FREQUENCY, "atomic left frequency")
    require(right_frequency == X_FREQUENCY + U_CLOCK_REL[2],
            "atomic right frequency")
    require(all(value % 7 != 0
                for value in U_CLOCK_ATOM_U + U_CLOCK_ATOM_V),
            "atomic septimal zero")
    return (
        U_CLOCK_ATOM_U,
        U_CLOCK_ATOM_V,
        address,
        left_frequency,
        M_DEEP,
        right_frequency,
        max(abs(value) for value in U_CLOCK_ATOM_U),
        max(abs(value) for value in U_CLOCK_ATOM_V),
    )


# Fast periodic-subtraction engine used for the full U_full bank.
def fast_in_comb(
    word: tuple[int, ...],
    t_den: int,
    index: int,
    shift: int,
) -> list[tuple[int, int]]:
    speed = word[index]
    unit = t_den // (182 * speed)
    low = (-13 - 14 * shift) % 182
    intervals = []
    for branch in range(speed):
        left = (low + 182 * branch) * unit
        right = left + 26 * unit
        if right <= t_den:
            intervals.append((left, right))
        else:
            intervals.append((left, t_den))
            intervals.append((0, right - t_den))
    intervals.sort()
    return intervals


def fast_subtract_comb(
    intervals: list[tuple[int, int]],
    t_den: int,
    speed: int,
    period_denominator: int,
    low: int,
    high: int,
) -> list[tuple[int, int]]:
    unit = t_den // (period_denominator * speed)
    step = period_denominator * unit
    window_length = (high - low) * unit
    base = (low % period_denominator) * unit
    output = []
    for left, right in intervals:
        first_possible = left - window_length + 1
        branch = -((base - first_possible) // step)
        window_left = base + branch * step
        cursor = left
        while window_left < right:
            window_right = window_left + window_length
            if window_right > cursor:
                if window_left > cursor:
                    output.append((cursor, window_left))
                cursor = window_right
                if cursor >= right:
                    break
            window_left += step
        if cursor < right:
            output.append((cursor, right))
    return output


def fast_build_set(
    word: tuple[int, ...],
    t_den: int,
    pattern: dict[int, str],
    shift: tuple[int, ...],
) -> list[tuple[int, int]]:
    positives = [index for index, mode in pattern.items() if mode == "in"]
    require(positives, "fast set has no positive comb")
    start = min(positives, key=lambda index: word[index])
    intervals = fast_in_comb(word, t_den, start, shift[start])
    for index, mode in pattern.items():
        if mode == "guard_safe":
            intervals = fast_subtract_comb(
                intervals, t_den, word[index], 91,
                -13 - 7 * shift[index], 13 - 7 * shift[index],
            )
    rest = sorted(
        (word[index], index)
        for index, mode in pattern.items()
        if index != start and mode in ("in", "out")
    )
    for _speed, index in rest:
        if pattern[index] == "out":
            low, high = -13 - 14 * shift[index], 13 - 14 * shift[index]
        else:
            low, high = 13 - 14 * shift[index], 169 - 14 * shift[index]
        intervals = fast_subtract_comb(
            intervals, t_den, word[index], 182, low, high,
        )
    return intervals


def validate_intervals(
    intervals: list[tuple[int, int]],
    denominator: int,
    label: str,
    merge_required: bool = False,
) -> int:
    previous_right = -1
    for left, right in intervals:
        require(0 <= left < right <= denominator,
                (label, "invalid endpoint", left, right))
        if merge_required:
            require(left > previous_right, (label, "adjacent or overlap"))
        else:
            require(left >= previous_right, (label, "overlap"))
        previous_right = right
    return sum(right - left for left, right in intervals)


def fast_make_tabs(
    q_intervals: list[tuple[int, int]],
    frequency: int,
    nn: int,
    embeddings: tuple[tuple[int, int], ...],
) -> tuple[tuple[list[int], list[int]], ...]:
    answer = []
    for prime, root in embeddings:
        left = [pow(root, (-frequency * a) % nn, prime)
                for a, _b in q_intervals]
        right = [pow(root, (-frequency * b) % nn, prime)
                 for _a, b in q_intervals]
        answer.append((left, right))
    return tuple(answer)


def fast_x_sweep(
    e_intervals: list[tuple[int, int]],
    q_intervals: list[tuple[int, int]],
    q_starts: list[int],
    frequency: int,
    t_den: int,
    nn: int,
    embeddings: tuple[tuple[int, int], ...],
    tabs: tuple[tuple[list[int], list[int]], ...],
    want_values: bool = True,
) -> tuple[tuple[int, ...], int]:
    totals = [0] * len(embeddings)
    overlap = 0
    wrap_forward = [
        pow(root, (-frequency * t_den) % nn, prime)
        for prime, root in embeddings
    ] if want_values else []
    wrap_backward = [
        pow(root, (frequency * t_den) % nn, prime)
        for prime, root in embeddings
    ] if want_values else []
    for e_left, e_right in e_intervals:
        scaled_left = R_DILATION * e_left
        start = scaled_left % t_den
        span = R_DILATION * (e_right - e_left)
        require(span < t_den, ("fast sweep span", span, t_den))
        stop = start + span
        index = bisect_right(q_starts, start) - 1
        offset = 0
        if index < 0:
            index = len(q_intervals) - 1
            offset = -t_den
        if want_values:
            base_exponent = (-frequency * (scaled_left - start)) % nn
            base_values = [
                pow(root, base_exponent, prime) for prime, root in embeddings
            ]
            accumulators = [0] * len(embeddings)
            wrap_values = [1] * len(embeddings)
            if offset:
                wrap_values = list(wrap_backward)
        while True:
            q_left_0, q_right_0 = q_intervals[index]
            q_left = q_left_0 + offset
            q_right = q_right_0 + offset
            if q_left >= stop:
                break
            if q_right > start:
                left = max(start, q_left)
                right = min(stop, q_right)
                if left < right:
                    overlap += right - left
                    if want_values:
                        for position, (prime, root) in enumerate(embeddings):
                            if left == q_left:
                                value_left = (
                                    tabs[position][0][index]
                                    * wrap_values[position]
                                ) % prime
                            else:
                                value_left = pow(root, (-frequency * left) % nn,
                                                 prime)
                            if right == q_right:
                                value_right = (
                                    tabs[position][1][index]
                                    * wrap_values[position]
                                ) % prime
                            else:
                                value_right = pow(
                                    root, (-frequency * right) % nn, prime
                                )
                            accumulators[position] = (
                                accumulators[position] + value_left - value_right
                            ) % prime
            index += 1
            if index == len(q_intervals):
                index = 0
                offset += t_den
                if want_values:
                    for position, (prime, _root) in enumerate(embeddings):
                        wrap_values[position] = (
                            wrap_values[position] * wrap_forward[position]
                        ) % prime
        if want_values:
            for position, (prime, _root) in enumerate(embeddings):
                totals[position] = (
                    totals[position]
                    + base_values[position] * accumulators[position]
                ) % prime
    return tuple(totals), overlap


def fast_endpoint_sum(
    intervals: list[tuple[int, int]],
    frequency: int,
    nn: int,
    embeddings: tuple[tuple[int, int], ...],
) -> tuple[int, ...]:
    totals = [0] * len(embeddings)
    for left, right in intervals:
        exponent_left = (-frequency * R_DILATION * left) % nn
        exponent_right = (-frequency * R_DILATION * right) % nn
        for position, (prime, root) in enumerate(embeddings):
            totals[position] = (
                totals[position]
                + pow(root, exponent_left, prime)
                - pow(root, exponent_right, prime)
            ) % prime
    return tuple(totals)


def target_representatives() -> tuple[
    tuple[int, ...], tuple[int, ...], tuple[tuple[int, ...], ...]
]:
    v1 = (0, 12, 0, 0, 0, 0, 0, 1, 0)
    v2 = (0, 0, 12, 0, 0, 0, 0, 0, 1)
    reps = tuple(
        tuple((alpha * v1[index] + beta * v2[index]) % 13
              for index in range(9))
        for alpha in range(13) for beta in range(13)
    )
    require(len(set(reps)) == 169, "target representatives collide")
    return v1, v2, reps


def full_bank_certificate() -> tuple[object, ...]:
    word = to_current(U_FULL_REL)
    t_den = 182 * lcm_tuple(word)
    nn = R_DILATION * t_den
    require(t_den == 483730250419703196, ("U_full T_DEN", t_den))
    require(nn == 81750412320929840124, ("U_full NN", nn))
    require(product_factorization(FULL_NN_FACTORS) == nn, "U_full NN factors")
    p1, h1 = FULL_EMBEDDINGS[0][0], FULL_EMBEDDINGS[0][1]
    embeddings = ((p1, h1),)
    zero = (0,) * 9
    q_intervals = fast_build_set(word, t_den, PATTERN_QA, zero)
    q_starts = [left for left, _right in q_intervals]
    q_measure = Fraction(
        validate_intervals(q_intervals, t_den, "U_full fast Q"), t_den
    )
    require(q_measure == Fraction(197889477091847, 6201669877175682),
            ("U_full fast Q measure", q_measure))
    tabs = fast_make_tabs(q_intervals, X_FREQUENCY, nn, embeddings)
    v1, v2, reps = target_representatives()
    packet = owner_packet_rows(word)
    word_mod = tuple(value % 13 for value in word)
    require(rank_mod_13(packet) == 6, "U_full packet rank")
    require(all(dot(row, vector) % 13 == 0
                for row in packet for vector in (word_mod, v1, v2)),
            "target generators leave L-perp")
    require(rank_mod_13((word_mod, v1, v2)) == 3, "target quotient rank")

    gammas = []
    interval_count = 0
    for ell in reps:
        e_intervals = fast_build_set(word, t_den, PATTERN_E, ell)
        interval_count += len(e_intervals)
        ax, _overlap = fast_x_sweep(
            e_intervals, q_intervals, q_starts, X_FREQUENCY,
            t_den, nn, embeddings, tabs,
        )
        by = fast_endpoint_sum(
            e_intervals, -(X_FREQUENCY + word[TARGET_B]), nn, embeddings
        )
        phase = pow(
            h1, (M_DEEP * ell[TARGET_B] * (nn // 13)) % nn, p1
        )
        gammas.append(phase * ax[0] % p1 * by[0] % p1)
    gamma_tuple = tuple(gammas)
    gamma_zero = gamma_tuple[0]
    gamma_v2 = gamma_tuple[1]
    require((gamma_zero, gamma_v2) == (
        248447851579556601771,
        510954897124935772821,
    ), ("U_full fast gamma controls", gamma_zero, gamma_v2))
    differing = sum(
        value != gamma_zero for value in gamma_tuple[1:]
    )
    require(differing == 168, ("U_full nonconstant twists", differing))

    z13 = pow(h1, nn // 13, p1)
    z13_powers = tuple(pow(z13, exponent, p1) for exponent in range(13))
    pair_map = tuple(
        (ell[TARGET_A] % 13, ell[TARGET_B] % 13) for ell in reps
    )
    require(len(set(pair_map)) == 169, "target pairing not bijective")
    targets = []
    for qa in range(13):
        for qb in range(13):
            value = 0
            for gamma, (ta, tb) in zip(gamma_tuple, pair_map):
                exponent = (-(ta * qa + tb * qb)) % 13
                value = (value + gamma * z13_powers[exponent]) % p1
            targets.append(value)
    target_tuple = tuple(targets)
    nonzero_targets = sum(value != 0 for value in target_tuple)
    require(nonzero_targets == 169,
            ("U_full target support", nonzero_targets))
    require((sum(target_tuple) - 169 * gamma_zero) % p1 == 0,
            "inverse DFT sum control")
    for index in (0, 31):
        ta, tb = pair_map[index]
        reconstruction = 0
        for target_index, value in enumerate(target_tuple):
            qa, qb = divmod(target_index, 13)
            reconstruction = (
                reconstruction
                + value * z13_powers[(ta * qa + tb * qb) % 13]
            ) % p1
        require((reconstruction - 169 * gamma_tuple[index]) % p1 == 0,
                ("inverse DFT reconstruction", index))
    return (
        t_den,
        nn,
        len(q_intervals),
        q_measure,
        interval_count,
        gamma_zero,
        gamma_v2,
        (gamma_v2 - gamma_zero) % p1,
        differing,
        nonzero_targets,
        digest_integers(gamma_tuple),
        digest_integers(target_tuple),
        gamma_tuple,
        target_tuple,
    )


# Independent boundary-refinement engine.
def push_interval(
    output: list[tuple[int, int]],
    left: int,
    right: int,
) -> None:
    if left >= right:
        return
    if output and output[-1][1] == left:
        output[-1] = (output[-1][0], right)
    else:
        output.append((left, right))


def boundary_points(
    left: int,
    right: int,
    speed: int,
    shift: int,
    width_denominator: int,
    t_den: int,
) -> list[int]:
    period_denominator = 13 * width_denominator
    require(t_den % (period_denominator * speed) == 0,
            "nonintegral comb boundary")
    unit = t_den // (period_denominator * speed)
    step = period_denominator * unit
    points = []
    for sign in (-1, 1):
        base = (sign * 13 - width_denominator * shift) * unit
        branch = (left - base) // step + 1
        point = base + branch * step
        while point < right:
            points.append(point)
            point += step
    return sorted(set(points))


def midpoint_inside_comb(
    left: int,
    right: int,
    speed: int,
    shift: int,
    width_denominator: int,
    t_den: int,
) -> bool:
    midpoint_sum = left + right
    total_denominator = 26 * t_den
    numerator = (
        13 * speed * midpoint_sum + 2 * t_den * shift
    ) % total_denominator
    distance = min(numerator, total_denominator - numerator)
    return width_denominator * distance < total_denominator


def refine(
    intervals: list[tuple[int, int]],
    speed: int,
    shift: int,
    width_denominator: int,
    want_inside: bool,
    t_den: int,
) -> list[tuple[int, int]]:
    output = []
    for left, right in intervals:
        cuts = [left]
        cuts.extend(
            boundary_points(
                left, right, speed, shift, width_denominator, t_den
            )
        )
        cuts.append(right)
        for a, b in zip(cuts, cuts[1:]):
            inside = midpoint_inside_comb(
                a, b, speed, shift, width_denominator, t_den
            )
            if inside == want_inside:
                push_interval(output, a, b)
    return output


def reference_build_set(
    word: tuple[int, ...],
    t_den: int,
    pattern: dict[int, str],
    shift: tuple[int, ...],
) -> list[tuple[int, int]]:
    intervals = [(0, t_den)]
    positives = sorted(
        (word[index], index)
        for index, mode in pattern.items() if mode == "in"
    )
    require(positives, "reference set has no positive comb")
    for _speed, index in positives:
        intervals = refine(
            intervals, word[index], shift[index], 14, True, t_den
        )
    if pattern.get(GUARD) == "guard_safe":
        intervals = refine(
            intervals, word[GUARD], shift[GUARD], 7, False, t_den
        )
    exclusions = sorted(
        (word[index], index)
        for index, mode in pattern.items() if mode == "out"
    )
    for _speed, index in exclusions:
        intervals = refine(
            intervals, word[index], shift[index], 14, False, t_den
        )
    validate_intervals(
        intervals, t_den, "reference Boolean set", merge_required=True
    )
    return intervals


def reference_endpoint_sum(
    intervals: list[tuple[int, int]],
    frequency: int,
    nn: int,
    embeddings: tuple[tuple[int, int], ...],
) -> tuple[int, ...]:
    totals = [0] * len(embeddings)
    for left, right in intervals:
        exponent_left = (-frequency * R_DILATION * left) % nn
        exponent_right = (-frequency * R_DILATION * right) % nn
        for position, (prime, root) in enumerate(embeddings):
            totals[position] = (
                totals[position]
                + pow(root, exponent_left, prime)
                - pow(root, exponent_right, prime)
            ) % prime
    return tuple(totals)


def reference_marked_sum(
    present: list[tuple[int, int]],
    word_intervals: list[tuple[int, int]],
    frequency: int,
    t_den: int,
    nn: int,
    embeddings: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, ...], int, int]:
    scaled_present = [
        (R_DILATION * left, R_DILATION * right)
        for left, right in present
    ]
    present_index = 0
    totals = [0] * len(embeddings)
    total_length = 0
    components = 0
    for branch in range(R_DILATION):
        offset = branch * t_den
        for q_left_0, q_right_0 in word_intervals:
            q_left = offset + q_left_0
            q_right = offset + q_right_0
            while (
                present_index < len(scaled_present)
                and scaled_present[present_index][1] <= q_left
            ):
                present_index += 1
            scan = present_index
            while (
                scan < len(scaled_present)
                and scaled_present[scan][0] < q_right
            ):
                e_left, e_right = scaled_present[scan]
                left = max(e_left, q_left)
                right = min(e_right, q_right)
                if left < right:
                    total_length += right - left
                    components += 1
                    exponent_left = (-frequency * left) % nn
                    exponent_right = (-frequency * right) % nn
                    for position, (prime, root) in enumerate(embeddings):
                        totals[position] = (
                            totals[position]
                            + pow(root, exponent_left, prime)
                            - pow(root, exponent_right, prime)
                        ) % prime
                scan += 1
    return tuple(totals), total_length, components


FULL_REFERENCE_EXPECTED = {
    "t_den": 483730250419703196,
    "nn": 81750412320929840124,
    "e_count": 33810,
    "ev2_count": 34560,
    "q_count": 28730,
    "e_measure": Fraction(13426314167250, 344537215398649),
    "q_measure": Fraction(197889477091847, 6201669877175682),
    "mass": Fraction(411318338170045, 524041104621345129),
    "components": (123752, 126429),
    "ax0": (263333472381374948713, 943470807284885900842),
    "axv2": (523949679655523736011, 930290002246125433456),
    "by0": (367849901592567656326, 478322449598622688750),
    "byv2": (87230428072980866510, 228380476345275410214),
    "gamma0": (248447851579556601771, 595260918637066437724),
    "gammav2": (510954897124935772821, 6877617630528951420),
    "diff": (262507045545379171050, 474372059165550435309),
}

CLOCK_REFERENCE_EXPECTED = {
    "t_den": 25517955731592084389899466157780,
    "nn": 4312534518639062261893009780664820,
    "e_count": 147372,
    "ev2_count": 147404,
    "q_count": 136158,
    "e_measure": Fraction(
        1077696235371926419093062388,
        28866465759719552477261839545,
    ),
    "q_measure": Fraction(
        1118165137485116246305534177,
        32715327861015492807563418151,
    ),
    "mass": Fraction(
        1397606991636352080199080533,
        1692517471993352536064760510465,
    ),
    "components": (558850, 572967),
    "ax0": (
        6191221409511609461312041450920889,
        27033469806232578724428009703504103,
    ),
    "axv2": (
        8625389282305270646556668770397976,
        32769511973016243540285242339670264,
    ),
    "by0": (
        898615200741111596208809819456216,
        53687201923418596111341183735189556,
    ),
    "byv2": (
        5143445290392984781545278954609013,
        10305622304044080397547915601054092,
    ),
    "gamma0": (
        56767723330345680038743661041266194,
        65234233976034532625816110096140982,
    ),
    "gammav2": (
        34870555972766792317398130208739733,
        74671298727704698408794173004769050,
    ),
    "diff": (
        55728453977924232992728645219440300,
        9437064751670165782978062908628068,
    ),
}


def reference_endpoint_certificate(
    label: str,
    relation_word: tuple[int, ...],
    embedding_data: tuple[
        tuple[int, int, tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]],
        ...,
    ],
    nn_factors: tuple[tuple[int, int], ...],
    expected: dict[str, object],
) -> tuple[object, ...]:
    word = to_current(relation_word)
    t_den = 182 * lcm_tuple(word)
    nn = R_DILATION * t_den
    require((t_den, nn) == (expected["t_den"], expected["nn"]),
            (label, "scale drift", t_den, nn))
    require(product_factorization(nn_factors) == nn, (label, "NN factors"))
    embeddings = tuple((row[0], row[1]) for row in embedding_data)
    for position, (prime, root, p_factors, witnesses) in enumerate(
        embedding_data, start=1
    ):
        verify_lucas_certificate(
            prime, p_factors, witnesses, f"{label}-p{position}"
        )
        verify_embedding(
            prime, root, nn, nn_factors, f"{label}-embedding-{position}"
        )

    zero = (0,) * 9
    _v1, v2, _reps = target_representatives()
    q_intervals = reference_build_set(word, t_den, PATTERN_QA, zero)
    e_zero = reference_build_set(word, t_den, PATTERN_E, zero)
    e_v2 = reference_build_set(word, t_den, PATTERN_E, v2)
    e_measure = Fraction(
        validate_intervals(
            e_zero, t_den, f"{label} E0", merge_required=True
        ),
        t_den,
    )
    q_measure = Fraction(
        validate_intervals(
            q_intervals, t_den, f"{label} Q", merge_required=True
        ),
        t_den,
    )
    require(
        (len(e_zero), len(e_v2), len(q_intervals))
        == (expected["e_count"], expected["ev2_count"], expected["q_count"]),
        (label, "interval counts"),
    )
    require((e_measure, q_measure)
            == (expected["e_measure"], expected["q_measure"]),
            (label, "set measures", e_measure, q_measure))

    ax_zero, mass_length, components_zero = reference_marked_sum(
        e_zero, q_intervals, X_FREQUENCY, t_den, nn, embeddings
    )
    ax_v2, _mass_v2, components_v2 = reference_marked_sum(
        e_v2, q_intervals, X_FREQUENCY, t_den, nn, embeddings
    )
    y_frequency = X_FREQUENCY + word[TARGET_B]
    by_zero = reference_endpoint_sum(e_zero, -y_frequency, nn, embeddings)
    by_v2 = reference_endpoint_sum(e_v2, -y_frequency, nn, embeddings)
    mass = Fraction(mass_length, nn)
    require(mass == expected["mass"], (label, "delayed mass", mass))
    require((components_zero, components_v2) == expected["components"],
            (label, "component counts"))
    require(ax_zero == expected["ax0"], (label, "AX0"))
    require(ax_v2 == expected["axv2"], (label, "AXv2"))
    require(by_zero == expected["by0"], (label, "BY0"))
    require(by_v2 == expected["byv2"], (label, "BYv2"))

    gamma_zero = []
    gamma_v2 = []
    differences = []
    for position, (prime, root) in enumerate(embeddings):
        phase = pow(
            root, (M_DEEP * v2[TARGET_B] * (nn // 13)) % nn, prime
        )
        zero_value = ax_zero[position] * by_zero[position] % prime
        v2_value = phase * ax_v2[position] % prime * by_v2[position] % prime
        difference = (v2_value - zero_value) % prime
        require(difference != 0, (label, "constant endpoint bank", position))
        gamma_zero.append(zero_value)
        gamma_v2.append(v2_value)
        differences.append(difference)
    require(tuple(gamma_zero) == expected["gamma0"], (label, "gamma0"))
    require(tuple(gamma_v2) == expected["gammav2"], (label, "gammav2"))
    require(tuple(differences) == expected["diff"], (label, "gamma diff"))
    return (
        t_den,
        nn,
        (len(e_zero), len(e_v2), len(q_intervals)),
        e_measure,
        q_measure,
        mass,
        (components_zero, components_v2),
        ax_zero,
        ax_v2,
        by_zero,
        by_v2,
        tuple(gamma_zero),
        tuple(gamma_v2),
        tuple(differences),
    )


def security_certificate(source: Path) -> tuple[object, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert statement found")
    forbidden_calls = {"eval", "exec", "compile", "__import__"}
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    require(not (called & forbidden_calls),
            ("forbidden calls", called & forbidden_calls))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "bisect", "fractions", "hashlib", "math",
        "pathlib",
    }
    require(imported <= allowed, ("unexpected imports", imported - allowed))
    return (
        "NO_ASSERT_AST",
        tuple(sorted(imported)),
        tuple(sorted(forbidden_calls)),
    )


def main() -> None:
    source = Path(__file__)
    dependency_hashes = tuple(
        (label, status, lf_sha256(path))
        for label, path, _expected, status in HASH_PINS
    )
    expected_hashes = tuple(
        (label, status, expected)
        for label, _path, expected, status in HASH_PINS
    )
    require(dependency_hashes == expected_hashes,
            ("dependency hash drift", dependency_hashes, expected_hashes))
    security = security_certificate(source)
    carrier_obstruction = carrier_map_obstruction_certificate()
    carrier_semantic_sha256 = sha256(
        repr(carrier_obstruction).encode("utf-8")
    ).hexdigest()
    if EXPECTED_CARRIER_SEMANTIC_SHA256 is not None:
        require(
            carrier_semantic_sha256 == EXPECTED_CARRIER_SEMANTIC_SHA256,
            ("carrier semantic digest", carrier_semantic_sha256,
             EXPECTED_CARRIER_SEMANTIC_SHA256),
        )

    relation_full = relation_packet_certificate(U_FULL_REL)
    relation_clock = relation_packet_certificate(U_CLOCK_REL)
    full_masks = mask_lift_certificate(
        "U_full", 11, U_FULL_REL, (1, 2, 3, 5, 7, 9),
        (13, 19, 20, 1, 7, 5, 21, 9, 5),
    )
    clock_masks = mask_lift_certificate(
        "U_clock", 11, U_CLOCK_REL, (1, 2, 3, 5, 7, 9),
        (21, 19, 20, 5, 9, 5, 7, 15, 5),
    )
    q27 = mask_lift_certificate(
        "U_q27", 27, U_Q27_REL, (3, 15, 18, 21),
        (1, 3, 3, 3, 18, 3, 15, 21, 3),
        (6, 6, 5, 1, 0, 6, 5, 4, 3),
        (1, 3, 5, 0, 0, 0, 0, 0, 0),
    )
    q51 = mask_lift_certificate(
        "U_q51", 51, U_Q51_REL, (1, 11, 12, 18, 23, 34, 35),
        (1, 1, 11, 11, 18, 12, 34, 35, 23),
        (6, 6, 5, 1, 0, 6, 5, 4, 3),
        (1, 3, 5, 0, 0, 0, 0, 0, 0),
    )
    q51_affine = q51_affine_certificate()
    clock_centre = clock_centre_certificate()
    atomic_term = atomic_term_certificate()

    for label, embedding_data, nn_factors, expected in (
        ("U_full", FULL_EMBEDDINGS, FULL_NN_FACTORS, FULL_REFERENCE_EXPECTED),
        ("U_clock", CLOCK_EMBEDDINGS, CLOCK_NN_FACTORS,
         CLOCK_REFERENCE_EXPECTED),
    ):
        for position, (prime, root, p_factors, witnesses) in enumerate(
            embedding_data, start=1
        ):
            verify_lucas_certificate(
                prime, p_factors, witnesses, f"{label}-preflight-p{position}"
            )
            verify_embedding(
                prime, root, int(expected["nn"]), nn_factors,
                f"{label}-preflight-embedding-{position}",
            )

    full_bank = full_bank_certificate()
    full_reference = reference_endpoint_certificate(
        "U_full", U_FULL_REL, FULL_EMBEDDINGS, FULL_NN_FACTORS,
        FULL_REFERENCE_EXPECTED,
    )
    clock_reference = reference_endpoint_certificate(
        "U_clock", U_CLOCK_REL, CLOCK_EMBEDDINGS, CLOCK_NN_FACTORS,
        CLOCK_REFERENCE_EXPECTED,
    )
    require(
        (full_bank[5], full_bank[6], full_bank[7])
        == (
            full_reference[11][0],
            full_reference[12][0],
            full_reference[13][0],
        ),
        "fast/reference U_full endpoint mismatch",
    )

    semantic_payload = (
        dependency_hashes,
        security,
        RELATION,
        relation_full,
        relation_clock,
        full_masks,
        clock_masks,
        q27,
        q51,
        q51_affine,
        clock_centre,
        atomic_term,
        full_bank,
        full_reference,
        clock_reference,
    )
    semantic_sha256 = sha256(
        repr(semantic_payload).encode("utf-8")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256,
                 EXPECTED_SEMANTIC_SHA256))

    print("THM-3479 half-twist relation-current two-transplant exact companion")
    print("status=PROVED STRUCTURAL + VERIFIED-EXACT + FINITE-EXACT + INDEPENDENTLY AUDITED; core companion")
    print(f"script={SCRIPT}")
    print(f"stored_output={OUTPUT}")
    print(f"dependency_hashes={dependency_hashes}")
    print(f"security_certificate={security}")
    print(f"relation_labels={REL_LABELS}")
    print(f"relation={RELATION}")
    print(f"relation_mod_91={relation_full[1]}")
    print(f"relation_mod_13={relation_full[2]}")
    print(f"target_axis_ranks=(U_full,U_clock)=({relation_full[3]},{relation_clock[3]})")
    print("")
    print("[U_full]")
    print(f"tuple_relation_order={U_FULL_REL}")
    print(f"tuple_current_order={to_current(U_FULL_REL)}")
    print(f"exact_relation={relation_full[0]}  primitive={gcd(*U_FULL_REL)}")
    print(f"q11_residues_mod22={full_masks[0]}")
    print(f"q11_normalized={full_masks[1]}")
    print(f"q11_masks={full_masks[2]}")
    print(f"q11_witness_profile={full_masks[4]}")
    print(f"delayed_mass={full_reference[5]}")
    print(f"full_bank_scales=(T_DEN,NN)={full_bank[:2]}")
    print(f"full_bank_gamma=(zero,v2,diff)={full_bank[5:8]}")
    print(f"full_bank_nonconstant_twists={full_bank[8]}/168")
    print(f"full_bank_nonzero_unrestricted_targets={full_bank[9]}/169")
    print(f"full_bank_gamma_sha256={full_bank[10]}")
    print(f"full_bank_target_sha256={full_bank[11]}")
    print(f"independent_intervals=(E0,Ev2,Qa)={full_reference[2]}")
    print(f"independent_endpoint_images=(AX0,AXv2,BY0,BYv2)={full_reference[7:11]}")
    print("")
    print("[U_clock]")
    print(f"tuple_relation_order={U_CLOCK_REL}")
    print(f"tuple_current_order={to_current(U_CLOCK_REL)}")
    print(f"exact_relation={relation_clock[0]}  primitive={gcd(*U_CLOCK_REL)}")
    print(f"valuation_profile={tuple(valuation(value,13) for value in U_CLOCK_REL)}")
    print(f"q11_residues_mod22={clock_masks[0]}")
    print(f"q11_normalized={clock_masks[1]}")
    print(f"q11_masks={clock_masks[2]}")
    print(f"q11_witness_profile={clock_masks[4]}")
    print(f"zero_centre={clock_centre[0]}  projective_wedge_values={clock_centre[5]}")
    print(f"same_clock_base_distances={clock_centre[1]}")
    print(f"same_clock_future_distances={clock_centre[2]}")
    print(f"same_clock_margins={clock_centre[6]}")
    print(f"delayed_mass={clock_reference[5]}")
    print(f"independent_intervals=(E0,Ev2,Qa)={clock_reference[2]}")
    print(f"endpoint_gamma=(zero,v2,diff)={clock_reference[11:14]}")
    print(f"atomic_term=(u,v,address,X,m,Y,max_u,max_v)={atomic_term}")
    print("")
    print("[q27/q51 congruence lifts]")
    print(f"q27_tuple={U_Q27_REL}")
    print(f"q27=(residues,normalized,witness,profile,mod7,v13)={(q27[0],q27[1],q27[3],q27[4],q27[5],q27[6])}")
    print(f"q51_tuple={U_Q51_REL}")
    print(f"q51=(residues,normalized,witness,profile,mod7,v13)={(q51[0],q51[1],q51[3],q51[4],q51[5],q51[6])}")
    print(f"q51_affine=(residues,bars_mod34,k_mod3,fibres)={q51_affine}")
    print("")
    print("[7x13 private-support carrier map audit]")
    print(f"carrier_source_hashes={carrier_obstruction[0]}")
    print(f"source_relation_labels={carrier_obstruction[1]}")
    print(f"target_private_owner_labels={carrier_obstruction[2]}")
    print(f"shared_label_dictionary={carrier_obstruction[3]}")
    print(f"target_packets_and_edges={(carrier_obstruction[4], carrier_obstruction[5])}")
    print(f"target_degrees={carrier_obstruction[6]}")
    print(f"target_automorphisms_and_orders={(carrier_obstruction[7], carrier_obstruction[8])}")
    print(f"target_edge_orbit_sizes={tuple(len(row) for row in carrier_obstruction[9])}")
    print(f"source_character_line_and_translation_order={(carrier_obstruction[10], carrier_obstruction[11])}")
    print(f"c13_equivariant_bijection_exists={carrier_obstruction[12]}")
    print(f"carrier_map_verdict={carrier_obstruction[13]}")
    print("SOURCE=U_full/U_clock nine-coordinate relation tuples and the F13^2 coordinate-twist bank")
    print("TARGET=thirteen edges of the proved THM-3473 private-support two-section")
    print("CANDIDATE_MAP=a bijection from one translation-stable 13-character fibre to the thirteen edges")
    print("PRESERVED_PREDICATE_SOUGHT=nonzero endpoint weights yielding nonzero bridge weight and both K4 tree sums")
    print("LOST_INFORMATION=relation address, clocks, q11 masks, endpoint phase, private-sheet counts, and k-mod-3 state")
    print("REQUIRED_SIDECAR=an explicit labelled role chart for the weighted determinant; phase-holonomy is still required for absolute H1 or physical realization")
    print("CHEAPEST_HOSTILES=no order-13 graph automorphism; unique bridge edge (u5,u7); two sixteen-term K4 tree sums")
    print("")
    print("PROVED_BY_COMPANION=explicit relation; literal masks; delayed masses; U_full all 169 unrestricted A(q); U_clock zero-centre same-clock delayed word and nonconstant endpoint bank; atomic term; q27/q51 lifts")
    print("OBSTRUCTION=the native 13-character C13 action has no equivariant bijection to the private-support edge carrier; later 72 labelled role charts are non-equivariant and noncanonical")
    print("OPEN=C(relation;X,m) for both tuples; U_clock A(e_c2) and refined bank; coupled all-91-unit B(q); source-native or C13-equivariant endpoint-to-carrier map; phase/holonomy; ancestry/bispectrum; lawful Boolean THM-2512 current; scalar-row exclusion; LRC(14)")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"carrier_semantic_sha256={carrier_semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")
    print("ALL DETERMINISTIC EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
