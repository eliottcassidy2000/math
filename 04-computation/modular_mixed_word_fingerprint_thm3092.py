#!/usr/bin/env python3
"""Exact marked C2*C3 quotient atlas for THM-3092."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter, defaultdict, deque
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / "04-computation/affine_projective_prime_power_handshake_thm3090.py"
DEPENDENCY_SHA256 = "ffc396b2bb731c5e89859f5ea3652d707c964e3d89a1ebf4dd1dbe47793ffcc9"
OUTPUT = ROOT / "05-knowledge/results/modular_mixed_word_fingerprint_thm3092.out"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3090 dependency changed")
spec = importlib.util.spec_from_file_location("thm3090", DEPENDENCY)
require(spec is not None and spec.loader is not None, "dependency import")
thm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(thm)


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def identity(degree):
    return tuple(range(degree))


def power(permutation, exponent):
    value = identity(len(permutation))
    base = permutation
    while exponent:
        if exponent & 1:
            value = compose(base, value)
        base = compose(base, base)
        exponent >>= 1
    return value


def permutation_order(permutation):
    for exponent in range(1, 1000):
        if power(permutation, exponent) == identity(len(permutation)):
            return exponent
    raise RuntimeError("order bound")


def cycle_type(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        current = min(unseen)
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def generated_subgroup(generators):
    one = identity(len(generators[0]))
    seen = {one}
    queue = deque([one])
    while queue:
        value = queue.popleft()
        for generator in generators:
            target = compose(generator, value)
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return frozenset(seen)


def point_orbits(group, degree):
    unseen = set(range(degree))
    sizes = []
    while unseen:
        seed = min(unseen)
        orbit = {element[seed] for element in group}
        require(orbit <= unseen, (degree, seed, "orbit overlap"))
        unseen -= orbit
        sizes.append(len(orbit))
    return tuple(sorted(sizes))


def translations(field):
    points = tuple(product(field.elements, repeat=2))
    lookup = {point: index for index, point in enumerate(points)}
    rows = []
    for shift in points:
        image = []
        for point in points:
            target = (
                field.add(point[0], shift[0]),
                field.add(point[1], shift[1]),
            )
            image.append(lookup[target])
        rows.append(tuple(image))
    return frozenset(rows)


def center(group):
    return tuple(
        element
        for element in group
        if all(compose(element, other) == compose(other, element) for other in group)
    )


def atlas(group):
    involutions = tuple(element for element in group if permutation_order(element) == 2)
    threes = tuple(element for element in group if permutation_order(element) == 3)
    joint = Counter()
    cycles = Counter()
    subgroup_orbits = Counter()
    word_subgroup_counts = defaultdict(Counter)
    generating_words = []
    for involution in involutions:
        for three in threes:
            word = compose(involution, three)
            word_order = permutation_order(word)
            subgroup = generated_subgroup((involution, three))
            subgroup_order = len(subgroup)
            joint[(word_order, subgroup_order)] += 1
            cycles[(word_order, cycle_type(word))] += 1
            subgroup_orbits[(word_order, subgroup_order, point_orbits(subgroup, len(word)))] += 1
            word_subgroup_counts[word][subgroup_order] += 1
            if subgroup_order == len(group):
                generating_words.append(word)
    return {
        "involutions": len(involutions),
        "threes": len(threes),
        "joint": joint,
        "cycles": cycles,
        "subgroup_orbits": subgroup_orbits,
        "word_subgroup_counts": word_subgroup_counts,
        "generating_words": tuple(generating_words),
    }


S4_AFFINE = atlas(thm.AFFINE2)
S4_PROJECTIVE = atlas(thm.PROJECTIVE3)
AFFINE_COUNTERFEIT = atlas(thm.AFFINE3)
PROJECTIVE_COUNTERFEIT = atlas(thm.PROJECTIVE8)

S4_JOINT = {(2, 6): 24, (3, 12): 24, (4, 24): 24}
AFFINE_JOINT = {
    (2, 6): 576,
    (6, 6): 144,
    (6, 18): 720,
    (6, 54): 864,
    (8, 48): 432,
    (8, 432): 864,
}
PROJECTIVE_JOINT = {(2, 6): 504, (7, 504): 1512, (9, 504): 1512}

require(S4_AFFINE["joint"] == S4_JOINT, S4_AFFINE["joint"])
require(S4_PROJECTIVE["joint"] == S4_JOINT, S4_PROJECTIVE["joint"])
require(AFFINE_COUNTERFEIT["joint"] == AFFINE_JOINT, AFFINE_COUNTERFEIT["joint"])
require(PROJECTIVE_COUNTERFEIT["joint"] == PROJECTIVE_JOINT, PROJECTIVE_COUNTERFEIT["joint"])

require(
    S4_AFFINE["cycles"]
    == {(2, (1, 1, 2)): 24, (3, (1, 3)): 24, (4, (4,)): 24},
    S4_AFFINE["cycles"],
)
require(
    AFFINE_COUNTERFEIT["cycles"][(8, (1, 8))] == 1296,
    AFFINE_COUNTERFEIT["cycles"],
)
require(
    AFFINE_COUNTERFEIT["subgroup_orbits"][(8, 48, (1, 8))] == 432
    and AFFINE_COUNTERFEIT["subgroup_orbits"][(8, 432, (9,))] == 864,
    AFFINE_COUNTERFEIT["subgroup_orbits"],
)
AFFINE_SPLIT_WORDS = frozenset(
    word
    for word, counts in AFFINE_COUNTERFEIT["word_subgroup_counts"].items()
    if counts[48]
)
AFFINE_FULL_WORDS = frozenset(
    word
    for word, counts in AFFINE_COUNTERFEIT["word_subgroup_counts"].items()
    if counts[432]
)
require(
    AFFINE_SPLIT_WORDS == AFFINE_FULL_WORDS and len(AFFINE_SPLIT_WORDS) == 108,
    (len(AFFINE_SPLIT_WORDS), len(AFFINE_FULL_WORDS)),
)
require(
    all(
        AFFINE_COUNTERFEIT["word_subgroup_counts"][word][48] == 4
        and AFFINE_COUNTERFEIT["word_subgroup_counts"][word][432] == 8
        for word in AFFINE_SPLIT_WORDS
    ),
    "AGL exact word multiplicities",
)
require(
    PROJECTIVE_COUNTERFEIT["cycles"]
    == {(2, (1, 2, 2, 2, 2)): 504, (7, (1, 1, 7)): 1512, (9, (9,)): 1512},
    PROJECTIVE_COUNTERFEIT["cycles"],
)

for group in (thm.AFFINE2, thm.PROJECTIVE3, thm.AFFINE3, thm.PROJECTIVE8):
    require(center(group) == (identity(len(group[0])),), (len(group), "center"))

BINARY_TRANSLATIONS = translations(thm.F2)
TERNARY_TRANSLATIONS = translations(thm.F3)
PROJECTIVE_V4 = frozenset(
    element
    for element in thm.PROJECTIVE3
    if element == identity(4) or cycle_type(element) == (2, 2)
)
require(len(BINARY_TRANSLATIONS) == 4, "binary translations")
require(len(TERNARY_TRANSLATIONS) == 9, "ternary translations")
require(len(PROJECTIVE_V4) == 4, "projective V4")
require(
    all(
        power(word, 2) in BINARY_TRANSLATIONS
        and power(word, 2) != identity(4)
        for word in S4_AFFINE["generating_words"]
    ),
    "S4 half-face",
)
require(
    all(
        power(word, 2) in PROJECTIVE_V4
        and power(word, 2) != identity(4)
        for word in S4_PROJECTIVE["generating_words"]
    ),
    "projective S4 half-face",
)
require(
    all(
        power(word, 4) not in TERNARY_TRANSLATIONS
        and cycle_type(power(word, 4)) == (1, 2, 2, 2, 2)
        for word in AFFINE_COUNTERFEIT["generating_words"]
    ),
    "AGL2(F3) half-face",
)

require(len(S4_AFFINE["generating_words"]) == 24, "S4 epimorphisms")
require(len(AFFINE_COUNTERFEIT["generating_words"]) == 864, "AGL epimorphisms")
require(len(PROJECTIVE_COUNTERFEIT["generating_words"]) == 3024, "PGL epimorphisms")
INNER_CLASSES = (24 // 24, 864 // 432, 3024 // 504)
require(INNER_CLASSES == (1, 2, 6), INNER_CLASSES)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    lines = [
        "Modular mixed-word fingerprint and septimal counterfeit separation",
        f"dependency_thm3090_sha256={DEPENDENCY_SHA256}",
        "S4;involutions=9;order3=8;joint=(2,6):24,(3,12):24,(4,24):24;generating_order=4;epimorphisms=24;inner_classes=1;word_cycle=4;half_face=nonzero_V4_translation",
        "AGL2F3;involutions=45;order3=80;joint=(2,6):576,(6,6):144,(6,18):720,(6,54):864,(8,48):432,(8,432):864;generating_order=8;epimorphisms=864;inner_classes=2;word_cycle=1+8;fourth_power=1+2+2+2+2_not_translation",
        "PGL2F8;involutions=63;order3=56;joint=(2,6):504,(7,504):1512,(9,504):1512;generating_orders=7,9;epimorphisms=3024;inner_classes=6;word_cycles=1+1+7_or_9",
        "affine_origin_hostile=split_and_full_rows_have_the_same_108_order8_word_permutations;pair_multiplicities=4_vs8",
        "septimal_split=order7_words_are_split_torus_cycles_1+1+7;order9_words_are_nonsplit_9-cycles;degree9_does_not_choose_between_them",
        "C7_scope=PGL_split_tori_are_only_abstractly_isomorphic_to_F8*;the_GL_scalar_fibre_dies_projectively",
        "mixed_word_gate=the_first_C2*C3_word_sc_separates_S4_from_both_degree9_counterfeits;only_S4_has_a_nonzero_translation_half-face",
        "scope=marked_finite_quotient_and_natural_permutation_actions_only;no_canonical_modular_marking,quartic_owner,Keller,LRC,or_tree_identification",
        "all_exact_pair_subgroup_cycle_center_and_half_face_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
