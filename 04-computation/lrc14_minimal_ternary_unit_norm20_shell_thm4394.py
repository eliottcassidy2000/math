#!/usr/bin/env python3
"""Exact producer verifier for the minimal ternary-unit l1=20 shell.

This verifier reuses the audited THM-4393 carrier primitives and freezes every
norm-twenty sector, proof window, maximum, and multiple-relation witness.
An independent verifier must not import this implementation.
"""

from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path


ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / "01-canon").is_dir() and
            (parent / "04-computation").is_dir())
SOURCE = ROOT / "04-computation" / "lrc14_minimal_ternary_unit_norm18_shell_thm4393.py"
SPEC = importlib.util.spec_from_file_location("lrc_l18", SOURCE)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

HEIGHT = 500
SHORT = M.shell_patterns(18)
SHELL = M.shell_patterns(20, exact=True)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


EXPECTED_SHELL = (
    (1, 2, 17), (1, 5, 14), (1, 8, 11), (2, 5, 13),
    (2, 7, 11), (4, 5, 11), (5, 7, 8),
)

EXPECTED = {
    (1, 2, 17): dict(all=1493, minimal=1476, rays=1477,
        i0=Fraction(9, 833), i3=Fraction(9, 833), bulk=Fraction(18, 833),
        maximum=Fraction(36, 1295), winners=((1, 11, 185),),
        masses=(Fraction(12, 1295),) * 3, threshold=Fraction(22015, 53),
        proof_height=415, proof_all=1029, proof_minimal=1012,
        proof_rays=1013, proof_positive=1012, proof_components=21652,
        proof_zeros=()),
    (1, 5, 14): dict(all=1818, minimal=1800, rays=1801,
        i0=Fraction(9, 686), i3=Fraction(9, 1372), bulk=Fraction(6, 343),
        maximum=Fraction(24, 1043), winners=((1, 11, 149),),
        masses=(Fraction(6, 1043), Fraction(12, 1043), Fraction(6, 1043)),
        threshold=Fraction(21903, 47), proof_height=466, proof_all=1567,
        proof_minimal=1549, proof_rays=1550, proof_positive=1549,
        proof_components=37274, proof_zeros=()),
    (1, 8, 11): dict(all=2310, minimal=2295, rays=2296,
        i0=Fraction(9, 539), i3=Fraction(45, 8624), bulk=Fraction(39, 2156),
        maximum=Fraction(412, 16093), winners=((11, 19, 121),),
        masses=(Fraction(92, 16093), Fraction(12, 847), Fraction(92, 16093)),
        threshold=Fraction(1158696, 3385), proof_height=342, proof_all=1071,
        proof_minimal=1056, proof_rays=1057, proof_positive=1056,
        proof_components=21946, proof_zeros=()),
    (2, 5, 13): dict(all=1956, minimal=1943, rays=1944,
        i0=Fraction(9, 637), i3=Fraction(18, 3185), bulk=Fraction(54, 3185),
        maximum=Fraction(166, 6853), winners=((7, 11, 89),),
        masses=(Fraction(50, 6853), Fraction(6, 623), Fraction(50, 6853)),
        threshold=Fraction(4009005, 11332), proof_height=353, proof_all=977,
        proof_minimal=964, proof_rays=965, proof_positive=963,
        proof_components=18460, proof_zeros=((1, 19, 41),)),
    (2, 7, 11): dict(all=2317, minimal=2304, rays=2304,
        i0=Fraction(9, 539), i3=Fraction(18, 3773), bulk=Fraction(6, 343),
        maximum=Fraction(608, 24563), winners=((11, 29, 121),),
        masses=(Fraction(130, 24563), Fraction(12, 847), Fraction(130, 24563)),
        threshold=Fraction(1547469, 4369), proof_height=354, proof_all=1160,
        proof_minimal=1147, proof_rays=1147, proof_positive=1147,
        proof_components=24512, proof_zeros=()),
    (4, 5, 11): dict(all=2312, minimal=2299, rays=2300,
        i0=Fraction(9, 539), i3=Fraction(81, 21560), bulk=Fraction(87, 5390),
        maximum=Fraction(894, 41503), winners=((11, 49, 121),),
        masses=(Fraction(153, 41503), Fraction(12, 847), Fraction(153, 41503)),
        threshold=Fraction(118580, 249), proof_height=476, proof_all=2092,
        proof_minimal=2079, proof_rays=2080, proof_positive=2078,
        proof_components=55166, proof_zeros=((1, 19, 41),)),
    (5, 7, 8): dict(all=2820, minimal=2809, rays=2809,
        i0=Fraction(279, 13720), i3=Fraction(81, 27440), bulk=Fraction(6, 343),
        maximum=Fraction(216, 8855), winners=((13, 23, 55),),
        masses=(Fraction(8, 1771), Fraction(136, 8855), Fraction(8, 1771)),
        threshold=Fraction(185955, 499), proof_height=372, proof_all=1549,
        proof_minimal=1538, proof_rays=1538, proof_positive=1538,
        proof_components=32424, proof_zeros=()),
}

EXPECTED_MULTIPLE_MEASURES = {
    (1, 19, 41): Fraction(0), (5, 13, 49): Fraction(32, 4459),
    (7, 11, 65): Fraction(6, 455), (7, 11, 97): Fraction(18, 679),
    (11, 29, 41): Fraction(138, 8323), (11, 29, 79): Fraction(360, 16037),
    (13, 23, 41): Fraction(38, 6601), (13, 23, 47): Fraction(102, 7567),
    (13, 23, 67): Fraction(282, 20033), (13, 31, 37): Fraction(142, 8029),
    (17, 23, 53): Fraction(2592, 145061), (17, 37, 53): Fraction(90, 6307),
    (25, 29, 43): Fraction(120, 8729),
}


def complete_multiple_rays():
    labelled = [(pattern, c) for pattern in SHELL for c in M.relation_vectors(pattern)]
    candidates = set()
    for (_, c), (_, d) in combinations(labelled, 2):
        raw = M.cross(c, d)
        if raw == (0, 0, 0):
            continue
        if not (all(x > 0 for x in raw) or all(x < 0 for x in raw)):
            continue
        w = tuple(sorted(M.primitive_direction(raw)))
        if any(x <= 0 or x % 2 == 0 or x % 3 == 0 for x in w):
            continue
        if len(set(w)) != 3 or M.gcd_many(w) != 1:
            continue
        if M.direct_relations(w, SHORT):
            continue
        relations = M.direct_relations(w, SHELL)
        if sum(len(cs) for cs in relations.values()) >= 2:
            candidates.add(w)
    return tuple(sorted(candidates))


def main():
    check(SHELL == EXPECTED_SHELL, "seven exact norm-twenty patterns")
    check(len(SHORT) == 26, "twenty-six earlier even patterns")
    short_union = set()
    for pattern in SHORT:
        short_union.update(M.generated_relation_map(pattern, HEIGHT))

    print("LRC14 MINIMAL TERNARY-UNIT L1=20 SHELL")
    print(f"height={HEIGHT} short_patterns={len(SHORT)} shell={SHELL}")
    global_rows = []
    all_by_pattern = {}
    minimal_by_pattern = {}
    physical_cache = {}
    for pattern in SHELL:
        all_rows = M.generated_relation_map(pattern, HEIGHT)
        rows = {w: cs for w, cs in all_rows.items() if w not in short_union}
        all_by_pattern[pattern] = all_rows
        minimal_by_pattern[pattern] = rows
        values = {}
        ray_count = 0
        for w, relations in rows.items():
            ray_count += len(relations)
            reference = None
            for relation in relations:
                components, _ = M.raw_relation_components(w, relation)
                value = sum(components.values(), Fraction(0))
                if reference is None:
                    reference = value
                elif value != reference:
                    raise RuntimeError(("presentation mismatch", pattern, w, reference, value))
            values[w] = reference
        maximum = max(values.values(), default=Fraction(0))
        winners = tuple(w for w, value in values.items() if value == maximum)
        i0 = M.slice_integral(pattern, 0)
        i3 = M.slice_integral(pattern, 3)
        bulk = Fraction(2, 3) * (i0 + 2 * i3)
        threshold = (Fraction(18, 7) / (maximum - bulk)
                     if maximum > bulk else None)
        closed = threshold is not None and threshold <= HEIGHT
        height = M.floor_fraction(threshold)
        proof = {w: cs for w, cs in rows.items() if max(w) <= height}
        proof_all = {w: cs for w, cs in all_rows.items() if max(w) <= height}
        proof_values = {}
        component_count = 0
        zero_triples = []
        for w, relations in proof.items():
            if w not in physical_cache:
                physical_cache[w] = M.literal_physical_components(w)
            literal, representatives = physical_cache[w]
            value = sum(literal.values(), Fraction(0))
            proof_values[w] = value
            component_count += len(literal)
            if not value:
                zero_triples.append(w)
            for relation in relations:
                raw, metadata = M.raw_relation_components(w, relation)
                if raw != literal:
                    raise RuntimeError(("raw/literal mismatch", pattern, w, relation))
                for C, n in representatives.items():
                    if metadata[C][0] != M.dot(relation, n):
                        raise RuntimeError(("defect mismatch", pattern, w, relation, C))
        if max(proof_values.values()) != maximum:
            raise RuntimeError(("proof maximum mismatch", pattern))
        winner = winners[0]
        relation = next(iter(proof[winner]))
        raw, metadata = M.raw_relation_components(winner, relation)
        masses = tuple(sum((raw[C] for C in raw if metadata[C][0] == delta), Fraction(0))
                       for delta in (-3, 0, 3))
        global_rows.append((maximum, pattern, winners))
        actual = dict(all=len(all_rows), minimal=len(rows), rays=ray_count,
            i0=i0, i3=i3, bulk=bulk, maximum=maximum, winners=winners,
            masses=masses, threshold=threshold, proof_height=height,
            proof_all=len(proof_all), proof_minimal=len(proof),
            proof_rays=sum(len(cs) for cs in proof.values()),
            proof_positive=sum(v > 0 for v in proof_values.values()),
            proof_components=component_count, proof_zeros=tuple(zero_triples))
        check(actual == EXPECTED[pattern], f"frozen sector row {pattern}")
        check(bulk + Fraction(18, 7 * (height + 1)) < maximum,
              f"strict analytic tail beyond proof window {pattern}")
        print(
            f"row={pattern} all={len(all_rows)} minimal={len(rows)} rays={ray_count} "
            f"I0={i0} I3={i3} bulk={bulk} maximum={maximum} "
            f"winners={winners} masses={masses} threshold={threshold} closed_at_H={closed} "
            f"proof_height={height} proof_all={len(proof_all)} proof_minimal={len(proof)} "
            f"proof_rays={sum(len(cs) for cs in proof.values())} "
            f"proof_positive={sum(v > 0 for v in proof_values.values())} "
            f"proof_components={component_count} proof_zeros={tuple(zero_triples)}"
        )
    physical_union = set().union(*(set(rows) for rows in minimal_by_pattern.values()))
    incidences = sum(len(rows) for rows in minimal_by_pattern.values())
    rays = sum(sum(len(cs) for cs in rows.values()) for rows in minimal_by_pattern.values())
    cross = {}
    within = {}
    positive = 0
    for w in physical_union:
        labels = tuple(pattern for pattern in SHELL if w in minimal_by_pattern[pattern])
        if len(labels) > 1:
            cross[w] = labels
        relations = M.direct_relations(w, SHELL)
        if sum(len(cs) for cs in relations.values()) > len(labels):
            within[w] = relations
        first = next(c for cs in relations.values() for c in cs)
        raw, _ = M.raw_relation_components(w, first)
        positive += bool(raw)
    complete = complete_multiple_rays()
    observed_multiple = tuple(sorted(set(cross) | set(within)))
    if complete != observed_multiple:
        raise RuntimeError(("multiple-ray completeness mismatch", complete, observed_multiple))
    check((incidences, rays, len(physical_union), positive) ==
          (14926, 14931, 14918, 14917), "height-500 incidence ledger")
    check(len(cross) == 8 and len(within) == 5, "8 cross plus 5 within overlaps")
    measured = {}
    print(f"H={HEIGHT} incidences={incidences} rays={rays} triples={len(physical_union)} positive={positive}")
    print(f"cross={cross}")
    print(f"within={within}")
    for w in complete:
        relations = M.direct_relations(w, SHELL)
        first = next(c for cs in relations.values() for c in cs)
        raw, _ = M.raw_relation_components(w, first)
        measure = sum(raw.values(), Fraction(0))
        measured[w] = measure
        print(f"multiple={w} relations={relations} measure={measure}")
    check(measured == EXPECTED_MULTIPLE_MEASURES, "complete multiple-ray measures")
    global_row = max(global_rows)
    check(global_row == (Fraction(36, 1295), (1, 2, 17), ((1, 11, 185),)),
          "unique global maximum")
    semantic = json.dumps({
        "shell": EXPECTED_SHELL,
        "rows": {str(k): {q: str(v) for q, v in row.items()}
                 for k, row in EXPECTED.items()},
        "multiple_measures": {str(k): str(v) for k, v in measured.items()},
        "global": str(global_row),
    }, sort_keys=True, separators=(",", ":"))
    print("global=" + repr(global_row))
    print(f"producer_checks={CHECKS}; inherited_live_checks={M.CHECKS}")
    print("semantic_sha256=" + sha256(semantic.encode()).hexdigest())
    print("verdict=PASS")


if __name__ == "__main__":
    main()

