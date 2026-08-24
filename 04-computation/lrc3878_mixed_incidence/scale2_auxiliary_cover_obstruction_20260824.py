#!/usr/bin/env python3
"""Exact auxiliary-comb covering obstruction for the (2,1,9) LRC14 seam."""

from fractions import Fraction as F
from hashlib import sha256
import json
import sys


sys.stdout.reconfigure(newline="\n")
I = (F(2, 21), F(8, 63))
C = (I, (F(55, 63), F(19, 21)))
CHECKS = 0
EXPECTED_SEMANTIC_SHA256 = "cb06eca35b81813fe2cb54aa96d1d160e44d827cd00cd5c7049cece194fb1d6a"


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def danger(a):
    r = F(1, 14 * a)
    out = []
    for k in range(a + 1):
        lo, hi = max(F(0), F(k, a) - r), min(F(1), F(k, a) + r)
        if lo < hi and any(lo < y and x < hi for x, y in C):
            out.append((lo, hi))
    return tuple(out)


def covers(targets, pieces):
    pieces = sorted(pieces)
    for left, right in targets:
        pos = left
        started = False
        for lo, hi in pieces:
            if hi <= pos:
                continue
            if not started:
                if lo > left:
                    return False
                started = True
                pos = hi
            else:
                # Open teeth that only touch leave their common point uncovered.
                if lo >= pos and pos < right:
                    return False
                pos = max(pos, hi)
            if pos >= right:
                break
        if not started or pos < right:
            return False
    return True


def complement_lengths(target, pieces):
    left, right = target
    clipped = sorted((max(left, a), min(right, b)) for a, b in pieces
                     if a < right and left < b)
    pos = left
    gaps = []
    for lo, hi in clipped:
        if lo > pos:
            gaps.append(lo - pos)
        pos = max(pos, hi)
    if pos < right:
        gaps.append(right - pos)
    return tuple(gaps)


def fmt(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def main():
    require(I[1] - I[0] == F(2, 63), "positive obstruction width")
    require(C[1] == (1-I[1], 1-I[0]), "reflection symmetry")
    ds = {a: danger(a) for a in range(1, 31)}

    # If a<=b and D_a U D_b covers I, tooth counting gives
    # |I| <= 2|I|/7+4/(7a), hence a<=floor(4/(5|I|))=25.
    require(F(4, 5 * (I[1]-I[0])) == F(126, 5), "measure cutoff exact")
    cutoff = 25
    table = []
    candidates = []
    for a in range(1, cutoff + 1):
        gaps = complement_lengths(I, ds[a])
        require(bool(gaps), f"single comb noncover a={a}")
        largest = max(gaps)
        # A connected D_a-complement gap must fit in one open b-tooth,
        # whose full length is 1/(7b).  Non-strict bmax over-includes equality.
        bmax = int(F(1, 7) / largest)
        table.append((a, fmt(largest), bmax))
        for b in range(a, bmax + 1):
            candidates.append((a, b))
            require(not covers(C, ds[a] + ds[b]), f"two-comb hostile {(a,b)}")

    require(len(candidates) == 20, "finite candidate count")
    require(not any(covers(C, ds[a] + ds[b]) for a, b in candidates),
            "no finite pair cover")

    triple = (8, 9, 10)
    require(covers(C, sum((ds[a] for a in triple), ())), "8910 positive control")
    expected_teeth = (
        (F(13, 140), F(3, 28)),
        (F(13, 126), F(5, 42)),
        (F(13, 112), F(15, 112)),
    )
    actual_teeth = tuple(next((x for x in ds[a] if x[0] < I[1] and I[0] < x[1]))
                         for a in (10, 9, 8))
    require(actual_teeth == expected_teeth, "explicit positive-arc teeth")
    require(expected_teeth[0][1] > expected_teeth[1][0], "10/9 strict overlap")
    require(expected_teeth[1][1] > expected_teeth[2][0], "9/8 strict overlap")
    require(expected_teeth[0][0] < I[0] and I[1] < expected_teeth[2][1],
            "target endpoints strictly covered")

    triples30 = []
    for c in range(1, 31):
        for a in range(1, c + 1):
            for b in range(a, c + 1):
                if covers(C, ds[a] + ds[b] + ds[c]):
                    triples30.append((a, b, c))
    require(triples30 == [triple], "unique nondecreasing triple max<=30")

    semantic = {
        "C": [[fmt(x), fmt(y)] for x, y in C],
        "pair_measure_cutoff": cutoff,
        "gap_table": table,
        "candidate_pairs": candidates,
        "covering_pairs": [],
        "triples_max_30": triples30,
        "positive_teeth": [[fmt(x), fmt(y)] for x, y in expected_teeth],
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("LRC14_SCALE2_AUXILIARY_COVER_OBSTRUCTION_20260824")
    print("scope=C_(1,9)=(2/21,8/63)U(55/63,19/21);D_a={||a*w||<1/14}")
    print("two_comb_measure_reduction=a<=b_implies_a<=25")
    print(f"finite_candidate_pairs={len(candidates)};covering_pairs=0")
    print("three_comb_cover=(8,9,10);unique_nondecreasing_cover_with_max<=30")
    print("positive_arc_teeth=(13/140,3/28)U(13/126,5/42)U(13/112,15/112)")
    print("LRC_meaning=settled_LRCUpTo13_allows_only_one_auxiliary_clock_beyond_an_11_speed_body")
    print("stopping_obstruction=three_clocks_are_needed;using_two_is_already_LRC14_circular")
    print("covector_join=reusing_8t_9t_or_10t_as_a_body_coordinate_creates_a_forbidden_small_crossing_relation_in_W=V_dec")
    print(f"semantic_sha256={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
