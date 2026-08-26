#!/usr/bin/env python3
"""Exact overlap census for the THM-4160 bodies and THM-4158 carriers.

This is a self-contained audit. It imports neither theorem computation and uses
only integer inequalities and itertools.combinations.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError("REQUIRE FAILED: " + message)


A = (120, 126, 143)
P = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
N = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)
O = tuple(h for h in P if h not in A)
MAX_RELEVANT_M = 120


def in_pm_band_formula(h, m):
    """Membership in the exact THM-4158 alphabet P_m, using its bands."""
    require(h > 0 and m > 0, "carrier inputs must be positive")
    if m == 1:
        return (
            1 <= h <= 10
            or 15 <= h <= 21
            or 29 <= h <= 33
            or 43 <= h <= 44
        )
    return (
        (m <= h and 16 * h <= 13 * (12 * m + 1))
        or (15 * m <= h and 16 * h <= 27 * (12 * m + 1))
        or (29 * m <= h and 16 * h <= 41 * (12 * m + 1))
    )


def in_pm_endpoint_rows(h, m):
    """Independent spelling of the endpoint inequalities in THM-4158(5)."""
    require(h > 0 and m > 0, "endpoint-row inputs must be positive")
    last_k = 3 if m == 1 else 2
    for k in range(last_k + 1):
        lower = m * (14 * k + 1)
        upper_numerator = (12 * m + 1) * (14 * k + 13)
        if lower <= h and 16 * h <= upper_numerator:
            return True
    return False


LABELS = tuple(sorted(set(P) | set(N)))
for h in LABELS:
    for m in range(1, MAX_RELEVANT_M + 1):
        require(
            in_pm_band_formula(h, m) == in_pm_endpoint_rows(h, m),
            "band/endpoint carrier mismatch at h=%d,m=%d" % (h, m),
        )


def carrier_mask(h):
    mask = 0
    for m in range(1, MAX_RELEVANT_M + 1):
        if in_pm_band_formula(h, m):
            mask |= 1 << (m - 1)
    return mask


LABEL_MASK = {h: carrier_mask(h) for h in LABELS}
ANCHOR_MASK = LABEL_MASK[A[0]] & LABEL_MASK[A[1]] & LABEL_MASK[A[2]]


require(len(A) == 3 and len(P) == 30 and len(O) == 27 and len(N) == 10,
        "THM-4160 source-set cardinality changed")
require(A == tuple(sorted(A)) and P == tuple(sorted(P))
        and O == tuple(sorted(O)) and N == tuple(sorted(N)),
        "source tuples must be strictly sorted for endpoint extraction")
require(set(A).issubset(P), "anchors are not contained in P")
require(set(O) == set(P) - set(A), "O is not exactly P minus A")
require(set(P).isdisjoint(N), "newcomers collide with P")
require(len(set(P)) == len(P) and len(set(N)) == len(N),
        "source labels are not distinct")
require(comb(len(O), 8) == 2_220_075, "old universe count changed")
require(len(N) * comb(len(O), 7) == 8_880_300,
        "new universe count changed")


CLASS_ORDER = ("old",) + N
EXPECTED_TOTAL = {"old": 2_220_075, **{q: 888_030 for q in N}}
EXPECTED_STABLE = {
    "old": 344_366,
    5: 0,
    66: 187_639,
    182: 182_920,
    298: 116_280,
    336: 116_280,
    340: 116_280,
    380: 116_280,
    386: 116_280,
    528: 50_388,
    572: 50_388,
}
EXPECTED_CANONICAL = {
    "old": 422_222,
    5: 0,
    66: 229_540,
    182: 182_920,
    298: 161_976,
    336: 128_656,
    340: 128_656,
    380: 128_656,
    386: 128_656,
    528: 104_652,
    572: 104_652,
}
EXPECTED_GLOBAL = {
    "old": 451_581,
    5: 0,
    66: 246_746,
    182: 190_291,
    298: 198_025,
    336: 152_329,
    340: 152_329,
    380: 153_408,
    386: 153_408,
    528: 182_920,
    572: 182_920,
}


def blank_stats():
    return {
        "total": 0,
        "stable": 0,
        "transition": 0,
        "beyond": 0,
        "canonical": 0,
        "canonical_stable": 0,
        "canonical_transition": 0,
        "canonical_beyond": 0,
        "global": 0,
        "global_stable": 0,
        "global_transition": 0,
        "global_beyond": 0,
        "outside_global_stable": 0,
        "outside_global_transition": 0,
        "outside_global_beyond": 0,
        "witness_incidence": 0,
        "multi_carrier": 0,
        "witness_histogram": Counter(),
    }


STATS = {key: blank_stats() for key in CLASS_ORDER}


def record(key, combo, newcomer):
    stats = STATS[key]
    body_min = min(A[0], combo[0], newcomer if newcomer is not None else A[0])
    body_max = max(A[-1], combo[-1], newcomer if newcomer is not None else A[-1])
    defect = 16 * body_max - 156 * body_min
    if defect <= 0:
        split = "stable"
    elif defect <= 13:
        split = "transition"
    else:
        split = "beyond"

    mask = ANCHOR_MASK
    if newcomer is not None:
        mask &= LABEL_MASK[newcomer]
    for h in combo:
        mask &= LABEL_MASK[h]
        if mask == 0:
            break

    canonical_bit = 1 << (body_min - 1)
    canonical = (mask & canonical_bit) != 0
    global_member = mask != 0
    require(not canonical or global_member,
            "canonical membership is not contained in global membership")

    witnesses = mask.bit_count()
    stats["total"] += 1
    stats[split] += 1
    if canonical:
        stats["canonical"] += 1
        stats["canonical_" + split] += 1
    if global_member:
        stats["global"] += 1
        stats["global_" + split] += 1
        stats["witness_incidence"] += witnesses
        stats["witness_histogram"][witnesses] += 1
        if witnesses > 1:
            stats["multi_carrier"] += 1
    else:
        stats["outside_global_" + split] += 1


for combo in combinations(O, 8):
    record("old", combo, None)

for q in N:
    for combo in combinations(O, 7):
        record(q, combo, q)


def summed(field):
    return sum(STATS[key][field] for key in CLASS_ORDER)


for key in CLASS_ORDER:
    stats = STATS[key]
    require(stats["total"] == EXPECTED_TOTAL[key],
            "class total changed for %s" % key)
    require(stats["stable"] == EXPECTED_STABLE[key],
            "stable count changed for %s" % key)
    require(stats["transition"] == 0,
            "a defect in [1,13] appeared for %s" % key)
    require(stats["canonical"] == EXPECTED_CANONICAL[key],
            "canonical overlap changed for %s" % key)
    require(stats["global"] == EXPECTED_GLOBAL[key],
            "global overlap changed for %s" % key)
    require(stats["canonical_stable"] == stats["stable"],
            "a stable body missed its minimum-indexed carrier for %s" % key)
    require(stats["global_stable"] == stats["stable"],
            "a stable body missed the global carrier union for %s" % key)
    require(
        stats["global"]
        == stats["global_stable"]
        + stats["global_transition"]
        + stats["global_beyond"],
        "global split does not conserve class %s" % key,
    )
    require(
        stats["total"]
        == stats["global"]
        + stats["outside_global_stable"]
        + stats["outside_global_transition"]
        + stats["outside_global_beyond"],
        "inside/outside split does not conserve class %s" % key,
    )

require(summed("total") == 11_100_375, "total body universe changed")
require(summed("stable") == 1_397_101, "total stable count changed")
require(summed("transition") == 0, "total transition count changed")
require(summed("beyond") == 9_703_274, "total beyond count changed")
require(summed("canonical") == 1_720_586,
        "minimum-indexed overlap total changed")
require(summed("global") == 2_063_957, "global overlap total changed")
require(summed("global_stable") == 1_397_101,
        "global stable overlap changed")
require(summed("global_beyond") == 666_856,
        "global beyond overlap changed")
require(summed("outside_global_stable") == 0,
        "a stable body lies outside all THM-4158 carriers")
require(summed("outside_global_beyond") == 9_036_418,
        "global-new beyond count changed")
require(sum(STATS[q]["global"] for q in N) == 1_612_376,
        "new-family global overlap changed")
require(sum(STATS[q]["total"] - STATS[q]["global"] for q in N) == 7_267_924,
        "new-family outside-global count changed")


all_witness_histogram = Counter()
for key in CLASS_ORDER:
    all_witness_histogram.update(STATS[key]["witness_histogram"])
require(sum(all_witness_histogram.values()) == summed("global"),
        "witness histogram does not count global bodies exactly once")
require(
    sum(mult * count for mult, count in all_witness_histogram.items())
    == summed("witness_incidence"),
    "witness incidence does not match its multiplicity histogram",
)


lines = []
lines.append("THM4160_VS_THM4158_EXACT_OVERLAP_AUDIT")
lines.append("status=PASS")
lines.append(
    "u4160=old:C(27,8);new:10*C(27,7);total=11100375;"
    "all bodies are sorted 11-element sets"
)
lines.append(
    "thm4158_relevant_carriers=m=1..120;P1=four exact bands;"
    "Pm(m>=2)=three exact bands"
)
lines.append(
    "global_membership=exists m in [1,min(H)] with H subset P_m;"
    "canonical_membership=H subset P_min(H)"
)
lines.append(
    "collision_rule=each H counted once by Boolean carrier-mask union;"
    "multiple witnessing m values are not added as bodies"
)
lines.append(
    "split=stable:D<=0;transition:1<=D<=13;beyond:D>=14;"
    "D=16*max(H)-156*min(H)"
)
lines.append(
    "class,total,stable,transition,beyond,canonical,canonical_beyond,"
    "global,global_beyond,outside_global_beyond,multi_carrier,witness_incidence"
)
for key in CLASS_ORDER:
    stats = STATS[key]
    lines.append(
        "%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"
        % (
            key,
            stats["total"],
            stats["stable"],
            stats["transition"],
            stats["beyond"],
            stats["canonical"],
            stats["canonical_beyond"],
            stats["global"],
            stats["global_beyond"],
            stats["outside_global_beyond"],
            stats["multi_carrier"],
            stats["witness_incidence"],
        )
    )
lines.append(
    "all,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"
    % (
        summed("total"),
        summed("stable"),
        summed("transition"),
        summed("beyond"),
        summed("canonical"),
        summed("canonical_beyond"),
        summed("global"),
        summed("global_beyond"),
        summed("outside_global_beyond"),
        summed("multi_carrier"),
        summed("witness_incidence"),
    )
)
lines.append(
    "global_witness_multiplicity="
    + ",".join(
        "%d:%d" % (multiplicity, all_witness_histogram[multiplicity])
        for multiplicity in sorted(all_witness_histogram)
    )
)
lines.append(
    "new_only=total:8880300;canonical:%d;global:%d;outside_global:%d"
    % (
        sum(STATS[q]["canonical"] for q in N),
        sum(STATS[q]["global"] for q in N),
        sum(STATS[q]["total"] - STATS[q]["global"] for q in N),
    )
)
lines.append(
    "stable_coverage=all_1397101_stable_bodies_are_canonical_and_global"
)
lines.append(
    "strictly_outside_global=9036418;all_are_in_the_9703274_beyond_bucket"
)
semantic_payload = "\n".join(lines) + "\n"
lines.append("semantic_sha256=" + sha256(semantic_payload.encode("utf-8")).hexdigest())
sys.stdout.write("\n".join(lines) + "\n")
