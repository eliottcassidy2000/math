#!/usr/bin/env python3
"""Clean-room content-one overlap audit for THM-4160 versus THM-4158.

This script imports no repository computation.  It uses two independent
counting routes:

1. literal enumeration of each THM-4160 core body, with THM-4158 carrier
   membership represented by a bit mask derived directly from the carrier
   endpoint inequalities;
2. inclusion-exclusion over the maximal optional-label sets of the displayed
   THM-4158 bands, plus a separate subset min/max dynamic program.

All overlap counts concern undilated/content-one eleven-label core bodies.
The final hostile demonstrates why they must not be promoted to an all-content
classification.
"""

from collections import Counter, defaultdict
from itertools import combinations
from math import comb
import sys


sys.stdout.reconfigure(newline="\n")


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(label for label in POOL if label not in ANCHORS)
NEWCOMERS = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


require(len(OPTIONAL) == 27, "optional pool size")
require(len(set(ANCHORS)) == 3, "distinct anchors")
require(len(set(OPTIONAL)) == 27, "distinct optional labels")
require(len(set(NEWCOMERS)) == 10, "distinct newcomer labels")
require(set(ANCHORS).isdisjoint(OPTIONAL), "anchors/optional disjoint")
require(set(NEWCOMERS).isdisjoint(POOL), "newcomers/base pool disjoint")


def endpoint_safe(m, h):
    """Whether h is safe on THM-4158's closed carrier for parameter m.

    This uses endpoint containment in a safe tooth rather than the displayed
    three/four-band formula.  A tooth index j is valid exactly when

        m(14j+1) <= h,
        16h <= (12m+1)(14j+13).
    """

    for j in range(h // (14 * m) + 1):
        if (m * (14 * j + 1) <= h
                and 16 * h <= (12 * m + 1) * (14 * j + 13)):
            return True
    return False


ALL_LABELS = set(ANCHORS + OPTIONAL + NEWCOMERS)


def content_one_masks():
    # Every core body contains 120.  Every label in P_m is at least m, so no
    # content-one body in this family can use a carrier parameter m>120.
    return {
        h: sum(1 << (m - 1) for m in range(1, 121) if endpoint_safe(m, h))
        for h in ALL_LABELS
    }


MASK = content_one_masks()
ANCHOR_MASK = (1 << 120) - 1
for anchor in ANCHORS:
    ANCHOR_MASK &= MASK[anchor]
require(ANCHOR_MASK != 0, "nonempty anchor carrier set")


direct = Counter()
per_newcomer = defaultdict(Counter)
carrier_multiplicity = Counter()
thm4148_gate_hits = 0


def direct_classify(variable_labels, family, newcomer=None):
    global thm4148_gate_hits

    carrier_mask = ANCHOR_MASK
    minimum = min(ANCHORS)
    maximum = max(ANCHORS)
    for label in variable_labels:
        carrier_mask &= MASK[label]
        minimum = min(minimum, label)
        maximum = max(maximum, label)

    defect = 16 * maximum - 156 * minimum
    gate = "stable" if defect <= 0 else ("affine" if defect <= 13 else "beyond")
    covered = bool(carrier_mask)
    direct[(family, gate, covered)] += 1
    if newcomer is not None:
        per_newcomer[newcomer][(gate, covered)] += 1
    if covered:
        carrier_multiplicity[carrier_mask.bit_count()] += 1

    # THM-4148's named first-window gate is
    # 27(13*minimum-maximum) >= 4*minimum*maximum.
    thm4148_gate_hits += (
        27 * (13 * minimum - maximum) >= 4 * minimum * maximum
    )


for choice in combinations(OPTIONAL, 8):
    direct_classify(choice, "old")
for newcomer in NEWCOMERS:
    for choice in combinations(OPTIONAL, 7):
        direct_classify((newcomer,) + choice, "new", newcomer)


def sum_direct(family=None, gate=None, covered=None):
    total = 0
    for (row_family, row_gate, row_covered), count in direct.items():
        if family is not None and family != row_family:
            continue
        if gate is not None and gate != row_gate:
            continue
        if covered is not None and covered != row_covered:
            continue
        total += count
    return total


OLD_TOTAL = comb(27, 8)
NEW_PER_LABEL = comb(27, 7)
NEW_TOTAL = len(NEWCOMERS) * NEW_PER_LABEL
TOTAL = OLD_TOTAL + NEW_TOTAL

require(sum_direct("old") == OLD_TOTAL, "old family enumeration")
require(sum_direct("new") == NEW_TOTAL, "new family enumeration")
require(sum_direct(gate="affine") == 0, "empty affine-defect stratum")
require(sum_direct(gate="stable", covered=False) == 0,
        "every content-stable body is in a wrapped carrier")
require(thm4148_gate_hits == 0, "all THM-4160 cores fail THM-4148 gate")


# Independent route: use the displayed carrier bands, remove dominated
# optional-label sets, and count their union by inclusion-exclusion.
def band_contains(m, h):
    if m == 1:
        bands = ((1, 10), (15, 21), (29, 33), (43, 44))
    else:
        s = 12 * m + 1
        bands = (
            (m, 13 * s // 16),
            (15 * m, 27 * s // 16),
            (29 * m, 41 * s // 16),
        )
    return any(lower <= h <= upper for lower, upper in bands)


ANCHOR_PARAMETERS = tuple(
    m for m in range(1, 121)
    if all(band_contains(m, anchor) for anchor in ANCHORS)
)


def maximal_optional_sets(newcomer=None):
    rows = []
    for m in ANCHOR_PARAMETERS:
        if newcomer is not None and not band_contains(m, newcomer):
            continue
        row = frozenset(label for label in OPTIONAL if band_contains(m, label))
        if row not in rows:
            rows.append(row)
    return tuple(row for row in rows if not any(row < other for other in rows))


def union_k_subsets(rows, k):
    total = 0
    for size in range(1, len(rows) + 1):
        sign = 1 if size % 2 else -1
        for selected in combinations(rows, size):
            intersection = set(selected[0])
            for row in selected[1:]:
                intersection.intersection_update(row)
            if len(intersection) >= k:
                total += sign * comb(len(intersection), k)
    return total


old_overlap_ie = union_k_subsets(maximal_optional_sets(), 8)
new_overlap_ie = {
    newcomer: union_k_subsets(maximal_optional_sets(newcomer), 7)
    for newcomer in NEWCOMERS
}


# A second independent min/max count does not inspect carrier masks at all.
def optional_minmax_profile(k):
    states = {(0, None, None): 1}
    for label in OPTIONAL:
        updated = dict(states)
        for (chosen, minimum, maximum), count in states.items():
            if chosen == k:
                continue
            key = (
                chosen + 1,
                label if minimum is None else min(minimum, label),
                label if maximum is None else max(maximum, label),
            )
            updated[key] = updated.get(key, 0) + count
        states = updated
    return {
        (minimum, maximum): count
        for (chosen, minimum, maximum), count in states.items()
        if chosen == k
    }


PROFILE_7 = optional_minmax_profile(7)
PROFILE_8 = optional_minmax_profile(8)


def gate_profile(profile, extra=()):
    base = ANCHORS + tuple(extra)
    base_minimum = min(base)
    base_maximum = max(base)
    result = Counter()
    for (minimum, maximum), count in profile.items():
        body_minimum = min(base_minimum, minimum)
        body_maximum = max(base_maximum, maximum)
        defect = 16 * body_maximum - 156 * body_minimum
        gate = "stable" if defect <= 0 else ("affine" if defect <= 13 else "beyond")
        result[gate] += count
    return result


old_gate_dp = gate_profile(PROFILE_8)
new_gate_dp = {q: gate_profile(PROFILE_7, (q,)) for q in NEWCOMERS}
all_gate_dp = old_gate_dp.copy()
for row in new_gate_dp.values():
    all_gate_dp.update(row)


OLD_OVERLAP = sum_direct("old", covered=True)
NEW_OVERLAP = sum_direct("new", covered=True)
OVERLAP = OLD_OVERLAP + NEW_OVERLAP
OUTSIDE = TOTAL - OVERLAP
STABLE_OVERLAP = sum_direct(gate="stable", covered=True)
BEYOND_OVERLAP = sum_direct(gate="beyond", covered=True)
NEW_OUTSIDE = sum_direct("new", covered=False)

require(old_overlap_ie == OLD_OVERLAP, "old inclusion-exclusion agreement")
require(sum(new_overlap_ie.values()) == NEW_OVERLAP,
        "new inclusion-exclusion agreement")
require(old_gate_dp == Counter({"stable": 344366, "beyond": 1875709}),
        "old gate DP")
require(all_gate_dp == Counter({"stable": 1397101, "beyond": 9703274}),
        "combined gate DP")
require(STABLE_OVERLAP == all_gate_dp["stable"],
        "all stable bodies overlap by the first carrier band")
require(OVERLAP == STABLE_OVERLAP + BEYOND_OVERLAP,
        "overlap gate partition")


# Scope hostile: carrier membership is not invariant under common dilation.
# This core lies outside every P_m at content one, but 2H lies in P_35.
SCALE_HOSTILE = (20, 30, 40, 42, 60, 63, 120, 126, 143, 168, 264)


def scaled_body_carriers(body, content):
    limit = content * min(body)
    result = 0
    for m in range(1, limit + 1):
        if all(endpoint_safe(m, content * label) for label in body):
            result |= 1 << (m - 1)
    return result


hostile_content_one = scaled_body_carriers(SCALE_HOSTILE, 1)
hostile_content_two = scaled_body_carriers(SCALE_HOSTILE, 2)
require(hostile_content_one == 0, "scale hostile outside at content one")
require(hostile_content_two != 0, "scale hostile enters after dilation")
require((hostile_content_two & -hostile_content_two).bit_length() == 35,
        "scale hostile least content-two carrier")


EXPECTED = {
    "overlap": 2_063_957,
    "outside": 9_036_418,
    "new_outside": 7_267_924,
    "stable_overlap": 1_397_101,
    "beyond_overlap": 666_856,
}
require(OVERLAP == EXPECTED["overlap"], "reported overlap")
require(OUTSIDE == EXPECTED["outside"], "reported outside count")
require(NEW_OUTSIDE == EXPECTED["new_outside"], "reported new outside count")
require(STABLE_OVERLAP == EXPECTED["stable_overlap"], "reported stable overlap")
require(BEYOND_OVERLAP == EXPECTED["beyond_overlap"], "reported beyond overlap")


print("THM4160_VS_THM4158_CONTENT_ONE_OVERLAP_INDEPENDENT_AUDIT")
print("status=PASS;scope=undilated_content_one_core_body_labels_only")
print("collision_model=old_has_zero_newcomers;new_has_exactly_one_disjoint_newcomer;representations_injective")
print(f"family_sizes=old:{OLD_TOTAL},new:{NEW_TOTAL},total:{TOTAL}")
print(f"carrier_parameter_bound=1..120;anchor_compatible={len(ANCHOR_PARAMETERS)}")
print("direct_endpoint_mask=one_body_one_count;multiple_carrier_parameters_ORed_before_counting")
print(f"old_split=stable_overlap:{sum_direct('old','stable',True)},beyond_overlap:{sum_direct('old','beyond',True)},beyond_outside:{sum_direct('old','beyond',False)}")
print(f"new_split=stable_overlap:{sum_direct('new','stable',True)},beyond_overlap:{sum_direct('new','beyond',True)},beyond_outside:{sum_direct('new','beyond',False)}")
print(f"total_overlap={OVERLAP}")
print(f"outside_all_THM4158_Pm_at_content_one={OUTSIDE}")
print(f"singleton_newcomers_outside_all_THM4158_Pm_at_content_one={NEW_OUTSIDE}")
print(f"overlap_gate_split=stable:{STABLE_OVERLAP},beyond:{BEYOND_OVERLAP},affine:0")
print(f"inclusion_exclusion=old:{old_overlap_ie},new:{sum(new_overlap_ie.values())},total:{old_overlap_ie + sum(new_overlap_ie.values())}")
print("newcomer_overlap=" + ",".join(f"{q}:{new_overlap_ie[q]}" for q in NEWCOMERS))
print(f"gate_DP=stable:{all_gate_dp['stable']},beyond:{all_gate_dp['beyond']},affine:{all_gate_dp['affine']}")
print(f"THM4148_named_gate_hits={thm4148_gate_hits}")
print(f"carrier_multiplicity=min:{min(carrier_multiplicity)},max:{max(carrier_multiplicity)},multi_body_count:{sum(v for k,v in carrier_multiplicity.items() if k>1)}")
print("scale_scope_hostile=(20,30,40,42,60,63,120,126,143,168,264);content1:none;content2:least_m35")
print("promotion_scope=FINITE_EXACT_content_one_overlap_only;no_all_content_overlap_claim")
print("result=PASS")
