#!/usr/bin/env python3
"""Grouped exact THM-4166 versus THM-4158 content-one overlap census.

The 1,032 THM-4166 qualifier labels are read from its frozen exact global
census. For each label, the union of all compatible THM-4158 carriers is
counted by inclusion-exclusion over at most five maximal optional-label sets.
Min/max subset counts are evaluated symbolically, so no 916-million-body scan
is performed.
"""

import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import re
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError("REQUIRE FAILED: " + message)


ROOT = Path(__file__).resolve().parents[1]
QUALIFIER_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_two_deletion_repair_graph_thm4166_global_census.out"
)
QUALIFIER_OUTPUT_SHA256 = (
    "baae4adc57739c852dec2a8c666328cea2c4031934547e2afa4044a314c9bff6"
)

ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(label for label in POOL if label not in ANCHORS)
ONE_DELETION = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)


raw_qualifier_output = QUALIFIER_OUTPUT.read_bytes()
require(
    sha256(raw_qualifier_output).hexdigest() == QUALIFIER_OUTPUT_SHA256,
    "THM-4166 qualifier source hash changed",
)
qualifier_text = raw_qualifier_output.decode("ascii")
match = re.search(r"^admitted_q (.*)$", qualifier_text, re.MULTILINE)
require(match is not None, "missing admitted_q row")
QUALIFIERS = tuple(ast.literal_eval(match.group(1)))

require(len(ANCHORS) == 3, "anchor count")
require(len(POOL) == 30 and len(OPTIONAL) == 27, "pool count")
require(len(QUALIFIERS) == 1_032 and QUALIFIERS[-1] == 8_265,
        "qualifier count/boundary")
require(tuple(sorted(set(QUALIFIERS))) == QUALIFIERS,
        "qualifiers sorted and unique")
require(set(QUALIFIERS).isdisjoint(POOL), "qualifiers collide with pool")
require(set(ONE_DELETION).issubset(QUALIFIERS),
        "THM-4160 newcomers not subsumed")
require(7 not in POOL and 7 not in QUALIFIERS, "named m=7 separator")


def band_contains(m, label):
    """Membership in the exact displayed THM-4158 alphabet P_m."""
    if m == 1:
        bands = ((1, 10), (15, 21), (29, 33), (43, 44))
    else:
        scale = 12 * m + 1
        bands = (
            (m, 13 * scale // 16),
            (15 * m, 27 * scale // 16),
            (29 * m, 41 * scale // 16),
        )
    return any(lower <= label <= upper for lower, upper in bands)


def endpoint_contains(m, label):
    """Independent carrier endpoint inequalities, including m=1's row 3."""
    last_tooth = 3 if m == 1 else 2
    return any(
        m * (14 * tooth + 1) <= label
        and 16 * label <= (12 * m + 1) * (14 * tooth + 13)
        for tooth in range(last_tooth + 1)
    )


relevant_labels = tuple(sorted(set(POOL) | set(QUALIFIERS)))
for label in relevant_labels:
    for m in range(1, 121):
        require(
            band_contains(m, label) == endpoint_contains(m, label),
            f"band/endpoint mismatch at m={m},label={label}",
        )


ANCHOR_PARAMETERS = tuple(
    m for m in range(1, 121)
    if all(band_contains(m, anchor) for anchor in ANCHORS)
)
require(len(ANCHOR_PARAMETERS) == 108, "anchor-compatible carrier count")


def maximal_optional_sets(newcomer=None):
    rows = []
    for m in ANCHOR_PARAMETERS:
        if newcomer is not None and not band_contains(m, newcomer):
            continue
        row = frozenset(
            label for label in OPTIONAL if band_contains(m, label)
        )
        if row not in rows:
            rows.append(row)
    return tuple(row for row in rows if not any(row < other for other in rows))


def gate_name(base, minimum, maximum):
    body_minimum = min((*base, minimum))
    body_maximum = max((*base, maximum))
    defect = 16 * body_maximum - 156 * body_minimum
    if defect <= 0:
        return "stable"
    if defect <= 13:
        return "affine"
    return "beyond"


def subset_gate_counts(labels, choose, base):
    """Count choose-subsets by body defect using their min/max labels."""
    ordered = tuple(sorted(labels))
    answer = Counter()
    if len(ordered) < choose:
        return answer
    for left in range(len(ordered)):
        for right in range(left + choose - 1, len(ordered)):
            count = comb(right - left - 1, choose - 2)
            answer[gate_name(base, ordered[left], ordered[right])] += count
    require(sum(answer.values()) == comb(len(ordered), choose),
            "min/max subset conservation")
    return answer


def union_gate_counts(rows, choose, base):
    """Inclusion-exclusion for the union of choose-subsets of carrier rows."""
    answer = Counter()
    for size in range(1, len(rows) + 1):
        sign = 1 if size % 2 else -1
        for selected in combinations(rows, size):
            intersection = set(selected[0])
            for row in selected[1:]:
                intersection.intersection_update(row)
            for gate, count in subset_gate_counts(
                    intersection, choose, base).items():
                answer[gate] += sign * count
    require(all(count >= 0 for count in answer.values()),
            "negative union stratum")
    return answer


def canonical_count(newcomer=None, choose=7):
    """Count H subset P_min(H) without enumerating optional subsets."""
    base = ANCHORS if newcomer is None else ANCHORS + (newcomer,)
    answer = 0

    # The newcomer itself can be the minimum; it is disjoint from the pool.
    if (newcomer is not None and newcomer < min(ANCHORS)
            and all(band_contains(newcomer, x) for x in ANCHORS)):
        available = sum(
            label > newcomer and band_contains(newcomer, label)
            for label in OPTIONAL
        )
        answer += comb(available, choose)

    # Otherwise a unique selected optional label is the minimum.
    for minimum in OPTIONAL:
        if minimum >= min(base):
            continue
        if all(band_contains(minimum, x) for x in base):
            available = sum(
                label > minimum and band_contains(minimum, label)
                for label in OPTIONAL
            )
            answer += comb(available, choose - 1)

    # Or the fixed anchor 120 is the minimum.
    if (all(x >= 120 for x in base)
            and all(band_contains(120, x) for x in base)):
        available = sum(
            label > 120 and band_contains(120, label)
            for label in OPTIONAL
        )
        answer += comb(available, choose)
    return answer


def thm4148_gate_count(base, choose):
    ordered = OPTIONAL
    answer = 0
    for left in range(len(ordered)):
        for right in range(left + choose - 1, len(ordered)):
            minimum = min((*base, ordered[left]))
            maximum = max((*base, ordered[right]))
            if 27 * (13 * minimum - maximum) >= 4 * minimum * maximum:
                answer += comb(right - left - 1, choose - 2)
    return answer


def family_row(newcomer=None, choose=7):
    base = ANCHORS if newcomer is None else ANCHORS + (newcomer,)
    total = subset_gate_counts(OPTIONAL, choose, base)
    rows = maximal_optional_sets(newcomer)
    overlap = union_gate_counts(rows, choose, base)
    canonical = canonical_count(newcomer, choose)
    require(sum(overlap.values()) <= sum(total.values()), "overlap exceeds family")
    require(canonical <= sum(overlap.values()), "canonical exceeds carrier union")
    return total, overlap, canonical, len(rows)


OLD_TOTAL, OLD_OVERLAP, OLD_CANONICAL, old_row_count = family_row(None, 8)
require(OLD_TOTAL == Counter(stable=344_366, beyond=1_875_709),
        "old gate census")
require(OLD_OVERLAP == Counter(stable=344_366, beyond=107_215),
        "old carrier census")
require(OLD_CANONICAL == 422_222, "old canonical census")

all_total = Counter()
all_overlap = Counter()
all_canonical = 0
strict_total = Counter()
strict_overlap = Counter()
strict_canonical = 0
one_total = Counter()
one_overlap = Counter()
one_canonical = 0
row_count_histogram = Counter()
nonzero_overlap_labels = 0
thm4148_hits = thm4148_gate_count(ANCHORS, 8)

for newcomer in QUALIFIERS:
    total, overlap, canonical, row_count = family_row(newcomer, 7)
    all_total.update(total)
    all_overlap.update(overlap)
    all_canonical += canonical
    row_count_histogram[row_count] += 1
    nonzero_overlap_labels += bool(sum(overlap.values()))
    thm4148_hits += thm4148_gate_count(ANCHORS + (newcomer,), 7)
    if newcomer in ONE_DELETION:
        one_total.update(total)
        one_overlap.update(overlap)
        one_canonical += canonical
    else:
        strict_total.update(total)
        strict_overlap.update(overlap)
        strict_canonical += canonical


EXPECTED = {
    "all_total": Counter(stable=24_915_898, affine=35_140,
                         beyond=891_495_922),
    "all_overlap": Counter(stable=24_915_898, affine=35_140,
                           beyond=57_935_164),
    "all_canonical": 49_558_619,
    "strict_total": Counter(stable=23_863_163, affine=35_140,
                            beyond=883_668_357),
    "strict_overlap": Counter(stable=23_863_163, affine=35_140,
                              beyond=57_375_523),
    "strict_canonical": 48_260_255,
    "one_total": Counter(stable=1_052_735, beyond=7_827_565),
    "one_overlap": Counter(stable=1_052_735, beyond=559_641),
    "one_canonical": 1_298_364,
    "row_histogram": Counter({0: 117, 1: 612, 2: 182, 3: 59, 4: 34, 5: 28}),
}
require(all_total == EXPECTED["all_total"], "all-q gate census")
require(all_overlap == EXPECTED["all_overlap"], "all-q overlap census")
require(all_canonical == EXPECTED["all_canonical"], "all-q canonical census")
require(strict_total == EXPECTED["strict_total"], "strict gate census")
require(strict_overlap == EXPECTED["strict_overlap"], "strict overlap census")
require(strict_canonical == EXPECTED["strict_canonical"],
        "strict canonical census")
require(one_total == EXPECTED["one_total"], "THM-4160 gate inheritance")
require(one_overlap == EXPECTED["one_overlap"], "THM-4160 overlap inheritance")
require(one_canonical == EXPECTED["one_canonical"],
        "THM-4160 canonical inheritance")
require(row_count_histogram == EXPECTED["row_histogram"],
        "maximal-row histogram")
require(thm4148_hits == 0, "THM-4148 named gate unexpectedly hit")


def scaled_carrier_mask(body, content):
    mask = 0
    for m in range(1, content * min(body) + 1):
        if all(band_contains(m, content * label) for label in body):
            mask |= 1 << (m - 1)
    return mask


SCALE_HOSTILE = (18, 20, 30, 40, 42, 60, 120, 126, 143, 168, 264)
hostile_one = scaled_carrier_mask(SCALE_HOSTILE, 1)
hostile_two = scaled_carrier_mask(SCALE_HOSTILE, 2)
require(18 in QUALIFIERS and 18 not in ONE_DELETION,
        "hostile not strict two-deletion newcomer")
require(hostile_one == 0, "hostile covered at content one")
require(hostile_two != 0, "hostile fails to enter after dilation")
require((hostile_two & -hostile_two).bit_length() == 35,
        "hostile least doubled carrier")


OLD_BODY_COUNT = comb(27, 8)
PER_NEWCOMER = comb(27, 7)
ALL_NEW_BODY_COUNT = len(QUALIFIERS) * PER_NEWCOMER
STRICT_BODY_COUNT = (len(QUALIFIERS) - len(ONE_DELETION)) * PER_NEWCOMER
COMBINED_BODY_COUNT = OLD_BODY_COUNT + ALL_NEW_BODY_COUNT

ALL_NEW_OVERLAP = sum(all_overlap.values())
ALL_NEW_OUTSIDE = ALL_NEW_BODY_COUNT - ALL_NEW_OVERLAP
STRICT_OVERLAP = sum(strict_overlap.values())
STRICT_OUTSIDE = STRICT_BODY_COUNT - STRICT_OVERLAP
COMBINED_OVERLAP = sum(OLD_OVERLAP.values()) + ALL_NEW_OVERLAP
COMBINED_OUTSIDE = COMBINED_BODY_COUNT - COMBINED_OVERLAP
COMBINED_CANONICAL = OLD_CANONICAL + all_canonical
COMBINED_TOTAL_GATES = OLD_TOTAL + all_total
COMBINED_OVERLAP_GATES = OLD_OVERLAP + all_overlap

require(ALL_NEW_BODY_COUNT == 916_446_960, "THM-4166 body count")
require(STRICT_BODY_COUNT == 907_566_660, "strict body count")
require(COMBINED_BODY_COUNT == 918_667_035, "combined body count")
require(ALL_NEW_OVERLAP == 82_886_202 and ALL_NEW_OUTSIDE == 833_560_758,
        "all-new overlap split")
require(STRICT_OVERLAP == 81_273_826 and STRICT_OUTSIDE == 826_292_834,
        "strict overlap split")
require(COMBINED_OVERLAP == 83_337_783
        and COMBINED_OUTSIDE == 835_329_252,
        "combined overlap split")
require(COMBINED_CANONICAL == 49_980_841, "combined canonical count")
require(COMBINED_OVERLAP_GATES["stable"] == COMBINED_TOTAL_GATES["stable"],
        "every stable body must be carrier-covered")
require(COMBINED_OVERLAP_GATES["affine"] == COMBINED_TOTAL_GATES["affine"],
        "every affine body must be carrier-covered")

FULL_M7_UNION = 81_572_506_886_508 + COMBINED_BODY_COUNT
STRUCTURED_M7_UNION = 38_620_298_376 + COMBINED_BODY_COUNT
require(FULL_M7_UNION == 81_573_425_553_543, "full m7 union")
require(STRUCTURED_M7_UNION == 39_538_965_411, "structured m7 union")


lines = [
    "THM4166_VS_THM4158_GROUPED_CONTENT_ONE_OVERLAP_AUDIT",
    "status=PASS;scope=FINITE_EXACT_undilated_core_body_labels_only",
    f"qualifier_input_sha256={QUALIFIER_OUTPUT_SHA256}",
    "algorithm=maximal_carrier_rows_plus_inclusion_exclusion;"
    "minmax_subset_formula;no_916M_body_scan",
    f"qualifiers={len(QUALIFIERS)};last={QUALIFIERS[-1]};"
    f"anchor_carriers={len(ANCHOR_PARAMETERS)};"
    f"maximal_row_histogram={tuple(sorted(row_count_histogram.items()))}",
    f"carrier_nonzero_qualifiers={nonzero_overlap_labels};"
    f"carrier_zero_qualifiers={len(QUALIFIERS)-nonzero_overlap_labels}",
    f"old=total:{OLD_BODY_COUNT};canonical:{OLD_CANONICAL};"
    f"global:{sum(OLD_OVERLAP.values())};"
    f"outside:{OLD_BODY_COUNT-sum(OLD_OVERLAP.values())}",
    f"all_new=total:{ALL_NEW_BODY_COUNT};canonical:{all_canonical};"
    f"global:{ALL_NEW_OVERLAP};outside:{ALL_NEW_OUTSIDE}",
    f"strict_beyond_THM4160=total:{STRICT_BODY_COUNT};"
    f"canonical:{strict_canonical};global:{STRICT_OVERLAP};"
    f"outside:{STRICT_OUTSIDE}",
    f"combined=total:{COMBINED_BODY_COUNT};canonical:{COMBINED_CANONICAL};"
    f"global:{COMBINED_OVERLAP};outside:{COMBINED_OUTSIDE}",
    "combined_gate_total="
    f"stable:{COMBINED_TOTAL_GATES['stable']};"
    f"affine:{COMBINED_TOTAL_GATES['affine']};"
    f"beyond:{COMBINED_TOTAL_GATES['beyond']}",
    "combined_gate_overlap="
    f"stable:{COMBINED_OVERLAP_GATES['stable']};"
    f"affine:{COMBINED_OVERLAP_GATES['affine']};"
    f"beyond:{COMBINED_OVERLAP_GATES['beyond']}",
    "strict_gate_overlap="
    f"stable:{strict_overlap['stable']};affine:{strict_overlap['affine']};"
    f"beyond:{strict_overlap['beyond']}",
    f"THM4148_named_gate_hits={thm4148_hits}",
    "named_m7_disjoint=q7_absent_from_pool_and_qualifiers;"
    f"full_union:{FULL_M7_UNION};structured_union:{STRUCTURED_M7_UNION}",
    "scale_hostile=(18,20,30,40,42,60,120,126,143,168,264);"
    "content1:none;content2:least_m35",
    "scope_warning=outside_carrier_counts_are_not_dilation_invariant;"
    "LRC14_OPEN",
]
semantic_payload = "\n".join(lines) + "\n"
lines.append("semantic_sha256=" + sha256(semantic_payload.encode()).hexdigest())
sys.stdout.write("\n".join(lines) + "\n")
