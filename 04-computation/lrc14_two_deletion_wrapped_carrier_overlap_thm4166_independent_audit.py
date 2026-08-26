#!/usr/bin/env python3
"""Clean-room THM-4166 versus THM-4158 content-one overlap census.

The script parses the 1,032 admitted newcomer labels from the frozen
THM-4166 global-census output.  It does not import either theorem's code.

The 916,446,960 body universe is never enumerated.  Instead, the 888,030
seven-subsets of the fixed 27-label optional pool are enumerated once and
compressed by

    (common THM-4158 carrier mask, optional minimum, optional maximum).

Only 376 groups survive.  Every admitted newcomer is classified against
those groups.  A second route forms the maximal THM-4158 optional-label rows
for that newcomer and uses inclusion-exclusion (at most five rows), with an
exact min/max count for every intersection.  A third direct minimum-indexed
count checks the canonical subfamily.  Thus carrier collisions, gate splits,
and the canonical/global distinction are all checked without a 916M sweep.
"""

from ast import literal_eval
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition, message):
    if not condition:
        raise RuntimeError("REQUIRE FAILED: " + message)


ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(label for label in POOL if label not in ANCHORS)
ONE_DELETION_QUALIFIERS = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)
MAX_CONTENT_ONE_M = 120
FULL_OPTIONAL_MASK = (1 << len(OPTIONAL)) - 1

REPO = Path(__file__).resolve().parents[1]
FROZEN = (
    REPO
    / "05-knowledge/results/"
    / "lrc14_two_deletion_repair_graph_thm4166_global_census.out"
)
FROZEN_SHA256 = "baae4adc57739c852dec2a8c666328cea2c4031934547e2afa4044a314c9bff6"


require(len(ANCHORS) == 3, "anchor count")
require(len(POOL) == 30 and len(OPTIONAL) == 27, "pool cardinalities")
require(set(ANCHORS).issubset(POOL), "anchors lie in pool")
require(len(set(POOL)) == len(POOL), "pool labels distinct")


frozen_bytes = FROZEN.read_bytes()
require(sha256(frozen_bytes).hexdigest() == FROZEN_SHA256,
        "frozen THM-4166 output hash")
require(b"\r" not in frozen_bytes, "frozen THM-4166 output is raw LF")
frozen_text = frozen_bytes.decode("utf-8")
admitted_lines = [
    line for line in frozen_text.splitlines() if line.startswith("admitted_q ")
]
require(len(admitted_lines) == 1, "unique admitted_q ledger")
QUALIFIERS = tuple(literal_eval(admitted_lines[0][len("admitted_q "):]))

require(len(QUALIFIERS) == 1_032, "qualifier count")
require(QUALIFIERS == tuple(sorted(set(QUALIFIERS))),
        "qualifiers strictly increasing and unique")
require(QUALIFIERS[-1] == 8_265, "last qualifier")
require(set(QUALIFIERS).isdisjoint(POOL), "every qualifier lies outside pool")
require("admitted_count 1032 last_admitted 8265" in frozen_text,
        "frozen admitted count/endpoint")
require("bodies_per_q 888030" in frozen_text, "frozen bodies-per-q count")
require("extension_bodies 916446960" in frozen_text,
        "frozen extension-body count")


def endpoint_safe(m, label):
    """Exact THM-4158 carrier membership from its endpoint inequalities."""

    require(m >= 1 and label >= 1, "positive endpoint inputs")
    last_tooth = 3 if m == 1 else 2
    return any(
        m * (14 * tooth + 1) <= label
        and 16 * label <= (12 * m + 1) * (14 * tooth + 13)
        for tooth in range(last_tooth + 1)
    )


def carrier_mask(label):
    return sum(
        1 << (m - 1)
        for m in range(1, MAX_CONTENT_ONE_M + 1)
        if endpoint_safe(m, label)
    )


ALL_LABELS = set(POOL) | set(QUALIFIERS)
LABEL_MASK = {label: carrier_mask(label) for label in ALL_LABELS}

ANCHOR_MASK = (1 << MAX_CONTENT_ONE_M) - 1
for anchor in ANCHORS:
    ANCHOR_MASK &= LABEL_MASK[anchor]
require(ANCHOR_MASK.bit_count() == 108, "anchor-compatible carrier count")


# Row m is the optional-label set P_m intersect OPTIONAL.  These rows power
# the independent inclusion-exclusion route.
OPTIONAL_ROW = {}
for m in range(1, MAX_CONTENT_ONE_M + 1):
    row = 0
    for index, label in enumerate(OPTIONAL):
        if endpoint_safe(m, label):
            row |= 1 << index
    OPTIONAL_ROW[m] = row


def gate_name(minimum, maximum):
    defect = 16 * maximum - 156 * minimum
    if defect <= 0:
        return "stable"
    if defect <= 13:
        return "transition"
    return "beyond"


# Enumerate the fixed 888,030 optional choices only once.  The carrier mask
# is the intersection of all label masks, already including the anchors.
GROUPS = Counter()
for choice in combinations(OPTIONAL, 7):
    mask = ANCHOR_MASK
    for label in choice:
        mask &= LABEL_MASK[label]
    GROUPS[(mask, choice[0], choice[-1])] += 1

require(sum(GROUPS.values()) == comb(27, 7) == 888_030,
        "single optional-universe enumeration")
require(len(GROUPS) == 376, "compressed group count")
require(len({key[0] for key in GROUPS}) == 160,
        "compressed carrier-mask count")


def optional_mask_gate_counts(optional_mask, q):
    """Count seven-subsets of one optional mask by THM-4151 defect gate."""

    labels = [
        label for index, label in enumerate(OPTIONAL)
        if optional_mask & (1 << index)
    ]
    result = Counter()
    for left in range(len(labels)):
        for right in range(left + 6, len(labels)):
            ways = comb(right - left - 1, 5)
            minimum = min(120, q, labels[left])
            maximum = max(143, q, labels[right])
            result[gate_name(minimum, maximum)] += ways
    require(sum(result.values()) == comb(len(labels), 7),
            "min/max gate count conserves optional mask")
    return result


def maximal_optional_rows(q):
    """Maximal optional rows among all carrier parameters for anchors+q."""

    base_mask = ANCHOR_MASK & LABEL_MASK[q]
    rows = []
    for m in range(1, MAX_CONTENT_ONE_M + 1):
        if not (base_mask & (1 << (m - 1))):
            continue
        row = OPTIONAL_ROW[m]
        if row.bit_count() >= 7 and row not in rows:
            rows.append(row)
    maximal = [
        row for row in rows
        if not any(row != other and row & ~other == 0 for other in rows)
    ]
    return tuple(sorted(maximal))


def inclusion_exclusion_global(q):
    """Independent global carrier-union count, including all gate splits."""

    rows = maximal_optional_rows(q)
    result = Counter()
    term_count = 0
    for size in range(1, len(rows) + 1):
        sign = 1 if size % 2 else -1
        for selected in combinations(rows, size):
            intersection = FULL_OPTIONAL_MASK
            for row in selected:
                intersection &= row
            profile = optional_mask_gate_counts(intersection, q)
            for gate, count in profile.items():
                result[gate] += sign * count
            term_count += 1
    require(all(result[gate] >= 0 for gate in ("stable", "transition", "beyond")),
            "nonnegative inclusion-exclusion result")
    return result, len(rows), term_count


def direct_canonical(q):
    """Direct minimum-indexed count, independent of common carrier masks."""

    base = ANCHORS + (q,)
    base_minimum = min(base)
    base_maximum = max(base)
    result = Counter()

    # No optional label below base_minimum is chosen, so the body's minimum is
    # the minimum already present in anchors+q.
    if all(endpoint_safe(base_minimum, label) for label in base):
        candidates = [
            label for label in OPTIONAL
            if label >= base_minimum and endpoint_safe(base_minimum, label)
        ]
        for index in range(6, len(candidates)):
            ways = comb(index, 6)
            maximum = max(base_maximum, candidates[index])
            result[gate_name(base_minimum, maximum)] += ways

    # Otherwise (or additionally), the unique least selected optional label r
    # is the body's minimum.  Choose the remaining six labels above r.
    for minimum in OPTIONAL:
        if minimum >= base_minimum:
            break
        if not all(endpoint_safe(minimum, label) for label in base):
            continue
        candidates = [
            label for label in OPTIONAL
            if label > minimum and endpoint_safe(minimum, label)
        ]
        for index in range(5, len(candidates)):
            ways = comb(index, 5)
            maximum = max(base_maximum, candidates[index])
            result[gate_name(minimum, maximum)] += ways

    return result


def grouped_classification(q):
    """Main compressed carrier-mask classification for one qualifier."""

    base_mask = ANCHOR_MASK & LABEL_MASK[q]
    total = Counter()
    canonical = Counter()
    global_union = Counter()
    outside = Counter()
    witness_histogram = Counter()

    for (choice_mask, choice_minimum, choice_maximum), multiplicity in GROUPS.items():
        minimum = min(120, q, choice_minimum)
        maximum = max(143, q, choice_maximum)
        gate = gate_name(minimum, maximum)
        common_mask = base_mask & choice_mask
        total[gate] += multiplicity
        if common_mask:
            global_union[gate] += multiplicity
            witness_histogram[common_mask.bit_count()] += multiplicity
            if common_mask & (1 << (minimum - 1)):
                canonical[gate] += multiplicity
        else:
            outside[gate] += multiplicity

    return total, canonical, global_union, outside, witness_histogram


GATES = ("stable", "transition", "beyond")
aggregate = {
    "total": Counter(),
    "canonical": Counter(),
    "global": Counter(),
    "outside": Counter(),
}
aggregate_witness_histogram = Counter()
maximal_row_histogram = Counter()
inclusion_exclusion_terms = 0
per_q_lines = []
per_q = {}

for q in QUALIFIERS:
    total, canonical, global_union, outside, witnesses = grouped_classification(q)
    ie_global, maximal_row_count, term_count = inclusion_exclusion_global(q)
    direct_minimum = direct_canonical(q)
    direct_total = optional_mask_gate_counts(FULL_OPTIONAL_MASK, q)

    require(total == direct_total, f"q={q} total min/max route")
    require(global_union == ie_global, f"q={q} global inclusion-exclusion route")
    require(canonical == direct_minimum, f"q={q} canonical direct route")
    require(total == global_union + outside, f"q={q} inside/outside partition")
    require(all(canonical[g] <= global_union[g] for g in GATES),
            f"q={q} canonical contained in global")

    for name, row in (
        ("total", total),
        ("canonical", canonical),
        ("global", global_union),
        ("outside", outside),
    ):
        aggregate[name].update(row)
    aggregate_witness_histogram.update(witnesses)
    maximal_row_histogram[maximal_row_count] += 1
    inclusion_exclusion_terms += term_count

    summary = {
        "total": sum(total.values()),
        "canonical": sum(canonical.values()),
        "global": sum(global_union.values()),
        "outside": sum(outside.values()),
        "multi": sum(count for witnesses_count, count in witnesses.items()
                     if witnesses_count > 1),
        "incidence": sum(witnesses_count * count
                         for witnesses_count, count in witnesses.items()),
        "rows": maximal_row_count,
    }
    per_q[q] = summary
    per_q_lines.append(
        ";".join((
            f"q={q}",
            "total=" + ",".join(str(total[g]) for g in GATES),
            "canonical=" + ",".join(str(canonical[g]) for g in GATES),
            "global=" + ",".join(str(global_union[g]) for g in GATES),
            "outside=" + ",".join(str(outside[g]) for g in GATES),
            f"multi={summary['multi']}",
            f"incidence={summary['incidence']}",
            f"maximal_rows={maximal_row_count}",
        ))
    )


EXPECTED = {
    "total": (24_915_898, 35_140, 891_495_922),
    "canonical": (24_915_898, 35_140, 24_607_581),
    "global": (24_915_898, 35_140, 57_935_164),
    "outside": (0, 0, 833_560_758),
}
for name, expected in EXPECTED.items():
    actual = tuple(aggregate[name][gate] for gate in GATES)
    require(actual == expected, f"frozen aggregate {name} gate split")

require(sum(aggregate["total"].values()) == 916_446_960,
        "full THM-4166 body count")
require(sum(aggregate["canonical"].values()) == 49_558_619,
        "canonical overlap count")
require(sum(aggregate["global"].values()) == 82_886_202,
        "global overlap count")
require(sum(aggregate["outside"].values()) == 833_560_758,
        "outside carrier-union count")

MULTI_BODY_COUNT = sum(
    count for witnesses_count, count in aggregate_witness_histogram.items()
    if witnesses_count > 1
)
WITNESS_INCIDENCE = sum(
    witnesses_count * count
    for witnesses_count, count in aggregate_witness_histogram.items()
)
require(MULTI_BODY_COUNT == 69_768_295, "multiple-carrier body count")
require(WITNESS_INCIDENCE == 830_026_747, "carrier witness incidence")
require(sum(aggregate_witness_histogram.values()) == 82_886_202,
        "witness histogram counts global bodies once")


# Recover the previously audited ten one-deletion rows as a positive control.
one_deletion = Counter()
for q in ONE_DELETION_QUALIFIERS:
    require(q in per_q, f"one-deletion qualifier {q} admitted by THM-4166")
    total, canonical, global_union, outside, _ = grouped_classification(q)
    for name, row in (
        ("total", total),
        ("canonical", canonical),
        ("global", global_union),
        ("outside", outside),
    ):
        one_deletion[name] += sum(row.values())
    for name, row in (
        ("total", total),
        ("canonical", canonical),
        ("global", global_union),
        ("outside", outside),
    ):
        for gate in GATES:
            one_deletion[f"{name}_{gate}"] += row[gate]

require(one_deletion["total"] == 8_880_300, "THM-4160 new-family total")
require(one_deletion["canonical"] == 1_298_364,
        "THM-4160 new-family canonical control")
require(one_deletion["global"] == 1_612_376,
        "THM-4160 new-family global control")
require(one_deletion["outside"] == 7_267_924,
        "THM-4160 new-family outside control")


# Representation injectivity is structural.  Every body contains exactly one
# label outside POOL, namely q, and its other variable labels form K subset O.
require(len(QUALIFIERS) * comb(27, 7) == 916_446_960,
        "injective qualifier/optional representation count")


def scaled_carriers(body, content):
    # P_m has minimum label m, so m <= content*min(body) is exhaustive.
    return tuple(
        m for m in range(1, content * min(body) + 1)
        if all(endpoint_safe(m, content * label) for label in body)
    )


# Scope hostile: one admitted new body is outside every carrier at content one
# but its common double lies in P_35.  Hence the outside count is not an
# all-content separation.
DILATION_HOSTILE = (20, 30, 40, 42, 60, 66, 120, 126, 143, 168, 264)
hostile_optional = set(DILATION_HOSTILE) - set(ANCHORS) - {66}
require(66 in QUALIFIERS, "hostile qualifier admitted")
require(hostile_optional.issubset(OPTIONAL) and len(hostile_optional) == 7,
        "hostile is a THM-4166 body")
hostile_content_one = scaled_carriers(DILATION_HOSTILE, 1)
hostile_content_two = scaled_carriers(DILATION_HOSTILE, 2)
require(hostile_content_one == (), "hostile outside at content one")
require(hostile_content_two == (35,), "hostile double has unique carrier P_35")


per_q_payload = "\n".join(per_q_lines) + "\n"
PER_Q_SHA256 = sha256(per_q_payload.encode("utf-8")).hexdigest()
transition_q = tuple(
    q for q in QUALIFIERS
    if grouped_classification(q)[0]["transition"] > 0
)
global_zero_q = tuple(q for q in QUALIFIERS if per_q[q]["global"] == 0)
canonical_zero_q = tuple(q for q in QUALIFIERS if per_q[q]["canonical"] == 0)
max_global = max(per_q[q]["global"] for q in QUALIFIERS)
max_global_q = tuple(q for q in QUALIFIERS if per_q[q]["global"] == max_global)
max_canonical = max(per_q[q]["canonical"] for q in QUALIFIERS)
max_canonical_q = tuple(
    q for q in QUALIFIERS if per_q[q]["canonical"] == max_canonical
)


def compact(counter):
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


semantic_lines = [
    "THM4166_VS_THM4158_CONTENT_ONE_GROUPED_OVERLAP_AUDIT",
    "status=PASS",
    f"frozen_input_sha256={FROZEN_SHA256}",
    "qualifiers=1032;bodies_per_q=888030;total=916446960",
    "carrier_parameter_bound=1..120;anchor_compatible=108",
    "grouped_optional_enumeration=888030_once;groups=376;carrier_masks=160",
    "gate_order=stable,transition,beyond",
]
for name in ("total", "canonical", "global", "outside"):
    semantic_lines.append(
        name + "=" + ",".join(str(aggregate[name][gate]) for gate in GATES)
        + ";sum=" + str(sum(aggregate[name].values()))
    )
semantic_lines.extend((
    f"multiple_carrier_bodies={MULTI_BODY_COUNT};witness_incidence={WITNESS_INCIDENCE}",
    "witness_histogram=" + compact(aggregate_witness_histogram),
    "maximal_row_histogram=" + compact(maximal_row_histogram),
    f"inclusion_exclusion_terms={inclusion_exclusion_terms}",
    f"per_q_ledger_sha256={PER_Q_SHA256}",
    "transition_q=" + repr(transition_q),
    f"global_zero_q_count={len(global_zero_q)};canonical_zero_q_count={len(canonical_zero_q)}",
    f"max_global={max_global};q={max_global_q}",
    f"max_canonical={max_canonical};q={max_canonical_q}",
    "thm4160_new_control="
    f"total:{one_deletion['total']},canonical:{one_deletion['canonical']},"
    f"global:{one_deletion['global']},outside:{one_deletion['outside']}",
    "representation_collisions=0;reason=unique_outside_POOL_qualifier",
    "dilation_hostile=" + repr(DILATION_HOSTILE)
    + ";content1:none;content2:P_35",
    "scope=FINITE_EXACT_undilated_content_one_new_THM4166_bodies_only",
))
semantic_payload = "\n".join(semantic_lines) + "\n"
semantic_lines.append(
    "semantic_sha256=" + sha256(semantic_payload.encode("utf-8")).hexdigest()
)
semantic_lines.append("result=PASS")
sys.stdout.write("\n".join(semantic_lines) + "\n")
