#!/usr/bin/env python3
"""Close the unique projected k=3, z1=380 frontier row of THM-2941.

THM-2941 leaves a unique projected k=3 row at its maximum first drift:

    E=(1,4,8,10,12,14), z1=380.

For this fixed row, the singleton-excess inequality requires

    delta(z1)+delta(z2)+delta(z3)+delta(z4) >= h * 3/91,

with 380=z1<z2<z3<z4 and at least one drift at the projected wall
z>=1159.  Exact suffix dynamic programming through the inherited horizon
7,000 finds all finite packets.  The inherited omitted-label estimate
delta(z)<=6r/(49z) proves that no packet using a label beyond the horizon
can reach the inequality.

For every resulting packet, reduce its four drift speeds on the body
ruler L to denominators d_i=L/gcd(L,z_i), and put D=lcm(d_i).  If q|D,
M=D/q, g_i=gcd(d_i,q), ell_i=ceil(d_i/7), and
H_i=D/lcm(d_i,q), one q-fibre of the four danger combs has capacity at
most

    sum_i H_i ceil(ell_i/g_i).

The exact projected body-safe support has a q-fibre exceeding this bound
for every packet.  Thus the unique z1=380 row is impossible, improving
the projected k=3 first-drift cap from 380 to 379.  This script makes no
claim about rows with z1<=379 and performs no wider census.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
SUFFIX_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)
FIBRE_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_frontier_fibre_closure_thm2941.out"
)

EXPECTED_SUFFIX_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_SUFFIX_OUTPUT_SHA256 = (
    "61e16aab8a368881c574047e576645e6b41837dc9f804f7a78d37230d843612b"
)
EXPECTED_FIBRE_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
EXPECTED_SEMANTIC_SHA256 = (
    "e9a42f99a855db99c4af95c383ac942683d57834a165f2f56733bebb6e0d9a26"
)

BODY = (1, 4, 8, 10, 12, 14)
FIRST = 380
K = 3
REMAINING = 3
EXPECTED_PACKETS = (
    (380, 410, 492, 1164),
    (380, 410, 492, 1220),
    (380, 410, 492, 1358),
    (380, 410, 492, 1500),
    (380, 410, 492, 1836),
)
EXPECTED_FRONTIER_SUMMARY = (
    "  mode=projected;surviving_first_rows=376020;"
    "max_first=380;frontier_rows=1"
)
EXPECTED_FRONTIER_ROW = (
    "    FRONTIER;E=1,4,8,10,12,14;h=1049/2940;r=34;L=11760;"
    "analytic_cap=1415;largest_floor=1159;z1=380;"
    "delta1=1531/391020;"
    "suffix=1500:237/171500,492:2981/843780,410:2687/843780;"
    "lower=1049/89180;upper=401287/33399625;"
    "gap=437649/1736780500"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(SUFFIX_PATH) == EXPECTED_SUFFIX_SHA256,
    "THM-2941 suffix dependency hash changed",
)
require(
    file_sha256(SUFFIX_OUTPUT_PATH) == EXPECTED_SUFFIX_OUTPUT_SHA256,
    "THM-2941 suffix output dependency hash changed",
)
require(
    file_sha256(FIBRE_PATH) == EXPECTED_FIBRE_SHA256,
    "THM-2928 fibre dependency hash changed",
)
suffix = load_module("thm2941_projected_suffix", SUFFIX_PATH)
fibre = load_module("thm2928_body_fibre", FIBRE_PATH)
support = fibre.support_module


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def fibre_capacity(D, d, q):
    ell = (d + 6) // 7
    common = gcd(d, q)
    height = D // lcm(d, q)
    return height * ((ell + common - 1) // common)


def main():
    inherited_output = SUFFIX_OUTPUT_PATH.read_text()
    require(
        EXPECTED_FRONTIER_SUMMARY in inherited_output,
        "unique projected k=3 frontier summary changed",
    )
    require(
        EXPECTED_FRONTIER_ROW in inherited_output,
        "unique projected k=3 frontier row changed",
    )
    require(
        suffix.HORIZON == 7_000
        and suffix.ETAS[K] == F(3, 91)
        and suffix.PROJECTED_RATIOS[K] == F(13, 132),
        "inherited k=3 constants changed",
    )

    carrier = suffix.A.carrier_for(BODY)
    h = F(sum(right - left for left, right in carrier), suffix.A.RULER)
    components = len(carrier)
    L = 14 * lcm(*BODY)
    lower = h * suffix.ETAS[K]
    wall = suffix.PROJECTED_RATIOS[K] * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    require(
        (h, components, L, lower, high_floor)
        == (F(1049, 2940), 34, 11760, F(1049, 89180), 1159),
        "frontier geometry changed",
    )

    delta = {
        label: suffix.A.singleton_coverage(carrier, label) - h / 7
        for label in range(FIRST, suffix.HORIZON + 1)
        if label % L
    }
    require(FIRST in delta, "first drift disappeared")
    labels = tuple(
        label
        for label in range(FIRST + 1, suffix.HORIZON + 1)
        if label % L
    )
    values = tuple(delta[label] for label in labels)
    count = len(labels)

    # best[index][slots][need_high] is the exact maximum suffix excess
    # using exactly ``slots`` labels at indices >= index.  The Boolean
    # coordinate requires at least one label at or beyond the wall.
    best = [
        [[None, None] for _slots in range(REMAINING + 1)]
        for _index in range(count + 1)
    ]
    best[count][0][0] = F(0)
    for index in range(count - 1, -1, -1):
        is_high = labels[index] >= high_floor
        best[index][0][0] = F(0)
        for slots in range(1, REMAINING + 1):
            for need_high in (0, 1):
                choices = []
                skipped = best[index + 1][slots][need_high]
                if skipped is not None:
                    choices.append(skipped)
                next_need = int(bool(need_high and not is_high))
                taken = best[index + 1][slots - 1][next_need]
                if taken is not None:
                    choices.append(values[index] + taken)
                best[index][slots][need_high] = max(choices) if choices else None

    nodes_by_depth = Counter()
    tests_by_depth = Counter()
    rejected_by_bound = Counter()
    packets = []

    def visit(start, slots, need_high, total, prefix):
        depth = len(prefix)
        nodes_by_depth[depth] += 1
        bound = best[start][slots][int(need_high)]
        require(bound is not None and total + bound >= lower, "invalid live node")
        if slots == 0:
            require(not need_high and total >= lower, "invalid packet leaf")
            packets.append((FIRST, *prefix))
            return

        for index in range(start, count - slots + 1):
            tests_by_depth[depth] += 1
            label = labels[index]
            next_high = bool(need_high and label < high_floor)
            tail_bound = best[index + 1][slots - 1][int(next_high)]
            next_total = total + values[index]
            if tail_bound is None or next_total + tail_bound < lower:
                rejected_by_bound[depth] += 1
                continue
            visit(
                index + 1,
                slots - 1,
                next_high,
                next_total,
                (*prefix, label),
            )

    visit(
        0,
        REMAINING,
        FIRST < high_floor,
        delta[FIRST],
        (),
    )
    require(tuple(packets) == EXPECTED_PACKETS, ("packet stream changed", packets))

    # Any tuple using an omitted label has at most the two best exact suffix
    # excesses plus one copy of the uniform omitted-label bound.  This also
    # bounds tuples with two or three omitted labels because the two best
    # exact excesses dominate that tail bound.
    ordinary_tail = F(6 * components, 49 * (suffix.HORIZON + 1))
    exact_top_two = sorted(
        ((value, label) for value, label in zip(values, labels)),
        reverse=True,
    )[:2]
    require(
        len(exact_top_two) == 2
        and exact_top_two[1][0] >= ordinary_tail,
        "omitted-label comparison changed",
    )
    omitted_upper = (
        delta[FIRST]
        + exact_top_two[0][0]
        + exact_top_two[1][0]
        + ordinary_tail
    )
    require(omitted_upper < lower, "an omitted suffix may still reach the wall")

    actual_L, ranges = support.safe_cell_ranges(BODY)
    require(actual_L == L, "body support ruler changed")
    packet_rows = []
    for packet in packets:
        ds = tuple(L // gcd(label, L) for label in packet)
        D = lcm(*ds)
        arcs = fibre.projected_support_arcs(D, ranges)
        target = sum(right - left for left, right in arcs)
        ambient = sum((D // d) * ((d + 6) // 7) for d in ds)
        require(target <= ambient, "packet unexpectedly fails scalar capacity")

        best_witness = None
        for q in support.divisors(D):
            histogram = fibre.residue_load_histogram(arcs, q)
            maximum = max(
                load for load, multiplicity in histogram if multiplicity
            )
            capacity = sum(fibre_capacity(D, d, q) for d in ds)
            candidate = (
                maximum - capacity,
                q,
                D // q,
                maximum,
                capacity,
            )
            if best_witness is None or candidate > best_witness:
                best_witness = candidate
        require(
            best_witness is not None and best_witness[0] > 0,
            f"packet survived all-divisor capacity: {packet}",
        )
        packet_rows.append((packet, ds, D, target, ambient, best_witness))

    require(len(packet_rows) == 5, "closed packet count changed")
    semantic_payload = (
        BODY,
        FIRST,
        h,
        components,
        L,
        lower,
        high_floor,
        tuple(packets),
        ordinary_tail,
        omitted_upper,
        tuple(packet_rows),
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    require(
        semantic_hash == EXPECTED_SEMANTIC_SHA256,
        "frontier-row closure semantic digest changed",
    )

    lines = [
        "LRC14 projected k=3 frontier-row fibre closure",
        f"suffix_source_sha256={file_sha256(SUFFIX_PATH)}",
        f"suffix_output_sha256={file_sha256(SUFFIX_OUTPUT_PATH)}",
        f"fibre_source_sha256={file_sha256(FIBRE_PATH)}",
        "scope=unique THM-2941 projected k=3 max-first row only;"
        "no claim for z1<=379 and no wider census",
        f"body={','.join(map(str, BODY))};h={ftext(h)};"
        f"components={components};L={L}",
        f"k={K};z1={FIRST};lower={ftext(lower)};"
        f"projected_high_floor={high_floor};exact_horizon={suffix.HORIZON}",
        f"suffix_labels={count};"
        f"stream_nodes_by_depth={tuple(nodes_by_depth[d] for d in range(4))};"
        f"candidate_tests_by_depth={tuple(tests_by_depth[d] for d in range(3))};"
        "bound_rejections_by_depth="
        f"{tuple(rejected_by_bound[d] for d in range(3))}",
        f"omitted_tail={ftext(ordinary_tail)};"
        f"omitted_upper={ftext(omitted_upper)};"
        f"omitted_gap={ftext(lower - omitted_upper)}",
        f"finite_packets={len(packet_rows)}",
    ]
    for packet, ds, D, target, ambient, witness in packet_rows:
        margin, q, M, maximum, capacity = witness
        lines.append(
            f"  packet={','.join(map(str, packet))};"
            f"denominators={','.join(map(str, ds))};D={D};"
            f"target={target};ambient_capacity={ambient};"
            f"q_witness={q};M={M};max_fibre_load={maximum};"
            f"fibre_capacity={capacity};margin={margin};"
            "all_divisor_screen=KILL"
        )
    lines.extend(
        (
            "closed_packets=5/5",
            "conclusion=projected k=3 first-drift cap improves 380->379",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    output = "\n".join(lines) + "\n"
    DEFAULT_OUTPUT.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
