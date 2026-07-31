#!/usr/bin/env python3
"""Independent exact referee for the fixed projected z1=1932 state.

This audit does not import the primary z1932 closure.  It independently
reconstructs the four all-label ray maxima, enumerates the 125 five-vertex
trees through Pruefer codes, recomputes all
32 Hunter status capacities, and checks the decisive rational Farkas
inequality cell by cell.
"""

from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1932_ray_status_independent_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46"
)
EXPECTED_SEMANTIC_SHA256 = (
    "ffcf43a7f1c256c4735b5ce8ccfbd2ca74b0f346db664059b1fab4f78ac50e51"
)
BODY = (1, 4, 8, 10, 12, 14)
FIRST = 1932
LABELS = (1932, 2060, 2142, 2172, 2534)
DS = (140, 280, 588, 840, 980)
D = 5880
Q = 840
M = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256,
    "uniform ray/status dependency changed",
)
spec = spec_from_file_location("k2_z2060_independent_support", UNIFORM)
require(spec is not None and spec.loader is not None, "cannot load uniform")
U = module_from_spec(spec)
spec.loader.exec_module(U)


def prufer_tree(code, vertex_count=5):
    """Decode one Pruefer word, independently of the primary tree search."""
    degree = [1] * vertex_count
    for vertex in code:
        degree[vertex] += 1
    edges = []
    for vertex in code:
        leaf = min(index for index, value in enumerate(degree) if value == 1)
        edges.append(tuple(sorted((leaf, vertex))))
        degree[leaf] -= 1
        degree[vertex] -= 1
    leaves = [index for index, value in enumerate(degree) if value == 1]
    require(len(leaves) == 2, ("Pruefer decoder failed", code, degree))
    edges.append(tuple(sorted(leaves)))
    return tuple(sorted(edges))


TREES5 = tuple(
    prufer_tree(code) for code in product(range(5), repeat=3)
)
require(len(TREES5) == 125 and len(set(TREES5)) == 125, "tree atlas changed")


def hunter_cap(M, ds, es):
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (left, right): U.pair_lower(
            M, ds[left], es[left], ds[right], es[right]
        )
        for left in range(5)
        for right in range(left + 1, 5)
    }
    invoice = max(
        sum(overlaps[edge] for edge in tree) for tree in TREES5
    )
    return min(M, sum(sizes) - invoice)


def main():
    carrier = U.suffix.A.carrier_for(BODY)
    h = F(
        sum(right - left for left, right in carrier),
        U.suffix.A.RULER,
    )
    L = 14 * lcm(*BODY)
    lower = h * U.suffix.ETAS[2]
    deltas = tuple(U.delta(carrier, h, label) for label in LABELS)
    require(
        (h, len(carrier), L, lower, deltas)
        == (
            F(1049, 2940),
            34,
            11760,
            F(1049, 267540),
            (
                F(101, 157780),
                F(1867, 2119740),
                F(821, 1049580),
                F(263, 310415),
                F(2911, 3724980),
            ),
        ),
        "fixed scalar data changed",
    )
    label_ds = tuple(sorted(L // gcd(L, label) for label in LABELS))
    upper = sum(deltas, F(0))
    gap = upper - lower
    require(
        (label_ds, upper, gap)
        == (
            DS,
            F(442380823, 112512089655),
            F(12804017, 1170125732412),
        ),
        "fixed scalar state changed",
    )

    # Independently maximize every required denominator ray beyond FIRST.
    ray_maxima = []
    for d in (280, 588, 840, 980):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = residue * U.delta(carrier, h, residue)
            if amplitude < 0:
                continue
            label = residue
            if label <= FIRST:
                label += ((FIRST + 1 - label + L - 1) // L) * L
            candidates.append((amplitude / label, label, residue))
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        ray_maxima.append((d, *candidates[0]))
    require(
        tuple(ray_maxima)
        == (
            (280, F(821, 1049580), 2142, 2142),
            (588, F(1867, 2119740), 2060, 2060),
            (840, F(2911, 3724980), 2534, 2534),
            (980, F(263, 310415), 2172, 2172),
        ),
        ("all-label ray maxima changed", ray_maxima),
    )

    inner_ds = []
    lows = []
    marginals = []
    for d in DS:
        common = gcd(d, Q)
        low, remainder = divmod((d + 6) // 7, common)
        inner_ds.append(d // common)
        lows.append(low)
        marginals.append((Q // common) * remainder)
    capacities = tuple(
        hunter_cap(
            M,
            tuple(inner_ds),
            tuple(
                lows[index] + ((pattern >> index) & 1)
                for index in range(5)
            ),
        )
        for pattern in range(32)
    )
    require(
        tuple(marginals) == (120, 120, 0, 120, 0)
        and tuple(sorted(set(capacities))) == (2, 3, 4, 7)
        and capacities[0] == 2,
        ("status atlas changed", marginals, capacities),
    )

    actual_L, ranges = U.support.safe_cell_ranges(BODY)
    require(actual_L == L, "safe-cell ruler changed")
    arcs = U.fibre.projected_support_arcs(D, ranges)
    histogram = U.fibre.residue_load_histogram(arcs, Q)
    require(
        histogram == ((0, 120), (3, 420), (4, 112), (5, 168), (6, 20)),
        ("target load histogram changed", histogram),
    )
    heavy_three = sum(count for load, count in histogram if load >= 3)
    exception_events = sum(marginals)

    # Farkas certificate for the threshold-three row.  With alpha=1/360 and
    # z_total=0,z_i=1/360, every one of the 32 status-cell slacks is exact.
    alpha = F(1, 360)
    z = (F(0),) + (F(1, 360),) * 5
    slacks = []
    for pattern, capacity in enumerate(capacities):
        good = int(capacity >= 3)
        equality_column = (1,) + tuple(
            (pattern >> index) & 1 for index in range(5)
        )
        slack = sum(
            z[index] * equality_column[index] for index in range(6)
        ) - alpha * good
        slacks.append(slack)
    contradiction = (
        z[0] * Q
        + sum(z[index + 1] * marginals[index] for index in range(5))
        - alpha * heavy_three
    )
    require(
        all(slack >= 0 for slack in slacks)
        and contradiction == -1
        and heavy_three == 720
        and exception_events == 360,
        ("exact Farkas audit changed", slacks, contradiction),
    )
    require(
        all(
            int(capacities[pattern] >= 3)
            <= sum((pattern >> index) & 1 for index in range(5))
            for pattern in range(32)
        ),
        "transparent status implication failed",
    )

    semantic_payload = (
        BODY,
        FIRST,
        LABELS,
        DS,
        h,
        len(carrier),
        L,
        lower,
        deltas,
        upper,
        gap,
        tuple(ray_maxima),
        tuple(inner_ds),
        tuple(lows),
        tuple(marginals),
        capacities,
        histogram,
        heavy_three,
        exception_events,
        tuple(slacks),
        contradiction,
        TREES5,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            "independent semantic digest changed",
        )

    lines = [
        "LRC14 projected k=2 z1=1932 independent fixed-state referee",
        f"uniform_ray_status_source_sha256={file_sha256(UNIFORM)}",
        "primary_z1932_source_imported=False;tree_engine=Pruefer;LP_used=False",
        (
            f"body={BODY};first={FIRST};labels={LABELS};ds={DS};"
            f"lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
        ),
        f"all_label_ray_maxima={tuple(ray_maxima)}",
        (
            f"D={D};q={Q};M={M};marginals={tuple(marginals)};"
            f"capacity_values={tuple(sorted(set(capacities)))};"
            f"histogram={histogram}"
        ),
        (
            f"transparent_threshold3=heavy_fibres:{heavy_three};"
            f"exception_events:{exception_events};720>360"
        ),
        (
            f"farkas=alpha:1/360,z:(0,1/360,1/360,1/360,1/360,1/360);"
            f"min_cell_slack={ftext(min(slacks))};"
            f"contradiction={ftext(contradiction)}"
        ),
        f"Pruefer_K5_trees={len(TREES5)};status_cells={len(capacities)}",
        "conclusion=fixed projected z1=1932 scalar state is impossible",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    OUTPUT_PATH.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
