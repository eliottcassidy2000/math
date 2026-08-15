#!/usr/bin/env python3
"""Exact all-649 coherent-shift census on the common no-wrap chamber.

For every active upper-median body, every c=1,...,155 keeps all shifted
residues e+c strictly between 0 and L.  This wrapper reuses the pinned exact
interval engine of the first-boundary census, retains its restricted K3,3
diagnostic, and requires an unrestricted K6 Hunter rescue at every restricted
failure.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from multiprocessing import get_context
from pathlib import Path
import argparse
import sys


ROOT = Path(__file__).resolve().parents[1]
BOUNDARY = ROOT / "04-computation/lrc14_coherent_residue_shift_boundary_census_20260814.py"
EXPECTED_BOUNDARY = "5c69fbd209c332d76cc51814bc1ef126d8cc3dea0e9c9d16f00ccf095d6b9a7e"
LAST_SHIFT = 155
EXPECTED_SHIFT_COUNTS = (
    (67, 1), (69, 1), (80, 1), (82, 1), (84, 1), (91, 1), (93, 1),
    (97, 1), (104, 1), (106, 1), (110, 1), (112, 1), (115, 1),
    (117, 1), (119, 1), (123, 1), (125, 1), (128, 1), (130, 1),
    (136, 1), (137, 1), (138, 2), (140, 1), (141, 1), (143, 1),
    (149, 1), (151, 2), (153, 1), (154, 1),
)
EXPECTED_BODY_HISTOGRAM = (
    ((1, 2, 3, 4, 6, 8), 2),
    ((1, 2, 3, 4, 6, 12), 26),
    ((1, 2, 4, 6, 8, 12), 2),
    ((1, 3, 4, 6, 8, 12), 1),
)
EXPECTED_WITHIN_EDGE_HISTOGRAM = ((3, 5), (4, 26))
EXPECTED_SEMANTIC = "16c6a68ae4f5aaba4d18416d07c20fc69b7776b14df27fb06c9e15d976bc8257"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


require(file_sha(BOUNDARY) == EXPECTED_BOUNDARY,
        (BOUNDARY, file_sha(BOUNDARY)))
B = load("coherent_shift_uniform_no_wrap_boundary", BOUNDARY)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=6)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    bodies = B.M.body_universe()
    require(len(bodies) == 649, len(bodies))
    common_last = min(L - max(body) - 1 for body, L in bodies)
    common_last_attainers = tuple(
        (body, L) for body, L in bodies
        if L - max(body) - 1 == common_last
    )
    require(common_last == LAST_SHIFT == 155, common_last)
    require(common_last_attainers == (((1, 2, 3, 4, 6, 12), 168),),
            common_last_attainers)
    require(all(0 < min(body) + 1 and max(body) + LAST_SHIFT < L
                for body, L in bodies), "common no-wrap chamber")
    shifts = tuple(range(1, LAST_SHIFT + 1))
    if args.processes == 1:
        results = tuple(map(B.scan_shift, shifts))
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=get_context("fork")
        ) as pool:
            results = tuple(pool.map(B.scan_shift, shifts))
    require(tuple(result[0] for result in results) == shifts, "shift order")
    packets = sum(result[1] for result in results)
    failures = tuple(
        (shift, failure)
        for shift, tested, minimum, shift_failures, shift_semantic in results
        for failure in shift_failures
    )
    shift_counts = tuple(sorted(Counter(shift for shift, _ in failures).items()))
    body_histogram = tuple(sorted(Counter(
        failure[0][1] for shift, failure in failures
    ).items()))
    within_histogram = Counter()
    full_failures = []
    union_failures = []
    for shift, failure in failures:
        row, full_credit, full_margin, full_tree, union, union_margin, all_edges = failure
        within = sum(1 for a, b, weight in full_tree if a[0] == b[0])
        within_histogram[within] += 1
        if full_margin <= 0:
            full_failures.append((shift, failure))
        if union_margin >= 0:
            union_failures.append((shift, failure))
    within_edge_histogram = tuple(sorted(within_histogram.items()))
    semantic_accumulator = sha256()
    for result in results:
        semantic_accumulator.update((repr(result) + "\n").encode())
    semantic = semantic_accumulator.hexdigest()
    require(packets == 2_011_900, packets)
    require(len(failures) == 31, len(failures))
    require(shift_counts == EXPECTED_SHIFT_COUNTS, shift_counts)
    require(body_histogram == EXPECTED_BODY_HISTOGRAM, body_histogram)
    require(within_edge_histogram == EXPECTED_WITHIN_EDGE_HISTOGRAM,
            within_edge_histogram)
    require(not full_failures, full_failures)
    require(not union_failures, union_failures)
    require(semantic == EXPECTED_SEMANTIC, semantic)
    weakest_cross = min(
        (failure[0][0], shift, failure[0][1], failure[0][5], failure[0][6])
        for shift, failure in failures
    )
    weakest_full = min(
        (failure[2], shift, failure[0][1], failure[0][5], failure[0][6])
        for shift, failure in failures
    )
    print("LRC14 ALL-649 COMMON NO-WRAP COHERENT-SHIFT CENSUS")
    print("dependency", BOUNDARY.relative_to(ROOT), EXPECTED_BOUNDARY)
    print("bodies", len(bodies), "assignments", len(B.ASSIGNMENTS),
          "shifts", (shifts[0], shifts[-1]), "packets", packets)
    print("common_no_wrap_last", common_last,
          "attainers", common_last_attainers)
    print("restricted_failure_count", len(failures),
          "shift_counts", shift_counts)
    print("restricted_failure_body_histogram", body_histogram)
    print("selected_full_tree_within_edge_histogram", within_edge_histogram)
    print("weakest_restricted_margin", weakest_cross)
    print("weakest_full_rescue_margin", weakest_full)
    print("full_tree_failures", len(full_failures),
          "literal_union_failures", len(union_failures))
    print("semantic_sha256", semantic)
    print("conclusion=every common-no-wrap packet closes by restricted cross-K3,3 or unrestricted full-K6 Hunter")
    print("status=FINITE-EXACT all-649 common no-wrap chamber; not an arbitrary-shift theorem")


if __name__ == "__main__":
    main()
