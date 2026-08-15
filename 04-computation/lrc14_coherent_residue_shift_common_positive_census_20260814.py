#!/usr/bin/env python3
"""Exact coherent-shift census beyond the common no-wrap chamber.

Scans every coherent shift for which every nominal q=3 drift remains
positive on every one of the 649 active bodies.  The pinned boundary engine
computes the restricted cross-K3,3 margin and, at each restricted failure,
the unrestricted K6 Hunter margin and literal union.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from multiprocessing import get_context
from pathlib import Path
import argparse
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_coherent_residue_shift_boundary_census_20260814.py"
EXPECTED_SOURCE = "5c69fbd209c332d76cc51814bc1ef126d8cc3dea0e9c9d16f00ccf095d6b9a7e"
FIRST_SHIFT = 156
EXPECTED_RESTRICTED_FAILURES = 1_696
EXPECTED_FULL_FAILURES = 269
EXPECTED_UNION_FAILURES = 187
EXPECTED_FIRST_FULL = (
    390,
    (1, 2, 3, 4, 6, 12),
    168,
    90,
    (393, 394, 402),
    (391, 392, 396),
    -717_819,
    1_739_713_360,
)
EXPECTED_FIRST_UNION = (
    392,
    (1, 2, 3, 4, 6, 12),
    168,
    90,
    (393, 394, 398),
    (395, 396, 404),
    107_944,
    12_562_795,
)
EXPECTED_SEMANTIC = "049acb7ee5d9ece50601943a2552c0644e8f42239b6dd16af28844401a7b0620"


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


require(file_sha(SOURCE) == EXPECTED_SOURCE, file_sha(SOURCE))
B = load("coherent_shift_common_positive_boundary", SOURCE)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    bodies = B.M.body_universe()
    common_last = min(3 * L - max(body) - 1 for body, L in bodies)
    require(common_last == 491, common_last)
    last = common_last
    shifts = tuple(range(FIRST_SHIFT, last + 1))
    if args.processes == 1:
        results = tuple(map(B.scan_shift, shifts))
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=get_context("fork")
        ) as pool:
            results = tuple(pool.map(B.scan_shift, shifts))
    failures = tuple(
        (shift, failure)
        for shift, tested, minimum, shift_failures, shift_semantic in results
        for failure in shift_failures
    )
    full_failures = tuple((shift, failure) for shift, failure in failures if failure[2] <= 0)
    union_failures = tuple((shift, failure) for shift, failure in failures if failure[5] >= 0)
    shift_counts = tuple(sorted(Counter(shift for shift, _ in failures).items()))
    body_counts = tuple(sorted(Counter(failure[0][1] for _, failure in failures).items()))
    semantic = sha256()
    for result in results:
        semantic.update((repr(result) + "\n").encode())
    semantic_hex = semantic.hexdigest()
    weakest_full = min(
        ((failure[2], shift, failure[0][1], failure[0][5], failure[0][6])
         for shift, failure in failures),
        default=None,
    )
    weakest_union = max(
        ((failure[5], shift, failure[0][1], failure[0][5], failure[0][6])
         for shift, failure in failures),
        default=None,
    )
    first_full = full_failures[0]
    first_union = union_failures[0]
    first_full_control = (
        first_full[0], first_full[1][0][1], first_full[1][0][2],
        first_full[1][0][3], first_full[1][0][5], first_full[1][0][6],
        first_full[1][2].numerator, first_full[1][2].denominator,
    )
    first_union_control = (
        first_union[0], first_union[1][0][1], first_union[1][0][2],
        first_union[1][0][3], first_union[1][0][5], first_union[1][0][6],
        first_union[1][5].numerator, first_union[1][5].denominator,
    )
    require(len(failures) == EXPECTED_RESTRICTED_FAILURES, len(failures))
    require(len(full_failures) == EXPECTED_FULL_FAILURES, len(full_failures))
    require(len(union_failures) == EXPECTED_UNION_FAILURES, len(union_failures))
    require(first_full_control == EXPECTED_FIRST_FULL, first_full_control)
    require(first_union_control == EXPECTED_FIRST_UNION, first_union_control)
    require(semantic_hex == EXPECTED_SEMANTIC, semantic_hex)
    print("LRC14 COHERENT SHIFT COMMON-POSITIVITY EXACT CENSUS")
    print("dependency", SOURCE.relative_to(ROOT), EXPECTED_SOURCE)
    print("shifts", (shifts[0], shifts[-1]), "packets", len(shifts) * 12_980,
          "common_positive_last", common_last)
    print("restricted_failures", len(failures), "shift_counts", shift_counts)
    print("body_counts", body_counts)
    print("full_failures", len(full_failures), "union_failures", len(union_failures))
    print("first_full_failure", full_failures[:1])
    print("first_union_failure", union_failures[:1])
    print("weakest_full", weakest_full)
    print("weakest_union", weakest_union)
    print("semantic_sha256", semantic_hex)
    print("status=FINITE-EXACT common-positive census; first full-tree and literal-union boundaries frozen")


if __name__ == "__main__":
    main()
