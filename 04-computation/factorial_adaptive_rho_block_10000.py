#!/usr/bin/env python3
"""Dual exact continuation of the factorial pair census, d=6001..10000."""

import ast
from collections import Counter, defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path


START, END = 6001, 10000
PRIMARY_NAME = "factorial_adaptive_rho_block_6000.py"
PRIMARY_SHA256 = "b65edcf2870714ca57456b8297afdd05284a09b302ec4b84d2e57829520c94d1"
INDEPENDENT_NAME = "factorial_adaptive_rho_block_6000_independent_audit.py"
INDEPENDENT_SHA256 = "d416cb2955fd745394cf1043ac8c2eba28a6a97beb264dd9cbe9919ed8c96724"
EXPECTED_COUNTS = (4000, 2524, 1476, 1364, 112, 0)
EXPECTED_EXIT_HISTOGRAM = (
    ("d_minus_1_prime_power", 453), ("d_minus_2_prime", 383),
    ("d_minus_3_prime", 382), ("d_minus_4_prime", 322),
    ("d_minus_5_prime", 322), ("d_minus_6_prime", 216),
    ("d_prime", 446),
)
EXPECTED_RHO_HISTOGRAM = (
    (3, 24), (5, 6), (7, 18), (11, 32), (13, 15),
    (17, 7), (19, 5), (23, 4), (29, 1),
)
EXPECTED_SEMANTIC_SHA256 = "d90179fdebd48dd82cd368b957c9602fbd287774287de0fecb73b4a84dca69f3"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_pinned(name, expected, module_name):
    path = Path(__file__).resolve().with_name(name)
    actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, (name, actual, expected))
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, ("bad import", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compact_row(row):
    if len(row) == 2:
        return row
    divisor = tuple((p, blocks, len(pre), len(post), post if len(post) <= 32 else ())
                    for p, blocks, pre, post in row[3])
    rho = tuple((p, status, blocks, pre, post)
                for p, status, _, _, blocks, pre, post in row[4])
    suffix = row[5:] if row[1] == "survivor" else (row[5],)
    return row[:3] + (divisor, rho) + suffix


def main():
    primary = load_pinned(PRIMARY_NAME, PRIMARY_SHA256, "rho10000_primary")
    independent = load_pinned(INDEPENDENT_NAME, INDEPENDENT_SHA256,
                              "rho10000_independent")
    primary.START, primary.END = START, END
    independent.START, independent.END = START, END
    divisor, rho_module, _, coefficient, digital = primary.load_engines()
    primary_semantic = primary.build_semantic(divisor, rho_module, coefficient, digital)
    independent_semantic, _ = independent.build_semantic()
    require(primary_semantic == independent_semantic,
            (digest(primary_semantic), digest(independent_semantic)))
    semantic = primary_semantic
    exits, divisor_rows, rho_rows, survivors, packets, killers, inadmissible = independent.summarize(semantic)
    counts = (END - START + 1, len(exits), len(semantic[2]), len(divisor_rows),
              len(rho_rows), len(survivors))
    exit_histogram = tuple(sorted(Counter(label for _, label in exits).items()))
    killer_histogram = tuple(sorted(Counter(p for _, p in killers).items()))
    first_by_killer = {}
    residue_histogram = defaultdict(Counter)
    for d, p in killers:
        first_by_killer.setdefault(p, d)
        residue_histogram[p][d % p] += 1
    first_residual = semantic[2][0] if semantic[2] else ()
    first_rho = rho_rows[0] if rho_rows else ()
    d6001 = ((6001, independent.exit_label(6001))
             if independent.exit_label(6001) else compact_row(semantic[2][0]))
    require(d6001[0] == 6001, d6001)
    require(counts == EXPECTED_COUNTS, counts)
    require(exit_histogram == EXPECTED_EXIT_HISTOGRAM, exit_histogram)
    require(killer_histogram == EXPECTED_RHO_HISTOGRAM, killer_histogram)
    require(len(inadmissible) == 61 and not survivors, (len(inadmissible), survivors))
    require(first_residual[:3] == (6001, "rho", 11), first_residual[:3])
    require(first_rho[:3] == (6001, "rho", 11), first_rho[:3])
    require(first_by_killer[29] == 6518, first_by_killer)
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256, digest(semantic))
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(Path(__file__).read_text()))),
            "assert node")

    print("FACTORIAL ADAPTIVE RHO BLOCK d=6001..10000 DUAL EXACT CENSUS")
    print("dependencies=%s:%s;%s:%s" %
          (PRIMARY_NAME, PRIMARY_SHA256, INDEPENDENT_NAME, INDEPENDENT_SHA256))
    print("implementation=pinned THM-3475/3483 engine route versus self-contained repeated-factor/hull/DP/rho route; exact semantic equality")
    print("universe=d=6001..10000; canonical seven exits; F=A_(d-2)^d,G=A_(d-1)^d; rho bank=%s" % (primary.PRIMES_101,))
    print("counts=%s exit_histogram=%s" % (counts, exit_histogram))
    print("rho_needed_divisor_packets=%s" % (packets,))
    print("rho_killers=%s histogram=%s first_by_killer=%s" %
          (killers, killer_histogram, tuple(sorted(first_by_killer.items()))))
    print("killer_residue_histogram=%s" %
          (tuple((p, tuple(sorted(hist.items()))) for p, hist in sorted(residue_histogram.items())),))
    print("inadmissible_count=%d survivors=%s" %
          (len(inadmissible), tuple(compact_row(r) for r in survivors)))
    print("first_residual=%s" % (compact_row(first_residual),))
    print("first_rho_row=%s" % (compact_row(first_rho),))
    print("d6001=%s" % (d6001,))
    print("semantic_sha256=%s cross_record_equal=True" % EXPECTED_SEMANTIC_SHA256)
    print("boundary_if_no_survivor=r<=9998; otherwise r<=first_survivor-3")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
