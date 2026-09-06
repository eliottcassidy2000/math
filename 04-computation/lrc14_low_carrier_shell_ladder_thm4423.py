#!/usr/bin/env python3
"""Deterministic exact low-carrier shell ladder for THM-4423.

One pass over each shell collects every nonautomatic even carrier count from
92 through MAX_COUNT. Each collected row is reconstructed independently by
the interval sweep and the direct relation-lattice box before its exact raw
projection envelopes are compared with 6/77. The only imported computation
and frozen dependency are the canonical THM-4423 height-499 atlas artifacts.
"""

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import json


ROOT = Path(__file__).resolve().parents[1]
SCANNER_PATH = ROOT / "04-computation" / "lrc14_height499_projection_atlas_thm4423.py"
FROZEN = ROOT / "05-knowledge" / "results" / "lrc14_height499_projection_atlas_thm4423.out"
FROZEN_SHA256 = "f47c6c89a4fbadf9df66b2d169855175bfda53b26ba2e0194904ba99a9f63a18"
ROW_SHA256 = "a6d482622d2c96a9981239c66aa40aeb21cd570e46fd1cd6d5239d831eb72144"
TARGET = Fraction(6, 77)
BASE_COUNT = 90
_SCANNER = None
CHECKS = 0


def require(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(payload)


def eligible(value):
    return value > 0 and value % 2 == 1 and value % 3 != 0


def nonautomatic_shells(count, census_height=499):
    """Eligible c>census_height satisfying count>2c/11."""
    stop = (11 * count + 1) // 2
    return tuple(c for c in range(census_height + 1, stop) if eligible(c))


def relation_type(w):
    a, b, c = w
    labels = []
    if c == 2 * a + b:
        labels.append("c=2a+b")
    if c == 2 * b - a:
        labels.append("c=2b-a")
    if c == a + 2 * b:
        labels.append("c=a+2b")
    return "+".join(labels) if labels else "other"


def primitive_ray(carrier):
    divisor = gcd(gcd(abs(carrier[0]), abs(carrier[1])), abs(carrier[2]))
    primitive = tuple(value // divisor for value in carrier)
    for value in primitive:
        if value:
            if value < 0:
                primitive = tuple(-entry for entry in primitive)
            break
    return primitive, divisor


def carrier_geometry(carriers):
    rays = {}
    for carrier in carriers:
        ray, multiplier = primitive_ray(carrier)
        rays.setdefault(ray, set()).add(multiplier)
    geometry = tuple((ray, tuple(sorted(multipliers)))
                     for ray, multipliers in sorted(rays.items()))
    minimum = min(sum(abs(value) for value in ray) for ray in rays)
    minimizers = tuple(ray for ray in sorted(rays)
                       if sum(abs(value) for value in ray) == minimum)
    return geometry, minimum, minimizers


def scanner_module():
    global _SCANNER
    if _SCANNER is None:
        spec = spec_from_file_location("lrc14_height499_projection_atlas_thm4423",
                                       SCANNER_PATH)
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        _SCANNER = module
    return _SCANNER


def fraction_pair(value):
    return value.numerator, value.denominator


def scan_shell(payload):
    global CHECKS
    c, max_count = payload
    scanner = scanner_module()
    local_checks_start = CHECKS
    scanner_checks_start = scanner.CHECKS
    values = tuple(v for v in range(1, c, 2) if v % 3)
    target_counts = set(range(BASE_COUNT + 2, max_count + 1, 2))
    triples = 0
    carrier_min = None
    carrier_max = 0
    rows = []

    for a, b in combinations(values, 2):
        w = (a, b, c)
        if gcd(gcd(a, b), c) != 1:
            continue
        triples += 1
        carriers = scanner.carrier_sweep(w)
        count = len(carriers)
        require(count % 2 == 0, ("odd-carrier-count", w, count))
        carrier_min = count if carrier_min is None else min(carrier_min, count)
        carrier_max = max(carrier_max, count)

        # Automatic rows need no finite projection computation. Retain only
        # target counts lying strictly above the deficit threshold.
        if count not in target_counts or 11 * count <= 2 * c:
            continue

        direct = scanner.carrier_lattice_box(w)
        require(set(carriers) == direct,
                ("sweep-vs-lattice", w,
                 sorted(set(carriers).symmetric_difference(direct))))
        numerators, _, denominator = scanner.edge_envelopes(w, carriers)
        sums = tuple(Fraction(numerator, denominator)
                     for numerator in numerators)
        best = min(sums)
        geometry, relation_minimum, relation_minimizers = carrier_geometry(carriers)
        rows.append({
            "count": count,
            "w": w,
            "best": fraction_pair(best),
            "sums": tuple(fraction_pair(value) for value in sums),
            "relation": relation_type(w),
            "geometry": geometry,
            "relation_minimum_l1": relation_minimum,
            "relation_minimizers": relation_minimizers,
            "hostile": best > TARGET,
        })

    rows.sort(key=lambda row: (row["count"], row["w"]))
    return {
        "c": c,
        "triples": triples,
        "carrier_range": (carrier_min, carrier_max),
        "rows": rows,
        "checks": ((CHECKS - local_checks_start)
                   + (scanner.CHECKS - scanner_checks_start)),
    }


def frozen_audit():
    raw = FROZEN.read_bytes()
    digest = sha256(raw).hexdigest()
    require(digest == FROZEN_SHA256, ("frozen-sha256", digest))
    text = raw.decode("utf-8")
    for needle in (
        "status=PASS FINITE_EXACT",
        "max_speed<=499",
        "triples=753853",
        "hostile_count=0",
        "rows_sha256=" + ROW_SHA256,
        "verdict=PASS",
    ):
        require(needle in text, ("frozen-field", needle))
    return digest


def as_fraction(pair):
    return Fraction(pair[0], pair[1])


def canonical_rows(shell_results):
    return [
        {
            "c": shell["c"],
            "count": row["count"],
            "w": list(row["w"]),
            "best": list(row["best"]),
            "sums": [list(pair) for pair in row["sums"]],
            "relation": row["relation"],
            "geometry": [[list(ray), list(multipliers)]
                         for ray, multipliers in row["geometry"]],
            "relation_minimum_l1": row["relation_minimum_l1"],
            "relation_minimizers": [list(ray)
                                    for ray in row["relation_minimizers"]],
            "hostile": row["hostile"],
        }
        for shell in shell_results
        for row in shell["rows"]
    ]


def render(max_count, shell_results, frozen_digest, checks):
    lines = []
    lines.append("LRC14 LOW-CARRIER EXACT SHELL LADDER THM-4423 VERIFIER")
    lines.append("status=PASS FINITE_EXACT_RELATIVE; universal_projection=OPEN; LRC14=OPEN")
    lines.append("frozen_h499_sha256=%s rows_sha256=%s" %
                 (frozen_digest, ROW_SHA256))
    lines.append("universe=eligible shells above 499; counts=92..%d even; both_kernels=YES" %
                 max_count)
    lines.append("shells=%s" % (nonautomatic_shells(max_count),))

    for shell in shell_results:
        by_count = {}
        for row in shell["rows"]:
            by_count[row["count"]] = by_count.get(row["count"], 0) + 1
        lines.append("shell=%d triples=%d carrier_range=%s relevant_counts=%s" %
                     (shell["c"], shell["triples"], shell["carrier_range"],
                      dict(sorted(by_count.items()))))

    all_rows = canonical_rows(shell_results)
    semantic = json.dumps(all_rows, sort_keys=True, separators=(",", ":"))
    lines.append("ladder_rows_sha256=%s row_count=%d" %
                 (sha256(semantic.encode("ascii")).hexdigest(), len(all_rows)))
    lines.append("checks=%d" % checks)

    first_hostile = None
    non_norm_four = []
    maximal_proved = BASE_COUNT
    for count in range(BASE_COUNT + 2, max_count + 1, 2):
        shells = set(nonautomatic_shells(count))
        rows = [row for shell in shell_results if shell["c"] in shells
                for row in shell["rows"] if row["count"] == count]
        hostiles = [row for row in rows if row["hostile"]]
        shell_counts = {
            c: sum(row["count"] == count for row in shell["rows"])
            for c in sorted(shells)
            for shell in shell_results if shell["c"] == c
        }
        relation_counts = {}
        for row in rows:
            relation_counts[row["relation"]] = relation_counts.get(row["relation"], 0) + 1
            if row["relation"] == "other":
                non_norm_four.append((count, row))

        if rows:
            leader = max(rows, key=lambda row: (as_fraction(row["best"]), row["w"]))
            leader_best = as_fraction(leader["best"])
            margin = TARGET - leader_best
            leader_text = "w:%s,best:%s,margin:%s,sums:%s,relation:%s" % (
                leader["w"], leader_best, margin,
                tuple(str(as_fraction(pair)) for pair in reversed(leader["sums"])),
                leader["relation"])
        else:
            leader_text = "NONE"

        lines.append("count=%d shell_counts=%s rows=%d relations=%s hostile=%d" %
                     (count, shell_counts, len(rows),
                      dict(sorted(relation_counts.items())), len(hostiles)))
        lines.append("  extremizer=%s" % leader_text)
        if hostiles and first_hostile is None:
            first_hostile = (count, hostiles[0])
            break
        maximal_proved = count

    non_norm_four.sort(key=lambda item: (item[0], item[1]["w"][2], item[1]["w"]))
    if non_norm_four:
        count, row = non_norm_four[0]
        lines.append("first_non_norm_four=count:%d,c:%d,w:%s" %
                     (count, row["w"][2], row["w"]))
        lines.append("  sums_S1_S2_S3=%s best=%s margin=%s" %
                     (tuple(str(as_fraction(pair))
                            for pair in reversed(row["sums"])),
                      as_fraction(row["best"]), TARGET - as_fraction(row["best"])))
        lines.append("  relation_minimum_l1=%d minimizers=%s primitive_ray_layers=%s" %
                     (row["relation_minimum_l1"], row["relation_minimizers"],
                      row["geometry"]))
        lines.append("maximal_norm_four_only_count=%d" % (count - 2))
    else:
        lines.append("first_non_norm_four=NONE")
        lines.append("maximal_norm_four_only_count=%d" % maximal_proved)

    if first_hostile is None:
        next_count = maximal_proved + 2
        next_shells = nonautomatic_shells(next_count)
        lines.append("maximal_proved_count=%d" % maximal_proved)
        lines.append("next_frontier=count:%d,shells:%s,first_slice:(%d,%d)" %
                     (next_count, next_shells, next_count, next_shells[0]))
        lines.append("first_hostile=NONE")
        lines.append("verdict=PASS")
    else:
        count, row = first_hostile
        lines.append("maximal_proved_count=%d" % maximal_proved)
        lines.append("first_hostile=count:%d,w:%s,best:%s,sums:%s" %
                     (count, row["w"], as_fraction(row["best"]),
                      tuple(str(as_fraction(pair)) for pair in row["sums"])))
        lines.append("verdict=HOSTILE_FOUND")
    return "\n".join(lines) + "\n", first_hostile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-count", type=int, default=112)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.max_count >= 92 and args.max_count % 2 == 0,
            ("max-count", args.max_count))
    require(args.workers >= 1, ("workers", args.workers))

    frozen_digest = frozen_audit()
    shells = nonautomatic_shells(args.max_count)
    payloads = tuple((c, args.max_count) for c in shells)
    if args.workers == 1:
        shell_results = [scan_shell(payload) for payload in payloads]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            shell_results = list(executor.map(scan_shell, payloads))
    shell_results.sort(key=lambda shell: shell["c"])

    checks = CHECKS + sum(shell["checks"] for shell in shell_results)
    output, first_hostile = render(args.max_count, shell_results, frozen_digest,
                                   checks)
    print(output, end="")
    if args.output is not None:
        args.output.write_text(output, encoding="utf-8", newline="\n")
    if first_hostile is not None:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
