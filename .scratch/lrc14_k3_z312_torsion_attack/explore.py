#!/usr/bin/env python3
"""Explore finite-low-pair / one-high-ray certificates at k3,z1=312."""

from __future__ import annotations

import importlib.util
import argparse
import math
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


ROOT = Path("/Users/e/Documents/GitHub/math")
FRONT = Path(
    "/tmp/math-wt-k3-z324-scout/04-computation/"
    "lrc14_j7_k3_z312_ray_status_frontier_thm2941.py"
)
PROJECTION = ROOT / "04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


front = load("front312_explore", FRONT)
ray = front.ray
proj = load("projection312_explore", PROJECTION)


def first_on_ray(r, L, threshold):
    return r if r >= threshold else r + ((threshold - r + L - 1) // L) * L


def label_rows(stream, d):
    low = []
    high = []
    step = stream.L // d
    for unit in range(1, d):
        if math.gcd(unit, d) != 1:
            continue
        residue = step * unit
        amplitude = residue * ray.local.delta(stream.carrier, stream.h, residue)
        for label in range(residue, stream.high_floor, stream.L):
            if label > stream.first:
                low.append((amplitude / label, label, unit))
        high_label = first_on_ray(residue, stream.L, stream.high_floor)
        high_upper = amplitude / high_label if amplitude > 0 else F(0)
        high.append((high_upper, high_label, unit, amplitude))
    return tuple(sorted(low, reverse=True)), tuple(sorted(high, reverse=True))


def low_pair_rows(first, second, *, same, threshold):
    if not same:
        for a in first:
            if not second or a[0] + second[0][0] < threshold:
                break
            for b in second:
                if a[0] + b[0] < threshold:
                    break
                if a[1] != b[1]:
                    yield (a, b)
        return
    for i, a in enumerate(first):
        if i + 1 >= len(first) or a[0] + first[i + 1][0] < threshold:
            break
        for b in first[i + 1 :]:
            if a[0] + b[0] < threshold:
                break
            if a[1] != b[1]:
                yield (a, b)


def safe_cells(cells, L, labels):
    out = []
    for cell in cells:
        okay = True
        for label in labels:
            residue = (label * cell) % L
            # The phase segment [residue,residue+label]/L misses every open
            # radius-1/14 tooth iff it lies in the closed safe middle arc.
            if 14 * residue < L or 14 * (residue + label) > 13 * L:
                okay = False
                break
        if okay:
            out.append(cell)
    return tuple(out)


def pair_certificate(safe, d, unit):
    if len(safe) < 2:
        return None
    j = safe[0]
    for k in safe[1:]:
        residue = (unit * (k - j)) % d
        if 7 * min(residue, d - residue) >= d:
            return j, k, residue
    # The first anchor could be exceptional; do a bounded/full quadratic
    # fallback only when needed.
    for index, j in enumerate(safe[1:], 1):
        for k in safe[index + 1 :]:
            residue = (unit * (k - j)) % d
            if 7 * min(residue, d - residue) >= d:
                return j, k, residue
    return None


def uniform_torsion_certificate(safe, d):
    """Pair whose phase gap is >=1/7 for every unit modulo d."""
    by_residue = defaultdict(list)
    for cell in safe:
        by_residue[cell % d].append(cell)
    for p in (2, 3, 5, 7):
        if d % p:
            continue
        shift = d // p
        for residue, js in by_residue.items():
            ks = by_residue.get((residue + shift) % d)
            if ks:
                return js[0], ks[0], shift, p
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--body-index", type=int)
    parser.add_argument("--dump", action="store_true")
    args = parser.parse_args()
    grand_families = 0
    grand_low_pairs = set()
    grand_fail = []
    body_rows = list(front.EXPECTED_RESIDUALS.items())
    if args.body_index is not None:
        body_rows = [body_rows[args.body_index]]
    for body, residuals in body_rows:
        stream = ray.Stream(body)
        required = stream.lower - stream.first_delta
        cells = proj.body_cells(stream.carrier, stream.L)
        tables = {}
        all_d = sorted({d for ds in residuals for d in ds})
        for d in all_d:
            tables[d] = label_rows(stream, d)

        # Deliberately duplicate-permitting upper bound for two high slots.
        two_high_upper = F(-10**9)
        two_high_arg = None
        for ds in residuals:
            suffix = list(ds)
            suffix.remove(stream.first_d)
            for high_indices in ((0, 1), (0, 2), (1, 2), (0, 1, 2)):
                value = F(0)
                for index, d in enumerate(suffix):
                    rows = tables[d][1 if index in high_indices else 0]
                    value += rows[0][0] if rows else F(-10**9)
                if value > two_high_upper:
                    two_high_upper = value
                    two_high_arg = ds, tuple(suffix), high_indices

        family_rows = []
        pair_keys = set()
        fail = []
        safe_cache = {}
        cert_cache = {}
        for ds in residuals:
            suffix = list(ds)
            suffix.remove(stream.first_d)
            for d_high in sorted(set(suffix)):
                low_ds = list(suffix)
                low_ds.remove(d_high)
                d_a, d_b = sorted(low_ds)
                high_max = tables[d_high][1][0][0]
                low_rows = low_pair_rows(
                    tables[d_a][0], tables[d_b][0], same=d_a == d_b,
                    threshold=required - high_max,
                )
                # A high unit ray is relevant only if some low pair reaches
                # the scalar wall against that ray's exact maximum/supremum.
                for low_a, low_b in low_rows:
                    base = low_a[0] + low_b[0]
                    low_labels = tuple(sorted((low_a[1], low_b[1])))
                    family = (body, ds, d_high, low_labels)
                    family_rows.append(family)
                    pair_key = (body, low_labels)
                    pair_keys.add(pair_key)
                    if pair_key not in safe_cache:
                        safe_cache[pair_key] = safe_cells(
                            cells, stream.L, (stream.first, *low_labels)
                        )
                    cert_key = (pair_key, d_high)
                    if cert_key not in cert_cache:
                        cert_cache[cert_key] = uniform_torsion_certificate(
                            safe_cache[pair_key], d_high
                        )
                    if cert_cache[cert_key] is None:
                        fail.append(family)

        grand_families += len(family_rows)
        grand_low_pairs.update(pair_keys)
        grand_fail.extend(fail)
        certs = [value for value in cert_cache.values() if value is not None]
        case_rows = tuple(
            sorted(
                (*family, *cert_cache[((body, family[-1]), family[2])])
                for family in set(family_rows)
            )
        )
        print(
            "body", body,
            "two_high_gap", required - two_high_upper,
            "two_high_arg", two_high_arg,
            "families", len(family_rows),
            "low_pairs", len(pair_keys),
            "failed", len(fail),
            "safe_cell_min", min(map(len, safe_cache.values()), default=None),
            "safe_cell_max", max(map(len, safe_cache.values()), default=None),
            "certificate_sample", certs[:3],
            flush=True,
        )
        print(
            "  d_high_histogram", dict(sorted(Counter(row[2] for row in family_rows).items())),
            "torsion_p_histogram", dict(sorted(Counter(row[-1] for row in case_rows).items())),
            "unique_cases", len(case_rows),
            flush=True,
        )
        if args.dump:
            for row in case_rows:
                print("  CASE", row, flush=True)
        if fail:
            print("  FAIL SAMPLE", fail[:20], flush=True)
    print(
        "TOTAL families", grand_families,
        "body_low_pairs", len(grand_low_pairs),
        "failed", len(grand_fail),
        flush=True,
    )


if __name__ == "__main__":
    main()
