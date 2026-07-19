#!/usr/bin/env python3
"""Discover/replay the separate THM-1197 BV extension certificate.

The original THM-1182 certificate for ``1 <= a <= 30`` is immutable.  This
script imports its exact clipped-tooth integrator and LP discovery routine,
but accepts and emits only the new carrier block ``31 <= a <= 40``.  Discovery
is resumable after every phase row; the floating LP merely proposes integer
weights, and every proposal is checked with ``Fraction`` arithmetic before it
is written.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "04-computation/lrc14_slow_gap_bv_needle_thm1182.py"
DEFAULT_CERTIFICATE = ROOT / (
    "05-knowledge/results/"
    "lrc14_slow_gap_bv_needle_extension_thm1197_certificate.json"
)
SCHEMA = "lrc14-slow-gap-bv-needle-extension-thm1197-v1"
MIN_A = 31
MAX_A = 40
MULTIPLIERS = (12, 16, 20, 28, 44, 60, 76, 108)


def load_base():
    spec = importlib.util.spec_from_file_location("thm1182_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load THM-1182 generator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def expected_keys() -> set[tuple[int, int]]:
    return {
        (a, k)
        for a in range(MIN_A, MAX_A + 1)
        for k in range((a + 1) // 2)
    }


def verify_row(base, row: dict) -> tuple[Fraction, Fraction, Fraction]:
    if (
        row["bins"] != 200
        or len(row["weights"]) != 200
        or row["multiplier"] not in MULTIPLIERS
    ):
        raise RuntimeError(("row shape", row.get("a"), row.get("k")))
    low, tail, variation = base.verify_data(
        row["a"],
        row["k"],
        row["multiplier"],
        row["bins"],
        row["weights"],
    )
    stored = (
        Fraction(row["low_num"], row["low_den"]),
        Fraction(row["tail_num"], row["tail_den"]),
        Fraction(row["variation_num"], row["variation_den"]),
    )
    if (low, tail, variation) != stored or sum(row["weights"]) != row["total"]:
        raise RuntimeError(("stored invariant", row["a"], row["k"]))
    return low, tail, variation


def verify_document(base, document: dict) -> None:
    if (
        document.get("schema") != SCHEMA
        or document.get("min_a") != MIN_A
        or document.get("max_a") != MAX_A
    ):
        raise RuntimeError("extension certificate schema/range mismatch")
    keys = [(row["a"], row["k"]) for row in document["rows"]]
    if set(keys) != expected_keys() or len(keys) != len(set(keys)):
        raise RuntimeError("extension phase coverage mismatch")
    for row in document["rows"]:
        verify_row(base, row)


def make_row(base, a: int, k: int, bins: int) -> dict:
    errors: list[str] = []
    for multiplier in MULTIPLIERS:
        try:
            rows, weights = base.discover(a, k, multiplier, bins)
            low, tail, variation = base.verify_data(
                a, k, multiplier, bins, weights, rows
            )
            return {
                "a": a,
                "k": k,
                "multiplier": multiplier,
                "bins": bins,
                "weights": weights,
                "total": sum(weights),
                "low_num": low.numerator,
                "low_den": low.denominator,
                "tail_num": tail.numerator,
                "tail_den": tail.denominator,
                "variation_num": variation.numerator,
                "variation_den": variation.denominator,
            }
        except (RuntimeError, ValueError) as error:
            errors.append(repr(error))
    raise RuntimeError((a, k, errors))


def write_document(path: Path, rows: list[dict]) -> None:
    document = {
        "schema": SCHEMA,
        "min_a": MIN_A,
        "max_a": MAX_A,
        "reflection": "k maps to a-1-k modulo a",
        "analytic_tail": "mu(D_b)<=1/7+3*TV(f)/(49b)",
        "rows": sorted(rows, key=lambda row: (row["a"], row["k"])),
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, separators=(",", ":")) + "\n")
    temporary.replace(path)


def print_summary(document: dict) -> None:
    rows = document["rows"]
    low = max(Fraction(row["low_num"], row["low_den"]) for row in rows)
    tail = max(Fraction(row["tail_num"], row["tail_den"]) for row in rows)
    count = Counter(row["multiplier"] for row in rows)
    print(
        "VERIFIED_EXACT",
        f"rows={len(rows)}",
        f"carriers={MIN_A}..{MAX_A}",
        "multipliers=" + ",".join(f"{key}:{count[key]}" for key in sorted(count)),
    )
    print(f"max_low={low} gap_to_1/6={Fraction(1, 6) - low}")
    print(f"max_tail={tail} gap_to_1/6={Fraction(1, 6) - tail}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--discover", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--bins", type=int, default=200)
    args = parser.parse_args()
    base = load_base()

    if not args.discover:
        document = json.loads(args.certificate.read_text())
        verify_document(base, document)
        print_summary(document)
        return

    if args.bins != 200:
        raise RuntimeError("THM-1197 certificate schema fixes bins=200")
    rows: list[dict] = []
    if args.resume and args.certificate.exists():
        partial = json.loads(args.certificate.read_text())
        rows = partial.get("rows", [])
    by_key = {(row["a"], row["k"]): row for row in rows}
    for a in range(MIN_A, MAX_A + 1):
        for k in range((a + 1) // 2):
            if (a, k) in by_key:
                continue
            row = make_row(base, a, k, args.bins)
            by_key[a, k] = row
            write_document(args.certificate, list(by_key.values()))
            print(
                f"ROW a={a} k={k} M={row['multiplier']} "
                f"low={row['low_num'] / row['low_den']:.15f} "
                f"tail={row['tail_num'] / row['tail_den']:.15f}",
                flush=True,
            )
    document = json.loads(args.certificate.read_text())
    verify_document(base, document)
    print_summary(document)


if __name__ == "__main__":
    main()
