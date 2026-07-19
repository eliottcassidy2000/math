#!/usr/bin/env python3
"""Discover or exactly check the THM-1182 BV needle certificates.

The floating-point LP is discovery-only.  Every emitted row is immediately
checked with ``Fraction`` arithmetic.  With a probability density ``f`` of
global total variation ``V``, the paper supplier used by THM-1182 is

    mu(D_b) <= 1/7 + 3 V/(49 b).

Thus speeds ``b > B`` are harmless if the right side at ``B+1`` is below
``1/6``; speeds ``a < b <= B`` are checked against the rational step density.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


DEFAULT_CERTIFICATE = Path(
    "05-knowledge/results/lrc14_slow_gap_bv_needle_thm1182_certificate.json"
)
SCHEMA = "lrc14-slow-gap-bv-needle-thm1182-v1"


def gap(a: int, k: int) -> tuple[Fraction, Fraction]:
    return Fraction(14 * k + 1, 14 * a), Fraction(14 * k + 13, 14 * a)


def teeth(a: int, k: int, b: int) -> list[tuple[Fraction, Fraction]]:
    left, right = gap(a, k)
    x, y = b * left, b * right
    out: list[tuple[Fraction, Fraction]] = []
    for q in range(x.numerator // x.denominator - 2,
                   y.numerator // y.denominator + 4):
        lo = max(left, Fraction(14 * q - 1, 14 * b))
        hi = min(right, Fraction(14 * q + 1, 14 * b))
        if lo < hi:
            out.append((lo, hi))
    return out


def overlap_row(a: int, k: int, b: int, bins: int) -> list[Fraction]:
    left, right = gap(a, k)
    width = (right - left) / bins
    row = [Fraction() for _ in range(bins)]
    for lo, hi in teeth(a, k, b):
        first = max(0, min(bins - 1, int((lo - left) / width)))
        last = max(0, min(bins - 1, int((hi - left) / width)))
        for j in range(first, last + 1):
            u, v = left + j * width, left + (j + 1) * width
            row[j] += max(Fraction(), min(v, hi) - max(u, lo)) / width
    return row


def exact_rows(a: int, k: int, cutoff: int, bins: int) -> list[list[Fraction]]:
    return [overlap_row(a, k, b, bins) for b in range(a + 1, cutoff + 1)]


def discover(a: int, k: int, multiplier: int, bins: int = 200,
             scale: int = 10**12):
    cutoff = multiplier * a
    rational_rows = exact_rows(a, k, cutoff, bins)
    matrix = np.array([[float(x) for x in row] for row in rational_rows])

    # Variables are bin masses m[0:N], variation majorants z[0:N+1], and
    # the maximum low-frequency load c.
    variables = 2 * bins + 2
    rows: list[np.ndarray] = []
    upper: list[float] = []
    for values in matrix:
        row = np.zeros(variables)
        row[:bins], row[-1] = values, -1
        rows.append(row)
        upper.append(0.0)

    row = np.zeros(variables)
    row[0], row[bins] = 1, -1
    rows.append(row)
    upper.append(0.0)
    for j in range(1, bins):
        row = np.zeros(variables)
        row[j], row[j - 1], row[bins + j] = 1, -1, -1
        rows.append(row)
        upper.append(0.0)
        row = np.zeros(variables)
        row[j - 1], row[j], row[bins + j] = 1, -1, -1
        rows.append(row)
        upper.append(0.0)
    row = np.zeros(variables)
    row[bins - 1], row[2 * bins] = 1, -1
    rows.append(row)
    upper.append(0.0)

    # Use 999/1000 of the exact tail threshold, leaving room to rationalize.
    bin_width = Fraction(6, 7 * a * bins)
    variation_limit = Fraction(999, 1000) * Fraction(7 * (cutoff + 1), 18)
    row = np.zeros(variables)
    row[bins:2 * bins + 1] = 1
    rows.append(row)
    upper.append(float(variation_limit * bin_width))

    objective = np.zeros(variables)
    objective[-1] = 1
    equality = np.zeros((1, variables))
    equality[0, :bins] = 1
    solution = linprog(
        objective,
        A_ub=np.array(rows),
        b_ub=np.array(upper),
        A_eq=equality,
        b_eq=[1.0],
        bounds=[(0, None)] * variables,
        method="highs",
        options={"time_limit": 120},
    )
    if not solution.success:
        raise RuntimeError((a, k, multiplier, bins, solution.message))
    if solution.fun >= 1 / 6 - 1e-8:
        raise ValueError((a, k, multiplier, bins, solution.fun))
    weights = np.floor(solution.x[:bins] * scale).astype(np.int64)
    if weights.sum() == 0:
        raise RuntimeError("rounding erased the density")
    return rational_rows, [int(x) for x in weights]


def verify_data(a: int, k: int, multiplier: int, bins: int,
                weights: list[int], rows=None):
    cutoff = multiplier * a
    rows = exact_rows(a, k, cutoff, bins) if rows is None else rows
    total = sum(weights)
    if not (total > 0 and all(x >= 0 for x in weights)):
        raise RuntimeError((a, k, multiplier, "invalid weights"))
    loads = [
        sum(Fraction(w) * q for w, q in zip(weights, row)) / total
        for row in rows
    ]
    low = max(loads, default=Fraction())
    variation_mass = (
        weights[0]
        + sum(abs(y - x) for x, y in zip(weights, weights[1:]))
        + weights[-1]
    )
    bin_width = Fraction(6, 7 * a * bins)
    variation = Fraction(variation_mass, total) / bin_width
    tail = Fraction(1, 7) + Fraction(3, 49 * (cutoff + 1)) * variation
    if low >= Fraction(1, 6):
        raise RuntimeError((a, k, multiplier, "low", low))
    if tail >= Fraction(1, 6):
        raise RuntimeError((a, k, multiplier, "tail", tail, variation))
    return low, tail, variation


def make_row(a: int, k: int, base: int, bins: int) -> dict:
    errors: list[str] = []
    for multiplier in (base, base + 4, base + 8, base + 16, base + 32):
        try:
            rows, weights = discover(a, k, multiplier, bins)
            low, tail, variation = verify_data(
                a, k, multiplier, bins, weights, rows
            )
            break
        except (ValueError, RuntimeError) as error:
            errors.append(repr(error))
    else:
        raise RuntimeError((a, k, errors))
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


def verify_row(row: dict) -> None:
    low, tail, variation = verify_data(
        row["a"], row["k"], row["multiplier"], row["bins"], row["weights"]
    )
    stored = (
        Fraction(row["low_num"], row["low_den"]),
        Fraction(row["tail_num"], row["tail_den"]),
        Fraction(row["variation_num"], row["variation_den"]),
    )
    if (low, tail, variation) != stored:
        raise RuntimeError((row["a"], row["k"], "stored invariant mismatch"))


def verify_document(document: dict) -> None:
    if document.get("schema") != SCHEMA:
        raise RuntimeError("certificate schema mismatch")
    expected = {
        (a, k)
        for a in range(1, document["max_a"] + 1)
        for k in range((a + 1) // 2)
    }
    actual = {(row["a"], row["k"]) for row in document["rows"]}
    if actual != expected or len(actual) != len(document["rows"]):
        raise RuntimeError("phase coverage mismatch")
    for row in document["rows"]:
        verify_row(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=30)
    parser.add_argument("--base", type=int, default=12)
    parser.add_argument("--bins", type=int, default=200)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--verify-certificate", action="store_true")
    args = parser.parse_args()

    if args.verify_certificate:
        document = json.loads(args.certificate.read_text())
        verify_document(document)
        print(f"VERIFIED rows={len(document['rows'])} max_a={document['max_a']}")
        return

    rows = []
    for a in range(1, args.max_a + 1):
        for k in range((a + 1) // 2):
            rows.append(make_row(a, k, args.base, args.bins))
        current = [row for row in rows if row["a"] == a]
        print(
            f"DONE a={a} reps={len(current)} "
            f"max_multiplier={max(row['multiplier'] for row in current)} "
            f"max_low={max(row['low_num'] / row['low_den'] for row in current)} "
            f"max_tail={max(row['tail_num'] / row['tail_den'] for row in current)}",
            flush=True,
        )
    document = {
        "schema": SCHEMA,
        "max_a": args.max_a,
        "reflection": "k maps to a-1-k modulo a",
        "analytic_tail": "mu(D_b)<=1/7+3*TV(f)/(49b)",
        "rows": rows,
    }
    verify_document(document)
    args.certificate.write_text(json.dumps(document, separators=(",", ":")) + "\n")
    print(f"WROTE {args.certificate} rows={len(rows)}")


if __name__ == "__main__":
    main()
