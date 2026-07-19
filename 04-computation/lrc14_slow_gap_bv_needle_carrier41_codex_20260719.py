#!/usr/bin/env python3
"""Discover or exactly replay THM-1255's carrier-41 BV needle certificate.

This is a separate continuation of THM-1182 (carriers 1..30) and THM-1197
(31..40).  A floating-point LP only proposes a nonnegative 200-bin density.
Every accepted row is then checked with ``fractions.Fraction`` against all
speeds ``42 <= d <= M*41`` and against the analytic tail

    mu_f(D_d) <= 1/7 + 3*TV(f)/(49*d),   d > M*41.

Strict load below 1/6 for both ranges forbids six danger combs from covering
the addressed carrier gap.  Reflection sends k to 40-k, so k=0..20 is a
complete phase-representative set.

Tournament/carrier audit
------------------------
Runner vertices are the wrong quotient: each row controls infinitely many
possible faster runners at once.  We challenge runners, teeth, bins, walls,
residues, phase gaps, and proof obligations as vertices.  The faithful object
is the 21 phase-obligation stalks with their full density vectors.  For the
required pairwise loss audit, orient two phase obligations by the exact
bottleneck max(low, tail), breaking ties by k.  This gives a transitive
tournament and records only computational difficulty; it destroys the density
vector and therefore does not prove the covering predicate.
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
    "lrc14_slow_gap_bv_needle_carrier41_codex_20260719_certificate.json"
)
SCHEMA = "lrc14-slow-gap-bv-needle-carrier41-codex-20260719-v1"
CARRIER = 41
BINS = 200
MULTIPLIERS = (12, 16, 20, 28, 44, 60, 76, 108)

# Successful discovery depths from the initial exact phase census.  Trying the
# known depth first avoids retaining several failed large LPs in one process.
PREFERRED = {
    0: 12, 1: 12, 2: 12, 3: 20, 4: 28, 5: 44, 6: 28,
    7: 28, 8: 16, 9: 20, 10: 12, 11: 44, 12: 28, 13: 12,
    14: 20, 15: 60, 16: 28, 17: 44, 18: 20, 19: 12, 20: 12,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"carrier-41 verification failed: {message}")


def optimization_safe_require_probe() -> None:
    caught = False
    try:
        require(False, "deliberate optimization-safety probe")
    except RuntimeError as error:
        caught = "deliberate optimization-safety probe" in str(error)
    require(caught, "require probe did not fire")


def load_base():
    spec = importlib.util.spec_from_file_location("thm1182_base_c41", BASE_SCRIPT)
    require(spec is not None and spec.loader is not None, "load base generator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def expected_keys() -> set[tuple[int, int]]:
    return {(CARRIER, k) for k in range((CARRIER + 1) // 2)}


def verify_row(base, row: dict) -> tuple[Fraction, Fraction, Fraction]:
    require(row.get("a") == CARRIER, "row carrier")
    require(row.get("k") in range(21), "row phase")
    require(row.get("bins") == BINS, "row bin count")
    require(row.get("multiplier") in MULTIPLIERS, "row multiplier")
    weights = row.get("weights", [])
    require(len(weights) == BINS and all(weight >= 0 for weight in weights),
            "row weights")
    require(sum(weights) > 0 and sum(weights) == row.get("total"), "row mass")
    low, tail, variation = base.verify_data(
        CARRIER, row["k"], row["multiplier"], BINS, weights
    )
    stored = (
        Fraction(row["low_num"], row["low_den"]),
        Fraction(row["tail_num"], row["tail_den"]),
        Fraction(row["variation_num"], row["variation_den"]),
    )
    require((low, tail, variation) == stored, "stored exact invariants")
    require(low < Fraction(1, 6), "finite-speed load is strict")
    require(tail < Fraction(1, 6), "infinite-tail load is strict")
    return low, tail, variation


def verify_document(base, document: dict) -> list[tuple[dict, Fraction, Fraction]]:
    require(document.get("schema") == SCHEMA, "certificate schema")
    require(document.get("carrier") == CARRIER, "certificate carrier")
    require(document.get("reflection") == "k maps to 40-k modulo 41",
            "reflection metadata")
    rows = document.get("rows", [])
    keys = [(row.get("a"), row.get("k")) for row in rows]
    require(set(keys) == expected_keys() and len(keys) == len(set(keys)),
            "complete phase representatives")
    checked = []
    for row in sorted(rows, key=lambda item: item["k"]):
        low, tail, _ = verify_row(base, row)
        checked.append((row, low, tail))
    return checked


def make_row(base, k: int) -> dict:
    order = (PREFERRED[k],) + tuple(
        multiplier for multiplier in MULTIPLIERS if multiplier > PREFERRED[k]
    )
    errors: list[str] = []
    for multiplier in order:
        try:
            rows, weights = base.discover(CARRIER, k, multiplier, BINS)
            low, tail, variation = base.verify_data(
                CARRIER, k, multiplier, BINS, weights, rows
            )
            return {
                "a": CARRIER,
                "k": k,
                "multiplier": multiplier,
                "bins": BINS,
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
    raise RuntimeError((CARRIER, k, errors))


def write_document(path: Path, rows: list[dict]) -> None:
    document = {
        "schema": SCHEMA,
        "carrier": CARRIER,
        "reflection": "k maps to 40-k modulo 41",
        "analytic_tail": "mu(D_d)<=1/7+3*TV(f)/(49*d)",
        "rows": sorted(rows, key=lambda row: row["k"]),
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, separators=(",", ":")) + "\n")
    temporary.replace(path)


def print_summary(checked: list[tuple[dict, Fraction, Fraction]]) -> None:
    low_row, max_low, _ = max(checked, key=lambda item: item[1])
    tail_row, _, max_tail = max(checked, key=lambda item: item[2])
    count = Counter(row["multiplier"] for row, _, _ in checked)
    print("CARRIER_41_BV_NEEDLE_CERTIFICATE")
    print("optimization_safe_require_probe=PASS")
    print(
        "VERIFIED_EXACT",
        f"rows={len(checked)}",
        "carrier=41",
        "phase_representatives=0..20",
        "multipliers=" + ",".join(
            f"{key}:{count[key]}" for key in sorted(count)
        ),
    )
    print(
        f"max_low={max_low} gap_to_1/6={Fraction(1, 6) - max_low} "
        f"at_k={low_row['k']} M={low_row['multiplier']}"
    )
    print(
        f"max_tail={max_tail} gap_to_1/6={Fraction(1, 6) - max_tail} "
        f"at_k={tail_row['k']} M={tail_row['multiplier']}"
    )
    print("reflection=t->1-t sends G_k(41) to G_(40-k)(41)")
    print("finite_layer=exact loads for every 42<=d<=M(k)*41")
    print("tail_layer=BV bound for every d>M(k)*41")
    print("conclusion=no six strict danger combs cover any complete 41-gap")

    # The bottleneck tournament is transitive by construction.  It is a loss
    # audit, not a proof input.
    path = sorted(
        ((max(low, tail), row["k"]) for row, low, tail in checked),
        key=lambda item: (item[0], item[1]),
    )
    require(len({k for _, k in path}) == 21, "tournament vertices")
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=phase-obligation stalks k=0..20")
    print("observable=max(finite_load,tail_load); gauge=k")
    print("score_histogram=" + ",".join(str(score) for score in range(21)))
    print("directed_3_cycles=0 scc_sizes=" + ",".join("1" for _ in range(21)))
    print("hamiltonian_path_count=1")
    print("tie_hamiltonian_path=" + ",".join(str(k) for _, k in path))
    print("faithful_vertices=density-bin stalks; tournament forgets weights")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--discover", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    optimization_safe_require_probe()
    base = load_base()

    if args.discover:
        rows: list[dict] = []
        if args.resume and args.certificate.exists():
            rows = json.loads(args.certificate.read_text()).get("rows", [])
        by_key = {(row["a"], row["k"]): row for row in rows}
        for k in range(21):
            if (CARRIER, k) in by_key:
                continue
            row = make_row(base, k)
            by_key[CARRIER, k] = row
            write_document(args.certificate, list(by_key.values()))
            print(
                f"ROW a={CARRIER} k={k} M={row['multiplier']} "
                f"low={row['low_num'] / row['low_den']:.15f} "
                f"tail={row['tail_num'] / row['tail_den']:.15f}",
                flush=True,
            )

    document = json.loads(args.certificate.read_text())
    checked = verify_document(base, document)
    print_summary(checked)


if __name__ == "__main__":
    main()
