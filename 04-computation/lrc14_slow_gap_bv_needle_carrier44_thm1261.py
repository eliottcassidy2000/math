#!/usr/bin/env python3
"""Discover or exactly replay THM-1261's carrier-44 BV certificate.

LP output is discovery-only.  Every stored 200-bin density is checked with
exact rationals through its finite cutoff and by the THM-1182 BV tail after
that.  Reflection k -> 43-k has no integral fixed phase, so k=0..21 is the
complete 22-row bank.

The tournament loss audit orders phase-obligation stalks by exact bottleneck
load.  The order is transitive and forgets the density vector that carries the
actual noncoverage proof.
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
    "05-knowledge/results/lrc14_slow_gap_bv_needle_carrier44_thm1261_certificate.json"
)
SCHEMA = "lrc14-slow-gap-bv-needle-carrier44-thm1261-v1"
CARRIER = 44
BINS = 200
REPRESENTATIVES = range(22)
MULTIPLIERS = (12, 16, 20, 28, 44, 60, 76, 108)
PREFERRED = {
    0: 12, 1: 12, 2: 12, 3: 28, 4: 28, 5: 44, 6: 76,
    7: 60, 8: 20, 9: 28, 10: 12, 11: 12, 12: 44, 13: 28,
    14: 12, 15: 28, 16: 44, 17: 28, 18: 28, 19: 44,
    20: 28, 21: 12,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1261 carrier-44 verification failed: {message}")


def optimization_safe_require_probe() -> None:
    caught = False
    try:
        require(False, "deliberate optimization-safety probe")
    except RuntimeError as error:
        caught = "deliberate optimization-safety probe" in str(error)
    require(caught, "require probe did not fire")


def load_base():
    spec = importlib.util.spec_from_file_location("thm1182_base_c44", BASE_SCRIPT)
    require(spec is not None and spec.loader is not None, "load base generator")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_row(base, row: dict) -> tuple[Fraction, Fraction, Fraction]:
    require(row.get("a") == CARRIER, "carrier")
    require(row.get("k") in REPRESENTATIVES, "phase")
    require(row.get("bins") == BINS, "bins")
    require(row.get("multiplier") in MULTIPLIERS, "multiplier")
    weights = row.get("weights", [])
    require(len(weights) == BINS and all(w >= 0 for w in weights), "weights")
    require(sum(weights) > 0 and sum(weights) == row.get("total"), "mass")
    low, tail, variation = base.verify_data(
        CARRIER, row["k"], row["multiplier"], BINS, weights
    )
    stored = (
        Fraction(row["low_num"], row["low_den"]),
        Fraction(row["tail_num"], row["tail_den"]),
        Fraction(row["variation_num"], row["variation_den"]),
    )
    require((low, tail, variation) == stored, "stored invariants")
    require(low < Fraction(1, 6) and tail < Fraction(1, 6), "strict loads")
    return low, tail, variation


def verify_document(base, document: dict) -> list[tuple[dict, Fraction, Fraction]]:
    require(document.get("schema") == SCHEMA, "schema")
    require(document.get("carrier") == CARRIER, "carrier document")
    require(document.get("reflection") == "k maps to 43-k modulo 44", "reflection")
    rows = document.get("rows", [])
    keys = [(row.get("a"), row.get("k")) for row in rows]
    expected = {(CARRIER, k) for k in REPRESENTATIVES}
    require(set(keys) == expected and len(keys) == len(set(keys)), "phase rows")
    checked = []
    for row in sorted(rows, key=lambda item: item["k"]):
        low, tail, _ = verify_row(base, row)
        checked.append((row, low, tail))
    return checked


def make_row(base, k: int) -> dict:
    order = (PREFERRED[k],) + tuple(
        multiplier for multiplier in MULTIPLIERS if multiplier > PREFERRED[k]
    )
    failures = []
    for multiplier in order:
        try:
            rows, weights = base.discover(CARRIER, k, multiplier, BINS)
            low, tail, variation = base.verify_data(
                CARRIER, k, multiplier, BINS, weights, rows
            )
            return {
                "a": CARRIER, "k": k, "multiplier": multiplier, "bins": BINS,
                "weights": weights, "total": sum(weights),
                "low_num": low.numerator, "low_den": low.denominator,
                "tail_num": tail.numerator, "tail_den": tail.denominator,
                "variation_num": variation.numerator,
                "variation_den": variation.denominator,
            }
        except (RuntimeError, ValueError) as error:
            failures.append(repr(error))
    raise RuntimeError((CARRIER, k, failures))


def write_document(path: Path, rows: list[dict]) -> None:
    document = {
        "schema": SCHEMA,
        "carrier": CARRIER,
        "reflection": "k maps to 43-k modulo 44",
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
    print("THM-1261 CARRIER-44 BV NEEDLE CERTIFICATE")
    print("optimization_safe_require_probe=PASS")
    print(
        "VERIFIED_EXACT", f"rows={len(checked)}", "carrier=44",
        "phase_representatives=0..21",
        "multipliers=" + ",".join(f"{m}:{count[m]}" for m in sorted(count)),
    )
    print(
        f"max_low={max_low} gap_to_1/6={Fraction(1, 6)-max_low} "
        f"at_k={low_row['k']} M={low_row['multiplier']}"
    )
    print(
        f"max_tail={max_tail} gap_to_1/6={Fraction(1, 6)-max_tail} "
        f"at_k={tail_row['k']} M={tail_row['multiplier']}"
    )
    print("reflection=t->1-t sends G_k(44) to G_(43-k)(44)")
    print("phase_orbits=22 two-element pairs; no fixed phase")
    print("finite_layer=exact loads for every 45<=d<=M(k)*44")
    print("tail_layer=BV bound for every d>M(k)*44")
    print("conclusion=no six strict danger combs cover any complete 44-gap")
    path = sorted(
        ((max(low, tail), row["k"]) for row, low, tail in checked),
        key=lambda item: (item[0], item[1]),
    )
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=phase-obligation stalks k=0..21")
    print("observable=max(finite_load,tail_load); gauge=k")
    print("score_histogram=" + ",".join(str(i) for i in range(22)))
    print("directed_3_cycles=0 scc_sizes=" + ",".join("1" for _ in range(22)))
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
        rows = []
        if args.resume and args.certificate.exists():
            rows = json.loads(args.certificate.read_text()).get("rows", [])
        by_key = {(row["a"], row["k"]): row for row in rows}
        for k in REPRESENTATIVES:
            if (CARRIER, k) in by_key:
                continue
            row = make_row(base, k)
            by_key[CARRIER, k] = row
            write_document(args.certificate, list(by_key.values()))
            print(f"ROW a={CARRIER} k={k} M={row['multiplier']}", flush=True)
    checked = verify_document(base, json.loads(args.certificate.read_text()))
    print_summary(checked)


if __name__ == "__main__":
    main()
