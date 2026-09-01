#!/usr/bin/env python3
"""Independent structural verifier for the THM-4311 candidate packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def same(left: Path, right: Path) -> None:
    need(left.read_bytes() == right.read_bytes(), f"byte mismatch: {left} {right}")


def fnv_values(values: list[int]) -> int:
    state = OFFSET
    for value in values:
        for byte in value.to_bytes(8, "little"):
            state = ((state ^ byte) * PRIME) & MASK64
    return state


def plain_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    need(rows == sorted(set(rows)), f"unordered/duplicate rows: {path}")
    return rows


def audit_rows(path: Path) -> dict[tuple[int, int], dict[str, str]]:
    result: dict[tuple[int, int], dict[str, str]] = {}
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["q"]), int(row["r"]))
            need(key not in result, f"duplicate audit key: {key}")
            result[key] = row
    return result


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    scratch = args.scratch.resolve()
    support = repo / (
        "05-knowledge/results/"
        "lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309"
    )
    layer594 = repo / (
        "05-knowledge/results/"
        "lrc14_endpoint594_complete_layer_closure_thm4310"
    )
    old = repo / (
        "05-knowledge/results/"
        "lrc14_mixed_rank_depth_recursive_signatures_thm4296"
    )

    top593 = plain_rows(scratch / "residual_top593.csv")
    need(len(top593) == 16 and all(r == 593 for _, r in top593),
         "top593 identity/count changed")
    need(fnv_values([value for row in top593 for value in row]) ==
         0x5424C07FA724011F, "top593 FNV changed")

    for stem, extension in (("endpoint593_raw", "out"),
                            ("endpoint593_pair_audit", "csv"),
                            ("endpoint593_failures", "csv")):
        same(scratch / f"{stem}_O3.{extension}",
             scratch / f"{stem}_O2.{extension}")
    failures = list(csv.DictReader(
        (scratch / "endpoint593_failures_O3.csv").open(
            newline="", encoding="ascii")))
    need(failures == [{"q": "96", "r": "593", "body_hex": "34087401"}],
         "hostile failure identity changed")

    for stem, extension in (("endpoint593_singleton_deletion_quotient", "out"),
                            ("singleton432", "csv"),
                            ("protected_nonjoint432", "txt"),
                            ("safe_delete1", "txt")):
        same(scratch / f"{stem}_O3.{extension}",
             scratch / f"{stem}_O2.{extension}")
    with (scratch / "singleton432_O3.csv").open(
            newline="", encoding="ascii") as handle:
        singletons = list(csv.DictReader(handle))
    protected = (scratch / "protected_nonjoint432_O3.txt").read_text(
        encoding="ascii").splitlines()
    need(len(singletons) == 1520 and len(set(protected)) == 412,
         "singleton/protected counts changed")
    addition_private = [row for row in singletons
                        if row["witness_hex"] == "0036092c"]
    need(addition_private == [{"q": "96", "r": "593",
                               "body_hex": "34087401",
                               "witness_hex": "0036092c",
                               "witness_rank": "9"}],
         "selected response is not the hostile obligation's unique witness")
    safe = (scratch / "safe_delete1_O3.txt").read_text(
        encoding="ascii").strip()
    need(safe == "0006e281" and safe not in set(protected),
         "selected deletion is not certified safe")

    for stem, extension in (("exchange_full432_raw", "out"),
                            ("exchange_full432_pair_audit", "csv"),
                            ("exchange_full432_failures", "csv")):
        same(scratch / f"{stem}_O3.{extension}",
             scratch / f"{stem}_O2.{extension}")
    prefix391 = set(audit_rows(
        support / "results/final_full_prefix_pair_audit.csv"))
    top594 = set(audit_rows(
        layer594 / "results/direct/endpoint594_pair_audit_O3.csv"))
    expected432 = prefix391 | top594 | set(top593)
    need(len(prefix391) == 391 and len(top594) == 25 and
         len(expected432) == 432, "full target construction changed")
    exchange = audit_rows(scratch / "exchange_full432_pair_audit_O3.csv")
    need(set(exchange) == expected432 and
         all(int(row["failures"]) == 0 for row in exchange.values()),
         "exchange does not close exact full432 target")
    with (scratch / "exchange_full432_failures_O3.csv").open(
            newline="", encoding="ascii") as handle:
        need(not list(csv.DictReader(handle)), "exchange failure ledger nonempty")
    target_fnv = fnv_values([value for row in sorted(expected432) for value in row])
    need(target_fnv == 0xA7ED492C64D1C0D8, "full432 FNV changed")

    universe = set(plain_rows(old / "inputs/current_residual22647.csv"))
    union = set(plain_rows(scratch / "typed/typed_union2052.csv"))
    residual = set(plain_rows(scratch / "typed/final_residual20595.csv"))
    top592 = plain_rows(scratch / "typed/residual_top592.csv")
    need(len(universe) == 22647 and len(union) == 2052 and
         len(residual) == 20595 and not union & residual and
         union | residual == universe, "typed partition changed")
    need(len(top592) == 35 and all(r == 592 for _, r in top592),
         "typed frontier changed")

    print("LRC14_ENDPOINT593_PACKET_STRUCTURAL_VERIFY_V1")
    print("TOP593 16 FNV 5424c07fa724011f HOSTILE_FAILURE 96,593,34087401")
    print("FULL_TARGET 432 FNV a7ed492c64d1c0d8 EXCHANGE_FAILURES 0")
    print("PRIVATE_OBLIGATIONS 1520 PROTECTED_NONJOINT 412 SAFE_NONJOINT 3093")
    print(f"EXCHANGE_PAIR_SHA256 {sha(scratch / 'exchange_full432_pair_audit_O3.csv')}")
    print(f"TYPED_UNION_SHA256 {sha(scratch / 'typed/typed_union2052.csv')}")
    print(f"TYPED_RESIDUAL_SHA256 {sha(scratch / 'typed/final_residual20595.csv')}")
    print("SCOPE STRUCTURAL_PACKET_CHECK_FIXED_INPUTS_NOT_INDEPENDENT_ARITHMETIC_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
