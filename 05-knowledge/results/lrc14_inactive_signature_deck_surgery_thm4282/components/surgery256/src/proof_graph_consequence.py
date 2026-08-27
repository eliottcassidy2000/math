#!/usr/bin/env python3
"""Detached exact set audit for the combined (256,663) proof graph."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


MASK64 = (1 << 64) - 1
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_edges(path: Path) -> list[tuple[int, int]]:
    raw = path.read_bytes()
    need(not raw or raw.endswith(b"\n"), f"missing terminal LF: {path}")
    need(b"\r" not in raw, f"non-LF newline: {path}")
    rows: list[tuple[int, int]] = []
    for line in raw.splitlines():
        fields = line.split(b",")
        need(len(fields) == 2, f"malformed row in {path}: {line!r}")
        need(all(field and field.isdigit() for field in fields),
             f"non-decimal row in {path}: {line!r}")
        edge = (int(fields[0]), int(fields[1]))
        need(0 < edge[0] < edge[1], f"invalid edge in {path}: {edge}")
        rows.append(edge)
    need(rows == sorted(set(rows)), f"not lexicographic/distinct: {path}")
    return rows


def serialize(edges: list[tuple[int, int]]) -> bytes:
    return b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)


def fnv_edges(edges: list[tuple[int, int]]) -> str:
    state = FNV_OFFSET
    for q, r in edges:
        for value in (q, r):
            for shift in range(0, 64, 8):
                state ^= (value >> shift) & 0xFF
                state = (state * FNV_PRIME) & MASK64
    return f"{state:016x}"


def sha_edges(edges: list[tuple[int, int]]) -> str:
    return hashlib.sha256(serialize(edges)).hexdigest()


def top_description(edges: list[tuple[int, int]]) -> tuple[int, str]:
    maximum = max((r for _, r in edges), default=0)
    top = [edge for edge in edges if edge[1] == maximum]
    return maximum, ";".join(f"{q},{r}" for q, r in top) or "-"


def describe(name: str, edges: list[tuple[int, int]]) -> str:
    maximum, top = top_description(edges)
    return (f"{name} COUNT {len(edges)} FNV {fnv_edges(edges)} "
            f"SHA256 {sha_edges(edges)} MAX_ENDPOINT {maximum} TOP {top}")


def emit(out: Path, name: str, edges: set[tuple[int, int]]) -> list[tuple[int, int]]:
    rows = sorted(edges)
    (out / name).write_bytes(serialize(rows))
    return rows


def check_identity(name: str, edges: list[tuple[int, int]], count: int,
                   fnv: str, sha: str) -> None:
    need(len(edges) == count, f"{name} count changed")
    need(fnv_edges(edges) == fnv, f"{name} FNV changed")
    need(sha_edges(edges) == sha, f"{name} SHA256 changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--surgery", type=Path, required=True)
    parser.add_argument("--gain", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    surgery_rows = read_edges(args.surgery)
    gain_rows = read_edges(args.gain)
    residual_rows = read_edges(args.prior_residual)
    check_identity("surgery", surgery_rows, 188, "6588121dbec57bcb",
                   "5946ff653c51a74eba09a14430c494074e53b5aba87c3159bd17bafbe9e605d5")
    check_identity("prior gain", gain_rows, 586, "9769742754535d84",
                   "04aed4c107b3244a5e488266c46d1cae3bfffb20cddd831fc2d474c7b8c16a0e")
    check_identity("prior residual", residual_rows, 23637, "e8b363d2b3d9ba6a",
                   "21ad2530da6adc7187b4b5829c9fddad07ab704ccab96adb90d546e43756e95d")

    surgery = set(surgery_rows)
    gain = set(gain_rows)
    prior_residual = set(residual_rows)
    need(not gain & prior_residual, "prior gain and residual overlap")
    universe = gain | prior_residual
    universe_rows = sorted(universe)
    check_identity("reconstructed post-THM-4281 universe", universe_rows, 24223,
                   "80ec0687d8c7dba7",
                   "75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552")
    need(surgery <= universe, "surgery ledger leaves reconstructed universe")

    carrier = {edge for edge in prior_residual if edge[1] >= 645}
    surgery_gain = surgery & gain
    surgery_carrier = surgery & carrier
    gain_carrier = gain & carrier
    triple = surgery & gain & carrier
    exclusive_surgery = surgery - gain - carrier
    exclusive_gain = gain - surgery - carrier
    exclusive_carrier = carrier - surgery - gain
    combined = surgery | gain | carrier
    final = universe - combined

    ledgers = [
        ("carrier90.csv", carrier),
        ("overlap_surgery_gain9.csv", surgery_gain),
        ("overlap_surgery_carrier5.csv", surgery_carrier),
        ("overlap_gain_carrier0.csv", gain_carrier),
        ("overlap_triple0.csv", triple),
        ("exclusive_surgery174.csv", exclusive_surgery),
        ("exclusive_gain577.csv", exclusive_gain),
        ("exclusive_carrier85.csv", exclusive_carrier),
        ("combined_union850.csv", combined),
        ("final_residual23373.csv", final),
    ]
    emitted = {name: emit(args.out, name, edges) for name, edges in ledgers}

    expected = {
        "carrier90.csv": (90, "942995bee7469430", "222afb7618d887f32847b4531ffedf5616f20c2196e92f52ca2c11b09e1eab1f"),
        "overlap_surgery_gain9.csv": (9, "1f049d3f1c65c5c9", "e46c6386a9591d96672ceaba38ac9474d3fe7a6ea1fd943ee35a2d8e0ef735d0"),
        "overlap_surgery_carrier5.csv": (5, "f3cb9446a58b5543", "71f7557bafc537088c30ebfc6265a13e596658b6f6d60552a90cbdbdf59d2685"),
        "overlap_gain_carrier0.csv": (0, "cbf29ce484222325", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        "overlap_triple0.csv": (0, "cbf29ce484222325", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
        "exclusive_surgery174.csv": (174, "07b6c9a8a9b4ce65", "1df3a568b48f744a2277a0a78dc0b477f2216ccb3b202e32ab24216f302ae335"),
        "exclusive_gain577.csv": (577, "9892130695fe8d20", "e0f0257015f76b54db21aa094c8b5260cf36aebccaf1ec77eb54ac2576135dbf"),
        "exclusive_carrier85.csv": (85, "ff2d28bc0a1ae156", "5a1d98b93fb5d9dccf775a03bbe378c5e2b64a820398331c5d893a7e6591a1c0"),
        "combined_union850.csv": (850, "8f595510210a5785", "7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b"),
        "final_residual23373.csv": (23373, "c6ab0ae49ee32273", "c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3"),
    }
    for name, (count, fnv, sha) in expected.items():
        check_identity(name, emitted[name], count, fnv, sha)

    need(len(combined) == 188 + 586 + 90 - 9 - 5,
         "inclusion-exclusion failed")
    need(final == prior_residual - surgery - carrier,
         "prior-residual subtraction identity failed")
    maximum, top = top_description(sorted(final))
    need(maximum == 644, "final maximum endpoint changed")
    need(top == "220,644;256,644;258,644;294,644;366,644;416,644;512,644",
         "final top boundary changed")

    lines = [
        "THM4281_SIG663_SOLVERFREE_OVERLAP_AUDIT_V1",
        describe("SURGERY", surgery_rows),
        describe("PRIOR_GAIN", gain_rows),
        describe("PRIOR_RESIDUAL", residual_rows),
        describe("RECONSTRUCTED_UNIVERSE", universe_rows),
    ]
    for name, _ in ledgers:
        lines.append(describe(name.removesuffix(".csv").upper(), emitted[name]))
    lines.extend([
        "INCLUSION_EXCLUSION 188+586+90-9-5-0+0=850",
        "PRIOR_RESIDUAL_DELETION 179+90-5=264",
        "OVERLAP_SURGERY_GAIN_ROWS " + ";".join(f"{q},{r}" for q, r in emitted["overlap_surgery_gain9.csv"]),
        "OVERLAP_SURGERY_CARRIER_ROWS " + ";".join(f"{q},{r}" for q, r in emitted["overlap_surgery_carrier5.csv"]),
        "VERDICT PASS EXACT_SET_UNION_AND_PROOF_GRAPH_RESIDUAL",
    ])
    summary = "\n".join(lines) + "\n"
    (args.out / "audit.out").write_text(summary, encoding="ascii", newline="\n")
    print(summary, end="")


if __name__ == "__main__":
    main()
