#!/usr/bin/env python3
"""Hardened verifier for the frozen endpoint-586 literal-closure packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
MASK128 = (1 << 128) - 1
GRID = 26723298545143200


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_word(state: int, word: int) -> int:
    word &= MASK64
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def fnv(words) -> int:
    state = OFFSET
    for word in words:
        state = add_word(state, word)
    return state


def verify_manifest(packet: Path) -> int:
    entries: dict[str, str] = {}
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing SHA256SUMS")
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and name not in entries,
                "bad/duplicate manifest row")
        entries[name] = digest
    actual = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and
        path.suffix.lower() != ".exe" and "__pycache__" not in path.parts
    }
    require(set(entries) == actual, "manifest closure mismatch")
    for name, digest in entries.items():
        require(sha(packet / name) == digest,
                f"manifest digest mismatch: {name}")
    return len(entries)


def verify_carrier(packet: Path) -> list[int]:
    for stem in ("pair", "failures", "raw"):
        suffix = "csv" if stem != "raw" else "out"
        o2 = packet / f"endpoint586_{stem}_O2.{suffix}"
        o3 = packet / f"endpoint586_{stem}_O3.{suffix}"
        require(o2.read_bytes() == o3.read_bytes(),
                f"O2/O3 carrier mismatch: {stem}")
    require(sha(packet / "endpoint586_pair_O2.csv") ==
            "8e80b5349180e06df0ef09d5600178fc7778e2dfea415ae3ab34e33d66d3835b",
            "pair SHA changed")
    require(sha(packet / "endpoint586_failures_O2.csv") ==
            "66127795e4a1ea97fc37540c2f53ea686bba76babb172ef3c240fb24971713ca",
            "failure SHA changed")
    with (packet / "endpoint586_pair_O2.csv").open(
            newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 12 and [int(row["q"]) for row in rows] ==
            [6, 25, 50, 72, 96, 100, 105, 192, 206, 210, 260, 294] and
            all(int(row["r"]) == 586 for row in rows),
            "pair frontier changed")
    words: list[int] = []
    observed_failures: dict[int, tuple[int, int]] = {}
    total_exposed = total_failures = 0
    for row in rows:
        values = [int(row["q"]), int(row["r"]), int(row["active"]),
                  int(row["active_fnv"], 16), int(row["active_joint"]),
                  int(row["active_nonjoint"]), int(row["exposed"]),
                  int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
                  int(row["maximum_hits"]), int(row["failures"]),
                  int(row["failure_fnv"], 16)]
        words.extend(values); total_exposed += values[6]
        total_failures += values[10]
        if values[10]: observed_failures[values[0]] = (values[10], values[11])
    require(total_exposed == 525048 and total_failures == 4090 and
            observed_failures == {50: (4090, 0x14CE094F4AB4BA94)} and
            fnv(words) == 0x46A2C17CAECC55DF,
            "carrier pair census/FNV changed")
    raw = (packet / "endpoint586_raw_O2.out").read_text(encoding="ascii")
    require("BODY_TESTS 171685800" in raw and
            "HIT_INCIDENCES 12125735" in raw and
            "FAILURE_FNV ffb884b2b17e6ef4" in raw and
            "VERDICT HOSTILE_FAIL" in raw,
            "carrier raw summary changed")

    bodies: list[int] = []
    global_words: list[int] = []
    with (packet / "endpoint586_failures_O2.csv").open(
            newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in reader:
            q, r = int(row["q"]), int(row["r"])
            body = int(row["body_hex"], 16)
            require(q == 50 and r == 586 and body.bit_count() == 9,
                    "failure body escaped universe")
            bodies.append(body); global_words.extend((q, r, body))
    require(len(bodies) == 4090 and bodies == sorted(set(bodies)) and
            fnv(bodies) == 0x14CE094F4AB4BA94 and
            fnv(global_words) == 0xFFB884B2B17E6EF4,
            "failure body ledger changed")
    return bodies


def verify_literal(packet: Path, bodies: list[int]) -> None:
    p2 = packet / "endpoint586_direct_primary_O2.csv"
    p3 = packet / "endpoint586_direct_primary_O3.csv"
    independent = packet / "endpoint586_direct_rawcell_independent_O3.csv"
    require(p2.read_bytes() == p3.read_bytes() == independent.read_bytes(),
            "literal implementations disagree")
    require(sha(p2) ==
            "6dbfc6178c7da5c844abc95e980f2696cf02c3845a83bdfff95f881c12c9c4fd",
            "literal detail SHA changed")
    state = OFFSET
    equality = 0
    minimum: tuple[int, int] | None = None
    count = 0
    with p2.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "ordinal", "body_hex",
                "truncated_mass", "full_mass", "truncated_scaled_ticks",
                "full_scaled_ticks"], "literal header changed")
        for row in reader:
            q, r, ordinal = int(row["q"]), int(row["r"]), int(row["ordinal"])
            body = int(row["body_hex"], 16)
            low, full = int(row["truncated_mass"]), int(row["full_mass"])
            low_ticks = int(row["truncated_scaled_ticks"])
            full_ticks = int(row["full_scaled_ticks"])
            require(q == 50 and r == 586 and ordinal == count and
                    body == bodies[ordinal] and 0 <= low <= full,
                    "literal alignment/bound changed")
            require(low_ticks == 63 * low - 4 * GRID and
                    full_ticks == 63 * full - 4 * GRID and low_ticks > 0,
                    "literal arithmetic/positivity changed")
            equality += low == full
            candidate = (low_ticks, body)
            minimum = candidate if minimum is None else min(minimum, candidate)
            for word in (q, r, ordinal, body, low, full,
                         low_ticks & MASK64, (low_ticks & MASK128) >> 64,
                         full_ticks & MASK64, (full_ticks & MASK128) >> 64):
                state = add_word(state, word)
            count += 1
    require(count == 4090 and equality == 3602 and
            minimum == (136976666519138544, 0x011C6405) and
            state == 0xE657AC566B502B62,
            "literal census/minimum/FNV changed")
    out2 = (packet / "endpoint586_direct_primary_O2.out").read_bytes()
    out3 = (packet / "endpoint586_direct_primary_O3.out").read_bytes()
    independent_out = (
        packet / "endpoint586_direct_rawcell_independent_O3.out").read_bytes()
    require(out2 == out3 and b"LOW_POSITIVE 4090 FULL_POSITIVE 4090" in out2 and
            b"VERDICT PASS" in out2 and
            b"LOW_POSITIVE 4090 FULL_POSITIVE 4090" in independent_out and
            b"NO_CLASS_AGGREGATION" in independent_out,
            "literal summaries changed")


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical rows: {path.name}")
    return rows


def rows_fnv(rows: list[tuple[int, int]]) -> int:
    return fnv(word for row in rows for word in row)


def verify_typed(packet: Path) -> None:
    normal = packet / "post586_typed"
    optimized = packet / "post586_typed_opt"
    names = ["typed_union2229.csv", "final_residual20418.csv",
             "residual_top585.csv"]
    for name in names:
        require((normal / name).read_bytes() == (optimized / name).read_bytes(),
                f"typed normal/-O mismatch: {name}")
    require((packet / "post586_typed.out").read_bytes() ==
            (packet / "post586_typed_opt.out").read_bytes(),
            "typed stdout normal/-O mismatch")
    universe = read_rows(packet / "inputs/typed/universe22647.csv")
    prior_union = read_rows(packet / "inputs/typed/typed_union2217.csv")
    prior_residual = read_rows(packet / "inputs/typed/final_residual20430.csv")
    top586 = read_rows(packet / "inputs/typed/residual_top586.csv")
    union = read_rows(normal / names[0]); residual = read_rows(normal / names[1])
    top585 = read_rows(normal / names[2])
    require([len(x) for x in (universe, prior_union, prior_residual, top586,
                              union, residual, top585)] ==
            [22647, 2217, 20430, 12, 2229, 20418, 23],
            "typed counts changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe) and
            union == sorted(set(prior_union) | set(top586)) and
            residual == sorted(set(prior_residual) - set(top586)) and
            set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "typed partition transition changed")
    require(rows_fnv(union) == 0x035EBF12F02ECC62 and
            rows_fnv(residual) == 0x89B73E31224821C4 and
            rows_fnv(top585) == 0x8F1B7C8DB8FD5E87 and
            top585 == [row for row in residual if row[1] == 585],
            "typed FNV/frontier changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet.resolve()
    require(packet.is_dir(), "packet directory missing")
    manifest_count = verify_manifest(packet)
    bodies = verify_carrier(packet)
    verify_literal(packet, bodies)
    verify_typed(packet)
    print("ENDPOINT586_LITERAL_CLOSURE_PACKET_VERIFY PASS")
    print(f"manifest_files={manifest_count}")
    print("endpoint586_rows=12 carrier_failures=4090 direct_positive=4090")
    print("typed_union=2229 residual=20418 next_endpoint=585 next_rows=23")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
