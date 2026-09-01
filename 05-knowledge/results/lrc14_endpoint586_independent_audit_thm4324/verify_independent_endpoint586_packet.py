#!/usr/bin/env python3
"""Hardened self-contained verifier for the endpoint-586 clean-room packet."""

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


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_word(state: int, word: int) -> int:
    word &= MASK64
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 255)) * PRIME) & MASK64
    return state


def fnv(words) -> int:
    state = OFFSET
    for word in words:
        state = add_word(state, word)
    return state


def verify_manifest(packet: Path) -> int:
    entries: dict[str, str] = {}
    for line in (packet / "SHA256SUMS").read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and name not in entries,
                "bad/duplicate manifest row")
        entries[name] = digest
    actual = {path.relative_to(packet).as_posix() for path in packet.rglob("*")
              if path.is_file() and path.name != "SHA256SUMS" and
              path.suffix.lower() != ".exe" and "__pycache__" not in path.parts}
    require(set(entries) == actual, "manifest closure mismatch")
    for name, digest in entries.items():
        require(sha(packet / name) == digest, f"manifest mismatch: {name}")
    return len(entries)


def read_masks(path: Path) -> list[int]:
    masks = [int(token, 16) for token in path.read_text(encoding="ascii").split()]
    require(all(0 <= mask < (1 << 30) for mask in masks),
            f"mask escaped universe: {path.name}")
    return masks


def verify_inherited(packet: Path) -> None:
    joint = read_masks(packet / "inputs/joint421_masks.txt")
    carrier = read_masks(packet / "inputs/inherited_carrier3925.txt")
    require(len(joint) == len(set(joint)) == 421 and
            all(mask.bit_count() == 8 for mask in joint) and
            fnv(joint) == 0x20D63DD42FE8150E,
            "joint deck changed")
    require(len(carrier) == len(set(carrier)) == 3925 and
            sum(mask.bit_count() == 8 for mask in carrier) == 3809 and
            sum(mask.bit_count() == 9 for mask in carrier) == 116 and
            fnv(carrier) == 0xEEAE5518D84CCAC5 and set(joint) <= set(carrier),
            "serialized THM4318 carrier changed")
    require(sha(packet / "inputs/inherited_carrier3925.txt") ==
            "75265312f0b1a52d92d7f741dd0596b398229c15bc0a1cdf873e3871840bd213",
            "carrier serialization SHA changed")
    provenance = (packet / "export_carrier.out").read_text(encoding="ascii")
    require("MASKS 3925 RANK8 3809 RANK9 116 FNV eeae5518d84ccac5" in provenance and
            "JOINT_RETAINED 421" in provenance,
            "carrier provenance summary changed")


def verify_carrier(packet: Path) -> list[int]:
    for stem in ("carrier_pair", "carrier_failures", "carrier"):
        suffix = "csv" if stem != "carrier" else "out"
        require((packet / f"{stem}_O2.{suffix}").read_bytes() ==
                (packet / f"{stem}_O3.{suffix}").read_bytes(),
                f"carrier O2/O3 mismatch: {stem}")
    pair_path = packet / "carrier_pair_O2.csv"
    failure_path = packet / "carrier_failures_O2.csv"
    require(sha(pair_path) ==
            "8e80b5349180e06df0ef09d5600178fc7778e2dfea415ae3ab34e33d66d3835b" and
            sha(failure_path) ==
            "66127795e4a1ea97fc37540c2f53ea686bba76babb172ef3c240fb24971713ca",
            "carrier ledger SHA changed")
    with pair_path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    expected_q = [6, 25, 50, 72, 96, 100, 105, 192, 206, 210, 260, 294]
    require(len(rows) == 12 and [int(row["q"]) for row in rows] == expected_q and
            all(int(row["r"]) == 586 for row in rows),
            "carrier frontier changed")
    words: list[int] = []
    total_exposed = total_failures = 0
    for row in rows:
        values = [int(row["q"]), int(row["r"]), int(row["active"]),
                  int(row["active_fnv"], 16), int(row["active_joint"]),
                  int(row["active_nonjoint"]), int(row["exposed"]),
                  int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
                  int(row["maximum_hits"]), int(row["failures"]),
                  int(row["failure_fnv"], 16)]
        words.extend(values); total_exposed += values[6]; total_failures += values[10]
    require(total_exposed == 525048 and total_failures == 4090 and
            fnv(words) == 0x46A2C17CAECC55DF,
            "carrier census/pair FNV changed")
    raw = (packet / "carrier_O2.out").read_text(encoding="ascii")
    require("BODY_TESTS 171685800" in raw and
            "HIT_INCIDENCES 12125735" in raw and
            "FAILURE_FNV ffb884b2b17e6ef4" in raw and
            "VERDICT HOSTILE_FAIL" in raw,
            "carrier summary changed")

    bodies: list[int] = []
    global_words: list[int] = []
    with failure_path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in reader:
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require(q == 50 and r == 586 and body.bit_count() == 9,
                    "failure escaped q50 endpoint586 rank9")
            bodies.append(body); global_words.extend((q, r, body))
    require(len(bodies) == 4090 and bodies == sorted(set(bodies)) and
            fnv(bodies) == 0x14CE094F4AB4BA94 and
            fnv(global_words) == 0xFFB884B2B17E6EF4,
            "failure census/FNV changed")
    return bodies


def verify_literal(packet: Path, bodies: list[int]) -> None:
    normal = packet / "literal_normal.csv"
    optimized = packet / "literal_opt.csv"
    require(normal.read_bytes() == optimized.read_bytes() and
            (packet / "literal_normal.out").read_bytes() ==
            (packet / "literal_opt.out").read_bytes(),
            "literal normal/-O mismatch")
    require(sha(normal) ==
            "6dbfc6178c7da5c844abc95e980f2696cf02c3845a83bdfff95f881c12c9c4fd",
            "literal detail SHA changed")
    count = equality = 0
    detail_state = OFFSET
    minimum: tuple[int, int] | None = None
    maximum: tuple[int, int] | None = None
    with normal.open(newline="", encoding="ascii") as handle:
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
            minimum = (low_ticks, body) if minimum is None else min(minimum, (low_ticks, body))
            if maximum is None or low_ticks > maximum[0] or (
                    low_ticks == maximum[0] and body < maximum[1]):
                maximum = (low_ticks, body)
            for word in (q, r, ordinal, body, low, full,
                         low_ticks & MASK64, (low_ticks & MASK128) >> 64,
                         full_ticks & MASK64, (full_ticks & MASK128) >> 64):
                detail_state = add_word(detail_state, word)
            count += 1
    require(count == 4090 and equality == 3602 and
            minimum == (136976666519138544, 0x011C6405) and
            maximum == (247691989153176624, 0x3340E400) and
            detail_state == 0xE657AC566B502B62,
            "literal census/extrema/FNV changed")
    out = (packet / "literal_normal.out").read_text(encoding="ascii")
    require("PAIR_SAFE_CELLS 5802" in out and "ALL_CLASSES 2424 LOW_CLASSES 2371" in out and
            "LOW_POSITIVE 4090 FULL_POSITIVE 4090 EQUALITY 3602" in out and
            "VERDICT PASS" in out,
            "literal summary changed")


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical rows: {path.name}")
    return rows


def rows_fnv(rows: list[tuple[int, int]]) -> int:
    return fnv(word for row in rows for word in row)


def verify_typed(packet: Path) -> None:
    normal, optimized = packet / "typed_normal", packet / "typed_opt"
    names = ["reconstructed_top586.csv", "typed_union2229.csv",
             "final_residual20418.csv", "residual_top585.csv"]
    for name in names:
        require((normal / name).read_bytes() == (optimized / name).read_bytes(),
                f"typed normal/-O mismatch: {name}")
    require((packet / "typed_normal.out").read_bytes() ==
            (packet / "typed_opt.out").read_bytes(),
            "typed stdout normal/-O mismatch")
    universe = read_rows(packet / "inputs/universe22647.csv")
    prior_union = read_rows(packet / "inputs/typed_union2217.csv")
    prior_residual = read_rows(packet / "inputs/final_residual20430.csv")
    top = read_rows(normal / names[0]); union = read_rows(normal / names[1])
    residual = read_rows(normal / names[2]); next_top = read_rows(normal / names[3])
    require([len(x) for x in (universe, prior_union, prior_residual, top,
                              union, residual, next_top)] ==
            [22647, 2217, 20430, 12, 2229, 20418, 23],
            "typed counts changed")
    require(top == [row for row in prior_residual if row[1] == 586] and
            union == sorted(set(prior_union) | set(top)) and
            residual == sorted(set(prior_residual) - set(top)) and
            next_top == [row for row in residual if row[1] == 585] and
            set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "typed transition changed")
    require(rows_fnv(top) == 0xA1B617FAA2E7F63F and
            rows_fnv(union) == 0x035EBF12F02ECC62 and
            rows_fnv(residual) == 0x89B73E31224821C4 and
            rows_fnv(next_top) == 0x8F1B7C8DB8FD5E87,
            "typed FNVs changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet.resolve()
    require(packet.is_dir(), "packet missing")
    files = verify_manifest(packet)
    verify_inherited(packet)
    bodies = verify_carrier(packet)
    verify_literal(packet, bodies)
    verify_typed(packet)
    print("ENDPOINT586_INDEPENDENT_PACKET_VERIFY PASS")
    print(f"manifest_files={files}")
    print("frontier_rows=12 carrier_failures=4090 direct_positive=4090")
    print("typed_union=2229 residual=20418 next_endpoint=585 next_rows=23")
    print("agreement=SCOUT_PAIR_FAILURE_AND_LITERAL_DETAIL_BYTE_IDENTICAL")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
