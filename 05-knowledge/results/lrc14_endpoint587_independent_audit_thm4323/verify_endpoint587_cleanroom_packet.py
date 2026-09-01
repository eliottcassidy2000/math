#!/usr/bin/env python3
"""Self-contained verifier for the endpoint-587 clean-room packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path

FNV_BASIS = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
EXPECTED_INPUT_SHA = {
    "inputs/final_residual20440.csv": "be2d63e98beefb062e9ae4436d490ee2f630352989bf309cf85b5a1ffc44278c",
    "inputs/inherited_carrier3925.txt": "75265312f0b1a52d92d7f741dd0596b398229c15bc0a1cdf873e3871840bd213",
    "inputs/joint421_masks.txt": "1aaa3865bfd76697d861574b40f99745429fca0f7643c60b1a574649b3af0a96",
    "inputs/typed_union2207.csv": "f03c84f15d9a149b0a083b50e922118e814d5644f5fa21e7011ae1c414ff3675",
    "inputs/universe22647.csv": "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
}


class Fnv:
    def __init__(self) -> None:
        self.value = FNV_BASIS

    def add(self, word: int) -> None:
        word &= MASK64
        for shift in range(0, 64, 8):
            self.value ^= (word >> shift) & 0xFF
            self.value = (self.value * FNV_PRIME) & MASK64

    def add_i128(self, value: int) -> None:
        bits = value % (1 << 128)
        self.add(bits)
        self.add(bits >> 64)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_candidates(packet: Path) -> list[str]:
    excluded = {
        "SHA256SUMS",
        "verify_endpoint587_cleanroom_packet.out",
        "verify_endpoint587_cleanroom_packet_opt.out",
    }
    return sorted(
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file()
        and path.name not in excluded
        and path.suffix.lower() != ".exe"
    )


def verify_manifest(packet: Path) -> int:
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing SHA256SUMS")
    entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, "malformed manifest line")
        digest, relative = match.groups()
        require(relative not in entries, "duplicate manifest entry")
        require(not Path(relative).is_absolute() and ".." not in Path(relative).parts,
                "unsafe manifest path")
        entries[relative] = digest
    candidates = manifest_candidates(packet)
    require(sorted(entries) == candidates, "manifest coverage changed")
    for relative, digest in entries.items():
        require(sha(packet / relative) == digest,
                f"manifest digest mismatch: {relative}")
    return len(entries)


def fnv_words(words: list[int] | tuple[int, ...]) -> int:
    ledger = Fnv()
    for word in words:
        ledger.add(word)
    return ledger.value


def read_pairs(path: Path) -> tuple[tuple[int, int], ...]:
    result: list[tuple[int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split(",")
        require(len(fields) == 2, f"malformed pair file: {path.name}")
        pair = (int(fields[0]), int(fields[1]))
        require(0 < pair[0] < pair[1], "invalid pair")
        require(not result or result[-1] < pair, "pair file not strictly sorted")
        result.append(pair)
    return tuple(result)


def pair_fnv(pairs: tuple[tuple[int, int], ...]) -> int:
    ledger = Fnv()
    for q, r in pairs:
        ledger.add(q); ledger.add(r)
    return ledger.value


def read_masks(path: Path) -> tuple[int, ...]:
    masks = tuple(int(token, 16) for token in path.read_text(encoding="ascii").split())
    require(len(set(masks)) == len(masks), "duplicate mask")
    return masks


def verify_inputs(packet: Path) -> None:
    for relative, expected in EXPECTED_INPUT_SHA.items():
        require(sha(packet / relative) == expected, f"input SHA changed: {relative}")
    carrier = read_masks(packet / "inputs/inherited_carrier3925.txt")
    joint = read_masks(packet / "inputs/joint421_masks.txt")
    require(len(carrier) == 3925 and fnv_words(carrier) == 0xEEAE5518D84CCAC5,
            "carrier identity changed")
    require(sum(mask.bit_count() == 8 for mask in carrier) == 3809 and
            sum(mask.bit_count() == 9 for mask in carrier) == 116,
            "carrier rank census changed")
    require(len(joint) == 421 and all(mask.bit_count() == 8 for mask in joint) and
            fnv_words(joint) == 0x20D63DD42FE8150E,
            "joint deck identity changed")
    require(set(joint) <= set(carrier), "joint deck not retained")


def parse_failures(path: Path) -> tuple[tuple[int, int, int], ...]:
    with path.open("r", encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        rows = tuple((int(row["q"]), int(row["r"]), int(row["body_hex"], 16))
                     for row in reader)
    require(all(body.bit_count() == 9 for _, _, body in rows),
            "failure body rank changed")
    return rows


def verify_carrier(packet: Path) -> tuple[int, ...]:
    require((packet / "carrier_pair_O2.csv").read_bytes() ==
            (packet / "carrier_pair_O3.csv").read_bytes(),
            "O2/O3 pair ledgers differ")
    require((packet / "carrier_failures_O2.csv").read_bytes() ==
            (packet / "carrier_failures_O3.csv").read_bytes(),
            "O2/O3 failure ledgers differ")
    require((packet / "carrier_O2.out").read_bytes() ==
            (packet / "carrier_O3.out").read_bytes(),
            "O2/O3 carrier transcripts differ")

    with (packet / "carrier_pair_O2.csv").open("r", encoding="ascii", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    expected_q = (24, 25, 50, 72, 96, 100, 105, 192, 210, 260)
    require(tuple(int(row["q"]) for row in rows) == expected_q and
            all(int(row["r"]) == 587 for row in rows),
            "carrier pair rows changed")
    pair_hash = Fnv()
    for row in rows:
        for name in ("q", "r", "active"):
            pair_hash.add(int(row[name]))
        pair_hash.add(int(row["active_fnv"], 16))
        for name in ("active_joint", "active_nonjoint", "exposed"):
            pair_hash.add(int(row[name]))
        pair_hash.add(int(row["exposed_fnv"], 16))
        for name in ("minimum_hits", "maximum_hits", "failures"):
            pair_hash.add(int(row[name]))
        pair_hash.add(int(row["failure_fnv"], 16))
    require(sum(int(row["exposed"]) for row in rows) == 53771,
            "exposed total changed")
    require(sum(int(row["failures"]) for row in rows) == 41,
            "failure count changed")
    require(pair_hash.value == 0x0062CC50BE726E54,
            "pair ledger FNV changed")

    failures = parse_failures(packet / "carrier_failures_O2.csv")
    require(len(failures) == 41 and all((q, r) == (50, 587)
                                       for q, r, _ in failures),
            "failure fibre changed")
    bodies = tuple(body for _, _, body in failures)
    require(fnv_words(bodies) == 0xC76719CED1D5C52B,
            "q50 body FNV changed")
    global_hash = Fnv()
    for q, r, body in failures:
        global_hash.add(q); global_hash.add(r); global_hash.add(body)
    require(global_hash.value == 0xBAE65DC3D3BD34D0,
            "global failure FNV changed")

    transcript = (packet / "carrier_O2.out").read_text(encoding="ascii")
    for needle in (
        "BODY_TESTS 143071500",
        "SUMMARY EXPOSED 53771 HIT_INCIDENCES 1587855 FAILURES 41 FAILED_ROWS 1",
        "FAILURE_FNV bae65dc3d3bd34d0 PAIR_LEDGER_FNV 62cc50be726e54",
        "VERDICT HOSTILE_FAIL",
    ):
        require(needle in transcript, f"carrier transcript lost: {needle}")
    return bodies


def verify_wall(packet: Path, bodies: tuple[int, ...]) -> None:
    require((packet / "wall_detail_normal.csv").read_bytes() ==
            (packet / "wall_detail_opt.csv").read_bytes(),
            "normal/-O wall details differ")
    require((packet / "wall_normal.out").read_bytes() ==
            (packet / "wall_opt.out").read_bytes(),
            "normal/-O wall transcripts differ")
    require(sha(packet / "wall_detail_normal.csv") ==
            "5a89b23b60fce1ec4be6d03bdd92c4a3cf3a861afe3d30226f8bfe831154c895",
            "wall detail SHA changed")

    with (packet / "wall_detail_normal.csv").open("r", encoding="ascii", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    require(len(rows) == 41, "wall detail row count changed")
    require(tuple(int(row["body_hex"], 16) for row in rows) == bodies,
            "wall detail body order changed")
    grid = 53537802887368800
    detail_hash = Fnv()
    equality = 0
    low_values: list[tuple[int, int]] = []
    full_values: list[tuple[int, int]] = []
    for ordinal, row in enumerate(rows):
        q = int(row["q"]); r = int(row["r"])
        body = int(row["body_hex"], 16)
        low_mass = int(row["truncated_mass"])
        full_mass = int(row["full_mass"])
        low_ticks = int(row["truncated_scaled_ticks"])
        full_ticks = int(row["full_scaled_ticks"])
        require((q, r, int(row["ordinal"])) == (50, 587, ordinal),
                "wall detail ordinal changed")
        require(low_mass <= full_mass and low_ticks == 63 * low_mass - 4 * grid and
                full_ticks == 63 * full_mass - 4 * grid,
                "wall mass/tick identity changed")
        require(low_ticks > 0 and full_ticks > 0,
                "wall positivity changed")
        equality += low_mass == full_mass
        low_values.append((low_ticks, body))
        full_values.append((full_ticks, body))
        for word in (q, r, ordinal, body, low_mass, full_mass):
            detail_hash.add(word)
        detail_hash.add_i128(low_ticks)
        detail_hash.add_i128(full_ticks)
    require(equality == 31, "wall equality count changed")
    require(min(low_values) == (283424219270897292, 0x11186405),
            "truncated minimum changed")
    require(min(full_values) == (283424219270897292, 0x11186405),
            "full minimum changed")
    require(detail_hash.value == 0x7F216C622B96B321,
            "wall detail FNV changed")
    transcript = (packet / "wall_normal.out").read_text(encoding="ascii")
    for needle in (
        "GRID 53537802887368800 CELLS 8385 PAIR_SAFE_CELLS 5792",
        "ALL_CLASSES 2420 LOW_CLASSES 2365",
        "ALL_CLASS_FNV b24c8782fd3e3c7c LOW_CLASS_FNV 4145ef7ae527e0e0",
        "LOW_POSITIVE 41 FULL_POSITIVE 41 EQUALITY 31",
        "SUMMARY_FNV d8d467d453d1f9f0",
        "VERDICT PASS",
    ):
        require(needle in transcript, f"wall transcript lost: {needle}")


def verify_typed(packet: Path) -> None:
    require((packet / "typed_normal.out").read_bytes() ==
            (packet / "typed_opt.out").read_bytes(),
            "normal/-O typed transcripts differ")
    names = (
        "residual_top587.csv", "typed_union2217.csv",
        "final_residual20430.csv", "residual_top586.csv",
    )
    for name in names:
        require((packet / "typed_normal" / name).read_bytes() ==
                (packet / "typed_opt" / name).read_bytes(),
                f"normal/-O typed file differs: {name}")

    universe = read_pairs(packet / "inputs/universe22647.csv")
    prior_union = read_pairs(packet / "inputs/typed_union2207.csv")
    prior_residual = read_pairs(packet / "inputs/final_residual20440.csv")
    top587 = read_pairs(packet / "typed_normal/residual_top587.csv")
    union = read_pairs(packet / "typed_normal/typed_union2217.csv")
    residual = read_pairs(packet / "typed_normal/final_residual20430.csv")
    top586 = read_pairs(packet / "typed_normal/residual_top586.csv")
    require(len(universe) == 22647 and pair_fnv(universe) == 0xDF5374D4ACA67677,
            "typed universe changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior typed partition changed")
    require(top587 == tuple(pair for pair in prior_residual if pair[1] == 587) and
            len(top587) == 10 and pair_fnv(top587) == 0xF48CA5F1904D6F52,
            "derived top587 changed")
    require(set(union) == set(prior_union) | set(top587) and
            set(residual) == set(prior_residual) - set(top587),
            "typed successor set identity changed")
    require(len(union) == 2217 and pair_fnv(union) == 0xE6592CBEF9B616D8 and
            sha(packet / "typed_normal/typed_union2217.csv") ==
            "d465a2f62c77ddaf921e7f3d7f32c674ea45e46fcdc348c90c711d7ba8a7a6b6",
            "typed union successor changed")
    require(len(residual) == 20430 and pair_fnv(residual) == 0x4710F750DFCF91EA and
            sha(packet / "typed_normal/final_residual20430.csv") ==
            "8203dcfd6cc26b67bfdb648c7d8b50f94d7af1ab7ddbd6af2ff68acb15941f0b",
            "typed residual successor changed")
    require(top586 == tuple(pair for pair in residual if pair[1] == 586) and
            len(top586) == 12 and pair_fnv(top586) == 0xA1B617FAA2E7F63F and
            sha(packet / "typed_normal/residual_top586.csv") ==
            "d38b7fd9ea2aa9afdd12446d646cc8a9466cc5d4429612f03c8dff3165edf7ea",
            "next top586 changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet.resolve()
    manifest_files = verify_manifest(packet)
    verify_inputs(packet)
    bodies = verify_carrier(packet)
    verify_wall(packet, bodies)
    verify_typed(packet)
    print("ENDPOINT587_CLEANROOM_PACKET_VERIFY PASS")
    print(f"manifest_files={manifest_files}")
    print("rows=10 body_tests=143071500 failures=41 failure_rows=1")
    print("q50_min_ticks=283424219270897292 equality=31")
    print("typed_union=2217 residual=20430 next_endpoint=586 next_rows=12")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
