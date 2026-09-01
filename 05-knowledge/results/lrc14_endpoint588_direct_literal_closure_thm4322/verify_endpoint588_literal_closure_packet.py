#!/usr/bin/env python3
"""Hardened verifier for the frozen endpoint-588 literal-closure packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
MASK128 = (1 << 128) - 1
EXPECTED_DETAIL_SHA = "8c2652b9cba7cfdacd4d6e908334220eb08445306fbc09db816d0b0b9c48664d"

EXPECTED_FAILURES = {
    25: (9, 0x8FE6C7431C54C6FF),
    50: (99836, 0x4A8D6E89F7AA5D75),
    96: (1130, 0x5766A653739CF191),
    100: (14, 0xF0FDDA5AD80A46B7),
    105: (60, 0xA6F76A994FD08192),
    206: (7, 0xDEF96E7FAE5DB299),
    210: (43799, 0xF756DC445E2790CD),
    256: (7, 0x392C4D67394B2E17),
    420: (1, 0x664E6DEF29193BD0),
    462: (4, 0x5F993DEACDDBF9D3),
}

GRIDS = {
    25: 638440579576800,
    50: 638440579576800,
    96: 255376231830720,
    100: 638440579576800,
    105: 127688115915360,
    206: 13151875939282080,
    210: 127688115915360,
    256: 2043009854645760,
    420: 127688115915360,
    462: 127688115915360,
}

EXPECTED_MIN = {
    25: (3845468183436180, 0x0D1CD000),
    50: (3359314066788540, 0x003C6405),
    96: (1528699623231060, 0x2D14C001),
    100: (3874870450697778, 0x36704001),
    105: (731279928922452, 0x0458D082),
    206: (81246765687627708, 0x14087C01),
    210: (706240375488648, 0x06B82090),
    256: (13754370514849146, 0x1F106001),
    420: (795987916404408, 0x07584088),
    462: (890738628297084, 0x1E582001),
}

EXPECTED_DETAIL_FNV = {
    25: 0x86377E6053C258A5,
    50: 0xA020B250B51F1D2C,
    96: 0x6529AFCA0B9F4CC7,
    100: 0x5922C0795B22C80F,
    105: 0x1288AFDB9761C138,
    206: 0xCA6C2268055BBF8D,
    210: 0x84B74FEE2DFDC6F0,
    256: 0x048A2B90E53E964D,
    420: 0xCCD2B716EEE3D861,
    462: 0x0AF422412400ED12,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_words(words) -> int:
    state = OFFSET
    for word in words:
        word &= MASK64
        for shift in range(0, 64, 8):
            state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def add_word(state: int, word: int) -> int:
    word &= MASK64
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def verify_manifest(packet: Path) -> int:
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing SHA256SUMS")
    entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and name not in entries,
                "malformed/duplicate manifest row")
        entries[name] = digest
    actual = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and
        path.suffix.lower() != ".exe" and "__pycache__" not in path.parts
    }
    require(set(entries) == actual, "manifest closure mismatch")
    for name, digest in entries.items():
        require(sha256(packet / name) == digest,
                f"manifest digest mismatch: {name}")
    return len(entries)


def read_pair_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 66 and [int(row["q"]) for row in rows] ==
            sorted(int(row["q"]) for row in rows), "pair row census/order changed")
    require(all(int(row["r"]) == 588 for row in rows),
            "pair row escaped endpoint588")
    return rows


def verify_pair_ledger(packet: Path) -> None:
    for suffix in ("pair_O2.csv", "failures_O2.csv", "raw_O2.out"):
        left = packet / f"endpoint588_{suffix}"
        right = packet / f"endpoint588_{suffix.replace('O2', 'O3')}"
        require(left.read_bytes() == right.read_bytes(),
                f"O2/O3 mismatch: {suffix}")
    require(sha256(packet / "endpoint588_pair_O2.csv") ==
            "93e05caa1e3a7502cf11be2850660bf111ed475c6b034f9b8aab9db014efa959",
            "pair CSV SHA changed")
    require(sha256(packet / "endpoint588_failures_O2.csv") ==
            "45e7cc8d3b7644a90b829ee251746a5c0a3d47747f82a0c165a971d789563af0",
            "failure CSV SHA changed")
    rows = read_pair_rows(packet / "endpoint588_pair_O2.csv")
    words: list[int] = []
    total_exposed = total_hits = total_failures = failed_rows = 0
    observed: dict[int, tuple[int, int]] = {}
    for row in rows:
        q = int(row["q"]); r = int(row["r"])
        values = [q, r, int(row["active"]), int(row["active_fnv"], 16),
                  int(row["active_joint"]), int(row["active_nonjoint"]),
                  int(row["exposed"]), int(row["exposed_fnv"], 16),
                  int(row["minimum_hits"]), int(row["maximum_hits"]),
                  int(row["failures"]), int(row["failure_fnv"], 16)]
        words.extend(values)
        total_exposed += values[6]
        total_hits_line = None
        # Hit incidences are intentionally an aggregate-only field in stdout.
        total_failures += values[10]
        failed_rows += values[10] != 0
        if values[10]: observed[q] = (values[10], values[11])
    require(total_exposed == 5001257 and total_failures == 144867 and
            failed_rows == 10 and observed == EXPECTED_FAILURES,
            "pair aggregate/failure census changed")
    require(fnv_words(words) == 0x204A72794170E186,
            "pair-ledger FNV changed")
    raw = (packet / "endpoint588_raw_O2.out").read_text(encoding="ascii")
    require("HIT_INCIDENCES 76551030" in raw and
            "FAILURE_FNV 6f51fa88f3b09cdc" in raw and
            "BODY_TESTS 944271900" in raw and "VERDICT HOSTILE_FAIL" in raw,
            "raw carrier summary changed")


def read_failures(path: Path) -> dict[int, list[int]]:
    by_q: dict[int, list[int]] = defaultdict(list)
    global_words: list[int] = []
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in reader:
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require(r == 588 and body.bit_count() == 9,
                    "failure row escaped universe")
            by_q[q].append(body)
            global_words.extend((q, r, body))
    observed = {q: (len(bodies), fnv_words(bodies)) for q, bodies in by_q.items()}
    require(observed == EXPECTED_FAILURES and
            fnv_words(global_words) == 0x6F51FA88F3B09CDC,
            "failure body FNV changed")
    require(all(bodies == sorted(set(bodies)) for bodies in by_q.values()),
            "failure body ordering/distinctness changed")
    return dict(by_q)


def verify_literal(packet: Path) -> None:
    primary2 = packet / "endpoint588_direct_primary_O2.csv"
    primary3 = packet / "endpoint588_direct_primary_O3.csv"
    independent = packet / "endpoint588_direct_rawcell_independent_O3.csv"
    require(primary2.read_bytes() == primary3.read_bytes() ==
            independent.read_bytes(), "literal implementations disagree")
    require(sha256(primary2) == EXPECTED_DETAIL_SHA, "literal detail SHA changed")
    failures = read_failures(packet / "endpoint588_failures_O2.csv")
    seen: dict[int, list[int]] = defaultdict(list)
    detail_states = {q: OFFSET for q in EXPECTED_FAILURES}
    minima: dict[int, tuple[int, int]] = {}
    positive = equality = 0
    with primary2.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["q", "r", "ordinal", "body_hex",
                "truncated_mass", "full_mass", "truncated_scaled_ticks",
                "full_scaled_ticks"], "literal detail header changed")
        for row in reader:
            q = int(row["q"]); r = int(row["r"]); ordinal = int(row["ordinal"])
            body = int(row["body_hex"], 16)
            low = int(row["truncated_mass"]); full = int(row["full_mass"])
            low_ticks = int(row["truncated_scaled_ticks"])
            full_ticks = int(row["full_scaled_ticks"])
            require(r == 588 and q in GRIDS and ordinal == len(seen[q]) and
                    body == failures[q][ordinal], "detail/failure alignment changed")
            require(0 <= low <= full and low_ticks == 63 * low - 4 * GRIDS[q] and
                    full_ticks == 63 * full - 4 * GRIDS[q] and low_ticks > 0,
                    "literal inequality/arithmetic changed")
            seen[q].append(body); positive += 1; equality += low == full
            candidate = (low_ticks, body)
            minima[q] = min(minima.get(q, candidate), candidate)
            state = detail_states[q]
            for word in (q, r, ordinal, body, low, full,
                         low_ticks & MASK64, (low_ticks & MASK128) >> 64,
                         full_ticks & MASK64, (full_ticks & MASK128) >> 64):
                state = add_word(state, word)
            detail_states[q] = state
    require(dict(seen) == failures and positive == 144867 and equality == 109452,
            "literal detail census changed")
    require(minima == EXPECTED_MIN and detail_states == EXPECTED_DETAIL_FNV,
            "literal minima/detail FNV changed")
    out2 = (packet / "endpoint588_direct_primary_O2.out").read_bytes()
    out3 = (packet / "endpoint588_direct_primary_O3.out").read_bytes()
    require(out2 == out3 and b"LOW_POSITIVE 144867 FULL_POSITIVE 144867" in out2 and
            b"VERDICT PASS" in out2, "primary literal summary changed")
    raw_out = (packet / "endpoint588_direct_rawcell_independent_O3.out").read_bytes()
    require(b"LOW_POSITIVE 144867 FULL_POSITIVE 144867" in raw_out and
            b"NO_CLASS_AGGREGATION" in raw_out and b"VERDICT PASS" in raw_out,
            "independent raw-cell summary changed")


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical row ledger: {path.name}")
    return rows


def rows_fnv(rows: list[tuple[int, int]]) -> int:
    return fnv_words(word for row in rows for word in row)


def verify_typed(packet: Path) -> None:
    normal = packet / "post588_typed"
    optimized = packet / "post588_typed_opt"
    files = ["typed_union2207.csv", "final_residual20440.csv",
             "residual_top587.csv"]
    for name in files:
        require((normal / name).read_bytes() == (optimized / name).read_bytes(),
                f"normal/-O typed mismatch: {name}")
    require((packet / "post588_typed.out").read_bytes() ==
            (packet / "post588_typed_opt.out").read_bytes(),
            "normal/-O typed stdout mismatch")
    universe = read_rows(packet / "inputs/typed/universe22647.csv")
    prior_union = read_rows(packet / "inputs/typed/typed_union2141.csv")
    prior_residual = read_rows(packet / "inputs/typed/final_residual20506.csv")
    top588 = read_rows(packet / "inputs/typed/residual_top588.csv")
    union = read_rows(normal / files[0]); residual = read_rows(normal / files[1])
    top587 = read_rows(normal / files[2])
    require(len(universe) == 22647 and len(prior_union) == 2141 and
            len(prior_residual) == 20506 and len(top588) == 66 and
            len(union) == 2207 and len(residual) == 20440 and len(top587) == 10,
            "typed counts changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe) and
            set(top588).issubset(prior_residual) and
            union == sorted(set(prior_union) | set(top588)) and
            residual == sorted(set(prior_residual) - set(top588)) and
            set(union).isdisjoint(residual) and set(union) | set(residual) == set(universe),
            "typed partition transition changed")
    require(rows_fnv(union) == 0x18D067B5614CF47F and
            rows_fnv(residual) == 0x794BD808E92E27CD and
            rows_fnv(top587) == 0xF48CA5F1904D6F52 and
            top587 == [row for row in residual if row[1] == 587],
            "typed FNV/top frontier changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet.resolve()
    require(packet.is_dir(), "packet directory missing")
    manifest_files = verify_manifest(packet)
    verify_pair_ledger(packet)
    verify_literal(packet)
    verify_typed(packet)
    print("ENDPOINT588_LITERAL_CLOSURE_PACKET_VERIFY PASS")
    print(f"manifest_files={manifest_files}")
    print("endpoint588_rows=66 carrier_failures=144867 direct_positive=144867")
    print("typed_union=2207 residual=20440 next_endpoint=587 next_rows=10")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
