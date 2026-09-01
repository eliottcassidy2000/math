#!/usr/bin/env python3
"""Hardened static verifier for the independent rank-two packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Python removes `assert` statements under -O.  Make optimized invocation a
# hardened launcher for a fresh non-optimized verifier rather than silently
# weakening the checks.
if not __debug__:
    environment = os.environ.copy()
    environment.pop("PYTHONOPTIMIZE", None)
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        env=environment,
    )
    raise SystemExit(completed.returncode)

EXPECTED_INPUT_SHA = {
    "inputs/universe22647.csv": "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    "inputs/thm4231_remainder181194.csv": "9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1",
}
SCREEN_HEADER = [
    "q", "r", "grid", "cells0", "cells1", "cells2", "mass0", "mass1",
    "mass2", "total_mass", "top9_degree_sum", "coarse_mass",
    "coarse_ticks", "coarse_positive", "top9_mask_hex",
]
EXACT_HEADER = [
    "q", "r", "grid", "coarse_ticks", "minimum_mass", "minimum_ticks",
    "minimum_body_hex", "direct_replay_mass", "identity_ok", "nodes", "prunes",
]
SELECTED_PAIR = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) COARSE_TICKS (-?\d+) "
    r"MIN_MASS (\d+) MIN_TICKS (-?\d+) MIN_BODY ([0-9a-f]{8}) "
    r"DIRECT_MASS (\d+) RATIO_CROSS_VS_50_70 (-?\d+) "
    r"NODES (\d+) PRUNES (\d+)$"
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_fnv(path: Path) -> str:
    value = 14_695_981_039_346_656_037
    lines = path.read_text(encoding="ascii").splitlines()
    assert lines
    for line in lines[1:]:
        for byte in (line + "\n").encode("ascii"):
            value ^= byte
            value = (value * 1_099_511_628_211) & ((1 << 64) - 1)
    return f"{value:016x}"


def verify_manifest() -> int:
    manifest = ROOT / "SHA256SUMS.txt"
    entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, relative = line.split("  ", 1)
        assert len(digest) == 64 and relative not in entries
        entries[relative] = digest
    actual = {
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    }
    assert actual == set(entries), (actual - set(entries), set(entries) - actual)
    for relative, digest in entries.items():
        assert sha(ROOT / relative) == digest, relative
    return len(entries)


def load_pairs(relative: str) -> list[tuple[int, int]]:
    pairs = []
    for line in (ROOT / relative).read_text(encoding="ascii").splitlines():
        fields = line.split(",")
        assert len(fields) == 2
        pair = int(fields[0]), int(fields[1])
        assert 0 < pair[0] < pair[1]
        pairs.append(pair)
    assert len(pairs) == len(set(pairs))
    return pairs


def load_screen(relative: str, input_pairs: list[tuple[int, int]]) -> list[dict[str, str]]:
    with (ROOT / relative).open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == SCREEN_HEADER
        rows = list(reader)
    assert [(int(row["q"]), int(row["r"])) for row in rows] == input_pairs
    for row in rows:
        grid = int(row["grid"])
        masses = [int(row[f"mass{rank}"]) for rank in range(3)]
        total = int(row["total_mass"])
        degree_sum = int(row["top9_degree_sum"])
        coarse_mass = int(row["coarse_mass"])
        ticks = int(row["coarse_ticks"])
        assert grid > 0 and all(value >= 0 for value in masses)
        assert sum(masses) == total
        assert total - degree_sum == coarse_mass
        assert 63 * coarse_mass - 4 * grid == ticks
        assert row["coarse_positive"] == ("1" if ticks > 0 else "0")
        assert int(row["top9_mask_hex"], 16).bit_count() == 9
    return rows


def load_exact(relative: str, bad_pairs: set[tuple[int, int]]) -> list[dict[str, str]]:
    with (ROOT / relative).open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == EXACT_HEADER
        rows = list(reader)
    assert {(int(row["q"]), int(row["r"])) for row in rows} == bad_pairs
    for row in rows:
        grid = int(row["grid"])
        mass = int(row["minimum_mass"])
        ticks = int(row["minimum_ticks"])
        body = int(row["minimum_body_hex"], 16)
        assert body.bit_count() == 9
        assert row["identity_ok"] == "1"
        assert int(row["direct_replay_mass"]) == mass
        assert 63 * mass - 4 * grid == ticks
        assert ticks > 0
        assert int(row["nodes"]) > 0 and int(row["prunes"]) >= 0
    return rows


def parse_selected(relative: str) -> dict[tuple[int, int], tuple[int, int, str]]:
    rows = {}
    text = (ROOT / relative).read_text(encoding="utf-8")
    assert "VERDICT PASS" in text
    for line in text.splitlines():
        match = SELECTED_PAIR.fullmatch(line.strip())
        if match:
            rows[(int(match.group(1)), int(match.group(2)))] = (
                int(match.group(3)), int(match.group(6)), match.group(7)
            )
    assert set(rows) == {(50, 70), (50, 212), (50, 274), (100, 110)}
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-rerun-rawcell", action="store_true")
    args = parser.parse_args()

    manifest_files = verify_manifest()
    for relative, digest in EXPECTED_INPUT_SHA.items():
        assert sha(ROOT / relative) == digest

    narrow_pairs = load_pairs("inputs/universe22647.csv")
    broad_pairs = load_pairs("inputs/thm4231_remainder181194.csv")
    assert len(narrow_pairs) == 22_647
    assert len(broad_pairs) == 181_194

    narrow = load_screen("results/coarse_screen_O3.csv", narrow_pairs)
    broad = load_screen("results/thm4231_coarse_screen_O3.csv", broad_pairs)
    assert sha(ROOT / "results/coarse_screen_O3.csv") == (
        "04dd34162cd27caf047b16294f46cb156a83b6f86805ed8493c880703fcd2847"
    )
    assert sha(ROOT / "results/thm4231_coarse_screen_O3.csv") == (
        "57826d842ee968251696be4757cf94b0ff57f55f3cc88f8020e2a27ccbd19258"
    )
    assert data_fnv(ROOT / "results/coarse_screen_O3.csv") == "1736553eb2add333"
    assert data_fnv(ROOT / "results/thm4231_coarse_screen_O3.csv") == "d7f639e550399b4c"
    assert (ROOT / "results/audit_O2.out").read_bytes() == (
        ROOT / "results/audit_O3.out"
    ).read_bytes()
    assert (ROOT / "results/thm4231_audit_O2.out").read_bytes() == (
        ROOT / "results/thm4231_audit_O3.out"
    ).read_bytes()
    broad_summary = (ROOT / "results/thm4231_audit_O3.out").read_text(
        encoding="utf-8"
    )
    assert "SCREEN_DATA_FNV d7f639e550399b4c" in broad_summary
    assert "INPUT_NORMALIZED_FNV 95425eabee50378c" in broad_summary
    parity = (ROOT / "results/optimization_parity.txt").read_text(
        encoding="ascii"
    ).lower()
    assert "broad_screen_o2_o3_sha256 57826d842ee968251696be4757cf94b0ff57f55f3cc88f8020e2a27ccbd19258" in parity
    assert "verdict pass" in parity
    narrow_bad = {
        (int(row["q"]), int(row["r"]))
        for row in narrow if row["coarse_positive"] == "0"
    }
    broad_bad = {
        (int(row["q"]), int(row["r"]))
        for row in broad if row["coarse_positive"] == "0"
    }
    assert len(narrow_bad) == len(broad_bad) == 107
    assert narrow_bad == broad_bad
    assert sum(row["coarse_positive"] == "1" for row in narrow) == 22_540
    assert sum(row["coarse_positive"] == "1" for row in broad) == 181_087
    wide = [row for row in broad if int(row["grid"]) > (1 << 63) - 1]
    assert len(wide) == 1
    assert (int(wide[0]["q"]), int(wide[0]["r"]), int(wide[0]["grid"])) == (
        713, 719, 9_351_275_651_380_222_560
    )

    narrow_exact = load_exact("results/hostile107_exact_O3.csv", narrow_bad)
    broad_exact = load_exact("results/thm4231_hostile_exact_O3.csv", broad_bad)
    assert (ROOT / "results/hostile107_exact_O3.csv").read_bytes() == (
        ROOT / "results/thm4231_hostile_exact_O3.csv"
    ).read_bytes()
    assert (ROOT / "results/hostile107_exact_O2.csv").read_bytes() == (
        ROOT / "results/hostile107_exact_O3.csv"
    ).read_bytes()
    assert (ROOT / "results/thm4231_hostile_exact_O2.csv").read_bytes() == (
        ROOT / "results/thm4231_hostile_exact_O3.csv"
    ).read_bytes()
    assert len(narrow_exact) == len(broad_exact) == 107

    # Narrow rows must be bit-for-bit numerically identical inside the broad scan.
    broad_map = {(int(row["q"]), int(row["r"])): row for row in broad}
    for row in narrow:
        other = broad_map[(int(row["q"]), int(row["r"]))]
        for field in SCREEN_HEADER:
            assert row[field] == other[field]

    selected = parse_selected("results/selected_ratio_exact_O3.out")
    assert (ROOT / "results/selected_ratio_exact_O2.out").read_bytes() == (
        ROOT / "results/selected_ratio_exact_O3.out"
    ).read_bytes()
    assert selected[(50, 70)] == (91_205_797_082_400, 245_428_469_244, "031c7400")
    assert selected[(100, 110)][1:] == (63_178_284_254_904, "04f06408")

    candidate_ticks, candidate_grid = 245_428_469_244, 91_205_797_082_400
    positive = [row for row in broad if row["coarse_positive"] == "1"]
    raw_filter = [
        (int(row["q"]), int(row["r"])) for row in positive
        if int(row["coarse_ticks"]) <= candidate_ticks
    ]
    ratio_filter = [
        (int(row["q"]), int(row["r"])) for row in positive
        if int(row["coarse_ticks"]) * candidate_grid
        <= candidate_ticks * int(row["grid"])
    ]
    assert raw_filter == [(100, 110)]
    assert ratio_filter == [(50, 212), (50, 274), (100, 110)]

    for name in (
        "results/audit_O2.out", "results/audit_O3.out",
        "results/thm4231_audit_O2.out", "results/thm4231_audit_O3.out",
        "results/global_minimum_filter.out", "results/rawcell_winner_replay.out",
        "results/primary_crosscheck.out",
    ):
        text = (ROOT / name).read_text(encoding="utf-8")
        assert "VERDICT PASS" in text, name
    assert "GLOBAL_MINIMUM_TICKS" not in (
        ROOT / "results/thm4231_audit_O3.out"
    ).read_text(encoding="utf-8")
    assert "HOSTILE107_EXACT_MINIMUM_TICKS" in (
        ROOT / "results/thm4231_audit_O3.out"
    ).read_text(encoding="utf-8")

    with (ROOT / "results/rawcell_winner_replay.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        replay = list(csv.DictReader(handle))
    assert len(replay) == 110 and all(row["match"] == "1" for row in replay)

    if not args.no_rerun_rawcell:
        with tempfile.TemporaryDirectory(prefix="lrc14-rank2-verify-") as temporary:
            output = Path(temporary) / "rawcell.csv"
            process = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "src/rawcell_winner_replay.py"),
                    str(ROOT / "results/thm4231_hostile_exact_O3.csv"),
                    str(ROOT / "results/selected_ratio_exact_O3.out"),
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            assert output.read_bytes() == (
                ROOT / "results/rawcell_winner_replay.csv"
            ).read_bytes()
            frozen = (ROOT / "results/rawcell_winner_replay.out").read_text(
                encoding="utf-8"
            ).replace("\r\n", "\n")
            assert process.stdout.replace("\r\n", "\n") == frozen

    print("LRC14_RANK2_INDEPENDENT_PACKET_VERIFY_V1")
    print(f"MANIFEST_FILES {manifest_files}")
    print("NARROW_PAIRS 22647 COARSE_POSITIVE 22540 HOSTILE 107 EXACT_POSITIVE 107")
    print("BROAD_PAIRS 181194 COARSE_POSITIVE 181087 HOSTILE 107 EXACT_POSITIVE 107")
    print("WIDE_GRID_CONTROL 713,719=9351275651380222560")
    print("RAW_FILTER 1 NORMALIZED_FILTER 3 RAWCELL_REPLAYS 110")
    print("GLOBAL_RANK2_MINIMUM 245428469244/91205797082400 AT 50,70 BODY 031c7400")
    print("VERDICT PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
