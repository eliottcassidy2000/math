"""Hardened verifier for the finite rank-two arbitrary-outsider packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from fractions import Fraction
from pathlib import Path

POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}
MASK64 = (1 << 64) - 1
MASK128 = (1 << 128) - 1


def require(value: bool, message: str) -> None:
    if not value:
        raise AssertionError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_add(state: int, word: int) -> int:
    word &= MASK64
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * 0x100000001B3) & MASK64
    return state


def fnv_add_i128(state: int, value: int) -> int:
    bits = value & MASK128
    state = fnv_add(state, bits)
    return fnv_add(state, bits >> 64)


def pair_fnv(rows: list[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for q, r in rows:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return state


def read_pairs(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split(",")
        require(len(fields) == 2, f"malformed pair row in {path.name}")
        rows.append((int(fields[0]), int(fields[1])))
    require(len(rows) == len(set(rows)), f"duplicate pair in {path.name}")
    return rows


def read_degree(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    expected = [
        "q", "r", "grid", "rank2_total", "degree_bound_mass",
        "degree_bound_ticks", "positive", "top9_hex",
    ]
    require(rows and list(rows[0]) == expected, f"degree header changed in {path.name}")
    return rows


def degree_semantics(rows: list[dict[str, str]]) -> tuple[dict[tuple[int, int], tuple[str, ...]], int]:
    keys: list[tuple[int, int]] = []
    mapping: dict[tuple[int, int], tuple[str, ...]] = {}
    state = 0xCBF29CE484222325
    for row in rows:
        q, r = int(row["q"]), int(row["r"])
        grid = int(row["grid"])
        total = int(row["rank2_total"])
        lower = int(row["degree_bound_mass"])
        ticks = int(row["degree_bound_ticks"])
        positive = int(row["positive"])
        top9 = int(row["top9_hex"], 16)
        key = (q, r)
        require(q > 0 and q < r and key not in mapping, "bad/duplicate degree pair")
        require(grid > 0 and total >= 0, "nonpositive grid/total")
        require(ticks == 63 * lower - 4 * grid, f"tick normalization failed at {key}")
        require(positive in (0, 1) and positive == (ticks > 0), f"sign flag failed at {key}")
        require(top9.bit_count() == 9 and top9 < 1 << 30, f"top-nine mask failed at {key}")
        keys.append(key)
        value = tuple(row[field] for field in row)
        mapping[key] = value
        state = fnv_add(state, q)
        state = fnv_add(state, r)
        state = fnv_add_i128(state, grid)
        state = fnv_add_i128(state, total)
        state = fnv_add_i128(state, lower)
        state = fnv_add_i128(state, ticks)
        state = fnv_add(state, positive)
        state = fnv_add(state, top9)
    require(keys == sorted(keys, key=lambda pair: (-pair[1], pair[0])), "degree row order changed")
    return mapping, state


ALLBODY = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) CELL_COUNTS (\d+),(\d+),(\d+) "
    r"RANK2_TOTAL (\d+) BODIES (\d+) POSITIVE (\d+) ZERO (\d+) "
    r"MIN_TICKS (-?\d+) MIN_BODY ([0-9a-f]{8}) MAX_TICKS (-?\d+) "
    r"MAX_BODY ([0-9a-f]{8}) DEGREE_BOUND_TICKS (-?\d+)$"
)


def read_allbody(path: Path) -> dict[tuple[int, int], dict[str, int]]:
    result: dict[tuple[int, int], dict[str, int]] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        match = ALLBODY.match(line)
        if not match:
            continue
        values = [int(value, 16) if index in (11, 13) else int(value)
                  for index, value in enumerate(match.groups())]
        q, r, grid = values[0], values[1], values[2]
        key = (q, r)
        require(key not in result, f"duplicate allbody pair {key}")
        result[key] = {
            "grid": grid,
            "bodies": values[7],
            "positive": values[8],
            "zero": values[9],
            "min_ticks": values[10],
            "min_body": values[11],
            "max_ticks": values[12],
            "max_body": values[13],
            "degree_ticks": values[14],
        }
    return result


BRANCH = re.compile(
    r"^PAIR (\d+),(\d+) GRID (\d+) BOUND_EXACT 1 MIN_TICKS (-?\d+) "
    r"MIN_BODY ([0-9a-f]{8}) NODES (\d+) PRUNES (\d+)$"
)


def read_branch(path: Path) -> dict[tuple[int, int], tuple[int, int, int]]:
    result: dict[tuple[int, int], tuple[int, int, int]] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        match = BRANCH.match(line)
        if not match:
            continue
        q, r, grid, ticks, body, nodes, prunes = match.groups()
        key = (int(q), int(r))
        require(key not in result, f"duplicate branch pair {key}")
        require(int(nodes) > 0 and int(prunes) > 0, f"empty branch proof {key}")
        result[key] = (int(grid), int(ticks), int(body, 16))
    return result


def verify_manifest(packet: Path) -> None:
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing manifest")
    listed: set[str] = set()
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(name not in listed, "duplicate manifest entry")
        listed.add(name)
        require(sha256(packet / name) == digest, f"manifest mismatch: {name}")
    actual = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and path.suffix != ".exe"
    }
    require(listed == actual, "manifest closure changed")


def verify_cleanroom_manifest(root: Path) -> None:
    manifest = root / "SHA256SUMS.txt"
    require(sha256(manifest) == "1fdd6025bae13f9bca08eb4fb97a82f149d845522e6dd6d9ac7d70629af79302",
            "clean-room manifest identity")
    listed: set[str] = set()
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(name not in listed, "duplicate clean-room manifest entry")
        listed.add(name)
        require(sha256(root / name) == digest, f"clean-room manifest mismatch: {name}")
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    }
    require(len(listed) == 29 and listed == actual, "clean-room manifest closure changed")
    require(sha256(root / "results/thm4231_coarse_screen_O3.csv") ==
            "57826d842ee968251696be4757cf94b0ff57f55f3cc88f8020e2a27ccbd19258",
            "clean-room broad screen identity")
    require(sha256(root / "results/thm4231_hostile_exact_O3.csv") ==
            "abf2ba0b0a0551ebf0ac89b77347d5e1ce8ea96dd5007bd5064cfd8c174e8039",
            "clean-room exact107 identity")
    crosscheck = (root / "results/primary_crosscheck.out").read_text(encoding="ascii")
    require("SCREEN_PAIRS 181194 FIELDWISE_MATCHES 181194" in crosscheck and
            "EXACT_PAIRS 107 MINIMUM_MATCHES 107" in crosscheck and
            "INDEPENDENT_GRAPH_PATH EVENT_STATE_MERGE" in crosscheck and
            "PRIMARY_GRAPH_PATH SORTED_WALL_MIDPOINT_CLASSIFICATION" in crosscheck and
            "VERDICT PASS" in crosscheck,
            "clean-room/primary crosscheck")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    packet = args.packet.resolve()
    verify_manifest(packet)
    verify_cleanroom_manifest(packet / "independent_cleanroom")

    remainder_path = packet / "thm4231_remainder181194.csv"
    universe_path = packet / "universe22647.csv"
    require(sha256(remainder_path) == "9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1", "remainder SHA")
    require(sha256(universe_path) == "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317", "universe SHA")
    remainder = read_pairs(remainder_path)
    universe = read_pairs(universe_path)
    require(len(remainder) == 181_194 and pair_fnv(remainder) == 0x3874FECAC4ECBD8A, "THM4231 remainder identity")
    require(remainder == sorted(remainder), "THM4231 remainder order")
    require(max(r for _, r in remainder) == 769, "THM4231 endpoint")
    require(all(q not in POOL and r not in POOL for q, r in remainder), "remainder entered pool")
    require(len(universe) == 22_647 and pair_fnv(universe) == 0xDF5374D4ACA67677, "THM4287 universe identity")
    require(universe == sorted(universe), "THM4287 universe order")
    require(set(universe) < set(remainder) and len(set(remainder) - set(universe)) == 158_547, "universe provenance arithmetic")

    allbody_source = (packet / "rank2_allbody_audit.cpp").read_text(encoding="ascii")
    require("using i128 = __int128_t;" in allbody_source and "using i64 = i128;" in allbody_source,
            "wide signed-i128 arithmetic domain")
    require("std::vector<i64> walls = {0, graph.grid};" in allbody_source and
            "walls.erase(std::unique(walls.begin(), walls.end()), walls.end());" in allbody_source and
            "const i64 width = right - left;" in allbody_source,
            "finite-null-wall cell convention")
    require(not (packet / "rank2_degree_scan_root.cpp").exists(), "obsolete signed-i64 scanner entered packet")

    full_o2 = packet / "thm4231_remainder181194_degree_O2.csv"
    full_o3 = packet / "thm4231_remainder181194_degree_O3.csv"
    require(full_o2.read_bytes() == full_o3.read_bytes(), "full O2/O3 degree mismatch")
    require(sha256(full_o3) == "5b008404482b6e23006c0c1cb97407c22ddf34eb7cf4e09fd5f604206d8a0356", "full degree SHA")
    full_rows = read_degree(full_o3)
    full_map, full_fnv = degree_semantics(full_rows)
    require(set(full_map) == set(remainder), "degree/remainder pair mismatch")
    require(full_fnv == 0xF8DA82C5A6D732ED, "wide full degree FNV")
    coarse_good = {key for key, row in full_map.items() if row[6] == "1"}
    hostile = set(full_map) - coarse_good
    require(len(coarse_good) == 181_087 and len(hostile) == 107, "coarse partition")
    highest = max(r for _, r in hostile)
    require(highest == 554 and {key for key in hostile if key[1] == highest} == {(50, 554)}, "highest coarse hostile")
    overflow = [key for key, row in full_map.items() if int(row[2]) > (1 << 63) - 1]
    require(overflow == [(713, 719)], "signed-i64 overflow control pair")
    require(int(full_map[(713, 719)][2]) == 9_351_275_651_380_222_560, "signed-i64 overflow control grid")

    u_degree = packet / "universe22647_degree_wide_O3.csv"
    require(sha256(u_degree) == "4e61fb2f0944891d95fc8e44e19813226d9a771bf786d7f78369fb78602564ad", "wide U degree SHA")
    u_map, u_fnv = degree_semantics(read_degree(u_degree))
    require(u_fnv == 0x611B1BA5C25594DD and u_fnv != 0xFE6111D297A72A3D, "wide U FNV / stale-FNV rejection")
    require(set(u_map) == set(universe), "wide U pair set")
    require(all(u_map[key] == full_map[key] for key in u_map), "U/full pair-local graph disagreement")
    require({key for key, row in u_map.items() if row[6] == "0"} == hostile, "all coarse hostiles must survive in U")

    primary_o2 = packet / "thm4231_bad107_allbody_O2.out"
    primary_o3 = packet / "thm4231_bad107_allbody_O3.out"
    require(primary_o2.read_bytes() == primary_o3.read_bytes(), "allbody O2/O3 mismatch")
    require(sha256(primary_o3) == "31d0d05988eeb7c789dc5fd3b23e5c54e2f806c026f407de2c5fa6d781193a81", "allbody SHA")
    primary = read_allbody(primary_o3)
    require(set(primary) == hostile, "allbody candidate set")
    for key, row in primary.items():
        require(row["grid"] == int(full_map[key][2]), f"allbody grid {key}")
        require(row["bodies"] == 14_307_150 and row["positive"] == 14_307_150 and row["zero"] == 0, f"allbody quantifier {key}")
        require(row["min_ticks"] > 0 and row["min_body"].bit_count() == 9, f"allbody strict/rank {key}")
        require(row["degree_ticks"] == int(full_map[key][5]), f"allbody degree echo {key}")

    branch_o2 = packet / "thm4231_bad107_branch_O2.out"
    branch_o3 = packet / "thm4231_bad107_branch_O3.out"
    require(branch_o2.read_bytes() == branch_o3.read_bytes(), "branch O2/O3 mismatch")
    require(sha256(branch_o3) == "ed19a104828e15386ded7579d234f08aa6b39b88eb3daef0ad3da25dae998246", "branch SHA")
    branch = read_branch(branch_o3)
    require(set(branch) == hostile, "branch candidate set")
    for key, (grid, ticks, body) in branch.items():
        require((grid, ticks, body) == (primary[key]["grid"], primary[key]["min_ticks"], primary[key]["min_body"]), f"branch/flat disagreement {key}")

    candidate_min = min(hostile, key=lambda key: Fraction(primary[key]["min_ticks"], primary[key]["grid"]))
    require(candidate_min == (50, 70), "candidate normalized minimizer")
    require(primary[candidate_min]["min_ticks"] == 245_428_469_244 and primary[candidate_min]["min_body"] == 0x031C7400, "candidate minimum data")
    candidate_ratio = Fraction(primary[candidate_min]["min_ticks"], primary[candidate_min]["grid"])
    require(candidate_ratio == Fraction(973_922_497, 361_927_766_200), "candidate ratio")
    require([key for key in hostile if Fraction(primary[key]["min_ticks"], primary[key]["grid"]) == candidate_ratio] == [(50, 70)],
            "candidate normalized minimizer uniqueness")

    threshold_o2 = packet / "normalized_threshold3_allbody_O2.out"
    threshold_o3 = packet / "normalized_threshold3_allbody_O3.out"
    require(threshold_o2.read_bytes() == threshold_o3.read_bytes(), "threshold O2/O3 mismatch")
    require(sha256(threshold_o3) == "b2361810d025589002c4795615fbe0076d750037f8fd3a1d43a1435e51443cee", "threshold SHA")
    threshold = read_allbody(threshold_o3)
    weak_coarse = {
        key for key in coarse_good
        if Fraction(int(full_map[key][5]), int(full_map[key][2])) <= candidate_ratio
    }
    require(weak_coarse == {(50, 212), (50, 274), (100, 110)} == set(threshold), "normalized threshold set")
    require(all(Fraction(row["min_ticks"], row["grid"]) > candidate_ratio for row in threshold.values()), "threshold exact minima")

    # Every other coarse-good row lies above candidate_ratio by its rigorous
    # degree bound; the three weaker bounds were exactly optimized above.
    # Hence (50,70), 031c7400 is the intrinsic normalized global L2 minimizer.
    require(Fraction(primary[(50, 70)]["min_ticks"], 63 * primary[(50, 70)]["grid"]) == Fraction(973_922_497, 22_801_449_270_600), "rank-two certified surplus fraction")

    full_out_o2 = (packet / "thm4231_remainder181194_degree_O2.out").read_text(encoding="ascii")
    full_out_o3 = (packet / "thm4231_remainder181194_degree_O3.out").read_text(encoding="ascii")
    require(full_out_o2 == full_out_o3, "full degree stdout mismatch")
    require("PAIRS 181194 POSITIVE 181087 NONPOSITIVE 107" in full_out_o3 and "LEDGER_FNV f8da82c5a6d732ed" in full_out_o3, "full stdout summary")
    u_out = (packet / "universe22647_degree_wide_O3.out").read_text(encoding="ascii")
    require("PAIRS 22647 POSITIVE 22540 NONPOSITIVE 107" in u_out and "LEDGER_FNV 611b1ba5c25594dd" in u_out and "fe6111d297a72a3d" not in u_out, "wide U stdout")

    print("RANK2_ARBITRARY_PAIR_PACKET_VERIFY PASS")
    print("thm4231_remainder=181194 coarse_positive=181087 exact_fallback=107")
    print("global_normalized_min_pair=50,70 body=031c7400 ticks=245428469244 grid=91205797082400")
    print("rank2_certified_surplus=973922497/22801449270600 actual_haar_surplus_at_least_same")
    print("overflow_control=713,719 grid=9351275651380222560")
    print("consequence=fixed_pool_arbitrary_two_outsider_haar_theorem_relative_to_THM4231")
    print("scope=no_physical_entry_no_LRC14")


if __name__ == "__main__":
    main()
