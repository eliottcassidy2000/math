#!/usr/bin/env python3
"""Hardened semantic verifier for the THM-4320 proof packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import Counter
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(str(message))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_words(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        for shift in range(0, 64, 8):
            state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def fnv_i128(state: int, word: int) -> int:
    bits = word % (1 << 128)
    return fnv_words_from(state, [bits & MASK64, bits >> 64])


def fnv_words_from(state: int, words: list[int]) -> int:
    for word in words:
        for shift in range(0, 64, 8):
            state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    return fnv_words([value for row in rows for value in row])


def read_rows(path: Path, expected: int) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)) and len(rows) == expected,
            f"noncanonical/count-changed row ledger: {path}")
    return rows


def same(left: Path, right: Path) -> None:
    require(left.read_bytes() == right.read_bytes(),
            f"byte mismatch: {left} {right}")


def verify_manifest(packet: Path) -> int:
    manifest = packet / "SHA256SUMS"
    require(manifest.is_file(), "missing SHA256SUMS")
    entries: list[tuple[str, str]] = []
    for line in manifest.read_text(encoding="ascii").splitlines():
        digest, separator, name = line.partition("  ")
        require(separator == "  " and len(digest) == 64 and name,
                f"malformed manifest line: {line}")
        require(name != "SHA256SUMS", "manifest cannot hash itself")
        path = packet / name
        require(path.is_file(), f"manifest path missing: {name}")
        require(sha256(path) == digest, f"manifest hash mismatch: {name}")
        entries.append((name, digest))
    actual = sorted(
        str(path.relative_to(packet)).replace("\\", "/")
        for path in packet.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS"
        and not path.name.endswith(".exe")
    )
    declared = sorted(name for name, _ in entries)
    require(declared == actual, "manifest closure changed")
    return len(entries)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    packet = args.packet.resolve()

    for left, right in [
        ("post590_typed.out", "post590_typed_opt.out"),
        ("post590_typed/typed_union2113.csv",
         "post590_typed_opt/typed_union2113.csv"),
        ("post590_typed/final_residual20534.csv",
         "post590_typed_opt/final_residual20534.csv"),
        ("post590_typed/residual_top589.csv",
         "post590_typed_opt/residual_top589.csv"),
        ("endpoint589_raw_O2.out", "endpoint589_raw_O3.out"),
        ("endpoint589_pair_O2.csv", "endpoint589_pair_O3.csv"),
        ("endpoint589_failures_O2.csv", "endpoint589_failures_O3.csv"),
        ("direct_primary_O2.out", "direct_primary_O3.out"),
        ("direct_primary_O2.csv", "direct_primary_O3.csv"),
        ("direct_primary_O2.csv", "direct_independent.csv"),
        ("direct_independent.csv", "direct_independent_opt.csv"),
        ("direct_independent.out", "direct_independent_opt.out"),
        ("literal_lower_bound_independent_O2.out",
         "literal_lower_bound_independent_O3.out"),
        ("literal_body_masses_independent_O2.csv",
         "literal_body_masses_independent_O3.csv"),
        ("petal_bridge.out", "petal_bridge_opt.out"),
        ("q50_active_carrier_O2.out", "q50_active_carrier_O3.out"),
        ("q50_active_carrier_O2.csv", "q50_active_carrier_O3.csv"),
        ("q50_structure.out", "q50_structure_opt.out"),
        ("q50_vertex_degrees.csv", "q50_vertex_degrees_opt.csv"),
        ("q50_hub_fibres.csv", "q50_hub_fibres_opt.csv"),
        ("q50_neither_bodies.csv", "q50_neither_bodies_opt.csv"),
        ("post589_typed.out", "post589_typed_opt.out"),
        ("post589_typed/typed_union2141.csv",
         "post589_typed_opt/typed_union2141.csv"),
        ("post589_typed/final_residual20506.csv",
         "post589_typed_opt/final_residual20506.csv"),
        ("post589_typed/residual_top588.csv",
         "post589_typed_opt/residual_top588.csv"),
    ]:
        same(packet / left, packet / right)

    old = repo / (
        "05-knowledge/results/"
        "lrc14_mixed_rank_depth_recursive_signatures_thm4296"
    )
    packet4314 = repo / (
        "05-knowledge/results/"
        "lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314"
    )
    universe = read_rows(old / "inputs/current_residual22647.csv", 22647)
    union2100 = read_rows(packet4314 / "typed/typed_union2100.csv", 2100)
    residual20547 = read_rows(
        packet4314 / "typed/final_residual20547.csv", 20547)
    top590 = read_rows(packet4314 / "typed/residual_top590.csv", 13)
    union2113 = read_rows(
        packet / "post590_typed/typed_union2113.csv", 2113)
    residual20534 = read_rows(
        packet / "post590_typed/final_residual20534.csv", 20534)
    top589 = read_rows(packet / "post590_typed/residual_top589.csv", 28)
    union2141 = read_rows(
        packet / "post589_typed/typed_union2141.csv", 2141)
    residual20506 = read_rows(
        packet / "post589_typed/final_residual20506.csv", 20506)
    top588 = read_rows(packet / "post589_typed/residual_top588.csv", 66)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677 and
            fnv_rows(union2100) == 0x3B2D991DA091A7DF and
            fnv_rows(residual20547) == 0x59CA49A11D140EC5 and
            fnv_rows(top590) == 0x44AA8A793D162CF9,
            "typed predecessor identity changed")
    require(union2113 == sorted(set(union2100) | set(top590)) and
            residual20534 == sorted(set(residual20547) - set(top590)) and
            set(union2113).isdisjoint(residual20534) and
            set(union2113) | set(residual20534) == set(universe),
            "post590 typed partition changed")
    require(fnv_rows(union2113) == 0xC806CCE6B836FDFF and
            fnv_rows(residual20534) == 0x11285B5A49F4150D and
            top589 == [row for row in residual20534 if row[1] == 589] and
            fnv_rows(top589) == 0x5D9429C9F9971322,
            "post590 frontier changed")
    require(union2141 == sorted(set(union2113) | set(top589)) and
            residual20506 == sorted(set(residual20534) - set(top589)) and
            set(union2141).isdisjoint(residual20506) and
            set(union2141) | set(residual20506) == set(universe),
            "post589 typed partition changed")
    require(fnv_rows(union2141) == 0xC84BB7F7EAA0F230 and
            fnv_rows(residual20506) == 0x3CD0863A93C7602E and
            max(r for _, r in residual20506) == 588 and
            top588 == [row for row in residual20506 if row[1] == 588] and
            fnv_rows(top588) == 0x18CF9A572CF9A5BE,
            "post589 successor identities changed")

    add9: list[int] = []
    with (packet / "endpoint590_add9.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            add9.append(int(row["mask_hex"], 16))
    delete9 = [int(line, 16) for line in
               (packet / "endpoint590_delete9.txt").read_text(
                   encoding="ascii").splitlines()]
    require(len(add9) == len(set(add9)) == 9 and
            all(mask.bit_count() == 9 for mask in add9) and
            fnv_words(add9) == 0xD1CF49E4B811B958,
            "endpoint590 add9 identity changed")
    require(len(delete9) == len(set(delete9)) == 9 and
            all(mask.bit_count() == 8 for mask in delete9) and
            fnv_words(delete9) == 0x3546EB56552B4CDE,
            "endpoint590 delete9 identity changed")

    with (packet / "endpoint589_pair_O2.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        pair_rows = list(csv.DictReader(handle))
    require([(int(row["q"]), int(row["r"])) for row in pair_rows] == top589,
            "endpoint589 pair target changed")
    require(sum(int(row["exposed"]) for row in pair_rows) == 1451292 and
            sum(int(row["failures"]) for row in pair_rows) == 20036 and
            [(int(row["q"]), int(row["failures"])) for row in pair_rows
             if int(row["failures"])] == [(50, 20025), (96, 11)],
            "endpoint589 hostile-row distribution changed")
    pair_words: list[int] = []
    for row in pair_rows:
        pair_words.extend([
            int(row["q"]), int(row["r"]), int(row["active"]),
            int(row["active_fnv"], 16), int(row["active_joint"]),
            int(row["active_nonjoint"]), int(row["exposed"]),
            int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
            int(row["maximum_hits"]), int(row["failures"]),
            int(row["failure_fnv"], 16),
        ])
    require(fnv_words(pair_words) == 0xED09598FB9F312A9,
            "endpoint589 pair-ledger FNV changed")

    failures: dict[int, list[int]] = {50: [], 96: []}
    global_failure_words: list[int] = []
    with (packet / "endpoint589_failures_O2.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == ["q", "r", "body_hex"],
                "failure header changed")
        for row in rows:
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require(q in failures and r == 589 and body.bit_count() == 9,
                    "malformed endpoint589 failure")
            failures[q].append(body)
            global_failure_words.extend((q, r, body))
    require(len(failures[50]) == len(set(failures[50])) == 20025 and
            len(failures[96]) == len(set(failures[96])) == 11 and
            fnv_words(failures[50]) == 0xFF421454F02D9099 and
            fnv_words(failures[96]) == 0x6F70A2D28FF6A957 and
            fnv_words(global_failure_words) == 0x2228C33B4EF3B01D,
            "endpoint589 failure ledger changed")

    raw = (packet / "endpoint589_raw_O2.out").read_text(encoding="ascii")
    for needle in [
        "CARRIER 3925 FNV eeae5518d84ccac5 RANK8 3809 RANK9 116",
        "ROWS 28 ENDPOINT 589 ROW_FNV 5d9429c9f9971322",
        "BODY_TESTS 400600200",
        "SUMMARY EXPOSED 1451292 HIT_INCIDENCES 30755521 FAILURES 20036 FAILED_ROWS 2",
        "FAILURE_FNV 2228c33b4ef3b01d PAIR_LEDGER_FNV ed09598fb9f312a9",
        "VERDICT HOSTILE_FAIL",
    ]:
        require(needle in raw, f"carrier transcript missing: {needle}")

    direct: dict[int, list[tuple[int, int, int]]] = {50: [], 96: []}
    direct_ledgers = {50: OFFSET, 96: OFFSET}
    grids = {50: 2827379709554400, 96: 1130951883821760}
    with (packet / "direct_primary_O2.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == ["q", "r", "ordinal", "body_hex",
                                     "truncated_mass", "scaled_ticks"],
                "direct detail header changed")
        for row in rows:
            q, r, ordinal = int(row["q"]), int(row["r"]), int(row["ordinal"])
            body = int(row["body_hex"], 16)
            mass, ticks = int(row["truncated_mass"]), int(row["scaled_ticks"])
            require(q in direct and r == 589 and ordinal == len(direct[q]),
                    "direct detail order changed")
            require(body == failures[q][ordinal] and body.bit_count() == 9,
                    "direct detail body mismatch")
            require(mass > 0 and ticks == 63 * mass - 4 * grids[q] and
                    ticks > 0,
                    "direct truncated-mass inequality failed")
            direct[q].append((body, mass, ticks))
            direct_ledgers[q] = fnv_words_from(
                direct_ledgers[q], [q, r, ordinal, body, mass])
            direct_ledgers[q] = fnv_i128(direct_ledgers[q], ticks)
    require([body for body, _, _ in direct[50]] == failures[50] and
            [body for body, _, _ in direct[96]] == failures[96],
            "direct audit does not exhaust carrier failures")
    require(min((ticks, body) for body, _, ticks in direct[50]) ==
            (14566818763788984, 0x013C6401) and
            max((ticks, -body) for body, _, ticks in direct[50]) ==
            (26685137010259728, -0x21884443) and
            min((ticks, body) for body, _, ticks in direct[96]) ==
            (7172391058639758, 0x0D0C6401) and
            max((ticks, -body) for body, _, ticks in direct[96]) ==
            (8558758749944214, -0x35126400) and
            direct_ledgers == {50: 0x9A0CD88B508499A2,
                               96: 0x313A0CDBA0E5AC5C},
            "direct detail extrema/FNV changed")

    for name, prefix in [
        ("direct_primary_O2.out", "LRC14_ENDPOINT589_DIRECT_LITERAL_PRIMARY_V1"),
        ("direct_independent.out",
         "LRC14_ENDPOINT589_DIRECT_LITERAL_INDEPENDENT_V1"),
    ]:
        transcript = (packet / name).read_text(encoding="ascii")
        for needle in [
            prefix,
            "ROW 50,589 BODIES 20025 GRID 2827379709554400 CELLS 8389 PAIR_TICKS 2077256602813392 LOW_CLASSES 2383 CLASS_FNV 88d3eb2d7a477232 POSITIVE 20025 MIN_TICKS 14566818763788984 MIN_BODY 013c6401",
            "ROW 96,589 BODIES 11 GRID 1130951883821760 CELLS 8501 PAIR_TICKS 830901383902590 LOW_CLASSES 2352 CLASS_FNV cb74e72091a68363 POSITIVE 11 MIN_TICKS 7172391058639758 MIN_BODY 0d0c6401",
            "TOTAL_BODIES 20036 ALL_STRICTLY_POSITIVE 1",
            "TRUNCATED_LOW_RANK_CELL_MASS_LE_LITERAL_SAFE_MASS",
            "VERDICT PASS",
        ]:
            require(needle in transcript, f"direct transcript missing: {needle}")

    full_counts = {50: Counter(), 96: Counter()}
    with (packet / "literal_body_masses_independent_O2.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == [
            "q", "r", "body_hex", "truncated_mass_ticks",
            "truncated_surplus", "full_mass_ticks", "full_surplus",
            "omitted_high_rank_ticks"],
            "independent literal header changed")
        ordinals = {50: 0, 96: 0}
        for row in rows:
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require(q in ordinals and r == 589,
                    "independent literal row escaped hostile rows")
            ordinal = ordinals[q]
            require(ordinal < len(direct[q]) and body == direct[q][ordinal][0],
                    "independent literal body order changed")
            low = int(row["truncated_mass_ticks"])
            low_ticks = int(row["truncated_surplus"])
            full = int(row["full_mass_ticks"])
            full_ticks = int(row["full_surplus"])
            omitted = int(row["omitted_high_rank_ticks"])
            require((body, low, low_ticks) == direct[q][ordinal] and
                    full >= low and omitted == full - low and
                    full_ticks == 63 * full - 4 * grids[q] and
                    full_ticks - low_ticks == 63 * omitted,
                    "independent full/truncated identity failed")
            full_counts[q]["equal" if omitted == 0 else "strict"] += 1
            ordinals[q] += 1
        require(ordinals == {50: 20025, 96: 11},
                "independent literal quantifier changed")
    require(full_counts == {
        50: Counter({"equal": 16788, "strict": 3237}),
        96: Counter({"strict": 7, "equal": 4})},
        "full/truncated dominance census changed")
    literal_independent = (
        packet / "literal_lower_bound_independent_O2.out"
    ).read_text(encoding="ascii")
    for needle in [
        "PAIR_SAFE_CELLS 5798 PAIR_TICKS 2077256602813392 CLASSES_ALL 2434 CLASSES_RANK_LE_9 2383",
        "PAIR_SAFE_CELLS 5780 PAIR_TICKS 830901383902590 CLASSES_ALL 2423 CLASSES_RANK_LE_9 2352",
        "FULL_EQUALS_TRUNCATED 16788 FULL_STRICTLY_DOMINATES 3237",
        "FULL_EQUALS_TRUNCATED 4 FULL_STRICTLY_DOMINATES 7",
        "QUANTIFIER ALL_20025_Q50_AND_ALL_11_Q96_EXACT_CARRIER_FAILURE_BODIES",
        "VERDICT PASS",
    ]:
        require(needle in literal_independent,
                f"independent literal transcript missing: {needle}")

    active: list[int] = []
    with (packet / "q50_active_carrier_O2.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == ["mask_hex", "rank", "joint", "margin_ticks"],
                "q50 active header changed")
        for row in rows:
            mask, rank = int(row["mask_hex"], 16), int(row["rank"])
            require(rank in (8, 9) and mask.bit_count() == rank and
                    int(row["margin_ticks"]) > 0,
                    "q50 active row changed")
            active.append(mask)
    require(len(active) == len(set(active)) == 1398 and
            Counter(mask.bit_count() for mask in active) == {8: 1347, 9: 51} and
            fnv_words(active) == 0xC075113890C7F5E1,
            "q50 active carrier identity changed")
    structure = (packet / "q50_structure.out").read_text(encoding="ascii")
    for needle in [
        "HUBS bit14=95 degree=18671 bit24=193 degree=18801 BOTH=17454 ONLY95=1217 ONLY193=1347 NEITHER=7",
        "NEITHER_COMMON 00100400 {80,168} SUPPORT 26fc3707",
        "AUTOMORPHISM_CERTIFICATE 30_DISTINCT_VERTEX_DEGREES_IMPLY_TRIVIAL_LABEL_PERMUTATION_GROUP_AND_20025_SINGLETON_ORBITS",
        "GCD_HOSTILE gcd(95,589)=gcd(190,589)=19 BUT_DEGREES 18671 3465",
        "VERDICT PASS",
    ]:
        require(needle in structure, f"q50 structure transcript missing: {needle}")

    forbidden_assert = "as" + "sert "
    for source in [
        "derive_post590_typed_successor.py",
        "derive_post589_typed_successor.py",
        "endpoint589_direct_literal_primary.cpp",
        "endpoint589_direct_literal_independent.py",
        "endpoint589_literal_lower_bound_independent.cpp",
        "q50_active_carrier.cpp",
        "q50_failure_structure.py",
        "analyze_fixed50_petal_bridge.py",
    ]:
        require(forbidden_assert not in
                (packet / source).read_text(encoding="utf-8"),
                f"optimization-sensitive Python assertion in {source}")

    bridge = (packet / "petal_bridge.out").read_text(encoding="ascii")
    for needle in [
        "PETALS 0 BODIES 870",
        "PETALS 1 BODIES 4519",
        "PETALS 2 BODIES 7803",
        "PETALS 3 BODIES 5502",
        "PETALS 4 BODIES 1293",
        "PETALS 5 BODIES 38",
        "THM4234_INHERITED_LAYERS_0_TO_3 18694",
        "DIRECT_NEW_LAYERS_4_TO_5 1331 MIN_TRUNCATED_TICKS 15777555364138176 MIN_BODY 20744601",
        "VERDICT PASS",
    ]:
        require(needle in bridge, f"petal bridge transcript missing: {needle}")

    manifest_count = verify_manifest(packet)
    print("ENDPOINT589_DIRECT_LITERAL_CLOSURE_PACKET_VERIFY PASS")
    print(f"manifest_files={manifest_count}")
    print("endpoint589_rows=28 carrier_failures=20036 direct_positive=20036")
    print("q50_min_ticks=14566818763788984 q96_min_ticks=7172391058639758")
    print("typed_union=2141 residual=20506 next_endpoint=588 next_rows=66")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
