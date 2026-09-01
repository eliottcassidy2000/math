#!/usr/bin/env python3
"""Hardened verifier for the complete THM-4318 endpoint-590 packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import ast
from collections import Counter
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
FULL100 = (1 << 100) - 1


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


def same(left: Path, right: Path) -> None:
    require(left.read_bytes() == right.read_bytes(),
            f"byte mismatch: {left} {right}")


def verify_manifest(packet: Path) -> None:
    manifest_path = packet / "SHA256SUMS"
    require(manifest_path.is_file(), "packet manifest is missing")
    manifest: dict[str, str] = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        digest, relative = line.split("  ", 1)
        require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
                f"malformed manifest digest: {line}")
        require(relative not in manifest, f"duplicate manifest path: {relative}")
        manifest[relative] = digest
    actual = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path != manifest_path
    }
    require(set(manifest) == actual, "manifest file set changed")
    for relative, digest in manifest.items():
        require(sha256(packet / relative) == digest,
                f"manifest mismatch: {relative}")


def verify_source_manifest(computation: Path, packet: Path) -> None:
    source_manifest = packet / "SOURCE-SHA256SUMS"
    require(source_manifest.is_file(), "source manifest is missing")
    recorded: dict[str, str] = {}
    for line in source_manifest.read_text(encoding="ascii").splitlines():
        digest, name = line.split("  ", 1)
        require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
                f"malformed source-manifest digest: {line}")
        require("/" not in name and "\\" not in name and name not in recorded,
                f"malformed/duplicate source name: {name}")
        recorded[name] = digest
    actual = {path.name for path in computation.iterdir() if path.is_file()}
    require(set(recorded) == actual, "computation source file set changed")
    for name, digest in recorded.items():
        require(sha256(computation / name) == digest,
                f"computation source changed: {name}")


def read_mask_lines(path: Path) -> list[int]:
    masks = [int(line, 16) for line in path.read_text(encoding="ascii").splitlines()
             if line]
    require(len(masks) == len(set(masks)), f"duplicate mask in {path}")
    return masks


def read_pair_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical row ledger: {path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--packet", type=Path, required=True)
    parser.add_argument("--bootstrap", action="store_true")
    args = parser.parse_args()
    repo = args.repo.resolve()
    packet = args.packet.resolve()
    response = packet / "response"
    exchange = packet / "exchange"
    computation = repo / (
        "04-computation/"
        "lrc14_endpoint590_exact_nine_response_exchange_thm4318"
    )
    if not args.bootstrap:
        verify_manifest(packet)
        verify_source_manifest(computation, packet)
    for path in computation.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
                f"assert statement forbidden in {path.name}")

    for left, right in [
        ("endpoint590_raw_O3.out", "endpoint590_raw_O2.out"),
        ("endpoint590_pair_O3.csv", "endpoint590_pair_O2.csv"),
        ("endpoint590_failures_O3.csv", "endpoint590_failures_O2.csv"),
        ("endpoint590_response_structure_O3.out",
         "endpoint590_response_structure_O2.out"),
        ("endpoint590_response_signatures_O3.csv",
         "endpoint590_response_signatures_O2.csv"),
        ("endpoint590_greedy_cover_O3.csv",
         "endpoint590_greedy_cover_O2.csv"),
        ("endpoint590_cover_direct_O3.out",
         "endpoint590_cover_direct_O2.out"),
        ("endpoint590_cover_incidences_O3.csv",
         "endpoint590_cover_incidences_O2.csv"),
        ("endpoint590_exact_no8_O3.out", "endpoint590_exact_no8_O2.out"),
    ]:
        same(response / left, response / right)

    packet4314 = repo / (
        "05-knowledge/results/"
        "lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314"
    )
    top590_path = packet4314 / "typed/residual_top590.csv"
    require(sha256(top590_path) ==
            "eefb84e9fec0bcd809830a3c5ed18b0bb2aea1a2f2cb8b4994fafeaf1f5d01ce",
            "canonical top590 SHA changed")
    top590 = read_pair_rows(top590_path)
    require(len(top590) == 13 and
            fnv_words([value for row in top590 for value in row]) ==
            0x44AA8A793D162CF9,
            "canonical top590 identity changed")

    with (response / "endpoint590_pair_O3.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        pair_rows = list(csv.DictReader(handle))
    require([(int(row["q"]), int(row["r"])) for row in pair_rows] == top590,
            "baseline pair ledger targets changed")
    require(sum(int(row["exposed"]) for row in pair_rows) == 154023,
            "baseline exposed total changed")
    require(sum(int(row["failures"]) for row in pair_rows) == 100,
            "baseline failure total changed")
    require([(int(row["q"]), int(row["failures"])) for row in pair_rows
             if int(row["failures"])] == [(210, 100)],
            "baseline failure-row distribution changed")
    require({int(row["q"]) for row in pair_rows
             if int(row["minimum_hits"]) == 1} == {100, 105, 210},
            "baseline minimum-hit rows changed")
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
    require(fnv_words(pair_words) == 0x0C0EB3343C5A35BF,
            "baseline pair ledger FNV changed")

    failures: list[int] = []
    failure_words: list[int] = []
    with (response / "endpoint590_failures_O3.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
            require((q, r) == (210, 590) and body.bit_count() == 9,
                    "malformed endpoint590 failure")
            failures.append(body)
            failure_words.extend([q, r, body])
    require(len(failures) == 100 and len(set(failures)) == 100 and
            fnv_words(failure_words) == 0x8D19CBA1E86E53B5,
            "endpoint590 failure identity changed")

    raw = (response / "endpoint590_raw_O3.out").read_text(encoding="ascii")
    for needle in [
        "CARRIER 3925 FNV a0d08a38c10bdab7 RANK8 3818 RANK9 107",
        "ROWS 13 ENDPOINT 590 ROW_FNV 44aa8a793d162cf9",
        "BODY_TESTS 185992950",
        "SUMMARY EXPOSED 154023 HIT_INCIDENCES 4256322 FAILURES 100 FAILED_ROWS 1",
        "FAILURE_FNV 8d19cba1e86e53b5 PAIR_LEDGER_FNV c0eb3343c5a35bf",
        "NO_PHYSICAL_ENTRY_NO_LRC14",
        "VERDICT HOSTILE_FAIL",
    ]:
        require(needle in raw, f"baseline transcript missing: {needle}")

    atlas: dict[int, dict[str, str]] = {}
    count8 = count9 = 0
    with (response / "endpoint590_response_signatures_O3.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            reply = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
            require(0 < reply <= FULL100 and reply not in atlas,
                    "malformed/duplicate response signature")
            require(reply.bit_count() == int(row["weight"]),
                    "response signature weight changed")
            atlas[reply] = row
            count8 += int(row["count8"])
            count9 += int(row["count9"])
    require(len(atlas) == 14368 and count8 == 36285 and count9 == 568812,
            "complete response atlas census changed")

    response_text = (response / "endpoint590_response_structure_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "RANK8 ACTIVE 1073813 RESPONDERS 36285 RESPONDER_FNV 12fda2562de3e9a6",
        "RANK9 ACTIVE 6446114 RESPONDERS 568812 RESPONDER_FNV b531e8d2059a473c",
        "SIGNATURES 14368 MULTIPLICITY_RANGE 5570..32922 MULTIPLICITY_FNV a8e79ecbc16f3c1c",
        "GREEDY_COVER_SIZE 10",
        "VERDICT PASS",
    ]:
        require(needle in response_text, f"response transcript missing: {needle}")

    expected_masks = [
        0x20490236, 0x22045017, 0x29224016, 0x0A439108, 0x0B220096,
        0x0120403F, 0x12844116, 0x10686016, 0x084A6016,
    ]
    cover_rows: list[dict[str, str]] = []
    with (response / "endpoint590_solver_cover.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        cover_rows = list(csv.DictReader(handle))
    cover_masks = [int(row["mask_hex"], 16) for row in cover_rows]
    require(cover_masks == expected_masks and
            all(int(row["rank"]) == 9 and
                int(row["mask_hex"], 16).bit_count() == 9
                for row in cover_rows),
            "frozen nine-mask cover changed")
    require(fnv_words(cover_masks) == 0xD1CF49E4B811B958,
            "frozen nine-mask cover FNV changed")
    union = 0
    hit_counts = [0] * 100
    for row in cover_rows:
        reply = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
        require(reply.bit_count() == int(row["weight"]),
                "frozen cover declared weight changed")
        require(reply in atlas and
                int(atlas[reply]["least9"], 16) == int(row["mask_hex"], 16),
                "cover response is absent from complete atlas")
        union |= reply
        for index in range(100):
            hit_counts[index] += (reply >> index) & 1
    require(union == FULL100, "frozen nine-mask response union is incomplete")
    require(sum(hit_counts) == 136 and Counter(hit_counts) == {1: 69, 2: 26, 3: 5},
            "frozen nine-mask overlap distribution changed")

    direct = (response / "endpoint590_cover_direct_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "COVER 9 RANK9 9 COVER_FNV d1cf49e4b811b958",
        "TOTAL_INCIDENCES 136 HIT_RANGE 1..3 HIT_DISTRIBUTION 1:69 2:26 3:5",
        "HIT_LEDGER_FNV 57b3865558b0767b INCIDENCE_LEDGER_FNV 921483aaf38df968",
        "UNCOVERED 0",
        "VERDICT PASS",
    ]:
        require(needle in direct, f"direct audit transcript missing: {needle}")

    weights = [0] * 100
    with (response / "endpoint590_cover_dual_weights.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        for row in csv.DictReader(handle):
            index, weight = int(row["failure_ordinal"]), int(row["weight"])
            require(0 <= index < 100 and weight > 0 and weights[index] == 0,
                    "malformed integer-dual support")
            weights[index] = weight
    require(sum(weights) == 22 and max(weights) == 2,
            "integer-dual weight identity changed")
    loads = []
    for reply in atlas:
        loads.append(sum(weights[index] for index in range(100)
                         if reply >> index & 1))
    require(max(loads) == 3,
            "integer-dual load exceeds/fails frozen denominator three")
    require((sum(weights) + 2) // 3 == 8,
            "integer-dual certificate no longer proves lower bound eight")

    no8 = (response / "endpoint590_exact_no8_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "SIGNATURES 14368 MAXIMAL 1165 DEPTH 8",
        "RAW_FNV 98a18b04fa42e318 MAXIMAL_FNV 40cc79e6322deab0 DUAL_FNV f34be558e560d67",
        "NODES 7163197 DEAD 6775810 SUM_PRUNES 943159 DUAL_PRUNES 5575606 TRACE_FNV c2354a48f73aef96",
        "RESULT UNSAT",
        "EXACT_COMPLETE_SIGNATURE_QUOTIENT_DOMINANCE_AND_DEPTH_FIRST_SEARCH",
        "VERDICT PASS",
    ]:
        require(needle in no8, f"exact no-eight transcript missing: {needle}")

    solver = json.loads(
        (response / "endpoint590_solver_cover.json").read_text(encoding="ascii")
    )
    require(solver["success"] is True and solver["status"] == 0 and
            solver["cover_size"] == 9 and solver["mip_gap"] == 0.0,
            "solver-reported optimum metadata changed")

    expected_hashes = {
        "endpoint590_raw_O3.out": "a4e2b0eedc92623eaf55d51c48d750352d63d21eeab94e24bcfb3ec133586f08",
        "endpoint590_pair_O3.csv": "2bed885ff1e037ee5e7bb02b0b67d4349feee3814ec6a72376ac338feb2e3931",
        "endpoint590_failures_O3.csv": "ad925dc3036cd3a504928366dd9a181320288f40000f006f5d9d77d9edae9729",
        "endpoint590_response_structure_O3.out": "16d4e4ad757faaf7c377397aa0f301b6c8cbd2071c5243192b34cd96e631c821",
        "endpoint590_response_signatures_O3.csv": "11018f816361c93d6dbc2eaa81477ea12b64ea7239944758d70c713a50348f74",
        "endpoint590_greedy_cover_O3.csv": "c7bb1d7623564bc3fb7cfbb31a6c474e16ee49eace3e9599b168c55978048455",
        "endpoint590_solver_cover.csv": "3cde9686a8f9c90f72a07cad05d0d9007f754ba61a0ba1bc5fe5006c2a16a5cf",
        "endpoint590_cover_direct_O3.out": "1eea3df2897ef77efd02da90cde4c8dc682076a40cc47dd0acf582c8c735dafc",
        "endpoint590_cover_incidences_O3.csv": "84927470cb35242366ad9926b7621b1397280767da65875cb590a63618f01d7b",
        "endpoint590_cover_dual_weights.csv": "87286f36a2cf30ab43648b22fc4579f4d952885f53c07d213600dd967bb4b990",
        "endpoint590_cover_dual.out": "11f93c488e8642fe8c8bb5cdfdaf6cf144dd6a185c6e3dce4fcf048d995fdaac",
        "endpoint590_exact_no8_search.cpp": "1442a1525b0cae8d612e853e7a3a9a04098d1947360b68571e3ccc5d940352c3",
        "endpoint590_exact_no8_O3.out": "9487f8769027e51e5ffbca45b3c8b15024d06713e1fc698920b597566789e062",
    }
    for name, expected in expected_hashes.items():
        target = (computation / name if name.endswith((".cpp", ".py"))
                  else response / name)
        require(sha256(target) == expected,
                f"frozen artifact SHA changed: {name}")

    for left, right in [
        ("endpoint590_low_witness_O3.out", "endpoint590_low_witness_O2.out"),
        ("endpoint590_all_low_O3.csv", "endpoint590_all_low_O2.csv"),
        ("endpoint590_protected_low_O3.csv", "endpoint590_protected_low_O2.csv"),
        ("endpoint590_deletion_quotient_capacity_cover_O3.out",
         "endpoint590_deletion_quotient_capacity_cover_O2.out"),
        ("endpoint590_singleton493_capacity_cover_O3.csv",
         "endpoint590_singleton493_capacity_cover_O2.csv"),
        ("endpoint590_protected493_capacity_cover_O3.csv",
         "endpoint590_protected493_capacity_cover_O2.csv"),
        ("endpoint590_safe_old493_capacity_cover_O3.csv",
         "endpoint590_safe_old493_capacity_cover_O2.csv"),
        ("endpoint590_delete9_capacity_cover_candidate_O3.txt",
         "endpoint590_delete9_capacity_cover_candidate_O2.txt"),
        ("endpoint590_exchange493_capacity_cover_O3.out",
         "endpoint590_exchange493_capacity_cover_O2.out"),
        ("endpoint590_exchange493_capacity_cover_pair_O3.csv",
         "endpoint590_exchange493_capacity_cover_pair_O2.csv"),
        ("endpoint590_exchange493_capacity_cover_failures_O3.csv",
         "endpoint590_exchange493_capacity_cover_failures_O2.csv"),
    ]:
        same(exchange / left, exchange / right)

    low_text = (exchange / "endpoint590_low_witness_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "ROWS 13 ROW_FNV 44aa8a793d162cf9 BODY_TESTS 185992950",
        "ALL_SUMMARY ZERO 100 ONE 612 TWO 2018 LOW_FNV a5fc0bfbe9f4cc13",
        "PROTECTED_SUMMARY EXPOSED 154023 ZERO 100 ONE 524 TWO 1433 LOW_FNV 9efa84c46993339f",
        "VERDICT PASS",
    ]:
        require(needle in low_text, f"low-witness transcript missing: {needle}")

    def audit_low(path: Path, expected_count: int, expected_fnv: int) -> None:
        words: list[int] = []
        with path.open(newline="", encoding="ascii") as handle:
            rows = list(csv.DictReader(handle))
        require(len(rows) == expected_count, f"low-witness count changed: {path}")
        for row in rows:
            values = [
                int(row["q"]), int(row["r"]), int(row["body_hex"], 16),
                int(row["hits"]),
                int(row["first_mask_hex"], 16) if row["first_mask_hex"] else 0,
                int(row["second_mask_hex"], 16) if row["second_mask_hex"] else 0,
            ]
            require(values[1] == 590 and values[2].bit_count() == 9 and
                    0 <= values[3] <= 2,
                    f"malformed endpoint590 low witness: {row}")
            words.extend(values)
        require(fnv_words(words) == expected_fnv,
                f"low-witness FNV changed: {path}")

    audit_low(exchange / "endpoint590_all_low_O3.csv", 2730,
              0xA5FC0BFBE9F4CC13)
    audit_low(exchange / "endpoint590_protected_low_O3.csv", 2057,
              0x9EFA84C46993339F)

    quotient = (exchange /
        "endpoint590_deletion_quotient_capacity_cover_O3.out").read_text(
            encoding="ascii")
    for needle in [
        "COVER9_FNV d1cf49e4b811b958 AUGMENTED 3934 AUGMENTED_FNV 35ae2a57a188e4c6",
        "INHERITED_ROWS 480 ROW_FNV 8c9eeb1443dfa2d5 TARGET_ROWS 493 ROW_FNV 1fef91ec25d074e5",
        "ENDPOINT590_FAILURES_COVERED 100 ONE_WITNESS 612 RESOLVED 386 RETAINED 226 RETAINED_JOINT 32 RETAINED_NONJOINT 194",
        "PRIVATE_NONJOINT_OBLIGATIONS 1600 SINGLETON_FNV e862a639d9536826 PROTECTED_OLD_MASKS 425 PROTECTED_FNV 470279b13b453834",
        "OLD_NONJOINT 3504 SINGLETON_SAFE_OLD 3079 SAFE_FNV 89be292ce50c2831",
        "NO_SIMULTANEOUS_DELETION_INFERENCE_NO_PHYSICAL_ENTRY_NO_LRC14",
        "VERDICT PASS",
    ]:
        require(needle in quotient, f"deletion quotient missing: {needle}")

    with (exchange / "endpoint590_singleton493_capacity_cover_O3.csv").open(
        newline="", encoding="ascii") as handle:
        singleton_rows = list(csv.DictReader(handle))
    singleton_words: list[int] = []
    for row in singleton_rows:
        body = int(row["body_hex"], 16)
        witness = int(row["witness_hex"], 16)
        require(body.bit_count() == 9 and witness.bit_count() ==
                int(row["witness_rank"]), "malformed private obligation")
        singleton_words.extend([int(row["q"]), int(row["r"]), body, witness])
    require(len(singleton_rows) == 1600 and
            fnv_words(singleton_words) == 0xE862A639D9536826,
            "private singleton ledger changed")

    with (exchange / "endpoint590_protected493_capacity_cover_O3.csv").open(
        newline="", encoding="ascii") as handle:
        protected_rows = list(csv.DictReader(handle))
    protected_masks = [int(row["mask_hex"], 16) for row in protected_rows]
    require(len(protected_masks) == 425 and len(set(protected_masks)) == 425 and
            fnv_words(protected_masks) == 0x470279B13B453834,
            "protected-mask ledger changed")

    with (exchange / "endpoint590_safe_old493_capacity_cover_O3.csv").open(
        newline="", encoding="ascii") as handle:
        safe_rows = list(csv.DictReader(handle))
    safe_masks = [int(row["mask_hex"], 16) for row in safe_rows]
    deletion_masks = read_mask_lines(
        exchange / "endpoint590_delete9_capacity_cover_candidate_O3.txt")
    expected_deletion = [
        0x06021829, 0x23222801, 0x12444083, 0x20827018, 0x29C04082,
        0x02916180, 0x13C00881, 0x070C4840, 0x2380408A,
    ]
    require(len(safe_masks) == 3079 and len(set(safe_masks)) == 3079 and
            fnv_words(safe_masks) == 0x89BE292CE50C2831,
            "singleton-safe old-mask ledger changed")
    require(deletion_masks == expected_deletion == safe_masks[:9] and
            fnv_words(deletion_masks) == 0x3546EB56552B4CDE,
            "selected nine deletions changed")

    exchange_text = (exchange /
        "endpoint590_exchange493_capacity_cover_O3.out").read_text(
            encoding="ascii")
    for needle in [
        "ADDITIONS 9 FNV d1cf49e4b811b958 DELETIONS 9 FNV 3546eb56552b4cde",
        "EXCHANGED_CARRIER 3925 FNV eeae5518d84ccac5 RANK8 3809 RANK9 116 JOINT_RETAINED 421",
        "ROWS 493 ROW_FNV 1fef91ec25d074e5 BODY_TESTS 7053424950",
        "SUMMARY EXPOSED 2401583 HIT_INCIDENCES 79803901 FAILURES 0 FAILED_ROWS 0 FAILURE_FNV cbf29ce484222325 PAIR_LEDGER_FNV 1092fd57a8581a34",
        "DIRECT_SIMULTANEOUS_FIXED_NINE_ADD_NINE_DELETE",
        "VERDICT PASS",
    ]:
        require(needle in exchange_text, f"exchange replay missing: {needle}")

    with (exchange /
          "endpoint590_exchange493_capacity_cover_pair_O3.csv").open(
              newline="", encoding="ascii") as handle:
        exchange_rows = list(csv.DictReader(handle))
    exchange_keys = [(int(row["q"]), int(row["r"])) for row in exchange_rows]
    require(exchange_keys == sorted(set(exchange_keys)) and len(exchange_keys) == 493,
            "exchange target row ledger changed")
    source4313 = repo / (
        "05-knowledge/results/"
        "lrc14_endpoint592_fortythree_response_exchange_thm4313"
    )
    with (source4313 / "exchange43_final_pair_O3.csv").open(
        newline="", encoding="ascii"
    ) as handle:
        inherited467 = [(int(row["q"]), int(row["r"]))
                        for row in csv.DictReader(handle)]
    top591 = read_pair_rows(source4313 / "typed/residual_top591.csv")
    expected_target = sorted(set(inherited467) | set(top591) | set(top590))
    require(exchange_keys == expected_target,
            "exchange target is not inherited467 plus endpoint591/590")
    exchange_words: list[int] = []
    for row in exchange_rows:
        require(int(row["failures"]) == 0 and int(row["minimum_hits"]) >= 1,
                f"exchange row is not closed: {row}")
        exchange_words.extend([
            int(row["q"]), int(row["r"]), int(row["active"]),
            int(row["active_fnv"], 16), int(row["active_joint"]),
            int(row["active_nonjoint"]), int(row["exposed"]),
            int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
            int(row["maximum_hits"]), int(row["failures"]),
            int(row["failure_fnv"], 16),
        ])
    require(sum(int(row["exposed"]) for row in exchange_rows) == 2401583 and
            fnv_words(exchange_words) == 0x1092FD57A8581A34,
            "exchange pair ledger summary changed")
    with (exchange /
          "endpoint590_exchange493_capacity_cover_failures_O3.csv").open(
              newline="", encoding="ascii") as handle:
        require(not list(csv.DictReader(handle)), "exchange failures reappeared")

    same(packet / "typed_endpoint590_consumer.out",
         packet / "typed_endpoint590_consumer_opt.out")
    for name in ["typed_union2113.csv", "final_residual20534.csv",
                 "residual_top589.csv"]:
        same(packet / "typed" / name, packet / "typed_opt" / name)
    prior_union = read_pair_rows(packet4314 / "typed/typed_union2100.csv")
    prior_residual = read_pair_rows(packet4314 / "typed/final_residual20547.csv")
    typed_union = read_pair_rows(packet / "typed/typed_union2113.csv")
    typed_residual = read_pair_rows(packet / "typed/final_residual20534.csv")
    top589 = read_pair_rows(packet / "typed/residual_top589.csv")
    require(typed_union == sorted(set(prior_union) | set(top590)) and
            typed_residual == sorted(set(prior_residual) - set(top590)),
            "typed successor is not exact endpoint590 consumption")
    require(len(typed_union) == 2113 and len(typed_residual) == 20534 and
            len(top589) == 28 and
            top589 == [row for row in typed_residual if row[1] == 589],
            "typed successor counts/top layer changed")
    require(fnv_words([value for row in typed_union for value in row]) ==
            0xC806CCE6B836FDFF and
            fnv_words([value for row in typed_residual for value in row]) ==
            0x11285B5A49F4150D and
            fnv_words([value for row in top589 for value in row]) ==
            0x5D9429C9F9971322,
            "typed successor FNV changed")
    require(sha256(packet / "typed/typed_union2113.csv") ==
            "0006156b27e794783a9ddc65a932ef005064b6bbd168193ddb5d092511881aa3" and
            sha256(packet / "typed/final_residual20534.csv") ==
            "664ffe4f24d281ccb6cbd0de250f50b008929f836696319e78fa345b591799bb" and
            sha256(packet / "typed/residual_top589.csv") ==
            "005306a28c756a93862fd1745414fc058032f4d99e32086c079c29607bfed0c6",
            "typed successor SHA changed")
    typed_text = (packet / "typed_endpoint590_consumer.out").read_text(
        encoding="ascii")
    for needle in [
        "NEW_UNION 2113 FNV c806cce6b836fdff",
        "NEW_RESIDUAL 20534 FNV 11285b5a49f4150d",
        "NEW_TOP 589 ROWS 28 FNV 5d9429c9f9971322",
        "NOT_PHYSICAL_ENTRY_NO_LRC14",
        "VERDICT PASS",
    ]:
        require(needle in typed_text, f"typed transcript missing: {needle}")

    print("ENDPOINT590_PACKET_VERIFY PASS")
    print("baseline_rows=13 body_tests=185992950 failures=100 failed_row=210,590")
    print("cover_size=9 rank9=9 uncovered=0 cover_fnv=d1cf49e4b811b958")
    print("proved_full_universe_lower=9 exact_cover_upper=9 optimum=9")
    print("exchange=9_for_9 target_rows=493 body_tests=7053424950 failures=0")
    print("typed_union=2113 residual=20534 next_endpoint=589 next_rows=28")
    print("solver_reported_optimum=9 role=CORROBORATIVE_NOT_PROOF_DEPENDENCY")


if __name__ == "__main__":
    main()
