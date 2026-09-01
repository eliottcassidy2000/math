#!/usr/bin/env python3
"""Independent consumer for the exact k=350 retention certificate."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MOD = (1 << 64) - 1


def need(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def fnv_words(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        need(0 <= word < 1 << 64, "FNV word outside u64")
        for byte in range(8):
            state ^= (word >> (8 * byte)) & 0xFF
            state = (state * PRIME) & MOD
    return state


def read_masks(path: Path) -> list[int]:
    tokens = path.read_text(encoding="ascii").split()
    need(all(re.fullmatch(r"[0-9a-f]{8}", token) for token in tokens),
         f"malformed mask ledger {path}")
    masks = [int(token, 16) for token in tokens]
    need(len(masks) == len(set(masks)), f"duplicate mask in {path}")
    need(all(mask.bit_count() in (8, 9) for mask in masks),
         f"rank escaped 8/9 in {path}")
    return masks


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    need(len(sys.argv) == 21,
         "usage: verify_k350_retention.py CERT_O3 CERT_O2 ATLAS_O3 ATLAS_O2 "
         "EDGES_O3 EDGES_O2 FULL_O3 FULL_O2 PAIRS_O3 PAIRS_O2 FAIL_O3 "
         "FAIL_O2 DELETE5141 RETAIN37 REPAIRED5104 PACKING37 FAILURES84 "
         "BASE_OUT BASE_PAIRS BASE_FAILURES")
    paths = [Path(word) for word in sys.argv[1:]]
    (cert_o3, cert_o2, atlas_o3, atlas_o2, edges_o3, edges_o2,
     full_o3, full_o2, pairs_o3, pairs_o2, fail_o3, fail_o2,
     delete_path, retain_path, repaired_path, packing_path,
     failures84_path, base_out, base_pairs, base_failures) = paths

    for left, right, label in (
        (cert_o3, cert_o2, "certificate transcript"),
        (atlas_o3, atlas_o2, "response atlas"),
        (edges_o3, edges_o2, "edge ledger"),
        (full_o3, full_o2, "full replay transcript"),
        (pairs_o3, pairs_o2, "full pair ledger"),
        (fail_o3, fail_o2, "final failure ledger"),
    ):
        need(left.read_bytes() == right.read_bytes(),
             f"O2/O3 differ for {label}")

    deleted = read_masks(delete_path)
    retained = read_masks(retain_path)
    repaired = read_masks(repaired_path)
    need((len(deleted), len(retained), len(repaired)) == (5141, 37, 5104),
         "deletion/retention counts changed")
    need(fnv_words(deleted) == 0x03921CF597EE9863,
         "D350 FNV changed")
    need(fnv_words(repaired) == 0xFF4C932F9A7ADAC8,
         "repaired-deletion FNV changed")
    retained_set = set(retained)
    need(retained_set <= set(deleted), "retention escaped D350")
    need(repaired == [mask for mask in deleted if mask not in retained_set],
         "repaired deletion is not ordered D350 minus H")

    packing = [int(word) for word in packing_path.read_text(
        encoding="ascii").split()]
    need(len(packing) == len(set(packing)) == 37 and
         all(0 <= index < 84 for index in packing),
         "packing37 ledger changed")
    packing_bits = sum(1 << index for index in packing)
    packing_fnv = fnv_words(packing)
    need(packing_fnv == 0xF62E44818F143DFE, "packing FNV changed")

    failures84 = read_csv(failures84_path)
    need(len(failures84) == 84, "failure obligation count changed")
    obligation_words: list[int] = []
    obligation_keys: list[tuple[int, int, int]] = []
    for row in failures84:
        q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
        need(0 < q < r and body.bit_count() == 9, "invalid failure obligation")
        obligation_keys.append((q, r, body))
        obligation_words.extend((q, r, body))
    need(len(obligation_keys) == len(set(obligation_keys)),
         "duplicate failure obligation")
    need(fnv_words(obligation_words) == 0x3A4207B0CB910C10,
         "failure obligation FNV changed")
    need(base_failures.read_bytes() == failures84_path.read_bytes(),
         "independent raw D350 replay did not regenerate the 84 failures")

    atlas = read_csv(atlas_o3)
    need(len(atlas) == 53, "response type count changed")
    responses: list[int] = []
    counts: list[int] = []
    least_to_response: dict[int, int] = {}
    for row in atlas:
        response = int(row["w0"], 16) | (int(row["w1"], 16) << 64)
        count = int(row["count"])
        least = int(row["least_mask"], 16)
        need(0 < response < 1 << 84 and count > 0,
             "invalid response type")
        need(least not in least_to_response, "duplicate least responder")
        least_to_response[least] = response
        responses.append(response)
        counts.append(count)
    need(len(responses) == len(set(responses)) == 53,
         "duplicate response type")
    need(sum(counts) == 70, "responder count changed")
    need(all((response & packing_bits).bit_count() <= 1
             for response in responses),
         "packing lower certificate failed")
    need(all(mask in least_to_response for mask in retained),
         "retention mask is not frozen least type representative")
    cover = 0
    for mask in retained:
        cover |= least_to_response[mask]
    need(cover == (1 << 84) - 1, "retention cover missed an obligation")

    edges = read_csv(edges_o3)
    need(len(edges) == 84, "residual edge count changed")
    edge_words: list[int] = []
    for index, row in enumerate(edges):
        got_index = int(row["index"])
        q, r, body = int(row["q"]), int(row["r"]), int(row["body_hex"], 16)
        size, witness_fnv = int(row["witness_count"]), int(
            row["witness_fnv"], 16)
        need(got_index == index and (q, r, body) == obligation_keys[index],
             "edge/obligation alignment changed")
        need(size > 0, "empty residual edge")
        edge_words.extend((q, r, body, size, witness_fnv))
    need(fnv_words(edge_words) == 0x2FAAB17FAD74C14A,
         "residual edge FNV changed")

    cert = cert_o3.read_text(encoding="ascii")
    exact_lines = (
        "SUPPORT_SIGN_CELLS 3526429 SIGN_FNV 853edf378527886 "
        "EQUALITIES 0 SUPPORT_LEDGER_FNV c93aae41a6892f37",
        "DELETE350 5141 FNV 3921cf597ee9863 EXACT_NONJOINT_THRESHOLD 1",
        "OBLIGATIONS 84 FNV 3a4207b0cb910c10 EDGE_FNV 2faab17fad74c14a",
        "RESPONDERS 70 FNV b18cf54b71b4a093 TYPES 53 "
        "TYPE_FNV 578ebd478ecee5a",
        "PACKING 37 FNV f62e44818f143dfe MAX_RESPONSE_HITS 1",
        "MIN_RETENTION 37 LOWER_PACKING 37 UPPER_COVER 37",
        "REPAIRED_DELETE 5104 FNV ff4c932f9a7adac8 FINAL_CARRIER 3925 "
        "FNV 6fbd0bffcf0ed78b RANK8 3858 RANK9 67 JOINT_RETAINED 421",
        "VERDICT PASS",
    )
    need(all(line in cert for line in exact_lines),
         "certificate transcript identity changed")

    pairs = read_csv(pairs_o3)
    need(len(pairs) == 391, "final pair row count changed")
    pair_words: list[int] = []
    exposed_sum = 0
    for row in pairs:
        q, r = int(row["q"]), int(row["r"])
        active, active_fnv = int(row["active"]), int(row["active_fnv"], 16)
        active_joint = int(row["active_joint"])
        active_nonjoint = int(row["active_nonjoint"])
        exposed, exposed_fnv = int(row["exposed"]), int(
            row["exposed_fnv"], 16)
        minimum, maximum = int(row["minimum_hits"]), int(row["maximum_hits"])
        failures, failure_fnv = int(row["failures"]), int(
            row["failure_fnv"], 16)
        need(0 < q < r and active == active_joint + active_nonjoint,
             "invalid final pair row")
        need(minimum >= 1 and maximum >= minimum and failures == 0 and
             failure_fnv == OFFSET, "final row failed closure")
        exposed_sum += exposed
        pair_words.extend((q, r, active, active_fnv, active_joint,
                           active_nonjoint, exposed, exposed_fnv, minimum,
                           maximum, failures, failure_fnv))
    need(exposed_sum == 1_418_344, "final exposed-body sum changed")
    need(fnv_words(pair_words) == 0xBB28F7E567C4A4B0,
         "final pair ledger FNV changed")
    need(fail_o3.read_text(encoding="ascii") == "q,r,body_hex\n",
         "final failure ledger is nonempty")

    raw_base_pairs = read_csv(base_pairs)
    need(len(raw_base_pairs) == 391 and
         sum(int(row["failures"]) for row in raw_base_pairs) == 84,
         "raw D350 pair audit failure census changed")
    raw_pair_words: list[int] = []
    for row in raw_base_pairs:
        raw_pair_words.extend((
            int(row["q"]), int(row["r"]), int(row["active"]),
            int(row["active_fnv"], 16), int(row["active_joint"]),
            int(row["active_nonjoint"]), int(row["exposed"]),
            int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
            int(row["maximum_hits"]), int(row["failures"]),
            int(row["failure_fnv"], 16)))
    need(fnv_words(raw_pair_words) == 0x8B4B347CCEA742D9,
         "raw D350 pair ledger FNV changed")
    base_text = base_out.read_text(encoding="ascii")
    need("SUMMARY TOTAL_EXPOSED 1418344 TOTAL_HIT_INCIDENCES 46919727 "
         "FAILURES 84 PAIR_LEDGER_FNV 8b4b347ccea742d9" in base_text and
         "VERDICT FAIL" in base_text,
         "raw D350 full replay identity changed")

    full = full_o3.read_text(encoding="ascii")
    full_lines = (
        "FINAL_CARRIER 3925 FNV 6fbd0bffcf0ed78b RANK8 3858 RANK9 67",
        "ROWS 391 ROW_FNV c732a1532c12b9f6 WORKERS 32 "
        "BODY_UNIVERSE_PER_ROW 14307150 TOTAL_BODY_TESTS 5594095650",
        "SUMMARY TOTAL_EXPOSED 1418344 TOTAL_HIT_INCIDENCES 47375188 "
        "FAILURES 0 PAIR_LEDGER_FNV bb28f7e567c4a4b0",
        "VERDICT PASS",
    )
    need(all(line in full for line in full_lines),
         "full replay identity changed")

    print("LRC14_K350_RETENTION_CONSUMER_V1")
    print("THRESHOLD NONJOINT_SUPPORT<=350 DELETE=5141 EXACT")
    print("HYPERGRAPH OBLIGATIONS=84 RESPONDERS=70 TYPES=53")
    print("MIN_RETENTION=37 PACKING=37 COVER=37")
    print("FINAL CARRIER=3925 RANKS=3858/67 FNV=6fbd0bffcf0ed78b")
    print("FULL391 BODY_TESTS=5594095650 FAILURES=0 PAIR_FNV=bb28f7e567c4a4b0")


if __name__ == "__main__":
    main()
