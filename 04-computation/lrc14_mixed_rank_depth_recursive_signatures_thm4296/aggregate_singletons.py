#!/usr/bin/env python3
"""Audit and summarize the exhaustive post-THM4287 singleton surgery census."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--targets", required=True, type=Path,
                    help="target ledger emitted by singleton_targets.py")
parser.add_argument("--runs-dir", required=True, type=Path,
                    help="directory containing one singleton run per target")
parser.add_argument("--live", required=True, type=Path,
                    help="sorted q,r current-residual CSV")
parser.add_argument("--output-dir", required=True, type=Path,
                    help="directory for the three aggregate CSV ledgers")
parser.add_argument(
    "--summary-output",
    type=Path,
    help="write the frozen ASCII/LF aggregate report here instead of stdout",
)
parser.add_argument(
    "--run-template",
    default="singleton_v2_{index}.out",
    help="filename template within --runs-dir (default: %(default)s)",
)
args = parser.parse_args()


def fnv_add(state: int, value: int) -> int:
    for byte in value.to_bytes(8, "little", signed=False):
        state ^= byte
        state = state * 0x100000001B3 & ((1 << 64) - 1)
    return state


def pair_fnv(rows: list[tuple[int, int]]) -> int:
    state = 0xCBF29CE484222325
    for q, r in rows:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return state


targets = []
for line in args.targets.read_text(encoding="ascii").splitlines():
    index, count, target = line.split()
    targets.append((int(index), int(count), target))
assert len(targets) == 192

with args.live.open(newline="", encoding="ascii") as h:
    live = [tuple(map(int, row)) for row in csv.reader(h)]
assert len(live) == 22647 and live == sorted(live)
live_set = set(live)

successes = []
union: set[tuple[int, int]] = set()
for index, expected_count, target in targets:
    path = args.runs_dir / args.run_template.format(index=index)
    text = path.read_text(encoding="ascii")
    ideal_match = re.search(r"^IDEAL_ROWS (\d+) FNV ([0-9a-f]+) ROWS (.*)$", text, re.M)
    old_match = re.search(r"^LIVE_ROWS \d+ TARGET \S+ J \d+ OLD_MASK ([0-9a-f]+)$", text, re.M)
    obligation_match = re.search(r"^OBLIGATIONS (\d+) FNV ([0-9a-f]+)", text, re.M)
    responder_match = re.search(r"^FULL_BODY_RESPONDERS (\d+) FNV ([0-9a-f]+).* COMMON_ACTIVE (\d+)", text, re.M)
    exact_match = re.search(r"^EXACT_MINIMUM (\S+).* NEW_LIVE_ROWS (\d+)$", text, re.M)
    assert ideal_match and old_match and obligation_match and responder_match and exact_match, path
    rows = [tuple(map(int, token.split(","))) for token in ideal_match.group(3).split()]
    assert len(rows) == expected_count == int(ideal_match.group(1))
    assert rows == sorted(rows) and set(rows) <= live_set
    assert pair_fnv(rows) == int(ideal_match.group(2), 16)
    minimum = exact_match.group(1)
    new_rows = int(exact_match.group(2))
    if minimum != "1":
        assert new_rows == 0 and int(responder_match.group(3)) == 0
        continue
    witness_match = re.search(
        r"^WITNESS ([0-9a-f]+) WEAKEST_PAIR (\S+) NORMALIZED_GAP (\S+) MARGIN_FNV ([0-9a-f]+)$",
        text,
        re.M,
    )
    rebuilt_match = re.search(
        r"^REBUILT_DECK 421 FNV ([0-9a-f]+) BODY_SCAN 14307150 CHECKS (\d+) FAILURES 0 BODY_FNV ([0-9a-f]+)$",
        text,
        re.M,
    )
    assert witness_match and rebuilt_match and new_rows == expected_count
    assert not (union & set(rows))
    union.update(rows)
    successes.append(
        {
            "index": index,
            "rows": rows,
            "ideal_fnv": ideal_match.group(2),
            "old_mask": old_match.group(1),
            "obligations": int(obligation_match.group(1)),
            "obligation_fnv": obligation_match.group(2),
            "responders": int(responder_match.group(1)),
            "responder_fnv": responder_match.group(2),
            "common": int(responder_match.group(3)),
            "witness": witness_match.group(1),
            "weakest": witness_match.group(2),
            "gap": witness_match.group(3),
            "margin_fnv": witness_match.group(4),
            "rebuilt_fnv": rebuilt_match.group(1),
            "body_checks": int(rebuilt_match.group(2)),
            "body_fnv": rebuilt_match.group(3),
        }
    )

ordered_union = sorted(union)
remaining = [pair for pair in live if pair not in union]
union_bytes = "".join(f"{q},{r}\n" for q, r in ordered_union).encode()
remaining_bytes = "".join(f"{q},{r}\n" for q, r in remaining).encode()

witness_lines = []
for item in sorted(successes, key=lambda x: x["index"]):
    witness_lines.append(
        f'{item["index"]},{item["old_mask"]},{item["witness"]},{len(item["rows"])},'
        f'{item["ideal_fnv"]},{item["obligations"]},{item["obligation_fnv"]},'
        f'{item["responders"]},{item["responder_fnv"]},{item["common"]},'
        f'{item["weakest"]},{item["gap"]},{item["margin_fnv"]},'
        f'{item["rebuilt_fnv"]},{item["body_checks"]},{item["body_fnv"]}'
    )
witness_bytes = ("\n".join(witness_lines) + "\n").encode()

summary = [
    "SINGLETON_IDEAL_AGGREGATE_V1",
    f"UNIVERSE LIVE {len(live)} SINGLETON_GROUPS {len(targets)} SINGLETON_ROWS {sum(c for _,c,_ in targets)}",
    f"SUCCESS_GROUPS {len(successes)} SUCCESS_ROWS {len(ordered_union)} DISJOINT YES",
    f"SUCCESS_UNION FNV {pair_fnv(ordered_union):016x} SHA256 {hashlib.sha256(union_bytes).hexdigest()} "
    f"MAX_ENDPOINT {max(r for _,r in ordered_union)}",
    f"REMAINING COUNT {len(remaining)} FNV {pair_fnv(remaining):016x} "
    f"SHA256 {hashlib.sha256(remaining_bytes).hexdigest()} MAX_ENDPOINT {max(r for _,r in remaining)}",
]
top_endpoint = max(r for _, r in remaining)
summary.append("REMAINING_TOP " + " ".join(
    f"{q},{r}" for q, r in remaining if r == top_endpoint
))
summary.append(
    f"WITNESS_LEDGER ROWS {len(witness_lines)} SHA256 {hashlib.sha256(witness_bytes).hexdigest()}"
)
for item in sorted(successes, key=lambda x: (-len(x["rows"]), x["index"])):
    summary.append(
        f'NODE J {item["index"]} ROWS {len(item["rows"])} IDEAL_FNV {item["ideal_fnv"]} '
        f'OLD {item["old_mask"]} OBLIGATIONS {item["obligations"]} OBLIGATION_FNV {item["obligation_fnv"]} '
        f'RESPONDERS {item["responders"]} COMMON {item["common"]} WITNESS {item["witness"]} '
        f'WEAKEST {item["weakest"]} GAP {item["gap"]} MARGIN_FNV {item["margin_fnv"]} '
        f'REBUILT_FNV {item["rebuilt_fnv"]}'
    )
summary.append("VERDICT PASS EXHAUSTIVE_192_SINGLETON_GROUP_ONE_REPLACEMENT_CENSUS")
summary_bytes = ("\n".join(summary) + "\n").encode("ascii")
if args.summary_output is None:
    print("\n".join(summary))
else:
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_bytes(summary_bytes)

args.output_dir.mkdir(parents=True, exist_ok=True)
(args.output_dir / "singleton_success_union.csv").write_bytes(union_bytes)
(args.output_dir / "singleton_success_remaining.csv").write_bytes(remaining_bytes)
(args.output_dir / "singleton_success_witnesses.csv").write_bytes(witness_bytes)
