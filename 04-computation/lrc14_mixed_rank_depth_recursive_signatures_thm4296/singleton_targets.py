#!/usr/bin/env python3
"""Emit live singleton signature groups, largest first.

All inputs are explicit so the script can be replayed against a frozen result
packet rather than resolving an older theorem packet implicitly.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", required=True, type=Path,
                        help="sorted q,r current-residual CSV")
    parser.add_argument("--signatures", required=True, type=Path,
                        help="full inactive-signature CSV")
    parser.add_argument("--output", type=Path,
                        help="write the target ledger here instead of stdout")
    args = parser.parse_args()

    with args.live.open(newline="", encoding="ascii") as handle:
        live_rows = [tuple(map(int, row)) for row in csv.reader(handle) if row]
    require(live_rows == sorted(live_rows), "live residual is not sorted")
    require(len(live_rows) == len(set(live_rows)),
            "live residual contains duplicate pairs")
    live = set(live_rows)

    groups: dict[int, list[tuple[int, int]]] = defaultdict(list)
    with args.signatures.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            pair = (int(row["q"]), int(row["r"]))
            if pair not in live or int(row["inactive_count"]) != 1:
                continue
            words = [int(row[f"w{i}"], 16) for i in range(7)]
            require(sum(word.bit_count() for word in words) == 1,
                    f"singleton row has non-singleton word at {pair}")
            index = next(64 * i + (word & -word).bit_length() - 1
                         for i, word in enumerate(words) if word)
            require(index < 421, f"signature index escapes deck at {pair}")
            groups[index].append(pair)

    lines = [
        f"{index} {len(rows)} {rows[-1][0]},{rows[-1][1]}"
        for index, rows in sorted(groups.items(),
                                  key=lambda item: (-len(item[1]), item[0]))
    ]
    # Keep the maintained packet byte-stable across host operating systems.
    payload = ("\n".join(lines) + "\n").encode("ascii")
    if args.output is None:
        import sys
        sys.stdout.buffer.write(payload)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(payload)


if __name__ == "__main__":
    main()
