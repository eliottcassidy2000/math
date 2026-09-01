#!/usr/bin/env python3
"""Exact petal-layer split of the q=50 endpoint-589 carrier failures."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
PETALS = {10, 15, 20, 30, 60, 63, 132, 170, 176, 190, 264, 290}
PETAL_MASK = sum(1 << bit for bit, speed in enumerate(POOL)
                 if speed in PETALS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("detail", type=Path)
    args = parser.parse_args()
    layers: dict[int, list[tuple[int, int, int, int]]] = {}
    with args.detail.open(newline="", encoding="ascii") as handle:
        rows = csv.DictReader(handle)
        require(rows.fieldnames == ["q", "r", "ordinal", "body_hex",
                                     "truncated_mass", "scaled_ticks"],
                "detail header changed")
        for row in rows:
            if int(row["q"]) != 50:
                continue
            body = int(row["body_hex"], 16)
            k = (body & PETAL_MASK).bit_count()
            layers.setdefault(k, []).append((
                int(row["scaled_ticks"]), body,
                int(row["truncated_mass"]), int(row["ordinal"])))
    expected = {0: 870, 1: 4519, 2: 7803, 3: 5502, 4: 1293, 5: 38}
    require(Counter({key: len(value) for key, value in layers.items()}) ==
            Counter(expected), "petal-layer census changed")
    require(sum(len(value) for value in layers.values()) == 20025,
            "q50 failure count changed")
    print("LRC14_ENDPOINT589_FIXED50_PETAL_BRIDGE_V1")
    for k in sorted(layers):
        minimum = min(layers[k], key=lambda item: (item[0], item[1]))
        print(f"PETALS {k} BODIES {len(layers[k])} MIN_TRUNCATED_TICKS "
              f"{minimum[0]} MIN_BODY {minimum[1]:08x}")
    inherited = sum(len(layers[k]) for k in range(4))
    new = sum(len(value) for k, value in layers.items() if k >= 4)
    new_minimum = min(
        (ticks, body) for k, value in layers.items() if k >= 4
        for ticks, body, _, _ in value)
    require(inherited == 18694 and new == 1331 and
            new_minimum == (15777555364138176, 0x20744601),
            "fixed50 bridge boundary changed")
    print(f"THM4234_INHERITED_LAYERS_0_TO_3 {inherited}")
    print(f"DIRECT_NEW_LAYERS_4_TO_5 {new} MIN_TRUNCATED_TICKS "
          f"{new_minimum[0]} MIN_BODY {new_minimum[1]:08x}")
    print("MAP FAILURE_BODY_TO_PETAL_COUNT PRESERVES_LABELLED_BODY "
          "DESTROYS_MARGIN_AND_CARRIER_MULTIPLICITY")
    print("SCOPE FINITE_EXACT_CLASSIFICATION_PLUS_CITED_THM4234_INHERITANCE_"
          "NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
