"""Round an LP dual down and verify it exactly over every exported column."""

from __future__ import annotations

import csv
import argparse
import json
from decimal import Decimal, ROUND_FLOOR
from pathlib import Path

SCALE = 1_000_000_000
COUNTS = (13, 3, 48, 1, 288, 2101, 14)
WORDS = (1, 1, 1, 1, 5, 33, 1)
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state = ((state ^ ((word >> shift) & 0xFF)) * FNV_PRIME) & MASK64
    return state

parser = argparse.ArgumentParser()
parser.add_argument("--pool", type=Path, required=True)
parser.add_argument("--dual", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--weights-output", type=Path, required=True)
parser.add_argument("--masks-output", type=Path, required=True)
args = parser.parse_args()

tokens = args.dual.read_text().split()
assert len(tokens) == sum(COUNTS)
weights = [int((Decimal(token) * SCALE).to_integral_value(rounding=ROUND_FLOOR)) for token in tokens]
assert min(weights) >= 0
args.weights_output.write_text(
    "\n".join(str(value) for value in weights) + "\n"
)
word_weights: list[list[int]] = []
ordinal = 0
for count, nw in zip(COUNTS, WORDS, strict=True):
    for local_word in range(nw):
        values = []
        for bit in range(64):
            local = 64 * local_word + bit
            values.append(weights[ordinal + local] if local < count else 0)
        word_weights.append(values)
    ordinal += count
assert ordinal == len(weights) and len(word_weights) == sum(WORDS)

maximum = -1
maximum_mask = None
columns = 0
violations = 0
mask_fnv = FNV_OFFSET
mask_rows: list[tuple[str, int]] = []
with args.pool.open(newline="") as handle:
    reader = csv.reader(handle); next(reader)
    for fields in reader:
        columns += 1
        mask = int(fields[0], 16)
        rank = int(fields[1])
        assert mask.bit_count() == rank in (8, 9)
        mask_fnv = fnv_add(mask_fnv, mask)
        mask_rows.append((fields[0], rank))
        score = 0
        for wi, token in enumerate(fields[3:]):
            word = int(token, 16)
            while word:
                bit = (word & -word).bit_length() - 1
                score += word_weights[wi][bit]
                word &= word - 1
        if score > maximum:
            maximum = score; maximum_mask = fields[0]
        if score > SCALE:
            violations += 1
assert columns == 37_497 and violations == 0
assert len({mask for mask, _ in mask_rows}) == columns
with args.masks_output.open("w", newline="", encoding="ascii") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(("mask_hex", "rank"))
    writer.writerows(mask_rows)
dual_numerator = sum(weights)
assert dual_numerator > 35 * SCALE
payload = {
    "scope": "all 37,497 exported retained-pool columns only",
    "pool_mask_count": columns,
    "pool_mask_fnv64": f"{mask_fnv:016x}",
    "scale": SCALE,
    "nonzero_weights": sum(value > 0 for value in weights),
    "dual_numerator": dual_numerator,
    "dual_value_decimal": f"{dual_numerator / SCALE:.9f}",
    "maximum_column_numerator": maximum,
    "maximum_column_value_decimal": f"{maximum / SCALE:.9f}",
    "maximum_column_mask_hex": maximum_mask,
    "column_violations": violations,
    "consequence": "Every integral cover from this pool has at least 36 masks by exact integer weak duality.",
    "limitation": "The certificate says nothing about rank-8/rank-9 responders omitted from the retained pool.",
}
args.output.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
