#!/usr/bin/env python3
"""Compare independently computed certificate outputs, without rerunning either."""

from hashlib import sha256
from pathlib import Path
import json

base = Path(__file__).resolve().parent.parent
production = json.loads((base / "negative-receipt.json").read_text())
independent = json.loads((base / "audit/independent-receipt.json").read_text())
comparisons = {
    "identical_cf_digits":
        production["continued_fraction_prefix"] == independent["continued_fraction_digits"],
    "identical_terminal_convergent":
        production["terminal_convergent"] == independent["terminal_convergent"],
    "identical_large_convergent_pairs":
        [(x["q"], x["p"]) for x in production["large_convergents"]]
        == [(x["q"], x["p"]) for x in independent["large_convergents"]],
    "identical_small_shape_rows":
        [(x["p"], x["q"], x["swaps"], x.get("maximum_content"))
         for x in production["small_shapes"]]
        == [(x["p"], x["q"], x["swaps"], x.get("maximum_content"))
            for x in independent["small_rows"]],
    "identical_small_hit_set":
        {(x["p"], x["q"], x["word"], x["fixed_point"])
         for x in production["integer_exceptions"]}
        == {(x["p"], x["q"], x["word"], x["orbit"][0])
            for x in independent["small_integer_hits"]},
    "identical_hostile_q120_count":
        production["broad_hostile_swap_occurrences"]
        == independent["hostile_census"]["through_q120_swaps"],
    "production_normal_optimized_identical":
        (base / "negative-receipt.json").read_bytes()
        == (base / "negative-receipt.optimized.json").read_bytes(),
    "independent_normal_optimized_identical":
        (base / "audit/independent-receipt.json").read_bytes()
        == (base / "audit/independent-receipt.optimized.json").read_bytes(),
}
if not all(comparisons.values()):
    raise RuntimeError(comparisons)
names = [
    "certify_negative_swap_exclusion.py",
    "negative-receipt.json",
    "NEGATIVE_ONE_SWAP_EXCLUSION.md",
    "audit/verify_negative_independent.py",
    "audit/independent-receipt.json",
    "audit/compare_negative_receipts.py",
]
print(json.dumps({
    "status": "PASS",
    "comparisons": comparisons,
    "audited_sha256": {
        name: sha256((base / name).read_bytes()).hexdigest() for name in names
    },
}, indent=2, sort_keys=True))
