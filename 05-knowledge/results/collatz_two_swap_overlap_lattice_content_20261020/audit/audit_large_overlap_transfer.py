#!/usr/bin/env python3
"""Exact threshold and replay audit for the p>=4096 overlap corollary."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, gcd
from pathlib import Path
import argparse
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--endpoint-dir", type=Path,
                    help="Frozen endpoint-exclusion packet; defaults to a staging or canonical sibling.")
args = parser.parse_args()
if args.endpoint_dir is not None:
    FROZEN = args.endpoint_dir.resolve()
else:
    candidates = [HERE.parent.parent / name for name in
                  ("two-edit-endpoint-closure", "collatz_two_swap_endpoint_exclusion_20261020")]
    FROZEN = next((p for p in candidates if p.is_dir()), None)
    if FROZEN is None:
        parser.error("Cannot find the frozen endpoint packet; pass --endpoint-dir.")
checks = 0


def need(ok, message):
    global checks
    checks += 1
    if not ok:
        raise ArithmeticError(message)


s = sum((F(23, 6) ** j / factorial(j) for j in range(11)), F(0))
need(s > 46, "log46 upper certificate")
need(2 * (F(1, 2) + F(1, 2) ** 3 / 3) == F(13, 12), "log3 lower certificate")
need(4096 == 16 ** 3, "threshold cube")
need(F(46, 3) / 16 - F(13, 12) == -F(1, 8), "exponential margin")
need(-F(1, 8) < -F(1, 9), "frozen estimate domination")
need(4096 >= 100 > 25, "inherited thresholds")

script = FROZEN / "audit" / "audit_endpoint_large_range.py"
frozen_receipt = FROZEN / "audit" / "large-range.receipt.json"
replay = subprocess.check_output([sys.executable, str(script)])
need(replay == frozen_receipt.read_bytes(), "frozen arithmetic replay differs")
data = json.loads(replay)
need(data["row_count"] == 18, "frozen row count")
rows = [r for r in data["rows"] if r["p"] >= 4096]
need(len(rows) == 16, "relevant row count")
for row in rows:
    p, q = row["p"], row["q"]
    need(4096 <= p < 10 ** 12 and p < q < 2 * p and gcd(p, q) == 1,
         "candidate range or primitive clock")
    need(F(row["gap_lower_bound"]) > F(row["coarse_upper_bound"]),
         "candidate gap contradiction")
need(data["terminal_convergent"][1] >= 10 ** 12, "terminal denominator")
print(json.dumps(dict(
    status="PASS_EXACT; excludes all integral legal k=1 overlaps with p>=4096",
    checks=checks,
    external_dependency=data["external_dependency"],
    threshold=4096,
    log46_exponential_partial_sum=str(s),
    exact_exponential_margin="-1/8 < -1/9",
    frozen_audit_script_sha256=sha256(script.read_bytes()).hexdigest(),
    frozen_audit_receipt_sha256=sha256(frozen_receipt.read_bytes()).hexdigest(),
    inherited_arithmetic_checks=data["checks"],
    relevant_convergents=[dict(q=r["q"], p=r["p"]) for r in rows],
    terminal_convergent=data["terminal_convergent"],
    unverified_residual="Every legal overlap, every coprime q=p+r with 3<=p<4096 and 1<=r<p; both clock signs.",
), indent=2))
