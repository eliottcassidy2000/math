#!/usr/bin/env python3
"""Review bindings and arithmetic controls, not a replacement for the proofs."""
from pathlib import Path
import hashlib
import json
from math import comb

BASE = Path(__file__).resolve().parents[1]
EXPECTED = "d8815510d96f84aa15de8c52e25139a61de14de6a8c3ca37571039f3d097ec06"
checks = []

def read(path):
    return (BASE / path).read_bytes()

def data(path):
    return json.loads(read(path))

def sha(path):
    return hashlib.sha256(read(path)).hexdigest()

def need(value, name):
    if not value:
        raise RuntimeError(name)
    checks.append(name)

report = "session23-root/SESSION23_SYNTHESIS.md"
need(sha(report) == EXPECTED, "exact final synthesis bytes")
bindings = [
    ("refined", "session23-refined/PROOF.md", data("session23-refined-review/receipt.json")["proof_sha256"]),
    ("renewal", "session23-root/RENEWAL_LIMIT_LAWS.md", data("session23-renewal-review/receipt.json")["proof_sha256"]),
    ("partial paths", "session23-symmetry/README.md", data("session23-symmetry-review/review-manifest.json")["sha256"]["session23-symmetry/README.md"]),
    ("slack compiler", "session23-connections/SYMMETRY_SLACK_RATIONALITY.md", data("session23-slack-review/final-review-receipt.json")["final_proof_sha256"]),
    ("slack renewal", "session23-root/SLACK_RENEWAL_TRANSFER.md", data("session23-slack-renewal-review/final-review-receipt.json")["final_proof_sha256"]),
    ("sharp slack", "session23-sharp-slack/SHARP_POLE_ORDER.md", data("session23-sharp-slack-review/source-pins.json")["session23-sharp-slack/SHARP_POLE_ORDER.md"]),
    ("slack spectrum", "session23-slack-spectrum/README.md", data("session23-slack-spectrum-review/source-pins.json")["reviewed_producer"]["session23-slack-spectrum/README.md"]),
]
for name, path, pin in bindings:
    need(sha(path) == pin, name + " current proof matches independent review pin")
need(b"Explicit sufficient recurrence onset" in read("session23-sharp-slack-review/AUDIT.md"), "sharp audit includes onset addition")

pairs = [
    ("refined independent", "session23-refined-review/normal.json", "session23-refined-review/optimized.json"),
    ("renewal independent", "session23-renewal-review/normal.json", "session23-renewal-review/optimized.json"),
    ("partial paths independent", "session23-symmetry-review/audit-normal.json", "session23-symmetry-review/audit-optimized.json"),
    ("slack compiler producer order seven", "session23-connections/VERIFICATION_ORDER7.txt", "session23-connections/VERIFICATION_ORDER7_OPTIMIZED.txt"),
    ("slack renewal independent", "session23-slack-renewal-review/normal.json", "session23-slack-renewal-review/optimized.json"),
    ("sharp slack independent", "session23-sharp-slack-review/receipt.json", "session23-sharp-slack-review/receipt-optimized.json"),
    ("slack spectrum independent", "session23-slack-spectrum-review/normal.json", "session23-slack-spectrum-review/optimized.json"),
]
for name, normal, optimized in pairs:
    need(read(normal) == read(optimized), name + " recorded normal/optimized outputs agree")

orders = [3 * (r + 1) + 6 * (2 if r == 1 else r // 2 + 1) for r in range(8)]
need(orders == [9, 18, 21, 24, 33, 36, 45, 48], "power-three recurrence-order prefix")

fib = [0, 1]
for _ in range(25):
    fib.append(sum(fib[-2:]))
need([fib[r+1] * fib[r+2] for r in range(7)] == [1,2,6,15,40,104,273], "inherited Fibonacci-product prefix")
for e in range(2,21):
    need(sum(comb(e-t,t) for t in range(e//2+1)) == fib[e+1], "strong fixed-spectrum coefficient sum e=" + str(e))
for r in range(11):
    need(1 + sum(fib[e+1]**2 for e in range(1,r+1)) == fib[r+1]*fib[r+2], "cumulative half-exponent budget r=" + str(r))

spectrum = data("session23-slack-spectrum-review/normal.json")
hostiles = [next(row for row in spectrum["rows"] if row["word"] == word) for word in ([1,3,3,1],[3,1,1,3])]
destinations = []
for row in hostiles:
    name = "".join(map(str,row["word"]))
    need(tuple(row[k] for k in ("n","H","aut","h","b")) == (10,585,9,65,9), "same-order hostile signature " + name)
    need(row["H"] == row["aut"] * row["h"], "raw/orbit quotient " + name)
    need((row["h"]-row["b"])//2 == 28, "local incidence degree " + name)
    keys = {tuple(port["target_merged_canonical_key"]) for port in row["literal_port_receipts"]}
    need(len(keys) == len(row["literal_port_receipts"]) == 28, "28 distinct recorded destinations " + name)
    destinations.append(keys)
need(not destinations[0].intersection(destinations[1]), "same-order recorded destination sets are disjoint")

publication = []
for file, commit in [
    ("SESSION23_REFINED_PUBLICATION_RECEIPT.json", "da7fc32234c75670a1a668fbb1db4f62a14c5021"),
    ("SESSION23_RENEWAL_PUBLICATION_RECEIPT.json", "9eea439e3ee7a1488b0269fba0d3ded20955d1cd"),
]:
    receipt = data(file)
    need(receipt["commit"] == commit and receipt["exact_readback"] is True, "recorded publication checkpoint " + commit[:8])
    publication.append({"receipt": file, "commit": commit, "exact_readback": True})

out = {
    "status": "PASS",
    "review_scope": "Synthesis fidelity, source bindings, arithmetic and recorded-output comparisons; not a rerun of the seven independent proof audits or fresh remote readback.",
    "report": report,
    "report_sha256": EXPECTED,
    "proof_bindings": [{"lane": n,"path": p,"sha256": h} for n,p,h in bindings],
    "checks": checks,
    "total_checks": len(checks),
    "power_three_recurrence_orders": orders,
    "same_order_hostile_destination_sizes": list(map(len,destinations)),
    "same_order_hostile_intersection_size": len(destinations[0] & destinations[1]),
    "recorded_publication_checkpoints": publication,
    "publication_scope": "Only these two checkpoints are verified in the supplied receipts. Later publication is not established by this synthesis review.",
    "repository_or_git_operations": False,
}
print(json.dumps(out, indent=2) + "\n", end="")
