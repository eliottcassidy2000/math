"""Bounded source/import audit; this is not a Lean build or an axiom audit.

Run: python openai_euler_static_audit_20260908.py --repo CHECKOUT [--pdf FILE]
Only Git-tracked Lean sources are examined; dependency packages are not scanned.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
from fractions import Fraction

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path, required=True)
parser.add_argument("--pdf", type=Path)
args = parser.parse_args()
root = args.repo.resolve()
tracked = subprocess.check_output(["git", "ls-files", "*.lean"], cwd=root, text=True).splitlines()
modules = {str(Path(p).with_suffix("")).replace("\\", ".").replace("/", "."): root / p for p in tracked}
contents = {m: p.read_text(encoding="utf-8") for m, p in modules.items()}
imports = {m: re.findall(r"^import\s+(\S+)", s, re.M) for m, s in contents.items()}

def reachable(start):
    seen, stack = set(), [start]
    while stack:
        m = stack.pop()
        if m in seen:
            continue
        seen.add(m)
        stack.extend(n for n in imports.get(m, []) if n in modules)
    return seen

result = {
    "status": "FINITE-EXACT static source audit, not a kernel or semantic proof audit",
    "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
    "roots": [],
}
for start in ["Euler", "Euler.Solution", "NavierStokes", "NavierStokes.ComparatorSolution"]:
    reach = reachable(start)
    text = "\n".join(contents[m] for m in sorted(reach))
    result["roots"].append({
        "module": start,
        "tracked_local_modules_including_root": len(reach),
        "reaches_Euler_Solution": "Euler.Solution" in reach,
        "comparator_challenge_modules": sorted(m for m in reach if m.startswith("ComparatorChallenges.")),
        "lexical_sorry_count": len(re.findall(r"\bsorry\b", text)),
        "lexical_admit_count": len(re.findall(r"\badmit\b", text)),
        "line_start_axiom_declaration_count": len(re.findall(r"^\s*axiom\b", text, re.M)),
    })

# Exact rational controls for the elementary profile expression; these finite
# values supplement the symbolic proof in the accompanying report.
def require(condition, message):
    if not condition:
        raise RuntimeError(message)

controls = []
for delta in [Fraction(1, 100), Fraction(1, 2), Fraction(1), Fraction(3)]:
    r = 1 / (1 + delta)
    for cosine in [Fraction(-1), Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(1)]:
        den = (1 + delta)**2 - 2 * (1 + delta) * cosine + 1
        derivative = ((1 + delta) * cosine - 1) / den
        poisson = (1 - r*r) / (1 - 2*r*cosine + r*r)
        require(derivative == (poisson - 1) / 2, "profile identity")
        require(-1/(2+delta) <= derivative <= 1/delta, "slope range")
        if cosine == -1:
            require(derivative == -1/(2+delta), "lower equality")
        if cosine == 1:
            require(derivative == 1/delta, "upper equality")
        controls.append({"delta": str(delta), "cosine": str(cosine), "derivative": str(derivative)})
result["profile_controls"] = controls
result["sign_hostile_control"] = {
    "description": "a=1, m=e1, w=e2, M(e2)=-e1, delta=1/100, slope=100",
    "flux": -1,
    "rank_one_positive_eigenvalue": 200,
    "proposed_sign_free_upper_bound_using_abs_flux": 2,
    "verdict": "fails: the nonnegative flux sidecar is necessary",
}
if args.pdf is not None:
    pdf = args.pdf
    result["pdf"] = {"sha256": hashlib.sha256(pdf.read_bytes()).hexdigest()}
print(json.dumps(result, indent=2))
