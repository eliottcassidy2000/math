#!/usr/bin/env python3
"""Split-outer-branch audit of the fixed R_6 finite-sheet value.

This companion realizes the complete 243-sheet universe differently from
the primary nested tower.  At p=71 the outer inverse cubic over the canonical
point q splits into three scalar branches.  The script enumerates those roots,
constructs and checks all three inverse points directly, evaluates R_5 on a
complete 81-sheet tower above each point, and multiplies the three branch
values.  Thus no 243-dimensional outer algebra and no expansion of R_6 is
used.  A direct nested-tower value at the same good prime is retained only as
an equality control after the split product has been formed.

Scope is one fixed map and one finite old-L sheet.  This is not the open
degree-243 separability/image gate or a general Jacobian-conjecture claim.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/keller_R5_finite_sheet_recursive_norm_probe_20260816.py"
SUPPORT_SHA256 = "a201191410e39d47fbf607191e8bd597453c697f134d3803466694b680d8c60d"
SUPPORT_SENTINEL = "\nrecords = [evaluate_prime(prime) for prime in (101, 103, 107)]\n"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


support_raw = SUPPORT_PATH.read_bytes()
require(hashlib.sha256(support_raw).hexdigest() == SUPPORT_SHA256,
        "THM-3521 recursive support script changed")
support_text = support_raw.decode("utf-8")
require(support_text.count(SUPPORT_SENTINEL) == 1,
        "THM-3521 definition boundary changed")
support_module = types.ModuleType("thm3521_split_branch_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(
    compile(support_text.split(SUPPORT_SENTINEL, 1)[0], str(SUPPORT_PATH), "exec"),
    support_module.__dict__,
)

ns = support_module.__dict__
PrimeField = ns["PrimeField"]
FINITE_Q = ns["FINITE_Q"]
fraction_mod = ns["fraction_mod"]
invariants = ns["invariants"]
fixed_map = ns["fixed_map"]
evaluate_from_L = ns["evaluate_from_L"]
TRANSITIONS = ns["TRANSITIONS"]

prime = 71
field = PrimeField(prime)
q = tuple(fraction_mod(value, prime) for value in FINITE_Q)
q_inv = invariants(field, q)
require(q_inv["L"] != 0 and q_inv["S"] != 0, "outer inverse denominator vanished")

outer_roots = [
    w
    for w in range(prime)
    if (q_inv["L"] * w**3 + q_inv["T"] * w - 2 * q[2]) % prime == 0
]
require(outer_roots == [10, 23, 38], "outer split-root ledger changed")

inverse_2s = pow(2 * q_inv["S"], prime - 2, prime)
inverse_8s = pow(8 * q_inv["S"], prime - 2, prime)
rows = []
norm_r5 = 1
all_ledgers = []
for w in outer_roots:
    y = (
        q_inv["Y0"]
        + 6 * q_inv["L"] * w
        - 3 * q_inv["K"] * q_inv["L"] * w * w
    ) * inverse_2s % prime
    z = (
        q_inv["Z0"]
        + 6 * q_inv["L"] * q_inv["A1"] * w
        - 9 * q_inv["L"] * q_inv["A2"] * w * w
    ) * inverse_8s % prime
    point = (w, y, z)
    require(fixed_map(field, point) == q, f"w={w}: split inverse graph changed")
    ledger: list[tuple] = []
    bottom: list[int] = []
    r5_value = evaluate_from_L(
        field, point, 4, ledger, constants=TRANSITIONS, bottom_capture=bottom
    )
    require(r5_value != 0 and len(bottom) == 1, f"w={w}: R5 branch vanished")
    require([entry[1] for entry in ledger] == [1, 3, 9, 27],
            f"w={w}: inner dimensions changed")
    require(all(entry[index] for entry in ledger for index in (2, 3, 4)),
            f"w={w}: inner unit/discriminant gate vanished")
    rows.append((w, y, z, r5_value, bottom[0]))
    all_ledgers.append(ledger)
    norm_r5 = norm_r5 * r5_value % prime

require(
    rows == [(10, 0, 2, 49, 6), (23, 1, 36, 22, 5), (38, 18, 60, 60, 54)],
    "split R5 branch ledger changed",
)

r6_split = pow(q_inv["L"], 1699, prime) * norm_r5 % prime
require(r6_split != 0, "split-branch R6 value vanished")

# Equality control: the inherited engine now keeps the outer cubic unsplit.
direct_ledger: list[tuple] = []
direct_bottom: list[int] = []
transitions_r6 = dict(TRANSITIONS)
transitions_r6[5] = (1, 1699)
r6_direct = evaluate_from_L(
    field, q, 5, direct_ledger, constants=transitions_r6,
    bottom_capture=direct_bottom,
)
require(r6_split == r6_direct, "split and nested R6 representations disagree")
require(r6_split == 9, "split R6 residue ledger changed")
require([entry[1] for entry in direct_ledger] == [1, 3, 9, 27, 81],
        "direct control dimensions changed")

# Hostile: deleting one named outer branch must not reproduce the norm.
hostile_omit_last = (
    pow(q_inv["L"], 1699, prime) * rows[0][3] * rows[1][3]
) % prime
require(hostile_omit_last != r6_split, "omitted-branch hostile did not fire")
require(hostile_omit_last == 25, "omitted-branch hostile ledger changed")
reverse_norm_r5 = 1
for row in reversed(rows):
    reverse_norm_r5 = reverse_norm_r5 * row[3] % prime
require(reverse_norm_r5 == norm_r5, "branch-order product control changed")

semantic_lines = [
    f"p={prime};outerL={q_inv['L']};roots={outer_roots};R6={r6_split};omit={hostile_omit_last}",
]
for row, ledger in zip(rows, all_ledgers):
    semantic_lines.append("row=" + ":".join(map(str, row)))
    semantic_lines.extend("gate=" + ":".join(map(str, entry)) for entry in ledger)
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
require(semantic_sha256 == "6384038b2f7e386289c7af89ac1c5a020aa18320e6957a59469ff013de396490",
        "split-branch semantic ledger changed")

print("== fixed Keller R6 split-outer-branch audit ==")
print(f"THM-3521 support sha256={SUPPORT_SHA256}")
print(f"p={prime}; outer L={q_inv['L']}; roots={outer_roots}")
print(f"branch rows (w,y,z,R5,bottom-N4L)={rows}")
for row, ledger in zip(rows, all_ledgers):
    compact = [(entry[1], entry[2], entry[3], entry[4], entry[5][:12]) for entry in ledger]
    print(f"  w={row[0]} gates (dim,Norm(L),Norm(S),Norm(disc),sha12)={compact}")
print(f"split product R6(q) mod {prime}={r6_split}; direct nested control={r6_direct}")
print(f"omit-w={rows[-1][0]} hostile={hostile_omit_last}")
print(f"semantic sha256={semantic_sha256}")
print("scope: fixed-map finite-sheet audit only; no image, degree-243 separability, all-level, or JC claim")
print("all exact checks passed")
