#!/usr/bin/env python3
"""Independent split-outer audit of the fixed R_7 finite-sheet value.

Over F_71 the first inverse cubic above q=(2,5/6,-7/8) splits into the three
named roots 10,23,38.  This script evaluates R_6 on a complete 243-sheet
tower above each scalar branch and then multiplies the three values with
L(q)^10663.  It therefore realizes all 729 sheets without constructing one
729-dimensional algebra and without importing the recursive R_7 candidate.

Each branch is represented twice: the old audited coefficient algebra descends
to L, while a frozen 361-term H polynomial supplies a disjoint bottom
representation whose 81-dimensional terminal norm is checked both
transitively and by a literal regular determinant.

Scope is one fixed map and one finite old-L sheet.  This is not an image,
irreducibility, all-level, arbitrary-map, or general-JC computation.
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
EXPECTED_SEMANTIC_SHA256 = (
    "d8fd9ddf1f8679b90434d1f6ffa6a717c1725e3dcb5703c040e1ab724081e72b"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
require(hashlib.sha256(support_raw).hexdigest() == SUPPORT_SHA256,
        "THM-3521 recursive support changed")
support_text = support_raw.decode("utf-8")
require(support_text.count(SUPPORT_SENTINEL) == 1,
        "THM-3521 definition boundary changed")
support_module = types.ModuleType("thm3521_R7_split_support")
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
evaluate_from_H = ns["evaluate_from_H"]
absolute_norm = ns["absolute_norm"]
flat_norm = ns["flat_norm"]
TRANSITIONS = ns["TRANSITIONS"]

require(
    TRANSITIONS == {1: (2**6, 1), 2: (2**35, 7), 3: (1, 43), 4: (1, 271)},
    "inherited transition table changed",
)
TRANSITIONS[5] = (1, 1699)

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
all_l_ledgers = []
all_h_ledgers = []
norm_r6 = 1
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
    require(fixed_map(field, point) == q, ("split inverse graph", w))

    l_ledger: list[tuple] = []
    h_ledger: list[tuple] = []
    l_bottom: list[int] = []
    h_capture: list[tuple] = []
    value_from_l = evaluate_from_L(
        field,
        point,
        5,
        l_ledger,
        constants=TRANSITIONS,
        bottom_capture=l_bottom,
    )
    value_from_h = evaluate_from_H(field, point, 5, h_ledger, h_capture)
    require(value_from_l == value_from_h, ("L/H branch disagreement", w))
    require(value_from_l != 0, ("R6 split branch vanished", w))
    require([entry[1] for entry in l_ledger] == [1, 3, 9, 27, 81],
            ("L-route dimensions", w))
    require([entry[1] for entry in h_ledger] == [1, 3, 9, 27],
            ("H-route dimensions", w))
    require(len(l_bottom) == 1 and len(h_capture) == 1, ("bottom capture", w))
    require(all(entry[index] for entry in l_ledger for index in (2, 3, 4)),
            ("L-route unit gate", w))
    require(all(entry[index] for entry in h_ledger for index in (2, 3, 4)),
            ("H-route unit gate", w))

    norm_orbit = tuple(entry[2] for entry in l_ledger) + (l_bottom[0],)
    unrolled = pow(2, 1431, prime)
    for orbit_value, exponent in zip(norm_orbit, (1699, 271, 43, 7, 1, 1)):
        unrolled = unrolled * pow(orbit_value, exponent, prime) % prime
    require(unrolled == value_from_l and all(norm_orbit), ("branch norm orbit", w))

    h_ring, h_value = h_capture[0]
    require(h_ring.dimension == 81, ("H bottom dimension", w))
    h_recursive = absolute_norm(h_ring, h_value)
    h_flat = flat_norm(h_ring, h_value)
    require(h_recursive == h_flat != 0, ("frozen-H determinant", w))

    rows.append((w, y, z, value_from_l, h_flat, norm_orbit))
    all_l_ledgers.append(tuple(l_ledger))
    all_h_ledgers.append(tuple(h_ledger))
    norm_r6 = norm_r6 * value_from_l % prime

require(tuple((row[0], row[1], row[2]) for row in rows)
        == ((10, 0, 2), (23, 1, 36), (38, 18, 60)),
        "split inverse-point ledger changed")

r7_split = pow(q_inv["L"], 10663, prime) * norm_r6 % prime
require(r7_split != 0, "split R7 value vanished")

# Hostiles retain named branches and exact H-leaf multiplicity.
omit_last = pow(q_inv["L"], 10663, prime) * rows[0][3] * rows[1][3] % prime
require(omit_last != r7_split, "omitted-branch hostile did not fire")
wrong_h_normalization = r7_split * pow(pow(64, 243, prime), prime - 2, prime) % prime
require(wrong_h_normalization != r7_split, "H-normalization hostile did not fire")
reverse_product = 1
for row in reversed(rows):
    reverse_product = reverse_product * row[3] % prime
require(reverse_product == norm_r6, "branch-order product changed")

semantic_lines = [
    f"support={SUPPORT_SHA256};p={prime};outerL={q_inv['L']};roots={outer_roots}",
    f"R7={r7_split};omit={omit_last};wrong64={wrong_h_normalization}",
]
for row, l_ledger, h_ledger in zip(rows, all_l_ledgers, all_h_ledgers):
    semantic_lines.append(f"row={row}")
    semantic_lines.extend("L=" + ":".join(map(str, entry)) for entry in l_ledger)
    semantic_lines.extend("H=" + ":".join(map(str, entry)) for entry in h_ledger)
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        ("split R7 semantic drift", semantic_sha256),
    )

print("== fixed Keller R7 split-outer independent audit ==")
print(f"THM-3521 support sha256={SUPPORT_SHA256}")
print(f"p={prime};outer L={q_inv['L']};roots={outer_roots}")
print("branch rows=(w,y,z,R6,frozen-H-flat-norm,(L,N(L),...,N^5(L)))")
for row in rows:
    print(f"  {row}")
for row, l_ledger, h_ledger in zip(rows, all_l_ledgers, all_h_ledgers):
    l_compact = tuple((entry[1], entry[2], entry[3], entry[4]) for entry in l_ledger)
    h_compact = tuple((entry[1], entry[2], entry[3], entry[4]) for entry in h_ledger)
    print(f"  w={row[0]} L-gates(dim,Norm(L),Norm(S),Norm(disc))={l_compact}")
    print(f"  w={row[0]} H-gates(dim,Norm(L),Norm(S),Norm(disc))={h_compact}")
print(f"split_product_R7(q)_mod_{prime}={r7_split};omit_w_{rows[-1][0]}={omit_last}")
print(f"wrong_H_normalization={wrong_h_normalization};semantic_sha256={semantic_sha256}")
print("independent verdict: complete split 3*243-sheet R7 finite value is nonzero")
print("scope: fixed-map old-L gate only; no image, irreducibility, all-level, arbitrary-map, or general JC claim")
print("all independent exact checks passed")
