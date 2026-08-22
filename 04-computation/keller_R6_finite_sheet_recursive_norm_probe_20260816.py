#!/usr/bin/env python3
"""Exact finite-sheet gate for the fixed Keller polynomial R_6.

The proved fixed-map normalizations are

    H   = 2^6 L N(L),
    J   = 2^35 L^7 N(H),
    G   = L^43 N(J),
    R_5 = L^271 N(G),
    R_6 = L^1699 N(R_5).

At the canonical finite inverse point q=(2,5/6,-7/8) above the old
boundary L=0, this probe evaluates R_6(q) without expanding R_6.  It
hash-pins and inherits the exact finite-algebra engine of THM-3521, then
extends its recursion by precisely one cubic rung.  The primary route
descends to L in a complete 243-sheet algebra.  The independent bottom
representation stops one rung earlier and evaluates the frozen polynomial
H in the complete 81-sheet algebra; its transitive norm is checked against
the determinant of the literal 81 by 81 multiplication matrix.

Every division is guarded by a full regular-representation unit test,
every cubic discriminant is a unit, and every universal inverse point is
checked by direct substitution into F.  Nonzero good reductions prove the
rational value R_6(q) is nonzero.  The omitted-64 run is a hostile control
for the nonmonic H normalization and must differ by exactly 64^(-81).

Scope: one fixed map, one named finite sheet, and the next old-L valuation.
This is not a sixth-image, degree-243 separability, irreducibility,
all-level, or general Jacobian-conjecture computation.  ``require`` remains
active under ``python -O``.
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

support_raw = SUPPORT_PATH.read_bytes()
if hashlib.sha256(support_raw).hexdigest() != SUPPORT_SHA256:
    raise RuntimeError("THM-3521 recursive support script changed")
support_text = support_raw.decode("utf-8")
if support_text.count(SUPPORT_SENTINEL) != 1:
    raise RuntimeError("THM-3521 definition boundary changed")
support_prefix = support_text.split(SUPPORT_SENTINEL, 1)[0]
support_module = types.ModuleType("thm3521_recursive_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
support_namespace = support_module.__dict__
exec(compile(support_prefix, str(SUPPORT_PATH), "exec"), support_namespace)

require = support_namespace["require"]
PrimeField = support_namespace["PrimeField"]
invariants = support_namespace["invariants"]
fixed_map = support_namespace["fixed_map"]
absolute_norm = support_namespace["absolute_norm"]
flat_norm = support_namespace["flat_norm"]
evaluate_from_L = support_namespace["evaluate_from_L"]
evaluate_from_H = support_namespace["evaluate_from_H"]
fraction_mod = support_namespace["fraction_mod"]
FINITE_Q = support_namespace["FINITE_Q"]
OLD_BOUNDARY_TARGET = support_namespace["OLD_BOUNDARY_TARGET"]
H_RAW_SHA256 = support_namespace["H_RAW_SHA256"]
H_terms = support_namespace["H_terms"]
H_poly = support_namespace["H_poly"]

# Extend the proved transition table by the single fixed-map rung under test.
# The dictionary is shared with evaluate_from_H's inherited global namespace.
TRANSITIONS = support_namespace["TRANSITIONS"]
require(TRANSITIONS == {1: (2**6, 1), 2: (2**35, 7), 3: (1, 43), 4: (1, 271)},
        "THM-3521 transition table changed")
TRANSITIONS[5] = (1, 1699)


def evaluate_prime_r6(prime: int):
    """Evaluate both complete finite-sheet representations at one good prime."""

    field = PrimeField(prime)
    q = tuple(fraction_mod(value, prime) for value in FINITE_Q)
    boundary_target = tuple(fraction_mod(value, prime) for value in OLD_BOUNDARY_TARGET)
    require(fixed_map(field, q) == boundary_target, f"p={prime}: canonical finite sheet changed")
    require(field.is_zero(invariants(field, boundary_target)["L"]),
            f"p={prime}: target left L=0")
    require(not field.is_zero(invariants(field, q)["L"]),
            f"p={prime}: finite point hit L=0")

    l_ledger: list[tuple] = []
    h_ledger: list[tuple] = []
    h_capture: list[tuple] = []
    l_bottom: list[int] = []
    value_from_l = evaluate_from_L(
        field, q, 5, l_ledger, constants=TRANSITIONS, bottom_capture=l_bottom
    )
    value_from_h = evaluate_from_H(field, q, 5, h_ledger, h_capture)
    require(value_from_l == value_from_h,
            f"p={prime}: L and explicit-H routes disagree")
    require(value_from_l != 0, f"p={prime}: R6 finite-sheet reduction vanished")
    require(len(l_ledger) == 5 and len(h_ledger) == 4,
            f"p={prime}: tower depth changed")
    require([row[1] for row in l_ledger] == [1, 3, 9, 27, 81],
            f"p={prime}: L dimensions changed")
    require([row[1] for row in h_ledger] == [1, 3, 9, 27],
            f"p={prime}: H dimensions changed")
    require(len(l_bottom) == 1 and len(h_capture) == 1,
            f"p={prime}: bottom capture changed")

    # Multiplicativity of the cubic norm gives the localized identity
    # R6=2^1431 L^1699 N(L)^271 N^2(L)^43 N^3(L)^7 N^4(L) N^5(L).
    l_norm_orbit = tuple(row[2] for row in l_ledger) + (l_bottom[0],)
    unrolled_value = pow(2, 1431, prime)
    for orbit_value, exponent in zip(l_norm_orbit, (1699, 271, 43, 7, 1, 1)):
        unrolled_value = unrolled_value * pow(orbit_value, exponent, prime) % prime
    require(unrolled_value == value_from_l,
            f"p={prime}: unrolled norm-orbit identity failed")
    require(all(l_norm_orbit), f"p={prime}: a norm-orbit factor vanished")

    h_ring, h_value = h_capture[0]
    require(h_ring.dimension == 81, f"p={prime}: H bottom dimension changed")
    h_norm_recursive = absolute_norm(h_ring, h_value)
    h_norm_flat = flat_norm(h_ring, h_value)
    require(h_norm_recursive == h_norm_flat,
            f"p={prime}: transitive/flat H norm mismatch")
    require(h_norm_flat != 0, f"p={prime}: explicit H became a zero divisor")

    hostile_constants = dict(TRANSITIONS)
    hostile_constants[1] = (1, 1)
    hostile_ledger: list[tuple] = []
    wrong_h_normalization = evaluate_from_L(
        field, q, 5, hostile_ledger, constants=hostile_constants
    )
    expected_wrong = value_from_l * pow(pow(64, 81, prime), prime - 2, prime) % prime
    require(wrong_h_normalization == expected_wrong,
            f"p={prime}: scalar norm exponent changed")
    require(wrong_h_normalization != value_from_l,
            f"p={prime}: normalization hostile did not fire")

    return {
        "prime": prime,
        "value": value_from_l,
        "wrong_h_normalization": wrong_h_normalization,
        "h_flat_norm": h_norm_flat,
        "l_norm_orbit": l_norm_orbit,
        "l_ledger": l_ledger,
        "h_ledger": h_ledger,
    }


PRIMES = (101, 103, 107)
records = [evaluate_prime_r6(prime) for prime in PRIMES]

R6_FACE_PAIR = (10663, 3867)
R7_FACE_PAIR = (
    7 * R6_FACE_PAIR[0] - 2 * R6_FACE_PAIR[1],
    3 * R6_FACE_PAIR[0] - 2 * R6_FACE_PAIR[1],
)
require(R7_FACE_PAIR == (66907, 24255), "THM-3522 next packet pair changed")

semantic_lines = [
    f"support={SUPPORT_SHA256}",
    "identity=2^1431;1699,271,43,7,1,1",
]
for record in records:
    semantic_lines.append(
        f"p={record['prime']};R6={record['value']};wrong64={record['wrong_h_normalization']};"
        f"Hflat={record['h_flat_norm']};Lorbit={record['l_norm_orbit']}"
    )
    for row in record["l_ledger"]:
        semantic_lines.append("L:" + ":".join(map(str, row)))
    for row in record["h_ledger"]:
        semantic_lines.append("H:" + ":".join(map(str, row)))
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
require(tuple(record["value"] for record in records) == (26, 70, 69),
        "R6 residue ledger changed")
require(tuple(record["h_flat_norm"] for record in records) == (37, 31, 13),
        "flat-H norm ledger changed")
require(tuple(record["wrong_h_normalization"] for record in records) == (67, 47, 3),
        "normalization-hostile ledger changed")
require(
    tuple(record["l_norm_orbit"] for record in records)
    == (
        (16, 12, 72, 9, 49, 97),
        (12, 53, 22, 85, 76, 94),
        (38, 45, 28, 3, 17, 17),
    ),
    "six-factor norm-orbit ledger changed",
)
require(semantic_sha256 == "a2ede01095e73ad727285743b83d6502ff37c0ab772a19e0a03fe9036ba5f7b8",
        "semantic ledger changed")

print("== fixed Keller R6 finite-sheet recursive norm gate ==")
print(f"THM-3521 support sha256={SUPPORT_SHA256}")
print(f"H raw sha256={H_RAW_SHA256}; terms={len(H_terms)}; degrees={H_poly.degree_list()}")
print("point q=(2,5/6,-7/8); F(q)=(2/27,1,1); L(F(q))=0; L(q)=241465/1728")
print("primary route: five cubic algebras, descending to L on 243 sheets")
print("independent bottom representation: four cubic algebras, evaluating frozen H on 81 sheets")
for record in records:
    print(
        f"p={record['prime']}: R6(q)={record['value']}; "
        f"explicit-H flat norm={record['h_flat_norm']}; "
        f"wrong-H-normalization control={record['wrong_h_normalization']}"
    )
    print(f"  (L,N(L),N^2(L),N^3(L),N^4(L),N^5(L))={record['l_norm_orbit']}")
    for route_name in ("l_ledger", "h_ledger"):
        compact = [
            (row[1], row[2], row[3], row[4], row[5][:12])
            for row in record[route_name]
        ]
        print(f"  {route_name}: (base_dim,Norm(L),Norm(S),Norm(disc),source_sha12)={compact}")
print(f"semantic sha256={semantic_sha256}")
print("finite-sheet verdict: R6(2,5/6,-7/8) != 0 over Q")
print("conditional geometric consequence with proved THM-3522 packet:")
print("  v_L(N(R6))=-10663 and R7=L^10663*N(R6) is polynomial and coprime to L")
print("  THM-3522 then gives the complete fixed-chart packet A(66907,24255) for R7")
print("scope: fixed map/old-L gate only; image, degree 243, all-level, and general JC remain open")
print("all exact checks passed")
