#!/usr/bin/env python3
"""Exact U_full guard-sheet address and joint-drift Fourier probe.

In the refined THM-3479 bank the third character parameter tau translates
the H guard window by 7*tau on the 91-circle.  Splitting a point as 7*a+r
exposes an ell-independent sheet a in F_13 and one of three r chambers.
For a left/right pair, common-tau dependence factors through the common
sheet a_L and relative drift d=a_R-a_L.  This is a lawful guard-factor
address law, not a THM-2471 ancestry or physical-current realization.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from pathlib import Path


P = 13
FIELD_PRIME = 547
XI = 475  # exact order 13 in F_547
EXPECTED_LEDGER_SHA256 = "f345d40c9b589910d83d2fd490ca9376b3cbb86d7aa6da4825677d06075bda7a"
BRIDGE_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


require(pow(XI, P, FIELD_PRIME) == 1 and XI != 1, "split-field root changed")

root = Path(__file__).resolve().parents[1]
bridge_path = root / "04-computation/lrc_half_twist_relation_current_bridge_thm3479.py"
bridge_bytes = bridge_path.read_bytes().replace(b"\r\n", b"\n")
bridge_hash = hashlib.sha256(bridge_bytes).hexdigest()
require(bridge_hash == BRIDGE_SHA256, "THM-3479 endpoint engine changed")
bridge_source = bridge_bytes.decode("utf-8")
require(
    "-13 - 7 * shift[index], 13 - 7 * shift[index]" in bridge_source,
    "H-guard 7*tau translation convention changed",
)
require(
    "U_FULL_REL = (13, 2197, 742586, 1, 183, 27, 131, 53, 313)" in bridge_source,
    "U_full relation tuple changed",
)


# For 7*a+r in [0,91), the forbidden cyclic interval [-13,13) has these
# sheet sets.  Boundaries r=1 are left-closed and r=6 right-open.
CHAMBERS = {
    "left": (Fraction(1, 2), frozenset((12, 0, 1))),
    "middle": (Fraction(7, 2), frozenset((11, 12, 0, 1))),
    "right": (Fraction(13, 2), frozenset((11, 12, 0))),
}


def direct_danger(sheet: int, residue: Fraction) -> bool:
    value = (Fraction(7 * sheet) + residue) % 91
    return value < 13 or value >= 78


def safe(chamber: str, sheet: int) -> int:
    return int(sheet % P not in CHAMBERS[chamber][1])


for chamber, (residue, danger) in CHAMBERS.items():
    direct = frozenset(sheet for sheet in range(P) if direct_danger(sheet, residue))
    require(direct == danger, f"{chamber} cyclic danger set changed")

# The common partition is half-open: r in [0,1), [1,6), [6,7).
# These endpoint controls are load-bearing because the endpoint engine sums
# boundary exponentials rather than only integrating open interiors.
for residue, chamber in ((Fraction(0), "left"), (Fraction(1), "middle"), (Fraction(6), "right")):
    direct = frozenset(sheet for sheet in range(P) if direct_danger(sheet, residue))
    require(direct == CHAMBERS[chamber][1], f"boundary r={residue} convention changed")


def fourier(values, frequency: int) -> int:
    return sum(
        value * pow(XI, -frequency * index % P, FIELD_PRIME)
        for index, value in enumerate(values)
    ) % FIELD_PRIME


# Single-point factorization:
# sum_tau safe_c(a+tau) xi^(-k tau) = xi^(k a) G_c(k).
single_kernels = {}
for chamber in CHAMBERS:
    base_values = tuple(safe(chamber, u) for u in range(P))
    for frequency in range(P):
        kernel = fourier(base_values, frequency)
        single_kernels[(chamber, frequency)] = kernel
        require(kernel != 0, f"{chamber} lost Fourier mode {frequency}")
        for address in range(P):
            direct = fourier(
                tuple(safe(chamber, address + tau) for tau in range(P)),
                frequency,
            )
            expected = pow(XI, frequency * address, FIELD_PRIME) * kernel % FIELD_PRIME
            require(direct == expected, "single guard-sheet covariance failed")

require(
    tuple(single_kernels[(chamber, 0)] for chamber in CHAMBERS) == (10, 9, 10),
    "safe-count profile changed",
)


# Pair factorization.  Put d=b-a.  Then
# sum_tau g_c(a+tau)g_e(b+tau)xi^(-k tau)
#   =xi^(ka) K_(c,e,d)(k).
pair_kernels = {}
pair_factorization_checks = 0
for left_chamber in CHAMBERS:
    for right_chamber in CHAMBERS:
        for drift in range(P):
            base_values = tuple(
                safe(left_chamber, u) * safe(right_chamber, u + drift)
                for u in range(P)
            )
            for frequency in range(P):
                kernel = fourier(base_values, frequency)
                pair_kernels[(left_chamber, right_chamber, drift, frequency)] = kernel
                for left_address in range(P):
                    right_address = (left_address + drift) % P
                    direct = fourier(
                        tuple(
                            safe(left_chamber, left_address + tau)
                            * safe(right_chamber, right_address + tau)
                            for tau in range(P)
                        ),
                        frequency,
                    )
                    expected = pow(XI, frequency * left_address, FIELD_PRIME) * kernel % FIELD_PRIME
                    require(direct == expected, "joint guard drift covariance failed")
                    pair_factorization_checks += 1


total_pair_modes = len(pair_kernels)
nonzero_pair_modes = sum(value != 0 for value in pair_kernels.values())
primitive_pair_modes = sum(
    value != 0
    for key, value in pair_kernels.items()
    if key[3] != 0
)
primitive_total = 3 * 3 * P * (P - 1)
frequency_one_nonzero = sum(
    pair_kernels[(left, right, drift, 1)] != 0
    for left in CHAMBERS
    for right in CHAMBERS
    for drift in range(P)
)


# The relative drift is genuinely joint.  Uniform diagonal and shifted
# couplings have identical sheet marginals, but different middle-chamber
# guard overlap.  This is a 13-state analogue of the 2x2 checkerboard debt.
middle_overlap_by_drift = tuple(
    sum(safe("middle", a) * safe("middle", a + drift) for a in range(P))
    for drift in range(P)
)
require(middle_overlap_by_drift[0] == 9, "diagonal middle overlap changed")
require(len(set(middle_overlap_by_drift)) > 1, "relative drift became invisible")

# Summing the common address before twisting annihilates every nonzero
# frequency, regardless of chamber pair and drift.
for left_chamber in CHAMBERS:
    for right_chamber in CHAMBERS:
        for drift in range(P):
            for frequency in range(1, P):
                total = sum(
                    pow(XI, frequency * address, FIELD_PRIME)
                    * pair_kernels[(left_chamber, right_chamber, drift, frequency)]
                    for address in range(P)
                ) % FIELD_PRIME
                require(total == 0, "early sheet marginalization retained a primitive mode")


ledger_rows = (
    bridge_hash,
    tuple((name, tuple(sorted(danger))) for name, (_residue, danger) in CHAMBERS.items()),
    tuple(sorted(single_kernels.items())),
    tuple(sorted(pair_kernels.items())),
    middle_overlap_by_drift,
)
ledger_hash = hashlib.sha256(repr(ledger_rows).encode("ascii")).hexdigest()
if EXPECTED_LEDGER_SHA256 != "TO_BE_PINNED":
    require(ledger_hash == EXPECTED_LEDGER_SHA256, "semantic ledger changed")


print("== U_full H-guard sheet and joint-drift address law ==")
print("common 91-circle partition: 13 sheets x 3 half-open residual chambers = 39 atoms")
print("boundary convention: r in [0,1), [1,6), [6,7); r=1 middle and r=6 right: PASS")
print("guard danger-sheet arcs: left=(12,0,1), middle=(11,12,0,1), right=(11,12,0)")
print("safe counts by chamber=(10,9,10); every 13 Fourier mode is nonzero")
print(f"single factorization checks={3*P*P}; joint factorization checks={pair_factorization_checks}")
print(
    f"joint kernels nonzero={nonzero_pair_modes}/{total_pair_modes}; "
    f"primitive nonzero={primitive_pair_modes}/{primitive_total}; frequency-one={frequency_one_nonzero}/{3*3*P}"
)
print(f"middle-chamber overlap by drift={middle_overlap_by_drift}")
print("same sheet marginals can carry different drift profiles: PASS_HOSTILE")
print("summing the common sheet before twisting kills every nonzero tau mode: PASS_HOSTILE")
print(f"semantic ledger sha256={ledger_hash}")
print("scope: actual U_full guard factor only; no common ancestry, endpoint coupling, current, row exclusion, or LRC(14)")
print("all exact checks passed")
