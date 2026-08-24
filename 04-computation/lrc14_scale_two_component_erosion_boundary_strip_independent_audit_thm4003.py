#!/usr/bin/env python3
"""Independent direct-grid audit for THM-4003.

Unlike the primary certificate, this path constructs literal Fraction wall
sets before using the closed residue formulas in the bounded census.  It also
checks the two-flank and top-balance threshold arithmetic.  LRC(14) remains
open.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def positive_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    require(residue != 0, "odd residue is nonzero")
    return residue


def formula_gaps(t: int, u: int) -> tuple[F, F, F, F]:
    common = gcd(t, u)
    return (
        F(positive_residue(3 * t - 4 * u, 42 * common), 42 * u),
        F(positive_residue(9 * t + 16 * u, 126 * common), 126 * u),
        F(positive_residue(9 * t - 110 * u, 126 * common), 126 * u),
        F(positive_residue(3 * t + 38 * u, 42 * common), 42 * u),
    )


def literal_gaps(t: int, u: int) -> tuple[F, F, F, F]:
    entering = {
        (F(t * (14 * index + 1), 14 * u)) % 1 for index in range(u)
    }
    exiting = {
        (F(t * (14 * index - 1), 14 * u)) % 1 for index in range(u)
    }
    boundaries = (F(2, 21), F(8, 63), F(55, 63), F(19, 21))
    return (
        min((event - boundaries[0]) % 1 for event in entering),
        min((boundaries[1] - event) % 1 for event in exiting),
        min((event - boundaries[2]) % 1 for event in entering),
        min((boundaries[3] - event) % 1 for event in exiting),
    )


literal_rows = 0
for t in range(11, 102, 2):
    for u in range(1, t + 1):
        direct = literal_gaps(t, u)
        formula = formula_gaps(t, u)
        require(direct == formula, f"literal directed grid t={t},u={u}")
        require(direct[0] == direct[3] and direct[1] == direct[2], "reflection")
        literal_rows += 1

correct = literal_gaps(11, 1)
wrong = (
    min((F(2, 21) - event) % 1 for event in {(F(11 * (14 * a + 1), 14)) % 1 for a in range(1)}),
    min((event - F(8, 63)) % 1 for event in {(F(11 * (14 * a - 1), 14)) % 1 for a in range(1)}),
    min((F(55, 63) - event) % 1 for event in {(F(11 * (14 * a + 1), 14)) % 1 for a in range(1)}),
    min((event - F(19, 21)) % 1 for event in {(F(11 * (14 * a - 1), 14)) % 1 for a in range(1)}),
)
require(correct == (F(29, 42), F(115, 126), F(115, 126), F(29, 42)), "sign control")
require(wrong == (F(13, 42), F(11, 126), F(11, 126), F(13, 42)), "wrong-way hostile")

require((3 * 12 - 4 * 9) % (42 * gcd(12, 9)) == 0, "even entering collision")
require((3 * 12 + 38 * 9) % (42 * gcd(12, 9)) == 0, "even exiting collision")

# Independently check the quadratic odd rounding and top-balance thresholds.
rounding_rows = 0
for t in range(11, 1002, 2):
    for U in range(max(11, (3 * t) // 4 + 1), t + 1):
        gate = 3 * t * (2 * U - 1) <= 8 * (U - 1) ** 2
        rounded = U >= (3 * t) // 4 + 2 + (t % 4 == 1)
        require(gate == rounded, f"quadratic rounding t={t},U={U}")
        rounding_rows += 1

def top_bound(ratio: F) -> F:
    if ratio < 22:
        return max(F(0), (3 * ratio - 11) / 154)
    return F(5, 14)


require(top_bound(F(11)) == F(1, 7), "all-type top-balance equality")
require(top_bound(F(1001, 189)) == F(2, 63), "scale-two top-balance equality")
require(top_bound(F(22)) == F(5, 14), "truncation equality")

# Cumulative owner minima and the strengthened body component, organized
# independently from the primary certificate.
cells = 0
symbolic_closed = []
extra_closed = []
survivors = []
old_residue_extra = set()
two_flank_extra = set()
combined_extra = set()
for t in range(11, 1002, 2):
    minima: list[F | None] = [None, None, None, None]
    for U in range(1, t + 1):
        row = formula_gaps(t, U)
        for index in range(4):
            minima[index] = row[index] if minima[index] is None else min(minima[index], row[index])
        if U < max(11, (3 * t) // 4 + 1):
            continue
        values = tuple(value for value in minima if value is not None)
        require(values[0] == values[3] and values[1] == values[2], "minimum reflection")
        residual = max(F(0), F(2, 63) - values[0] - values[1])
        body_image = t * F(2 * U - 1, 84 * U * (U - 1))
        old_simple = 3 * t > 4 * (U - 1)
        old_residue = F(t, 42 * U) > residual
        symbolic = 3 * t * (2 * U - 1) > 8 * (U - 1) ** 2
        exact = body_image > residual
        require(not symbolic or exact, "exact gate contains symbolic gate")
        if symbolic:
            symbolic_closed.append((t, U))
        elif exact:
            extra_closed.append((t, U))
        else:
            survivors.append((t, U))
        if not old_simple and old_residue:
            old_residue_extra.add((t, U))
        if not old_simple and symbolic:
            two_flank_extra.add((t, U))
        if not old_simple and exact:
            combined_extra.add((t, U))
        cells += 1

require((cells, len(symbolic_closed), len(extra_closed), len(survivors)) == (62989, 742, 77, 62170), "finite partition")
require(extra_closed[0] == (11, 11) and extra_closed[-1] == (185, 141), "extra endpoints")
separate_union = old_residue_extra | two_flank_extra
synergy = combined_extra - separate_union
require(
    (
        len(old_residue_extra),
        len(two_flank_extra),
        len(old_residue_extra & two_flank_extra),
        len(separate_union),
        len(combined_extra),
        len(synergy),
    )
    == (92, 248, 34, 306, 325, 19),
    "synergy partition",
)
serialization = ";".join(f"{t},{U}" for t, U in extra_closed).encode("ascii")
pair_hash = sha256(serialization).hexdigest()

print("THM4003_INDEPENDENT_DIRECT_GRID_AUDIT")
print("scope=conditional_scale2_(2,1,9);LRC14=OPEN")
print(f"literal_fraction_rows={literal_rows}")
print(f"quadratic_rounding_rows={rounding_rows}")
print("reflection_identity=gap1=gap4_and_gap2=gap3")
print("top_balance_equalities=1/7_at_11;2/63_at_1001/189;5/14_at_22")
print(f"finite_cells={cells};symbolic_closed={len(symbolic_closed)};exact_extra={len(extra_closed)};survivors={len(survivors)}")
print("synergy=old_residue_92;two_flank_248;overlap_34;union_306;combined_325;new_19")
print(f"synergy_pairs={tuple(sorted(synergy))}")
print(f"exact_extra_pair_sha256={pair_hash}")
print(f"exact_extra_pairs={tuple(extra_closed)}")
print(f"checks={CHECKS}")
print("RESULT=PASS")
