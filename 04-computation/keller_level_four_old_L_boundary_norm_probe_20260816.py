#!/usr/bin/env python3
"""Exact old-boundary valuation gate for the fourth Keller norm rung.

This companion reconstructs THM-3495's primitive 66,146-term polynomial J,
extracts its Newton face for the two inverse sheets diverging over L=0, and
checks a finite-sheet hostile value.  It proves only

    v_L(Norm(J)) = -43.

It does not construct Norm(J), identify its numerator, prove a level-four
image divisor, or certify a degree-81 eliminant.
"""

from __future__ import annotations

import contextlib
import io
import runpy
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
GLOBAL_PROBE = ROOT / "04-computation/keller_level_three_global_norm_probe_20260816.py"

# Reconstruct J through the already frozen exact route, but suppress its
# progress report so this companion has a small deterministic output.
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(str(GLOBAL_PROBE))
require(
    captured.getvalue().rstrip().endswith("all exact checks passed"),
    "the imported global J reconstruction did not finish cleanly",
)

J = namespace["J"]
a, b, c = namespace["a"], namespace["b"], namespace["c"]
ZERO = namespace["ZERO"]
J_LEDGER_SHA256 = "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe"
require(namespace["coefficient_hash"](J) == J_LEDGER_SHA256, "global J ledger changed")

# On a divergent root of L*w^3+T*w-2c, with u=1/w, the inverse chart has
# x=w~u^-1, y~D/S, z~-3(D/S)u.  Thus the relevant monomial weight on J is
# i-k for x^i y^j z^k.
terms = {
    (int(i), int(j), int(k)): int(coefficient)
    for (i, j, k), coefficient in J.to_dict().items()
}
face_weight = max(i - k for i, _j, k in terms)
face = ZERO
face_terms = 0
for (i, j, k), coefficient in terms.items():
    if i - k == face_weight:
        face += coefficient * a**i * b**j * c**k
        face_terms += 1

FACE_CONSTANT = -(2**58) * (3**51) * (13**8) * (79**4) * (313**2)
expected_face = FACE_CONSTANT * a**43 * (3 * a * c - 2 * b) ** 15
require(face_weight == 43, "the J Newton-face weight changed")
require(face_terms == 16, "the J Newton-face support changed")
require(face == expected_face, "the J Newton face no longer has the pinned factorization")

# The leading divergent relation x*z=-3y makes 3*x*z-2*y=-11y, a unit at
# the generic point of L.  The exact reduced coefficient is retained as a
# compact hostile check against accidental face cancellation.
reduced_face_coefficient = FACE_CONSTANT * (-11) ** 15
require(reduced_face_coefficient > 0, "reduced face sign changed")

# The same target used in THM-3495 has one finite inverse sheet
# q=(2,5/6,-7/8).  Evaluate J there exactly to prove that sheet is a unit.
finite_q = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
finite_value = sum(
    Fraction(coefficient)
    * finite_q[0] ** i
    * finite_q[1] ** j
    * finite_q[2] ** k
    for (i, j, k), coefficient in terms.items()
)
require(finite_value != 0, "J vanished on the finite hostile inverse sheet")

# Simultaneous generic-DVR unit witness at the target (2/27,1,1).
ta, tb, tc = Fraction(2, 27), Fraction(1), Fraction(1)
L_value = 27 * ta**2 * tc**2 - 18 * ta * tb * tc + 16 * ta + tb**3 * tc - tb**2
T_value = 4 - 3 * tb * tc
S_value = 27 * ta * tc**2 - 9 * tb * tc + 8
D_value = 18 * ta * tc - 3 * tb**2 * tc + 2 * tb
require(
    (L_value, tc, T_value, S_value, D_value)
    == (0, 1, 1, 1, Fraction(1, 3)),
    "generic-DVR unit witness changed",
)

print("== exact level-four old-L boundary norm gate ==")
print(f"J coefficient-ledger sha256={J_LEDGER_SHA256}")
print(f"J face: max(i-k)={face_weight}; terms={face_terms}")
print(
    "in(J)=-2^58*3^51*13^8*79^4*313^2"
    "*x^43*(3*x*z-2*y)^15"
)
print(f"divergent reduction: coefficient={reduced_face_coefficient}; monomial=x^43*y^15")
print("divergent sheets: v_L(J(q))=-43/2 each")
print(f"finite sheet: q={finite_q}; J(q)={finite_value} != 0")
print("therefore v_L(Norm(J))=-43")
print(
    "scope: old-L boundary valuation only; Norm(J), its numerator/image, "
    "degree-81 separability, and level-four JC consequences remain open"
)
print("all exact checks passed")
