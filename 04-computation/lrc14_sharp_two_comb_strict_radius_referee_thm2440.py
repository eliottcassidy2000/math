#!/usr/bin/env python3
"""Independent strict/closed exact referee for THM-2440.

The script reconstructs the centred component of two radius-1/14
integer combs with exact Fraction endpoints.  It distinguishes:

* the closed (equivalently, open almost-everywhere) component; and
* the literal connected component of the union of the open combs.

It also checks the finite branches left by the analytic fragmentation
bound, the endpoint-seam witness, scaled equality controls, and a
deliberately redundant coprime scan.  All truth-bearing checks use
explicit exceptions and therefore survive ``python3 -O``.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HALF = Fraction(1, 2)
LAMBDA = Fraction(1, 14)
RADIUS_AE = Fraction(15, 182)
RADIUS_OPEN = Fraction(15, 196)


def positive_teeth(speed: int) -> list[tuple[Fraction, Fraction]]:
    """Closed endpoint hulls of the teeth meeting [0,1/2]."""

    require(speed > 0, "speed must be positive")
    teeth: list[tuple[Fraction, Fraction]] = []
    for centre in range(speed + 1):
        left = Fraction(14 * centre - 1, 14 * speed)
        right = Fraction(14 * centre + 1, 14 * speed)
        if right < 0 or left > HALF:
            continue
        teeth.append((max(Fraction(0), left), min(HALF, right)))
    return teeth


def centred_radius(first: int, second: int, *, strict: bool) -> Fraction:
    """Radius of the centred component for two integer danger combs.

    For ``strict=False``, touching closed teeth join.  For ``strict=True``,
    two noncentral open teeth join only when their endpoint hulls overlap
    strictly.  The special first tooth contains zero in the actual open
    circle comb, so it is allowed to start the merge.
    """

    intervals = sorted(positive_teeth(first) + positive_teeth(second))
    radius = Fraction(0)
    for left, right in intervals:
        separated = left > radius
        strict_seam = strict and radius > 0 and left >= radius
        if separated or strict_seam:
            break
        radius = max(radius, right)
    return radius


def distance_to_integer(value: Fraction) -> Fraction:
    floor = value.numerator // value.denominator
    residue = value - floor
    return min(residue, 1 - residue)


# ---------------------------------------------------------------------------
# 1. The strict seam that refutes the reserved literal statement.
# ---------------------------------------------------------------------------

seam = Fraction(1, 14)
require(distance_to_integer(seam) == LAMBDA, "first comb seam changed")
require(distance_to_integer(13 * seam) == LAMBDA, "13-comb seam changed")
require(distance_to_integer(14 * seam) == 0, "14-comb does not bridge seam")
require(centred_radius(1, 13, strict=False) == RADIUS_AE,
        "closed {1,13} radius mismatch")
require(centred_radius(1, 13, strict=True) == seam,
        "open {1,13} must stop at the shared seam")
require(centred_radius(1, 14, strict=False) == RADIUS_OPEN,
        "closed {1,14} radius mismatch")
require(centred_radius(1, 14, strict=True) == RADIUS_OPEN,
        "open {1,14} radius mismatch")


# ---------------------------------------------------------------------------
# 2. Arithmetic left by the sharp fragmentation inequality.
#
# If [0,R] is covered a.e. by D_p union D_q, then
#
#   R <= 2R/7 + 6/(49p) + 6/(49q),
#
# hence 1/p+1/q >= 35R/6.
# ---------------------------------------------------------------------------

threshold_ae = Fraction(35, 6) * RADIUS_AE
threshold_open = Fraction(35, 6) * RADIUS_OPEN
require(threshold_ae == Fraction(175, 364), "wrong a.e. threshold")
require(threshold_open == Fraction(175, 392), "wrong open threshold")

# Closed/a.e. target: p>=5 and p=4 are impossible by reciprocals;
# for p=3 only q=4,5 remain.
require(Fraction(2, 5) < threshold_ae, "p>=5 reduction failed")
require(Fraction(1, 4) + Fraction(1, 5) < threshold_ae,
        "p=4 reduction failed")
require(Fraction(1, 3) + Fraction(1, 7) < threshold_ae,
        "p=3 tail reduction failed")
for q in (4, 5):
    require(centred_radius(3, q, strict=False) == Fraction(1, 42),
            f"closed p=3 candidate q={q} changed")

# For p=2 the interval after its central tooth has length 17/364.
# One q-tooth can cover it only for q<=3, and q=3 misses it.
tail_p2_ae = RADIUS_AE - Fraction(1, 28)
require(tail_p2_ae == Fraction(17, 364), "wrong p=2 a.e. tail")
require(Fraction(1, 7 * 4) < tail_p2_ae <= Fraction(1, 7 * 3),
        "wrong p=2 a.e. tooth cutoff")
require(centred_radius(2, 3, strict=False) == Fraction(1, 28),
        "closed {2,3} control changed")

# For p=1 the interval after 1/14 has length 1/91.  Speeds q<=12
# do not reach the seam; speeds q>=14 have teeth too short.  q=13
# is the unique exact closed/a.e. handoff.
tail_p1_ae = RADIUS_AE - seam
require(tail_p1_ae == Fraction(1, 91), "wrong p=1 a.e. tail")
for q in range(1, 13):
    require(centred_radius(1, q, strict=False) == seam,
            f"q={q} unexpectedly reaches past the closed seam")
require(Fraction(13, 14 * 13) == seam, "q=13 left endpoint mismatch")
require(Fraction(15, 14 * 13) == RADIUS_AE,
        "q=13 right endpoint mismatch")
require(Fraction(1, 7 * 14) < tail_p1_ae,
        "q>=14 tooth-length exclusion failed")

# Literal-open target: p>=5 is impossible; p=4 leaves q=5; p=3
# leaves q=4,5,7,8.  All direct candidates stop at their central tooth.
require(Fraction(2, 5) < threshold_open, "strict p>=5 reduction failed")
require(Fraction(1, 4) + Fraction(1, 6) < threshold_open,
        "strict p=4 tail reduction failed")
require(Fraction(1, 3) + Fraction(1, 9) < threshold_open,
        "strict p=3 tail reduction failed")
require(centred_radius(4, 5, strict=True) == Fraction(1, 56),
        "strict {4,5} control changed")
for q in (4, 5, 7, 8):
    require(centred_radius(3, q, strict=True) == Fraction(1, 42),
            f"strict p=3 candidate q={q} changed")

# For p=2 the remaining interval has length 2/49, so q<=3; q=3
# again misses it.
tail_p2_open = RADIUS_OPEN - Fraction(1, 28)
require(tail_p2_open == Fraction(2, 49), "wrong p=2 strict tail")
require(Fraction(1, 7 * 4) < tail_p2_open <= Fraction(1, 7 * 3),
        "wrong p=2 strict tooth cutoff")
require(centred_radius(2, 3, strict=True) == Fraction(1, 28),
        "strict {2,3} control changed")

# For p=1 the remaining interval has length 1/196, hence q<=28.
# Literal coverage of x=1/14 requires q divisible by 14.  The only
# candidates are q=14,28; q=28 stops one cell too early.
tail_p1_open = RADIUS_OPEN - seam
require(tail_p1_open == Fraction(1, 196), "wrong p=1 strict tail")
require(Fraction(1, 7 * 29) < tail_p1_open <= Fraction(1, 7 * 28),
        "wrong p=1 strict tooth cutoff")
seam_bridgers = [
    q for q in range(1, 29)
    if distance_to_integer(Fraction(q, 14)) < LAMBDA
]
require(seam_bridgers == [14, 28], "wrong strict seam-bridger list")
require(centred_radius(1, 28, strict=True) == Fraction(29, 392),
        "strict {1,28} control changed")
require(Fraction(29, 392) < RADIUS_OPEN,
        "q=28 should stop before the sharp strict radius")


# ---------------------------------------------------------------------------
# 3. Redundant exact scan and scale covariance controls.
# ---------------------------------------------------------------------------

closed_rows: list[tuple[Fraction, int, int]] = []
strict_rows: list[tuple[Fraction, int, int]] = []
for p in range(1, 41):
    for q in range(p, 401):
        if gcd(p, q) != 1:
            continue
        closed_rows.append((centred_radius(p, q, strict=False), p, q))
        strict_rows.append((centred_radius(p, q, strict=True), p, q))

closed_rows.sort(reverse=True)
strict_rows.sort(reverse=True)
require(closed_rows[0] == (RADIUS_AE, 1, 13),
        "bounded closed scan has wrong unique maximum")
require(strict_rows[0] == (RADIUS_OPEN, 1, 14),
        "bounded strict scan has wrong unique maximum")

for scale in range(1, 51):
    require(
        scale * centred_radius(scale, 13 * scale, strict=False) == RADIUS_AE,
        f"closed scale covariance failed at {scale}",
    )
    require(
        scale * centred_radius(scale, 14 * scale, strict=True) == RADIUS_OPEN,
        f"strict scale covariance failed at {scale}",
    )


def row_text(row: tuple[Fraction, int, int]) -> str:
    radius, p, q = row
    return f"{{{p},{q}}}:{radius}"


print("THM-2440 independent strict/closed two-comb referee")
print(f"fragmentation thresholds: closed/ae={threshold_ae} strict={threshold_open}")
print(
    "reserved strict witness: "
    f"x={seam}, dist(x,Z)={distance_to_integer(seam)}, "
    f"dist(13x,Z)={distance_to_integer(13 * seam)}"
)
print(
    "sharp normalized radii: "
    f"closed/ae {{1,13}}={RADIUS_AE}; "
    f"strict {{1,14}}={RADIUS_OPEN}"
)
print(
    "cross-convention controls: "
    f"strict {{1,13}}={centred_radius(1, 13, strict=True)}; "
    f"closed {{1,14}}={centred_radius(1, 14, strict=False)}"
)
print(
    "finite branch controls: "
    f"p2_ae_tail={tail_p2_ae}, p1_ae_tail={tail_p1_ae}, "
    f"p2_strict_tail={tail_p2_open}, p1_strict_tail={tail_p1_open}"
)
print(f"strict seam bridgers q<=28: {seam_bridgers}")
print(
    f"redundant coprime scan rows: closed={len(closed_rows)} "
    f"strict={len(strict_rows)}"
)
print("top closed rows: " + ", ".join(row_text(row) for row in closed_rows[:5]))
print("top strict rows: " + ", ".join(row_text(row) for row in strict_rows[:5]))
print("scaled equality controls: n=1..50 PASS")
print("RESULT: PASS")
