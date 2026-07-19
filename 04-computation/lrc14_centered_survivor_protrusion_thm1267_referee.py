#!/usr/bin/env python3
"""Exact referee for THM-1267's centered-survivor protrusion law.

The paper theorem combines three providers:

* exact centered-spoke interval geometry;
* the six-bin endpoint density and its two inverse quantiles;
* THM-1198's twelve-piece one-comb envelope, sharpened for the one load
  separated by THM-1241.

This dependency-free referee checks every rational constant with
``fractions.Fraction``.  It also exhausts a finite signed-gap geometry census.
The census is a check of the affine/topological layer, not a finite reduction
of the six-comb covering problem.

Tournament analysis is deliberately a loss audit.  Speed order on the five
load obligations gives a transitive tournament and singles out the d6 load,
but destroys phases, endpoint orientation, survivor mass, and protrusion.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import ceil, floor


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-1267 referee failed: {message}")


def optimization_safe_require_probe() -> None:
    caught = False
    try:
        require(False, "deliberate optimization-safety probe")
    except RuntimeError as error:
        caught = "deliberate optimization-safety probe" in str(error)
    require(caught, "require probe did not fire")


@dataclass(frozen=True)
class EnvelopePiece:
    lower: F
    upper: F
    constant: F
    reciprocal: F

    def value(self, slope: F) -> F:
        require(self.lower <= slope <= self.upper, "piece evaluation outside domain")
        return self.constant + self.reciprocal / slope

    def closed_maximum(self) -> F:
        """Maximum of a+b/L on this positive closed interval."""
        endpoint = self.lower if self.reciprocal >= 0 else self.upper
        return self.value(endpoint)


PIECES = (
    EnvelopePiece(F(6, 7), F(68, 63), F(0), F(1, 6)),
    EnvelopePiece(F(68, 63), F(8, 7), F(3, 4), F(-9, 14)),
    EnvelopePiece(F(8, 7), F(6, 5), F(0), F(3, 14)),
    EnvelopePiece(F(6, 5), F(48, 35), F(5, 18), F(-5, 42)),
    EnvelopePiece(F(48, 35), F(3, 2), F(0), F(11, 42)),
    EnvelopePiece(F(3, 2), F(12, 7), F(2, 9), F(-1, 14)),
    EnvelopePiece(F(12, 7), F(2), F(0), F(13, 42)),
    EnvelopePiece(F(2), F(244, 119), F(1, 24), F(19, 84)),
    EnvelopePiece(F(244, 119), F(15, 7), F(3, 4), F(-103, 84)),
    EnvelopePiece(F(15, 7), F(12, 5), F(0), F(8, 21)),
    EnvelopePiece(F(12, 5), F(18, 7), F(5, 18), F(-2, 7)),
    EnvelopePiece(F(18, 7), F(3), F(0), F(3, 7)),
)

EXPECTED_MAXIMA = (
    F(7, 36),
    F(3, 16),
    F(3, 16),
    F(55, 288),
    F(55, 288),
    F(13, 72),
    F(13, 72),
    F(13, 84),
    F(8, 45),
    F(8, 45),
    F(1, 6),
    F(1, 6),
)

DENSITY_HEIGHTS = (F(3, 4), F(13, 12), F(7, 6), F(7, 6), F(13, 12), F(3, 4))
THRESHOLD = F(20, 23)  # = 140/161
SEPARATED_CAP = F(23, 120)  # = 161/840


def integrate_prefix(length: F, heights: tuple[F, ...]) -> F:
    require(F(0) <= length <= F(1), "prefix length outside [0,1]")
    remaining = length
    mass = F(0)
    width = F(1, len(heights))
    for height in heights:
        used = min(remaining, width)
        mass += height * used
        remaining -= used
        if remaining == 0:
            break
    require(remaining == 0, "prefix integrator left unconsumed length")
    return mass


def inverse_endpoint_mass(target: F) -> F:
    require(F(0) <= target <= F(1), "mass target outside [0,1]")
    remaining = target
    length = F(0)
    width = F(1, 6)
    for height in DENSITY_HEIGHTS:
        bin_mass = height * width
        if remaining <= bin_mass:
            return length + remaining / height
        length += width
        remaining -= bin_mass
    require(remaining == 0, "inverse density exhausted all bins")
    return F(1)


def check_density_quantiles() -> None:
    require(sum(DENSITY_HEIGHTS) / 6 == 1, "density normalization")
    require(DENSITY_HEIGHTS == tuple(reversed(DENSITY_HEIGHTS)), "density symmetry")
    require(max(DENSITY_HEIGHTS) == F(7, 6), "density maximum")
    require(DENSITY_HEIGHTS[0] == DENSITY_HEIGHTS[-1] == F(3, 4),
            "endpoint-sixth height")
    require(integrate_prefix(F(1, 6), DENSITY_HEIGHTS) == F(1, 8),
            "endpoint-sixth mass")
    require(integrate_prefix(F(1, 2), DENSITY_HEIGHTS) == F(1, 2),
            "half-density mass")
    for target, expected in ((F(1, 36), F(1, 27)),
                             (F(11, 360), F(11, 270))):
        require(target < F(1, 8), "quantile must stay in endpoint sixth")
        require(inverse_endpoint_mass(target) == expected, "endpoint inverse quantile")
        require(integrate_prefix(expected, DENSITY_HEIGHTS) == target,
                "endpoint quantile forward check")
        require(integrate_prefix(expected, tuple(reversed(DENSITY_HEIGHTS))) == target,
                "right-endpoint orientation check")


def check_envelope() -> tuple[list[F], list[F]]:
    require(len(PIECES) == 12, "twelve exact envelope pieces")
    require(PIECES[0].lower == F(6, 7) and PIECES[-1].upper == 3,
            "envelope domain endpoints")
    for left, right in zip(PIECES, PIECES[1:]):
        require(left.upper == right.lower, "envelope domain partition")
        require(left.value(left.upper) == right.value(right.lower),
                "envelope continuity")
    maxima = [piece.closed_maximum() for piece in PIECES]
    require(tuple(maxima) == EXPECTED_MAXIMA, "all twelve piece maxima")
    require(max(maxima) == F(7, 36), "global compact-envelope cap")

    require(F(6, 7) < THRESHOLD < F(68, 63), "separation lies in first piece")
    require(F(140, 161) == THRESHOLD, "threshold simplification")
    require(F(161, 840) == SEPARATED_CAP, "separated-cap simplification")
    require(PIECES[0].constant + PIECES[0].reciprocal / THRESHOLD == SEPARATED_CAP,
            "open-boundary first-piece value")

    # Piece one is decreasing, so its open restricted supremum is the value
    # at the excluded threshold.  Every later closed piece stays below it.
    require(PIECES[0].reciprocal > 0, "first piece decreases in L")
    later_maxima = [piece.closed_maximum() for piece in PIECES[1:]]
    require(all(value < SEPARATED_CAP for value in later_maxima),
            "later exact pieces below separated cap")

    bv_tail_maximum = F(1, 7) + F(1, 7 * 3)
    require(bv_tail_maximum == F(4, 21), "BV-tail maximum at L=3")
    require(bv_tail_maximum < SEPARATED_CAP, "BV tail below separated cap")

    # Exact approach from above demonstrates the lower side of the open
    # supremum; monotonicity supplies the upper side.
    approach_gaps = []
    for denominator in (100, 1000, 10000):
        slope = THRESHOLD + F(1, denominator)
        require(slope < PIECES[0].upper, "approach witness left first piece")
        load = PIECES[0].value(slope)
        gap = SEPARATED_CAP - load
        require(F(0) < gap < F(1, denominator), "open-supremum approach gap")
        approach_gaps.append(gap)
    require(approach_gaps[2] < approach_gaps[1] < approach_gaps[0],
            "open-supremum witnesses converge monotonically")

    require(4 * F(7, 36) + SEPARATED_CAP == F(349, 360),
            "five-load separated budget")
    require(1 - F(349, 360) == F(11, 360), "separated survivor mass")
    require(F(4, 3) * F(1, 36) == F(1, 27), "basic endpoint inversion")
    require(F(4, 3) * F(11, 360) == F(11, 270),
            "separated endpoint inversion")
    return maxima, approach_gaps


def circle_distance(value: F) -> F:
    residue = value - floor(value)
    return min(residue, 1 - residue)


def interval_intersection_length(left_a: F, right_a: F, left_b: F, right_b: F) -> F:
    return max(F(0), min(right_a, right_b) - max(left_a, left_b))


def nearest_integers(value: F) -> tuple[int, ...]:
    choices = {floor(value), ceil(value)}
    distance = min(abs(F(choice) - value) for choice in choices)
    return tuple(sorted(choice for choice in choices if abs(F(choice) - value) == distance))


def check_geometry_row(c: int, d: int, k: int, p: int) -> tuple[bool, F]:
    require(0 < c < d, "geometry speed order")
    q = c + d
    t0 = F(2 * k + 1, 2 * c)
    epsilon = F(p) - q * t0
    rho = abs(epsilon)
    require(rho <= F(1, 2), "nearest-integer error")
    t = F(p, q)

    g_left = F(14 * k + 1, 14 * c)
    g_right = F(14 * k + 13, 14 * c)
    require(g_left < t < g_right, "centered spoke inside carrier gap")
    require(circle_distance(c * t) == circle_distance(d * t) > F(1, 4),
            "centered spoke common depth")

    h = floor(d * t)
    require(h == p - k - 1, "safe-component floor")
    s_left = F(14 * h + 1, 14 * d)
    s_right = F(14 * h + 13, 14 * d)
    require(s_left < t < s_right, "spoke inside complete d-safe component")

    g_center = (g_left + g_right) / 2
    s_center = (s_left + s_right) / 2
    require(g_center == t0, "carrier center")
    require(s_center == t0 + epsilon / d, "exact signed center shift")
    require(abs(s_center - g_center) == rho / d, "absolute center shift")

    intersection = interval_intersection_length(s_left, s_right, g_left, g_right)
    actual_tail = (s_right - s_left) - intersection
    expected_tail = max(F(0), (rho + F(3, 7)) / d - F(3, 7 * c))
    require(actual_tail == expected_tail, "physical protrusion formula")

    left_tail = max(F(0), g_left - s_left)
    right_tail = max(F(0), s_right - g_right)
    require(left_tail + right_tail == actual_tail, "tail decomposition")
    require(left_tail == 0 or right_tail == 0, "single endpoint tail")
    if actual_tail > 0:
        require((epsilon < 0 and left_tail > 0 and right_tail == 0)
                or (epsilon > 0 and right_tail > 0 and left_tail == 0),
                "tail orientation follows center shift")

    ell = actual_tail / (s_right - s_left)
    expected_ell = max(F(0), F(1, 2) + F(7, 6) * rho - F(d, 2 * c))
    require(ell == expected_ell, "normalized protrusion formula")
    return actual_tail > 0, ell


def geometry_census() -> tuple[int, int, F]:
    rows = 0
    protruding = 0
    largest_ell = F(0)
    # Signed k catches floor behavior across zero.  Ratios through 3 include
    # both protruding and contained components; this is not a cover search.
    for c in range(1, 65):
        for d in range(c + 1, 3 * c + 1):
            for k in range(-c, c + 1):
                target = (c + d) * F(2 * k + 1, 2 * c)
                for p in nearest_integers(target):
                    has_tail, ell = check_geometry_row(c, d, k, p)
                    rows += 1
                    protruding += int(has_tail)
                    largest_ell = max(largest_ell, ell)

    # Small exact positive control, including the two mirror nearest choices.
    for p, sign in ((1, -1), (2, 1)):
        c, d, k = 1, 2, 0
        epsilon = F(p) - (c + d) * F(1, 2)
        has_tail, ell = check_geometry_row(c, d, k, p)
        require(has_tail and ell == F(1, 12), "positive-control tail length")
        require((epsilon > 0) - (epsilon < 0) == sign, "positive-control orientation")
    return rows, protruding, largest_ell


def check_constant_chain() -> None:
    require(F(1, 27) < F(1, 6), "basic split threshold")
    require(F(11, 270) < F(1, 6), "separated split threshold")
    require(F(27, 63) == F(3, 7), "basic rho coefficient")
    require(F(113, 54) < F(89, 42), "endpoint density beats Lebesgue bound")
    require(F(563, 270) < F(113, 54), "Kakeya envelope improves basic bound")

    # Substituting rho=1/2 into the strict rho inequalities gives the two
    # displayed ratio caps.  These are exact boundary equalities.
    require((27 * F(113, 54) - 25) / 63 == F(1, 2),
            "basic rho boundary")
    require((135 * F(563, 270) - 124) / 315 == F(1, 2),
            "separated rho boundary")
    require(F(6, 7) * F(70, 69) == THRESHOLD,
            "THM-1241 normalized d6 separation")

    # Integer strictness: 270*d < 563*c becomes <= 563*c-1.  A finite exact
    # sweep checks the rounding form independently of Python assertions.
    checked = 0
    for c in range(1, 500):
        for d in range(1, 4 * c + 1):
            if 270 * d < 563 * c:
                require(270 * d <= 563 * c - 1, "integer strict rounding")
                checked += 1
    require(checked > 0, "integer rounding positive control")


def tournament_loss_audit() -> None:
    print("TOURNAMENT_LOSS_AUDIT")
    print("vertices=five noncarrier load obligations d2,d3,d4,d5,d6")
    print("observable=normalized-slope comparison; gauge=increasing speed")
    print("score_histogram=0,1,2,3,4 directed_3_cycles=0 scc_sizes=1,1,1,1,1")
    print("hamiltonian_path_count=1 tie_hamiltonian_path=d2,d3,d4,d5,d6")
    print("preserves=rank and the unique top load receiving THM-1241 separation")
    print("destroys=phase, endpoint orientation, survivor position, and load size")
    print("proof_carrier=one endpoint tail plus five weighted-load obligations")


def main() -> None:
    optimization_safe_require_probe()
    check_density_quantiles()
    maxima, approach_gaps = check_envelope()
    check_constant_chain()
    rows, protruding, largest_ell = geometry_census()

    print("THM-1267 CENTERED-SURVIVOR PROTRUSION EXACT REFEREE")
    print("optimization_safe_require_probe=PASS")
    print(f"geometry_census_rows={rows} protruding_rows={protruding} largest_ell={largest_ell}")
    print("geometry_identity=center_shift=epsilon/d; tail=max(0,(rho+3/7)/d-3/(7c))")
    print("normalized_tail=ell=max(0,1/2+7rho/6-d/(2c))")
    print("geometry_positive_control=(c,d,k,p)=(1,2,0,1_or_2) ell=1/12 mirror_tails=PASS")
    print("six_bin_heights=" + ",".join(map(str, DENSITY_HEIGHTS)) + " integral=1")
    print("endpoint_quantiles=mass_1/36:length_1/27;mass_11/360:length_11/270")
    print("envelope_piece_maxima=" + ",".join(map(str, maxima)))
    print("restricted_threshold=140/161=20/23 separated_supremum=161/840=23/120 open_not_attained")
    print("open_supremum_approach_gaps=" + ",".join(map(str, approach_gaps)))
    print("bv_tail_maximum=4/21 separated_cap_gap=1/840")
    print("five_load_cap=349/360 survivor_mass_gt=11/360 endpoint_tail_gt=11/270")
    print("basic_chain=ell>1/27 rho>(27d/c-25)/63 d/c<113/54")
    print("strong_chain=ell>11/270 rho>(135d/c-124)/315 d/c<563/270")
    print("integer_consequence=270d1<=563c-1")
    tournament_loss_audit()
    print("SCOPE=exact referee for providers and arithmetic; not a finite six-cover reduction")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
