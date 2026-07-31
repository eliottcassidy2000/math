#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2599.

The theorem itself is analytic.  This companion exhausts the translated-tooth
identity, the labelled wall-count mechanism, and a bounded bank of complete
root-cell paired-slice measures using only integer and Fraction arithmetic.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_q(value):
    return value.numerator // value.denominator


def fractional_part(value):
    return value - floor_q(value)


def circle_norm(value):
    residue = fractional_part(value)
    return min(residue, 1 - residue)


def danger(value, level=1):
    return circle_norm(value) < Fraction(level, 14)


def factor_walls(coefficient, shift_sign, shift, level, left, right):
    """Strict danger walls for c*y + shift_sign*shift/13 on (left,right)."""
    translated_left = coefficient * left + shift_sign * Fraction(shift, P)
    translated_right = coefficient * right + shift_sign * Fraction(shift, P)
    walls = set()
    for integer in range(floor_q(translated_left) - 2,
                         floor_q(translated_right) + 4):
        for side in (-1, 1):
            wall = (
                Fraction(integer)
                + side * Fraction(level, 14)
                - shift_sign * Fraction(shift, P)
            ) / coefficient
            if left < wall < right:
                walls.add(wall)
    return walls


def vector_walls(coefficient, shift_sign, level, root):
    """All labelled walls (shift,point) in one open root cell."""
    left = Fraction(root, P)
    right = Fraction(root + 1, P)
    labelled = []
    for shift in range(P):
        for wall in factor_walls(
            coefficient, shift_sign, shift, level, left, right
        ):
            labelled.append((shift, wall))
    return labelled


def paired_measures(target, blocker, level, root):
    """Exact measures of all thirteen opposite-shift paired slices."""
    left = Fraction(root, P)
    right = Fraction(root + 1, P)
    walls = {left, right}
    for shift in range(P):
        walls |= factor_walls(target, +1, shift, level, left, right)
        walls |= factor_walls(blocker, -1, shift, 1, left, right)
    walls = sorted(walls)

    target_masses = [Fraction(0) for _ in range(P)]
    blocker_masses = [Fraction(0) for _ in range(P)]
    paired_masses = [Fraction(0) for _ in range(P)]
    pointwise_minimum = P
    for chamber_left, chamber_right in zip(walls, walls[1:]):
        midpoint = (chamber_left + chamber_right) / 2
        width = chamber_right - chamber_left
        live = 0
        for shift in range(P):
            target_bit = danger(
                target * midpoint + Fraction(shift, P), level
            )
            blocker_bit = danger(
                blocker * midpoint - Fraction(shift, P), 1
            )
            if target_bit:
                target_masses[shift] += width
            if blocker_bit:
                blocker_masses[shift] += width
            if target_bit and not blocker_bit:
                paired_masses[shift] += width
                live += 1
        pointwise_minimum = min(pointwise_minimum, live)

    return target_masses, blocker_masses, paired_masses, pointwise_minimum


def paired_chambers(target, blocker, level, root):
    """Open common wall chambers and their full paired shift profiles."""
    left = Fraction(root, P)
    right = Fraction(root + 1, P)
    walls = {left, right}
    for shift in range(P):
        walls |= factor_walls(target, +1, shift, level, left, right)
        walls |= factor_walls(blocker, -1, shift, 1, left, right)
    walls = sorted(walls)
    result = []
    for chamber_left, chamber_right in zip(walls, walls[1:]):
        midpoint = (chamber_left + chamber_right) / 2
        profile = tuple(
            shift
            for shift in range(P)
            if danger(target * midpoint + Fraction(shift, P), level)
            and not danger(blocker * midpoint - Fraction(shift, P), 1)
        )
        result.append((chamber_left, chamber_right, profile))
    return result


print("== translated-tooth identity ==")
count_values = {1: set(), 2: set()}
for level in (1, 2):
    # Every relevant wall lies on the 1/182 grid; these are all open chambers.
    for chamber in range(182):
        y = Fraction(2 * chamber + 1, 364)
        for sign in (-1, 1):
            count = sum(
                danger(y + sign * Fraction(shift, P), level)
                for shift in range(P)
            )
            expected = 2 * level - danger(P * y, level)
            require(count == expected, "translated-tooth identity failed")
            count_values[level].add(count)
require(count_values == {1: {1, 2}, 2: {3, 4}},
        "translated-tooth count range changed")
print("  ordinary translated danger count: 1 or 2")
print("  guard translated danger count: 3 or 4")


print("\n== labelled ordinary wall census ==")
wall_cases = 0
for coefficient in range(1, 81):
    for root in range(P):
        for shift_sign in (-1, 1):
            labelled = vector_walls(coefficient, shift_sign, 1, root)
            require(len(labelled) == 2 * coefficient,
                    "ordinary vector wall count changed")
            require(len(set(labelled)) == len(labelled),
                    "a labelled ordinary wall was duplicated")
            require(len({wall for _, wall in labelled}) == len(labelled),
                    "two translated components shared an ordinary wall")
            wall_cases += 1
print(f"  exact coefficient/root/sign cases: {wall_cases}")
print("  vector wall count on every root: 2*coefficient")
print("  no wall is a root endpoint or belongs to two shift components")


print("\n== prime Boolean-profile Fourier support ==")


def cyclotomic_reduce(coefficients):
    """Reduce a length-13 rational coefficient vector modulo Phi_13."""
    require(len(coefficients) == P, "cyclotomic vector has wrong length")
    top = coefficients[-1]
    return tuple(value - top for value in coefficients[:-1])


profile_frequency_checks = 0
for mask in range(1, (1 << P) - 1):
    profile = [shift for shift in range(P) if mask & (1 << shift)]
    for frequency in range(1, P):
        coefficients = [0 for _ in range(P)]
        for shift in profile:
            coefficients[(frequency * shift) % P] += 1
        require(any(cyclotomic_reduce(coefficients)),
                "a nonempty proper Boolean profile lost a primitive colour")
        profile_frequency_checks += 1
require(profile_frequency_checks == 8190 * 12,
        "Boolean-profile Fourier check count changed")
print("  nonempty proper Boolean profiles: 8190")
print(f"  exact nonzero-frequency checks: {profile_frequency_checks}")
print("  every proper Boolean profile has all 12 primitive colours")


print("\n== complete opposite-shift root-cell bank ==")
targets = [k for k in range(1, 21) if k % P]
blockers = [P * b for b in range(1, 5)]
case_count = 0
ordinary_minimum = None
guard_minimum = None
ordinary_argmin = None
guard_argmin = None
for target in targets:
    for blocker in blockers:
        require(target != blocker, "typed target/blocker unexpectedly equal")
        for root in range(P):
            for level in (1, 2):
                target_masses, blocker_masses, paired, pointwise_minimum = (
                    paired_measures(target, blocker, level, root)
                )
                require(sum(target_masses, Fraction(0)) == Fraction(level, 7),
                        "target shift-total integral changed")
                require(sum(blocker_masses, Fraction(0)) == Fraction(1, 7),
                        "blocker shift-total integral changed")
                largest = max(paired)
                require(largest > 0, "a typed root cell lost all paired slices")
                if level == 1:
                    floor = Fraction(1, 182 * target * blocker)
                    require(largest >= floor,
                            "ordinary rational-grid floor failed")
                    record = (largest, target, blocker, root, paired.index(largest))
                    if ordinary_minimum is None or largest < ordinary_minimum:
                        ordinary_minimum = largest
                        ordinary_argmin = record
                else:
                    require(sum(paired, Fraction(0)) >= Fraction(1, P),
                            "guard integrated capacity floor failed")
                    require(largest >= Fraction(1, P * P),
                            "guard pigeonhole floor failed")
                    require(pointwise_minimum >= 1,
                            "guard pointwise cover failed")
                    record = (largest, target, blocker, root, paired.index(largest))
                    if guard_minimum is None or largest < guard_minimum:
                        guard_minimum = largest
                        guard_argmin = record
                case_count += 1

print(f"  exact typed (target,blocker,root,level) cases: {case_count}")
print("  every ordinary root has a positive slice and obeys 1/(182ak)")
print("  every guard root is pointwise covered; some slice has mass >=1/169")
print(
    "  finite-bank smallest largest ordinary slice: "
    f"{ordinary_minimum} at (k,a,h,s)={ordinary_argmin[1:]}"
)
print(
    "  finite-bank smallest largest guard slice: "
    f"{guard_minimum} at (k,a,h,s)={guard_argmin[1:]}"
)


print("\n== valuation-free unequal-slope controls ==")
controls = ((2, 3), (7, 11), (13, 14), (26, 39))
control_cases = 0
for target, blocker in controls:
    for root in range(P):
        _, _, paired, _ = paired_measures(target, blocker, 1, root)
        require(max(paired) >= Fraction(1, 182 * target * blocker),
                "unequal-slope ordinary control failed")
        control_cases += 1
print(f"  non-LRC unequal-slope root controls: {control_cases}")
print("  the ordinary mechanism uses k!=a; 13-adic typing guarantees it")


print("\n== sharp ordinary pointwise hostile ==")
hostile_target = 14
hostile_blocker = 2197
hostile_y = Fraction(733, 737)
target_profile = tuple(
    shift
    for shift in range(P)
    if danger(hostile_target * hostile_y + Fraction(shift, P), 1)
)
blocker_profile = tuple(
    shift
    for shift in range(P)
    if danger(hostile_blocker * hostile_y - Fraction(shift, P), 1)
)
paired_profile = tuple(
    shift for shift in target_profile if shift not in blocker_profile
)
require(target_profile == (1,), "hostile target profile changed")
require(blocker_profile == (1,), "hostile blocker profile changed")
require(not paired_profile, "ordinary pointwise hostile acquired a repair")
hostile_root = floor_q(P * hostile_y)
require(hostile_root == 12, "hostile physical root changed")
require(
    any(
        profile
        for _, _, profile in paired_chambers(
            hostile_target, hostile_blocker, 1, hostile_root
        )
    ),
    "hostile root lost the theorem's positive chamber",
)
print(f"  (k,a,y,h)=({hostile_target},{hostile_blocker},{hostile_y},{hostile_root})")
print(f"  target danger shifts = blocker danger shifts = {target_profile}")
print("  all paired slices vanish at y, but another chamber survives in I_12")


print("\n== toothpick full-shift cylinder control ==")
sample_target = 12
sample_blocker = 13
sample_level = 1
root_chambers = []
for root in range(P):
    positive = [
        chamber
        for chamber in paired_chambers(
            sample_target, sample_blocker, sample_level, root
        )
        if chamber[2]
    ]
    require(positive, "sample root has no positive paired chamber")
    root_chambers.append(positive)


def strict_subcylinder(chambers, depth):
    scale = P**depth
    for left, right, profile in chambers:
        numerator = floor_q(left * scale) + 1
        lower = Fraction(numerator, scale)
        upper = Fraction(numerator + 1, scale)
        if left < lower and upper < right:
            return numerator, profile, (left, right)
    return None


selected_cylinders = None
for common_depth in range(1, 13):
    candidates = [
        strict_subcylinder(chambers, common_depth)
        for chambers in root_chambers
    ]
    if all(candidate is not None for candidate in candidates):
        selected_cylinders = candidates
        break
require(selected_cylinders is not None, "common-depth cylinder search failed")

scale = P**common_depth
profile_sizes = []
for root, (numerator, profile, chamber) in enumerate(selected_cylinders):
    lower = Fraction(numerator, scale)
    upper = Fraction(numerator + 1, scale)
    require(chamber[0] < lower < upper < chamber[1],
            "selected cylinder escaped its wall chamber")
    require(numerator // (P ** (common_depth - 1)) == root,
            "selected cylinder lost its physical first digit")
    require(1 <= len(profile) <= 2 * sample_level,
            "selected cylinder profile has wrong size")
    profile_sizes.append(len(profile))

affine_path = tuple((3 - 2 * step) % P for step in range(8))
require(len(set(affine_path)) == 8, "affine control path repeated early")
path_measure = Fraction(1, P ** (8 * common_depth))
require(path_measure > 0, "specification cylinder vanished")

print(
    f"  control (k,a,L)=({sample_target},{sample_blocker},{sample_level}); "
    f"common cylinder depth: {common_depth}"
)
print(f"  paired profile sizes by root: {tuple(profile_sizes)}")
print(f"  affine 8-word h_j=3-2j: {affine_path}")
print(f"  exact disjoint-block intersection mass: {path_measure}")
print("  concatenating the 13 codewords embeds the full one-sided 13-shift")

print("\nall exact checks passed")
