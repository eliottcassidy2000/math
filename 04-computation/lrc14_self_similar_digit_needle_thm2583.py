#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2583.

Checks the uniform-grid description of the 13-target tooth boundaries,
old/future equal-root piercing by sufficiently deep base-13 preimages,
digit-cylinder isolation inside a prescribed rational carrier interval, and
the translated deepest-probe cover.  All geometry uses Fraction arithmetic.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def circle_distance(x, y):
    difference = (x - y) % 1
    return min(difference, 1 - difference)


def boundary(speed, ell, tooth, epsilon, target):
    return (
        (Fraction(tooth, 1) + Fraction(epsilon * ell, 14)
         - Fraction(target, P))
        / speed
    ) % 1


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def cylinder_containing(x, depth):
    scale = P**depth
    index = (x * scale).numerator // (x * scale).denominator
    return Fraction(index, scale), Fraction(index + 1, scale)


print("== THM-2583: self-similar digit needles ==")


print("\n== base boundary grid and future-root coverage ==")
grid_packets = 0
grid_points = 0
digit_cells = 0
for speed in range(1, 201):
    for ell in (1, 2):
        for epsilon in (-1, 1):
            points = sorted(
                boundary(speed, ell, tooth, epsilon, target)
                for tooth in range(speed)
                for target in range(P)
            )
            require(len(set(points)) == P * speed,
                    "fixed-sign target/tooth grid collided")
            gaps = [
                (points[(index + 1) % len(points)] - points[index]) % 1
                for index in range(len(points))
            ]
            require(set(gaps) == {Fraction(1, P * speed)},
                    "fixed-sign boundary set is not a uniform 13k-grid")
            counts = [0] * P
            for point in points:
                require(point.denominator % 7 == 0,
                        "boundary lost its septimal denominator obstruction")
                counts[int(P * point)] += 1
            require(counts == [speed] * P,
                    "future digit cells do not each contain k boundaries")
            grid_packets += 1
            grid_points += len(points)
            digit_cells += P

print(f"  exact (k,L,sign) grids: {grid_packets}")
print(f"  exact base boundary points: {grid_points}")
print(f"  future digit-cell counts checked: {digit_cells}")
print("  each sign is a shifted uniform 13k-grid with k points per root digit")
print("  every boundary denominator retains a factor 7")


print("\n== equal-root carrier piercing and digit-cylinder isolation ==")
piercing_cases = 0
atlas_trace_checks = 0
root_checks = 0
maximum_cylinder_depth = 0
for speed in range(1, 21):
    for ell in (1, 2):
        old_root = (7 * speed + 3 * ell) % P

        # A hostile-looking small carrier interval strictly inside I_h.
        carrier_left = Fraction(old_root, P) + Fraction(1, 52)
        carrier_right = Fraction(old_root, P) + Fraction(1, 26)
        require(carrier_right < Fraction(old_root + 1, P),
                "test carrier escaped its old-root digit")

        # Choose one base boundary whose immediate future digit is h.
        candidates = []
        for target in range(P):
            for tooth in range(speed):
                for epsilon in (-1, 1):
                    y = boundary(speed, ell, tooth, epsilon, target)
                    if int(P * y) == old_root:
                        candidates.append((y, tooth, epsilon, target))
        require(candidates, "no base boundary in the prescribed future digit")
        y, _, _, target0 = candidates[(speed + ell) % len(candidates)]
        require(int(P * y) == old_root,
                "chosen future boundary has the wrong root")

        # The interval has length 1/52, so N=2 is already beyond the
        # all-N grid-piercing threshold 13^N*length>1.
        delay = 2
        scale = P**delay
        require(scale * (carrier_right - carrier_left) > 1,
                "preimage mesh is not fine enough to pierce the carrier")
        lower = scale * carrier_left - y
        upper = scale * carrier_right - y
        preimage_index = ceil_fraction(lower)
        require(Fraction(preimage_index, 1) < upper,
                "uniform preimage grid missed an interval longer than its mesh")
        x0 = (Fraction(preimage_index, 1) + y) / scale
        require(carrier_left < x0 < carrier_right,
                "selected preimage is not inside the carrier")
        require(int(P * x0) == old_root,
                "old physical root changed")
        require((scale * x0) % 1 == y and int(P * y) == old_root,
                "old/future root diagonal failed")

        high_speed = scale * speed
        atlas = {}
        for target in range(P):
            for tooth in range(high_speed):
                for epsilon in (-1, 1):
                    point = boundary(high_speed, ell, tooth, epsilon, target)
                    require(point not in atlas,
                            "complete high-speed boundary atlas collided")
                    atlas[point] = (tooth, epsilon, target)
        require(x0 in atlas and atlas[x0][2] == target0,
                "piercing point is not the intended target boundary")

        # Refine to a deep physical base-13 digit cylinder.  Besides lying
        # in the source carrier, it must stay in the same old digit and in
        # the same future digit under T^N, and contain no other atlas point.
        cylinder = None
        for depth in range(delay + 1, 40):
            left, right = cylinder_containing(x0, depth)
            if not (carrier_left < left < x0 < right < carrier_right):
                continue
            branch_index = (scale * x0).numerator // (scale * x0).denominator
            if not (Fraction(branch_index, scale) <= left
                    and right <= Fraction(branch_index + 1, scale)):
                continue
            future_left = scale * left - branch_index
            future_right = scale * right - branch_index
            if not (Fraction(old_root, P) < future_left < y
                    < future_right < Fraction(old_root + 1, P)):
                continue
            traces = [point for point in atlas if left < point < right]
            if traces == [x0]:
                cylinder = (depth, left, right)
                break
        require(cylinder is not None,
                "no internal digit cylinder isolated the chosen boundary")
        depth, left, right = cylinder
        maximum_cylinder_depth = max(maximum_cylinder_depth, depth)

        require((left * P**depth).denominator == 1
                and (right * P**depth).denominator == 1
                and right - left == Fraction(1, P**depth),
                "selected interval is not a base-13 cylinder")
        require(x0.denominator % 7 == 0,
                "piercing boundary could lie on a digit wall")
        for point in atlas:
            atlas_trace_checks += 1
            require((left < point < right) == (point == x0),
                    "digit cylinder has the wrong complete-atlas trace")
        root_checks += 2
        piercing_cases += 1

print(f"  exact carrier-piercing cases: {piercing_cases}")
print(f"  complete-atlas cylinder traces: {atlas_trace_checks}")
print(f"  old/future equal-root checks: {root_checks}")
print(f"  maximum selected digit depth: {maximum_cylinder_depth}")
print("  one physical digit cylinder lies in the carrier and isolates one target boundary")


print("\n== translated deepest-probe cover ==")
probe_checks = 0
single_counts = 0
double_counts = 0
for numerator in range(997):
    theta = Fraction(numerator, 997)
    count = sum(
        circle_distance(theta, Fraction(residue, P)) < Fraction(1, 14)
        for residue in range(P)
    )
    require(count in (1, 2),
            "thirteen translated ordinary danger arcs lost the 1/2 cover law")
    single_counts += count == 1
    double_counts += count == 2
    probe_checks += 1

require(single_counts > 0 and double_counts > 0,
        "translated-probe cover did not exercise both multiplicities")
print(f"  exact rational circle points checked: {probe_checks}")
print(f"  multiplicity-one / multiplicity-two: {single_counts} / {double_counts}")
print("  every point lies in one or two of the thirteen translated deep probes")


print("\nabstract consequence:")
print("  positive rational neutral carrier U + sufficiently late 13^N k gate")
print("    -> internal base-13 cylinder H <= U with one absolute boundary trace")
print("    -> fixed old/future root and every target Abel character nonzero")
print("  inside THM-2559 W: owner/word/deep base incidence only; shifted transport unproved")
print("\nscope: nonneutral live provenance cannot be combined with shifted colours")
print("target-informed selector covariance and paired endpoint action remain open")
print("all exact checks passed")
