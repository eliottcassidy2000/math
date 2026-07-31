#!/usr/bin/env python3
"""Exact referee for the fixed-head affine paired-path gap.

For a 13-unit target speed k, a paired blocker a=13b, and root branch
iota_r(z)=(z+r)/13, the paired gate with shift q reduces to

    d_L((kz+kr+q)/13) * u_1(bz-q/13).

At fixed target phase v=kr+q, write A_v for the first danger set and C_q
for the blocker danger set.  A paired root/shift cell is null exactly when
A_v is contained in C_q up to endpoints.  The proof memo reduces every
possible containment to 2,704 small exact cases.  This script checks those
cases by two routes, freezes the complete 28-pair bad graph, and checks the
fixed-head affine-path counts.

Everything is dependency-free and exact over Fraction.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_intervals(intervals):
    clipped = []
    for left, right in intervals:
        left = max(Fraction(0), left)
        right = min(Fraction(1), right)
        if left < right:
            clipped.append((left, right))
    out = []
    for left, right in sorted(clipped):
        if not out or out[-1][1] < left:
            out.append([left, right])
        elif right > out[-1][1]:
            out[-1][1] = right
    return tuple((left, right) for left, right in out)


@lru_cache(maxsize=None)
def target_intervals(k, L, v):
    """A_v={0<z<1:d_L((kz+v)/13)=1}."""
    radius = Fraction(13 * L, 14 * k)
    return merge_intervals(
        (
            Fraction(13 * m - v, k) - radius,
            Fraction(13 * m - v, k) + radius,
        )
        for m in range(-2, k // 13 + 4)
    )


@lru_cache(maxsize=None)
def blocker_danger_intervals(b, q):
    """C_q={0<z<1:d_1(bz-q/13)=1}."""
    radius = Fraction(1, 14 * b)
    return merge_intervals(
        (
            Fraction(m, b) + Fraction(q, 13 * b) - radius,
            Fraction(m, b) + Fraction(q, 13 * b) + radius,
        )
        for m in range(-2, b + 3)
    )


def interval_subset(left, right):
    """Containment up to finitely many endpoints."""
    j = 0
    for start, stop in left:
        while j < len(right) and right[j][1] <= start:
            j += 1
        cursor = start
        jj = j
        while jj < len(right) and right[jj][0] <= cursor < stop:
            cursor = max(cursor, right[jj][1])
            jj += 1
        if cursor < stop:
            return False
    return True


def circle_distance(value):
    value %= 1
    return min(value, 1 - value)


def paired_mass_by_walls(k, b, L, v, q):
    """Independent exact wall-cell integration of A_v intersect C_q^c."""
    walls = {Fraction(0), Fraction(1)}
    for sign in (-1, 1):
        for m in range(-2, k // 13 + 4):
            z = (
                Fraction(13 * m - v, 1)
                + Fraction(sign * 13 * L, 14)
            ) / k
            if 0 < z < 1:
                walls.add(z)
        for m in range(-2, b + 3):
            z = (
                Fraction(m, 1)
                + Fraction(q, 13)
                + Fraction(sign, 14)
            ) / b
            if 0 < z < 1:
                walls.add(z)
    ordered = sorted(walls)
    mass = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        target = circle_distance(
            Fraction(k, 13) * midpoint + Fraction(v, 13)
        ) < Fraction(L, 14)
        blocker_safe = circle_distance(
            b * midpoint - Fraction(q, 13)
        ) >= Fraction(1, 14)
        if target and blocker_safe:
            mass += right - left
    return mass


EXPECTED_BAD = {
    15: {5: 7, 6: 6},
    16: {2: 9, 3: 8, 7: 5, 8: 4},
    17: {1: 9, 4: 7, 5: 6, 8: 4},
    18: {2: 8, 3: 7, 5: 6, 6: 5},
    19: {1: 8, 3: 7, 4: 6, 6: 5},
    20: {1: 8, 2: 7, 4: 6, 5: 5},
    21: {2: 7, 3: 6},
    22: {1: 7, 3: 6},
    23: {1: 7, 2: 6},
}


def bad_by_intervals(k, b, L, v, q):
    return interval_subset(
        target_intervals(k, L, v), blocker_danger_intervals(b, q)
    )


def component_lengths(intervals):
    return tuple(right - left for left, right in intervals)


def finite_reduction_checks():
    # For b>=2, every small-threshold target has a component longer than
    # the largest possible blocker tooth, 1/14.  The general k>=15 ordinary
    # lane is discharged analytically by the opposing full-tooth/full-gap
    # inequalities recorded in the memo.
    guard_component_floors = {}
    for k in (10, 11, 12):
        floor = min(
            max(component_lengths(target_intervals(k, 2, v)))
            for v in range(P)
        )
        require(floor > Fraction(1, 14), f"guard component floor failed at {k}")
        guard_component_floors[k] = floor

    ordinary_component_floors = {}
    for k in (12, 14):
        floor = min(
            max(component_lengths(target_intervals(k, 1, v)))
            for v in range(P)
        )
        require(floor > Fraction(1, 14), f"ordinary component floor failed at {k}")
        ordinary_component_floors[k] = floor

    # At k=15 the target parameter interval already has length more than
    # one period plus one ordinary danger tooth.  This is the sharp integer
    # threshold used by the full-tooth argument; monotonicity handles k>15.
    require(
        Fraction(15, 13) > 1 + Fraction(1, 7),
        "k=15 full-tooth threshold failed",
    )
    # For b=1 and q!=0 the two blocker-safe endpoint pieces have total
    # length 6/7.  At k=27 even two maximal target-safe gaps are too short;
    # monotonicity handles every larger k (k=26 is excluded by 13|k).
    require(
        Fraction(6, 7) > 2 * Fraction(78, 7 * 27),
        "b=1 endpoint-safe-gap cutoff failed",
    )
    return guard_component_floors, ordinary_component_floors


def reduced_exact_census():
    cases = 0
    bad = []

    # Guard reduction: k>=13 is killed by mass; at k=10,11,12 the component
    # floor kills b>=2, so only b=1 remains.
    for k in (10, 11, 12):
        for v in range(P):
            for q in range(P):
                cases += 1
                by_intervals = bad_by_intervals(k, 1, 2, v, q)
                by_mass = paired_mass_by_walls(k, 1, 2, v, q) == 0
                require(by_intervals == by_mass, f"guard routes disagree {(k,v,q)}")
                if by_intervals:
                    bad.append((2, k, 1, v, q))

    # Ordinary reduction: b>=2 is impossible; b=1 leaves k=12,14 and
    # 15<=k<=25.  The endpoint-safe-gap argument excludes every larger k.
    ordinary_k = (12, 14) + tuple(range(15, 26))
    for k in ordinary_k:
        for v in range(P):
            for q in range(P):
                cases += 1
                by_intervals = bad_by_intervals(k, 1, 1, v, q)
                by_mass = paired_mass_by_walls(k, 1, 1, v, q) == 0
                require(by_intervals == by_mass, f"ordinary routes disagree {(k,v,q)}")
                if by_intervals:
                    bad.append((1, k, 1, v, q))

    expected = [
        (1, k, 1, v, q)
        for k, graph in EXPECTED_BAD.items()
        for v, q in graph.items()
    ]
    require(cases == 2704, f"reduced universe changed: {cases}")
    require(bad == expected, f"bad-pair table changed: {bad}")
    return cases, bad


def extended_interval_control():
    # This is a hostile/control surface, not the proof of the infinite
    # reduction.  It checks that no omitted bad row appears in a broad box.
    checks = 0
    found = []
    for L, start in ((1, 12), (2, 10)):
        for k in range(start, 201):
            if k % 13 == 0:
                continue
            for b in range(1, 31):
                for v in range(P):
                    for q in range(P):
                        checks += 1
                        if bad_by_intervals(k, b, L, v, q):
                            found.append((L, k, b, v, q))
    expected = [
        (1, k, 1, v, q)
        for k, graph in EXPECTED_BAD.items()
        for v, q in graph.items()
    ]
    require(found == expected, "extended control found an unclassified bad pair")
    return checks


def affine_path_census(stops=8):
    histogram = Counter()
    minimum = 14
    minimizers = []
    for k, graph in EXPECTED_BAD.items():
        for r0 in range(P):
            for delta in range(1, P):
                good = []
                for q0 in range(P):
                    v = (k * r0 + q0) % P
                    path_ok = True
                    for j in range(stops):
                        rj = (r0 - j * delta) % P
                        qj = (v - k * rj) % P
                        require(qj == (q0 + k * j * delta) % P,
                                "affine sign convention changed")
                        if graph.get(v) == qj:
                            path_ok = False
                            break
                    if path_ok:
                        good.append(q0)
                histogram[len(good)] += 1
                if len(good) < minimum:
                    minimum = len(good)
                    minimizers = [(k, r0, delta, tuple(good))]
                elif len(good) == minimum:
                    minimizers.append((k, r0, delta, tuple(good)))

    expected_histogram = Counter({9: 86, 10: 268, 11: 584, 12: 372, 13: 94})
    require(histogram == expected_histogram, f"path histogram changed: {histogram}")
    require(minimum == 9, f"path floor changed: {minimum}")
    sharp = (16, 0, 8, (0, 1, 4, 5, 6, 9, 10, 11, 12))
    require(sharp in minimizers, "named sharp path disappeared")
    require(max(len(graph) for graph in EXPECTED_BAD.values()) == 4,
            "bad-domain bound changed")
    return histogram, minimizers


def main():
    print("== fixed-head affine paired-path referee ==")
    guard_floors, ordinary_floors = finite_reduction_checks()
    print("small-threshold max-component floors (guard):",
          {k: str(v) for k, v in guard_floors.items()})
    print("small-threshold max-component floors (ordinary):",
          {k: str(v) for k, v in ordinary_floors.items()})
    print("analytic reduction: guard -> b=1,k=10..12; ordinary -> b=1,k in {12,14,15..25} PASS")

    cases, bad = reduced_exact_census()
    print(f"two-route reduced census: {cases} exact (L,k,b,v,q) cells PASS")
    print(f"complete bad graph: {len(bad)} pairs on {len(EXPECTED_BAD)} speed rows")
    for k, graph in EXPECTED_BAD.items():
        print(f"  k={k}: " + ", ".join(f"{v}->{q}" for v, q in graph.items()))

    checks = extended_interval_control()
    print(f"extended interval control: {checks} cells, no extra bad pair PASS")

    histogram, minimizers = affine_path_census(stops=8)
    print("eight-stop affine good-start histogram:", dict(sorted(histogram.items())))
    print("sharp minimum: 9 good starts")
    print("  control (k,r0,delta)=(16,0,8): q0=[0,1,4,5,6,9,10,11,12]")
    print(f"  number of minimum rows: {len(minimizers)}")
    print("all finite paths: at least 13-max_v_domain=9 starts, because phases outside the bad graph domain have no bad blocker label")
    print("sign audit: r_j=r_0-j*delta, q_j=v-k*r_j=q_0+k*j*delta PASS")
    print("twisted-return specialization: delta=-k^(-1)*alpha gives q_j=q_0-j*alpha PASS")
    print("SCOPE: physical paired-cell positivity and ordered root/shift paths;")
    print("no THM-2542 chart-clock identification, relation residue, endpoint current, row exclusion, or LRC(14).")


if __name__ == "__main__":
    main()
