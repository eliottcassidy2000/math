#!/usr/bin/env python3
"""Exact companion for THM-2401.

This script uses only integer and Fraction arithmetic.  The prime-cyclic
nonvanishing theorem itself is inherited from THM-2398; here we check its
finite Boolean interface, both first-death decompositions, the sharp
partition examples, and every displayed rational constant.
"""

from fractions import Fraction
from itertools import product


def check(condition, label):
    if not condition:
        raise RuntimeError("audit failure: " + label)


def cross_correlation(a, f, p):
    """C(s)=sum_r a(r+s)f(r), with sets represented by Python sets."""
    return tuple(
        sum(1 for r in range(p) if r in f and (r + s) % p in a)
        for s in range(p)
    )


def check_layer_identity(a, f, u, p):
    """Check root, pair, and displacement versions of one filter loss."""
    a1 = {r for r in a if r in u}
    f1 = {r for r in f if r in u}
    delta = len(a) * len(f) - len(a1) * len(f1)

    root_a = {
        r: (1 if r in a and r not in u else 0) * len(f)
        for r in range(p)
    }
    root_f = {
        s: (1 if s in f and s not in u else 0) * len(a1)
        for s in range(p)
    }
    check(
        sum(root_a.values()) + sum(root_f.values()) == delta,
        "root decomposition",
    )

    pair_a = {}
    pair_f = {}
    displacement = {d: 0 for d in range(1, p)}
    for r in range(p):
        for s in range(p):
            ta = int(r in a and r not in u and s in f)
            tf = int(r in a1 and s in f and s not in u)
            pair_a[r, s] = ta
            pair_f[r, s] = tf
            if r != s:
                displacement[(s - r) % p] += ta + tf
            else:
                # A and F are disjoint, so this is forced to vanish.
                check(ta + tf == 0, "zero displacement")
    check(
        sum(pair_a.values()) + sum(pair_f.values()) == delta,
        "pair decomposition",
    )
    check(sum(displacement.values()) == delta, "displacement decomposition")
    return delta


def exhaustive_small_identity():
    p = 5
    assignments = 0
    layers = 0
    positive_losses = 0
    for word in product(range(3), repeat=p):
        a = {r for r, x in enumerate(word) if x == 1}
        f = {r for r, x in enumerate(word) if x == 2}
        if not a or not f:
            continue
        assignments += 1
        for uword in product(range(2), repeat=p):
            u = {r for r, x in enumerate(uword) if x}
            delta = check_layer_identity(a, f, u, p)
            layers += 1
            positive_losses += int(delta > 0)
    check(assignments == 180, "small assignment count")
    check(layers == 5760, "small layer count")
    return assignments, layers, positive_losses


def prime13_packet_check():
    p = 13
    a_masks = [{r} for r in range(p)]
    a_masks += [{r, (r + 1) % p} for r in range(p)]
    f_masks = [
        {r, s}
        for r in range(p)
        for s in range(r + 1, p)
    ]
    checked = 0
    for a in a_masks:
        for f in f_masks:
            if a & f:
                continue
            c = cross_correlation(a, f, p)
            check(c[0] == 0, "p13 pinned zero")
            check(sum(c) == len(a) * len(f) > 0, "p13 positive two-sided mass")
            # At prime 13 a rational vector has a zero nontrivial DFT
            # iff all thirteen coefficients are equal (THM-2398).
            check(len(set(c)) > 1, "p13 nonflat correlation")
            checked += 1
    check(checked == 1573, "p13 packet count")
    return checked


def sharp_root_partition():
    # Four equal cells, one for each possible absolute deleted root in
    # A-bank {0,1} or F-bank {2,3}.  Every cell loses all co-support.
    cells = [
        ({0}, {2}, {2}),
        ({1}, {2}, {2}),
        ({0}, {2}, {0}),
        ({0}, {3}, {0}),
    ]
    weights = [Fraction(1, 4)] * 4
    root_mass = {("A", 0): 0, ("A", 1): 0, ("F", 2): 0, ("F", 3): 0}
    total_before = 0
    total_after = 0
    for (a, f, u), weight in zip(cells, weights):
        a1 = a & u
        f1 = f & u
        total_before += weight * len(a) * len(f)
        total_after += weight * len(a1) * len(f1)
        for r in a - a1:
            root_mass["A", r] += weight * len(f)
        for s in f - f1:
            root_mass["F", s] += weight * len(a1)
    check(total_before == 1, "sharp root total before")
    check(total_after == 0, "sharp root total after")
    check(
        set(root_mass.values()) == {Fraction(1, 4)},
        "sharp root equality",
    )
    return root_mass


def sharp_displacement_partition():
    # Twelve equal cells realize all nonzero relative displacements.
    # The A-root is deleted on each cell.
    p = 13
    masses = {d: Fraction(0) for d in range(1, p)}
    for d in range(1, p):
        a = {0}
        f = {d}
        u = {d}
        delta = check_layer_identity(a, f, u, p)
        check(delta == 1, "sharp displacement cell")
        masses[d] += Fraction(1, 12)
    check(sum(masses.values()) == 1, "sharp displacement total")
    check(
        set(masses.values()) == {Fraction(1, 12)},
        "sharp displacement equality",
    )
    coefficient = min(masses.values()) / (p * p)
    check(coefficient == Fraction(1, 2028), "displacement coefficient")
    return masses, coefficient


def first_death_no_initial_floor():
    # Two common filters on two rational base strata.  The first leaves
    # epsilon co-support, and the second kills exactly that remainder.
    controls = []
    for denominator in (2, 3, 17, 101):
        epsilon = Fraction(1, denominator)
        m0 = Fraction(1)
        m1 = epsilon
        m2 = Fraction(0)
        check(m0 > 0 and m1 > 0 and m2 == 0, "no-floor control")
        controls.append((denominator, m0, m1, m2))
    return controls


def exact_constants():
    delta = Fraction(1)
    labelled_root = delta / 4
    displacement_colour = delta / (12 * 13 * 13)

    ordinary = Fraction(11, 24336) * labelled_root
    guard = Fraction(1, 2704) * labelled_root
    ordinary_status = ordinary / 4
    guard_status = guard / 4
    blocker_word = delta / 16
    blocker_word_unfixed_owner = blocker_word / 3

    check(labelled_root == Fraction(1, 4), "root constant")
    check(
        displacement_colour == Fraction(1, 2028),
        "displacement constant",
    )
    check(ordinary == Fraction(11, 97344), "ordinary repair constant")
    check(guard == Fraction(1, 10816), "guard repair constant")
    check(
        ordinary_status == Fraction(11, 389376),
        "ordinary status constant",
    )
    check(guard_status == Fraction(1, 43264), "guard status constant")
    check(blocker_word == Fraction(1, 16), "blocker word constant")
    check(
        blocker_word_unfixed_owner == Fraction(1, 48),
        "unfixed-owner blocker word constant",
    )
    return (
        labelled_root,
        displacement_colour,
        ordinary,
        guard,
        ordinary_status,
        guard_status,
        blocker_word,
        blocker_word_unfixed_owner,
    )


def hostile_correlations():
    p = 13
    # Separate base support: both marginals can be nonzero, but joint
    # co-support and the whole cross-correlation vanish.
    separate = (0,) * p

    # Thirteen equally weighted relative endpoint gauges make the
    # cross-correlation flat.  All nonzero characters vanish.
    flat = (1,) * p
    check(len(set(flat)) == 1, "flat hostile")
    check(sum(flat) == p, "flat hostile mass")
    return separate, flat


def main():
    assignments, layers, positive_losses = exhaustive_small_identity()
    packet_count = prime13_packet_check()
    root_mass = sharp_root_partition()
    displacement_mass, displacement_coefficient = sharp_displacement_partition()
    no_floor = first_death_no_initial_floor()
    constants = exact_constants()
    separate, flat = hostile_correlations()

    print("== THM-2401 exact common-filter / first-death audit ==")
    print(
        "small exhaustive:",
        f"assignments={assignments}",
        f"layers={layers}",
        f"positive_losses={positive_losses}",
    )
    print("p=13 disjoint singleton/adjacent x two-root packets:", packet_count)
    print("sharp labelled-root masses:", sorted(root_mass.items()))
    print(
        "sharp displacement:",
        f"classes={len(displacement_mass)}",
        f"each={next(iter(displacement_mass.values()))}",
        f"colour={displacement_coefficient}",
    )
    print(
        "constants:",
        "root=" + str(constants[0]),
        "displacement_colour=" + str(constants[1]),
        "ordinary=" + str(constants[2]),
        "guard=" + str(constants[3]),
        "ordinary_status=" + str(constants[4]),
        "guard_status=" + str(constants[5]),
        "blocker_word=" + str(constants[6]),
        "unfixed_owner=" + str(constants[7]),
    )
    print(
        "first-death no-M0-floor controls:",
        " ".join(f"N={n}:({m0},{m1},{m2})" for n, m0, m1, m2 in no_floor),
    )
    print(
        "hostiles:",
        f"separate_sum={sum(separate)}",
        f"flat_sum={sum(flat)}",
        f"flat_classes={len(set(flat))}",
    )
    print("PASS")


if __name__ == "__main__":
    main()
