#!/usr/bin/env python3
"""Direct exact referee for the mixed D/slope-seven support witness.

This companion does not reconstruct the component search.  It freezes the
rational midpoint and checks every retained factor there by direct strict
membership, then exhausts all twelve nonzero slope labels.
"""

from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_dilation_reversed_clock_fibre_product_probe as d
import lrc14_predecessor_carry_private_root_atlas_thm2640 as m


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(q):
    return q - q.numerator // q.denominator


def in_intervals(q, intervals):
    return any(Fraction(a) < q < Fraction(b) for a, b in intervals)


def in_weighted(q, pieces):
    return any(weight and Fraction(a) < q < Fraction(b)
               for a, b, weight in pieces)


def atom_membership(atom, q, prefixes):
    base = frac(q) * m.T
    if not in_weighted(base, atom["pieces"]):
        return False
    y = frac(m.R * q) * m.T
    prefix = prefixes[atom["future"]][atom["h"]]
    return in_intervals(y, tuple(
        (a, a + length) for a, length in zip(prefix[0], prefix[1])
    ))


def private_membership(module, pair_prefixes, rails, present, q,
                       carry, root):
    rail, sector, edge, kappa, h, shallow = 0, 0, 0, 1, 6, 1
    if not in_weighted(frac(q) * m.T, rails[rail][3]):
        return False
    if not in_intervals(frac(q) * m.T, present[shallow, (-h) % m.P]):
        return False
    deep = frac(module.C3 * q) * 182
    if not (Fraction(14 * root - 13) < deep < Fraction(14 * root)):
        return False
    actual_carry = ((m.R * q).numerator // (m.R * q).denominator) % m.P
    if actual_carry != carry:
        return False
    y = frac(m.R * q) * m.T
    prefix = pair_prefixes[sector][shallow][h][kappa]
    return in_intervals(y, tuple(
        (a, a + length) for a, length in zip(prefix[0], prefix[1])
    ))


def main():
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    component = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    require(component[0] < x < component[1],
            "frozen point left the exact D-edge component")
    require(frac(13 * x) == z, "frozen following coordinate is not D(x)")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    current = d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require((current["j"], current["h"], following["j"], following["h"])
            == (5, 2, 2, 6), "frozen D labels changed")
    require(current["h"] == following["j"],
            "D digit covariance failed at the labelled witness")
    require(atom_membership(current, x, prefixes),
            "frozen point left the current D atom")
    require(atom_membership(following, z, prefixes),
            "frozen point left the following D atom")

    pair_prefixes = m.build_pair_prefixes(module)
    shard = m.shard((0, 1))
    require(shard[1] == 26 and shard[5] == ((1, 0, 12),),
            "private-rail normalization changed")
    rows = shard[6][0]
    sector, edge, kappa, h = 0, 0, 1, 6
    carry, root = 2, 6
    require(root == (2 * carry + (2 * h + kappa) // m.P + 1) % m.P,
            "base private-root law failed")
    require(m.is_unit(rows[sector][edge][carry][kappa][h], root, 26),
            "base private row stopped being a unit")
    require(private_membership(
        module, pair_prefixes, rails, present, z, carry, root
    ), "frozen following point left the base private packet")

    successes = []
    failures = []
    for delta in range(1, m.P):
        shifted = frac(z + Fraction(7 * delta, m.R))
        carry2 = (carry + 7 * delta) % m.P
        root2 = (root + delta) % m.P
        formula_root = (2 * carry2 + (2 * h + kappa) // m.P + 1) % m.P
        require(root2 == formula_root,
                "slope-seven carry/root covariance failed")
        if root2 == 0:
            failures.append((delta, "root_zero", carry2, root2))
            continue
        unit = m.is_unit(
            rows[sector][edge][carry2][kappa][h], root2, 26
        )
        physical = private_membership(
            module, pair_prefixes, rails, present,
            shifted, carry2, root2,
        )
        if unit and physical:
            successes.append((delta, carry2, root2, shifted))
        else:
            failures.append((delta, "unit" if not unit else "physical",
                             carry2, root2))

    expected = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12)
    require(tuple(row[0] for row in successes) == expected,
            "mixed slope-success bank changed")
    require(tuple(failures) == ((7, "root_zero", 12, 0),
                                (8, "unit", 6, 1)),
            "mixed slope-failure boundary changed")

    print("LRC14 MIXED D/SLOPE-SEVEN DIRECT EXACT REFEREE")
    print(f"component={component}")
    print(f"x={x} z=D(x)={z}")
    print("D_clock_path=(3,1,0) D_digit_match=(h_current,j_following)=(2,2)")
    print("D_current=(source1,rail2,j5,h2,epsilon0,kappa0)")
    print("D_following=(source1,rail0,j2,h6,epsilon1,kappa1)")
    print("private_base=(rail0,sector0,edge0,local_kappa1,h6,shallow1,carry2,root6)")
    print(f"successful_slope_switches={tuple(successes)}")
    print(f"sharp_failures={tuple(failures)}")
    print("verdict=PASS: every printed factor holds strictly at one rational point")


if __name__ == "__main__":
    main()
