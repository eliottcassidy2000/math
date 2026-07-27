#!/usr/bin/env python3
"""Exact arithmetic audit for the owner-side acute-current lemma.

The nonvanishing proof is the open-half-plane argument in the companion
note.  This script checks all threshold inequalities and the exact rational
antipodal hostile controls at the first failing multipliers.
"""

from fractions import Fraction


def phase_is_antipodal(multiplier, x, y):
    value = multiplier * (x - y)
    return (value - Fraction(1, 2)).denominator == 1


def main():
    # In units of pi, the unsplit owner cone has width 2m/7 and one
    # owner side has width m/7.
    unsplit_nonzero = [m for m in range(1, 20) if Fraction(2 * m, 7) < 1]
    side_nonzero = [m for m in range(1, 20) if Fraction(m, 7) <= 1]
    side_autocorrelation = [
        m for m in range(1, 20) if Fraction(m, 7) < Fraction(1, 2)
    ]
    cross_pairs = [
        (m, n)
        for m in range(1, 20)
        for n in range(1, 20)
        if Fraction(m + n, 7) < 1
    ]

    assert unsplit_nonzero == [1, 2, 3]
    assert side_nonzero == [1, 2, 3, 4, 5, 6, 7]
    assert side_autocorrelation == [1, 2, 3]
    assert cross_pairs == [
        (m, n)
        for m in range(1, 6)
        for n in range(1, 7 - m)
    ]
    outer_thresholds = {
        m: Fraction(7 - m, 14 * m)
        for m in (4, 5, 6)
    }
    assert outer_thresholds == {
        4: Fraction(3, 56),
        5: Fraction(1, 35),
        6: Fraction(1, 84),
    }

    # First unsplit failure: two equal intervals centered at +/-1/(4m).
    m4 = 4
    epsilon4 = Fraction(1, 1000)
    centre4 = Fraction(1, 4 * m4)
    left4 = (centre4 - epsilon4, -centre4 - epsilon4)
    right4 = (centre4 + epsilon4, -centre4 + epsilon4)
    assert all(abs(x) < Fraction(1, 14) for x in left4 + right4)
    assert phase_is_antipodal(m4, left4[0], left4[1])
    assert phase_is_antipodal(m4, right4[0], right4[1])

    # First one-sided phase-cone failure: two equal intervals in (0,1/14)
    # whose centers differ by 1/(2m).
    m8 = 8
    interval1 = (Fraction(15, 4096), Fraction(17, 4096))
    interval2 = (Fraction(271, 4096), Fraction(273, 4096))
    assert Fraction(0) < interval1[0] < interval1[1] < Fraction(1, 14)
    assert Fraction(0) < interval2[0] < interval2[1] < Fraction(1, 14)
    assert interval2[0] - interval1[0] == Fraction(1, 2 * m8)
    assert interval2[1] - interval1[1] == Fraction(1, 2 * m8)
    assert phase_is_antipodal(m8, interval2[0], interval1[0])
    assert phase_is_antipodal(m8, interval2[1], interval1[1])

    # The same antipodal geometry occurs on one actual deepest-comb side
    # whenever the scale ratio A=c_*/c_j is divisible by 16.  This checks
    # geometry only; activation by a canonical lower Boolean gate is not
    # asserted.
    scale_ratio = 32
    theta1 = Fraction(1, 14 * scale_ratio)
    theta2 = Fraction(1, 16) + theta1
    assert Fraction(0) < theta1 < theta2 < Fraction(1, 14)
    assert theta2 - theta1 == Fraction(1, 2 * m8)
    assert phase_is_antipodal(m8, theta2, theta1)

    print("LRC14 OWNER-SIDE ACUTE-CURRENT AUDIT")
    print("unsplit owner-cone guaranteed multipliers:", unsplit_nonzero)
    print("one-sided open-cone guaranteed multipliers:", side_nonzero)
    print("one-sided THM-2355 autocorrelation multipliers:",
          side_autocorrelation)
    print("strict two-cell cross-cone pairs m+n<=6:", len(cross_pairs))
    print("two-side cancellation outer-strip thresholds:",
          tuple(sorted(outer_thresholds.items())))
    print("first unsplit hostile m=4 centers:", -centre4, centre4)
    print("  both endpoint orientations are exact antipodal pairs: PASS")
    print("first one-sided cone hostile m=8 intervals:",
          interval1, interval2)
    print("  both endpoint orientations are exact antipodal pairs: PASS")
    print("deepest-comb geometric m=8 control theta:", theta1, theta2)
    print("  Boolean activation is intentionally not asserted")
    print("m=7 boundary uses strict open half-plane, no uniform margin")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
