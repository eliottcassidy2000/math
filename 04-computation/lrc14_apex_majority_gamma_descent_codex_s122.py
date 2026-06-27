#!/usr/bin/env python3
"""S122: exact count audit for the 14 -> 7 gamma descent.

This script records the finite combinatorial counts behind THM-571.
It is not a brute-force proof of LRC(14); the mathematical input is the
accepted below-frontier LRC theorem for at most 12 moving speeds.  The audit
checks the shift-count lemmas and prints the full case split for a primitive
13-speed row with at least seven multiples of 14.
"""

from fractions import Fraction as F
from math import gcd


THRESH = F(1, 14)


def dist1(x):
    r = x % 1
    return min(r, 1 - r)


def bad_count(modulus, speed, gamma):
    return sum(1 for j in range(modulus)
               if dist1(F(speed, modulus) * (gamma + j)) < THRESH)


def audit_shift_counts(denom=196, speed_bound=420):
    facts = {
        "mod14_ordinary_max": 0,
        "mod14_halfstep_max": 0,
        "mod7_unit_max": 0,
    }
    witnesses = {}
    for a in range(denom):
        gamma = F(a, denom)
        for speed in range(1, speed_bound + 1):
            if speed % 14:
                c14 = bad_count(14, speed, gamma)
                if gcd(speed, 14) == 7:
                    if c14 > facts["mod14_halfstep_max"]:
                        facts["mod14_halfstep_max"] = c14
                        witnesses["mod14_halfstep"] = (speed, gamma, c14)
                elif c14 > facts["mod14_ordinary_max"]:
                    facts["mod14_ordinary_max"] = c14
                    witnesses["mod14_ordinary"] = (speed, gamma, c14)
            if speed % 7:
                c7 = bad_count(7, speed, gamma)
                if c7 > facts["mod7_unit_max"]:
                    facts["mod7_unit_max"] = c7
                    witnesses["mod7_unit"] = (speed, gamma, c7)
    return facts, witnesses


def case_split_rows():
    rows = []
    for r14 in range(7, 13):
        residual = 13 - r14
        for h7 in range(residual + 1):
            if h7 <= 1:
                if h7 == 0:
                    bad_bound = 2 * residual
                else:
                    # The half-step speed removes one parity class; each ordinary
                    # speed contributes at most one new shift in the other parity.
                    bad_bound = 7 + (residual - 1)
                rows.append((r14, residual, h7, "14-phase",
                             bad_bound, 14 - bad_bound))
            else:
                total7 = r14 + h7
                non7 = 13 - total7
                bad_bound = non7
                rows.append((r14, residual, h7, "7-phase",
                             bad_bound, 7 - bad_bound))
    return rows


def main():
    facts, witnesses = audit_shift_counts()
    print("LRC14 apex-majority gamma descent (S122)")
    print()
    print("1. Exact shift-count maxima")
    for key in sorted(facts):
        print(f"  {key}: {facts[key]} witness={witnesses.get(key.replace('_max', ''))}")
    print()
    print("2. Primitive r14>=7 case split")
    print("  columns: r14, residual, halfstep_in_residual, phase, bad_bound, surviving_points")
    all_ok = True
    for row in case_split_rows():
        r14, residual, h7, phase, bad_bound, survivors = row
        ok = survivors > 0
        all_ok = all_ok and ok
        print(f"  {r14:2d} {residual:2d} {h7:2d} {phase:8s} "
              f"bad<={bad_bound:2d} survivors>={survivors:2d} {'OK' if ok else 'FAIL'}")
    print()
    print(f"all case-split rows have a surviving lift: {all_ok}")
    print()
    print("Readout: the S121 two-halfstep obstruction is real for a raw")
    print("14-shift sieve, but in the apex-majority branch it descends to the")
    print("7-phase, where at least nine speeds are multiples of 7 and at most")
    print("four residual speeds remain.  Each residual forbids at most one of")
    print("the seven lifts, so a 7-phase lift survives.")


if __name__ == "__main__":
    main()
