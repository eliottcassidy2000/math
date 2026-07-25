#!/usr/bin/env python3
"""Exact arithmetic audit for THM-2271."""

from fractions import Fraction


def main():
    strict_profiles = sum(c - 2 for c in range(5, 20))
    repeated_profiles = len(range(5, 20))
    expansion = Fraction(169, 20)
    one_comb = Fraction(1, 7)

    strict_owner = Fraction(15041431, 593783190)
    strict_image = strict_owner * expansion
    strict_cut = strict_image - one_comb

    repeated_owner = Fraction(5229541, 593783190)
    repeated_image = repeated_owner * expansion
    repeated_gap = repeated_image - one_comb

    assert strict_profiles == 150
    assert repeated_profiles == 15
    assert strict_image == Fraction(15041431, 70270200)
    assert strict_cut == Fraction(5002831, 70270200)
    assert strict_cut > 0
    assert repeated_image == Fraction(5229541, 70270200)
    assert repeated_image > Fraction(1, 14)
    assert repeated_gap == -Fraction(4809059, 70270200)
    assert repeated_gap < 0

    print("THM-2271 exact arithmetic audit")
    print("strict profiles:", strict_profiles)
    print("strict labelled owner floor:", strict_owner)
    print("strict expiration-image floor:", strict_image)
    print("one-comb capacity:", one_comb)
    print("strict owner-absorber cut floor:", strict_cut)
    print("repeated-first profiles:", repeated_profiles)
    print("repeated expiration-image floor:", repeated_image)
    print("repeated image minus one-comb capacity:", repeated_gap)
    print("all exact assertions passed")


if __name__ == "__main__":
    main()
