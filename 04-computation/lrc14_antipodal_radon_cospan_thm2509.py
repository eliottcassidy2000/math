#!/usr/bin/env python3
"""Exact consequence checks for the THM-2509 antipodal Radon cospan."""

from fractions import Fraction
from itertools import combinations


P = 13
R = 7
SLOPES = tuple(range(1, P))
PAIRS = tuple((tau, (-tau) % P) for tau in range(1, 7))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def pair_count(slopes: set[int]) -> int:
    return sum(tau in slopes and minus in slopes for tau, minus in PAIRS)


def pi(tau: int, h: int, r: int) -> tuple[int, int]:
    return ((h + tau * r) % P, (h - tau * r) % P)


def inverse_pi(tau: int, u: int, v: int) -> tuple[int, int]:
    half = pow(2, -1, P)
    h = half * (u + v) % P
    representative = half * pow(tau, -1, P) * (u - v) % P
    require(0 <= representative < R, "point is outside the seven-strip image")
    return h, representative


def main() -> None:
    # Every seven-set among six antipodal pairs contains a complete pair.
    counts = {1: 0, 2: 0, 3: 0}
    total = 0
    for subset in combinations(SLOPES, 7):
        complete = pair_count(set(subset))
        require(complete >= 1, "seven-set missed every antipodal pair")
        counts[complete] += 1
        total += 1
    require(total == 792 and counts == {1: 192, 2: 480, 3: 120},
            "seven-set antipodal census drifted")

    # Pi_tau is a bijection from F_13 x {0,...,6} to a 91-point strip.
    injection_rows = 0
    translation_rows = 0
    for tau in SLOPES:
        image = {pi(tau, h, r) for h in range(P) for r in range(R)}
        require(len(image) == P * R, "antipodal chart is not injective")
        for h in range(P):
            for r in range(R):
                u, v = pi(tau, h, r)
                require(inverse_pi(tau, u, v) == (h, r), "inverse chart failed")
                require(pi((-tau) % P, h, r) == (v, u), "sign is not leg swap")
                injection_rows += 1
                for A in range(P):
                    for C in range(R):
                        moved = pi(tau, (h + A) % P, (r + C) % R)
                        # Conjugating physical CRT translation by Pi retains
                        # the representative carry exactly and is a permutation
                        # of the same 91-point strip.
                        h2, r2 = inverse_pi(tau, *moved)
                        require((h2, r2) == ((h + A) % P, (r + C) % R),
                                "conjugated CRT action lost the carry")
                        translation_rows += 1

    # Imported exact zero-slope histogram from the complete THM-2436 atlas
    # companion lrc14_truncated_radon_atlas_thm2507.cpp/.out.
    bad_counts = {
        1: 112, 2: 28, 3: 16, 4: 77, 5: 51, 6: 93,
        7: 106, 8: 49, 9: 42, 10: 85, 11: 71, 12: 103,
    }
    nonflat = 39643
    require(sum(bad_counts.values()) == 833, "atlas bad-slope total drifted")
    require(nonflat - sum(bad_counts.values()) == 38810,
            "atlas all-slope assignment count drifted")
    fixed_pair_counts = {
        (tau, minus): nonflat - bad_counts[tau] - bad_counts[minus]
        for tau, minus in PAIRS
    }
    require(min(fixed_pair_counts.values()) == 39428,
            "atlas fixed-pair minimum drifted")
    require(max(fixed_pair_counts.values()) == 39544,
            "atlas fixed-pair maximum drifted")

    abstract_pair_fraction = Fraction(1, 6)
    atlas_pair_fraction = Fraction(5, 6)
    k1_parent_floor = atlas_pair_fraction * Fraction(2, 7)
    k2_parent_floor = atlas_pair_fraction * Fraction(3, 7)
    require(k1_parent_floor == Fraction(5, 21), "k=1 floor drifted")
    require(k2_parent_floor == Fraction(5, 14), "k=2 floor drifted")

    print("THM-2509 antipodal Radon cospan exact companion: PASS")
    print(f"seven_slope_sets={total}; complete_pair_histogram={counts}")
    print(f"lossless_chart_rows={injection_rows}; crt_conjugacy_rows={translation_rows}")
    print("fixed_atlas_pair_assignment_counts=" + ",".join(
        f"{tau}/{minus}:{fixed_pair_counts[(tau, minus)]}"
        for tau, minus in PAIRS
    ))
    print(f"abstract_fixed_pair_fraction={abstract_pair_fraction}")
    print(f"atlas_fixed_pair_fraction={atlas_pair_fraction}")
    print(f"atlas_parent_mass_floors=k1:{k1_parent_floor},k2:{k2_parent_floor}")
    print("tau sign is only a leg-swap gauge; the unordered cospan is intrinsic")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
