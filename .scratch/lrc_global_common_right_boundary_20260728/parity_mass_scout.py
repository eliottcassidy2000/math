#!/usr/bin/env python3
"""Parity-class mass scout for the THM-2825 bank."""

from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
H = copies.T // (2 * P**5)
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cell_objects(details, full_module, e3, clocks, clock, sigma, target):
    section = copies.physical.target.source_present_section(
        full_module, e3, clock, sigma, target, clocks
    )
    source_base, target_base = details[clock]
    source = copies.weighted_intersection(source_base, section)
    target = copies.weighted_intersection(target_base, section)
    pulled = copies.physical.overlap.shift_weighted(target, -copies.SHIFT)
    aligned = copies.physical.overlap.intersect_weighted_profiles(
        source, pulled
    )
    require(
        all(a == b for _left, _right, a, b in aligned),
        "unequal common weight",
    )
    common = tuple((left, right, a) for left, right, a, _b in aligned)
    right = copies.physical.subtract_weighted(pulled, common)
    return common, right


def parity_masses(pieces, origin):
    masses = Counter()
    for left, _right, _weight in pieces:
        difference = left - origin
        require(difference % H == 0, "piece left the h-lattice")
        masses[(difference // H) % 2] += 1
    return masses[0], masses[1]


def main():
    (
        _module, _rails, _present, details, full_module, e3, clocks,
        _q_pairs, _delayed, _sw, _tw, _rail,
    ) = copies.physical_setup()

    pair_types = Counter()
    zero_mod_p = []
    signed_zero_mod_p = []
    cells = 0
    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                common, right = cell_objects(
                    details, full_module, e3, clocks, clock, sigma, target
                )
                require(bool(common) == bool(right), "support shadow split")
                if not right:
                    continue
                cells += 1
                origin = min(
                    left for left, _stop, _weight in (*common, *right)
                )
                m_mass = parity_masses(common, origin)
                r_mass = parity_masses(right, origin)
                pair_types[(clock, m_mass, r_mass)] += 1
                for name, masses in (("M", m_mass), ("R", r_mass)):
                    for parity, mass in enumerate(masses):
                        if mass % P == 0:
                            zero_mod_p.append(
                                ((clock, sigma, target), name, parity, mass)
                            )
                    if (masses[0] - masses[1]) % P == 0:
                        signed_zero_mod_p.append(
                            ((clock, sigma, target), name, masses)
                        )

    print("LRC THM-2825 PARITY-CLASS MASS SCOUT")
    print(f"cells={cells}")
    print(f"pair_types={tuple(sorted(pair_types.items()))}")
    print(f"parity_class_zero_mod13={tuple(zero_mod_p)}")
    print(
        "signed_augmentation_zero_mod13="
        f"{tuple(signed_zero_mod_p)}"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
