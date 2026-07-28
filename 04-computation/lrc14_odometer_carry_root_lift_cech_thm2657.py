#!/usr/bin/env python3
"""Exact finite referee for THM-2657's nonsplit odometer clutch."""

from collections import Counter
from math import gcd


P = 13
R = P ** 6
S = R // P
C3 = 2 * S
W = (1, 14, 27, 40, 53, 66, 13, P ** 3, C3)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def additive_order(k):
    return R // gcd(k, R)


def main():
    lift_counts = Counter()
    order_counts = Counter()
    carry_root_checks = 0
    for k in range(R):
        delta = 2 * k % P
        carry = k % P
        require(carry == 7 * delta % P,
                "carry/root slope-seven equivalence failed")
        lift_counts[delta] += 1
        order_counts[(delta == 0, additive_order(k))] += 1
        carry_root_checks += 1

    require(set(lift_counts.values()) == {S}
            and sum(lift_counts.values()) == R,
            "lift fibres changed size")
    require(all(additive_order(k) == R
                for k in range(R) if k % P),
            "a nonzero carry lift lost full odometer order")

    torsion13 = tuple(k for k in range(R) if 13 * k % R == 0)
    require(torsion13 == tuple(j * S for j in range(P)),
            "13-torsion subgroup changed")
    require(all(k % P == 0 for k in torsion13),
            "the extension acquired a splitting element")

    cocycle = []
    section_checks = 0
    for left in range(P):
        for right in range(P):
            total = left + right
            reduced = total % P
            omega_k = 7 * left + 7 * right - 7 * reduced
            expected = 91 * (total // P)
            require(omega_k == expected and omega_k % P == 0,
                    "minimal-section wrap cocycle changed")
            cocycle.append(omega_k)
            section_checks += 1
    cocycle_hist = tuple(sorted(Counter(cocycle).items()))
    require(cocycle_hist == ((0, 91), (91, 78)),
            "wrap-cocycle census changed")
    require((91 // P) % P == 7,
            "extension class changed in kernel/13-kernel")

    phase_grid = tuple(v for v in W if (7 * v) % S == 0)
    require(phase_grid == (C3,),
            "canonical lawful-phase speed census changed")

    print("THM-2657 exact odometer carry/root lift certificate")
    print(f"p={P} R={R} kernel={S} total_lifts={R}")
    print(f"lift_fibre_size={S} nonzero_lifts={(P-1)*S}")
    print(f"carry_root_checks={carry_root_checks}; all nonzero lifts order={R}")
    print(f"13_torsion={torsion13}; quotient_nonzero_intersection=0")
    print(f"minimal_section=7*delta; cocycle_k_hist={cocycle_hist}")
    print("physical_wrap=7/13^5; H2_class_mod13=7")
    print(f"section_checks={section_checks}")
    print(f"canonical_phase_grid_speeds={phase_grid}")
    print("nonsplit odometer extension: PASS")


if __name__ == "__main__":
    main()
