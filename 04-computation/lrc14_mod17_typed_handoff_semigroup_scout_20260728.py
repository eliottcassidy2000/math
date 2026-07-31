#!/usr/bin/env python3
"""Exact mod-17 audit of the currently typed LRC14 affine handoffs.

This scout separates chronological circle maps from the THM-2640/2657
carry-root deck translations.  It also identifies the least denominator-17
affine phase which makes the old THM-789 point 4/17 recurrent for the full
raw THM-2693 delayed word.

The calculation is standalone and uses only exact integer/Fraction arithmetic.
It proves no packet typing, endpoint transport, row exclusion, or LRC(14).
"""

from fractions import Fraction
from itertools import product
from math import gcd


P = 13
R = P**6
MOD = 17
RHO = Fraction(1, 14)
TARGET = P**3
GUARDS = (14, 27, 40, 53, 66, 2 * P**5)
SPEEDS = (TARGET,) + GUARDS


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def centered(value):
    """Representative in [-1/2,1/2)."""
    return frac(value + Fraction(1, 2)) - Fraction(1, 2)


def danger(speed, y):
    z = centered(speed * y)
    return -RHO <= z < RHO


def in_raw_word(y):
    return danger(TARGET, y) and all(not danger(v, y) for v in GUARDS)


def nearest_word_walls(y):
    """Exact component walls around a strict point of the raw word.

    Every factor changes status only at v*y = n+1/14 or n+13/14.  Five
    neighbouring n values around floor(v*y) contain the nearest walls.
    """
    lower_walls = []
    upper_walls = []
    for speed in SPEEDS:
        n0 = (speed * y).numerator // (speed * y).denominator
        walls = []
        for n in range(n0 - 2, n0 + 3):
            walls.append(Fraction(14 * n + 1, 14 * speed))
            walls.append(Fraction(14 * n + 13, 14 * speed))
        lower_walls.append(max(w for w in walls if w < y))
        upper_walls.append(min(w for w in walls if w > y))
    return max(lower_walls), min(upper_walls)


def delayed_affine(epsilon, phase, y):
    return frac(epsilon * P * y + phase)


def chamber_arrows(epsilon, phase_numerator, chamber):
    return tuple(
        (a, (epsilon * P * a + phase_numerator) % MOD)
        for a in chamber
        if (epsilon * P * a + phase_numerator) % MOD in chamber
    )


def main():
    # At denominator 17, the half-open danger tooth contains exactly 0,+/-1.
    danger_residues = tuple(
        a for a in range(MOD) if danger(1, Fraction(a, MOD))
    )
    require(danger_residues == (0, 1, 16),
            "denominator-17 danger residues changed")

    chamber = tuple(
        a for a in range(MOD) if in_raw_word(Fraction(a, MOD))
    )
    require(chamber == (4, 13), "raw W denominator-17 chamber changed")

    target_candidates = tuple(
        a for a in range(MOD)
        if (TARGET * a) % MOD in danger_residues
    )
    require(target_candidates == (0, 4, 13),
            "target residue candidates changed")

    factor_residues = {}
    for a in chamber:
        factor_residues[a] = tuple(centered(Fraction(v * a, MOD))
                                    for v in SPEEDS)
    require(tuple(abs(z) for z in factor_residues[4])
            == tuple(Fraction(a, MOD) for a in (1, 5, 6, 7, 8, 8, 2)),
            "4/17 factor clearance word changed")
    require(factor_residues[13]
            == tuple(-z for z in factor_residues[4]),
            "reflection stopped exchanging the two chamber points")

    signed_images = {
        epsilon: tuple((epsilon * P * a) % MOD for a in chamber)
        for epsilon in (1, -1)
    }
    require(signed_images == {1: (1, 16), -1: (16, 1)},
            "signed one-step residue images changed")
    outside = tuple(sorted(set(signed_images[1])))
    require(outside == (1, 16) and not set(outside).intersection(chamber),
            "one-step parity separation failed")
    require(all(not danger(TARGET, Fraction(a, MOD)) for a in outside),
            "the intermediate residues no longer fail the target tooth")

    powers = tuple(pow(P, e, MOD) for e in range(4))
    require(powers == (1, 13, 16, 4), "order-four residue orbit changed")
    require((P * -P) % MOD == 1 and pow(P, 4, MOD) == 1,
            "signed-two/unsigned-four return identities changed")

    phase_arrows = {}
    for epsilon, t in product((1, -1), range(MOD)):
        arrows = chamber_arrows(epsilon, t, chamber)
        if arrows:
            phase_arrows[epsilon, t] = arrows
    phase_set = tuple(sorted({t for _, t in phase_arrows}))
    require(phase_set == (3, 5, 12, 14),
            "internal-edge phase set changed")
    expected_arrows = {
        (1, 3): ((4, 4),),
        (1, 5): ((13, 4),),
        (1, 12): ((4, 13),),
        (1, 14): ((13, 13),),
        (-1, 3): ((13, 4),),
        (-1, 5): ((4, 4),),
        (-1, 12): ((13, 13),),
        (-1, 14): ((4, 13),),
    }
    require(phase_arrows == expected_arrows,
            "signed chamber-arrow table changed")

    # The symmetric slope-seven stalk derived in the companion probe selects
    # the quarter-torsion phase y=-1/4 modulo an odd transverse modulus q.
    # The first such modulus accepted by the raw delayed word is exactly 17.
    accepted_quarter_moduli = []
    for modulus in range(3, 100, 2):
        if gcd(modulus, P) != 1:
            continue
        numerator = (-pow(4, -1, modulus)) % modulus
        if in_raw_word(Fraction(numerator, modulus)):
            accepted_quarter_moduli.append((modulus, numerator))
    require(tuple(accepted_quarter_moduli)
            == ((17, 4), (43, 32), (47, 35), (51, 38),
                (63, 47), (71, 53), (95, 71)),
            "quarter-torsion transverse-modulus scan changed")

    # The smallest transverse phase fixes the old THM-789 hostile point.
    y0 = Fraction(4, MOD)
    theta = Fraction(3, MOD)
    require(delayed_affine(1, theta, y0) == y0,
            "the +13,3/17 delayed fixed point was lost")
    left, right = nearest_word_walls(y0)
    require(left < y0 < right and in_raw_word(y0),
            "4/17 is no longer strict in its raw-word component")
    left_margin = y0 - left
    right_margin = right - y0
    require((left, right)
            == (Fraction(2446165, 10396204),
                Fraction(2446177, 10396204)),
            "strict component around 4/17 changed")
    require((left_margin, right_margin)
            == (Fraction(11, 176735468),
                Fraction(193, 176735468)),
            "strict component margins changed")
    component_width = right - left
    require(component_width == Fraction(3, 2599051),
            "strict component width changed")

    # For every symbolic horizon H>=1, this is a positive raw-W corridor:
    # y=y0+u, -L/13^(H-1)<u<Rr/13^(H-1).  Every iterate is
    # y0+13^j*u and stays in the strict component for 0<=j<H.
    for horizon in range(1, 13):
        scale = P ** (horizon - 1)
        corridor = (y0 - left_margin / scale,
                    y0 + right_margin / scale)
        require(corridor[1] - corridor[0]
                == component_width / scale,
                "all-H corridor width formula failed")
        probes = (
            (corridor[0] + y0) / 2,
            y0,
            (y0 + corridor[1]) / 2,
        )
        for point in probes:
            orbit = point
            for _ in range(horizon):
                require(in_raw_word(orbit),
                        "finite-horizon corridor left the raw word")
                orbit = delayed_affine(1, theta, orbit)

    # Lift the delayed phase to an actual affine circle map.  It has an exact
    # fixed x over y0, but its translation is not a THM-2657 carry/root lift.
    beta = Fraction(3, MOD * R)
    inv12 = pow(12, -1, R)
    predecessor_integer = (-3 * inv12) % R
    x0 = Fraction(MOD * predecessor_integer + 4, MOD * R)
    require(frac(R * x0) == y0, "circle fixed point has wrong delayed phase")
    require(frac(P * x0 + beta) == x0,
            "transverse affine circle fixed point failed")
    carry = predecessor_integer % P
    future_digit = (P * y0).numerator // (P * y0).denominator
    require((carry, future_digit) == (3, 3),
            "local digit compatibility at the fixed point changed")
    require(R * beta == Fraction(3, MOD)
            and (R * beta).denominator == MOD,
            "transverse delayed phase denominator changed")
    require((2 * P**5) * beta == Fraction(6, P * MOD),
            "private-root probe phase changed")
    require(((2 * P**5) * beta * P).denominator == MOD,
            "transverse phase accidentally entered the F13 root grid")

    # A pure translation by beta has a nonconstant predecessor-integer jump:
    # zero below y=14/17 and one above it.  Hence THM-2657's global lift
    # classification cannot type this denominator-17 phase.
    jump_low = (Fraction(1, 2) + R * beta).numerator // (
        Fraction(1, 2) + R * beta).denominator
    jump_high = (Fraction(16, 17) + R * beta).numerator // (
        Fraction(16, 17) + R * beta).denominator
    require((jump_low, jump_high) == (0, 1),
            "pure transverse translation jump should be nonconstant")

    print("LRC14 mod-17 typed handoff semigroup scout")
    print("status=PROVED-ELEMENTARY projection/residue laws + FINITE-EXACT corridor")
    print(f"parameters=R:{R},target:{TARGET},guards:{GUARDS}")
    print(f"denominator17_danger_residues={danger_residues}")
    print(f"target_candidates={target_candidates}; raw_W17={chamber}")
    print(f"factor_centered_residues_at_4={factor_residues[4]}")
    print("typing=G_k:x->x+k/R gives y->y; slope-seven k=7delta is a deck/carry-root sidecar")
    print("typing=D_(eps,k):x->eps*13x+k/R gives y->eps*13y for integer k")
    print("typing=mixed one-D plus any G-word has delayed multiplier eps*13, not 7 or 18")
    print(f"signed_images_of_W17={signed_images}; intermediate_target_failure={outside}")
    print(f"mod17_powers_of_13={powers}; opposite_signed_two_step=1; unsigned_four_step=1")
    for key in sorted(phase_arrows):
        print(f"phase_arrows_eps{key[0]:+d}_t{key[1]}={phase_arrows[key]}")
    print(f"internal_edge_phase_numerators={phase_set}")
    print(f"quarter_torsion_raw_W_moduli_below100={tuple(accepted_quarter_moduli)}; first=17")
    print(f"fixed_delayed_map=y->13y+3/17; fixed_y={y0}")
    print(f"fixed_component=({left},{right}); margins=({left_margin},{right_margin})")
    print("all_H_corridor=(4/17-(11/176735468)/13^(H-1),4/17+(193/176735468)/13^(H-1))")
    print(f"all_H_width=({component_width})/13^(H-1)")
    print(f"circle_lift_beta={beta}; fixed_x={x0}; local_(carry,h)=({carry},{future_digit})")
    print(f"typing_debt=R*beta:{R * beta}; c3*beta:{(2 * P**5) * beta}; pure_jump_examples:{(jump_low, jump_high)}")
    print("scope=raw delayed word and typed-operation projection only; no globally typed carry/root lift, packet, endpoint, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
