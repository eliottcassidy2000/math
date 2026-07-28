#!/usr/bin/env python3
"""Exact companion for THM-2701's semantic-word nilpotence.

This script is dependency-free.  It verifies the two factor-sparse
implications used in the proof, classifies all 64 length-six target words by
their first obstruction, checks a strict five-state ``bbbbb`` witness, and
exhausts the period-four lattice relevant to the nearby phase-scheduled BABA
carrier.

The literal terminal words are THM-2305's singleton strata on the canonical
typed row.  In particular they retain the *unshifted* source-safe factor and
the guard-safe factor.  They are not THM-2623 guard-cospan sectors and are not
source-deleted delayed tooth words.
"""

from fractions import Fraction
from itertools import product


P = 13
GUARD = 1
UNITS = (14, 27, 40, 53, 66)
C1 = P
CA = P**3
CB = 2 * P**5
ROW = (GUARD,) + UNITS + (C1, CA, CB)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered(value):
    """Half-open representative in [-1/2,1/2)."""
    residue = value % 1
    return residue if residue < Fraction(1, 2) else residue - 1


def danger(speed, value, denominator=14):
    phase = centered(speed * value)
    radius = Fraction(1, denominator)
    return -radius <= phase < radius


def literal_word(value, target):
    """Membership in literal Q_a=E_a or Q_b=E_b."""
    require(target in ("a", "b"), "unknown target letter")
    if danger(GUARD, value, 7):
        return False
    if any(danger(speed, value) for speed in UNITS):
        return False
    if danger(C1, value):
        return False
    a_hit = danger(CA, value)
    b_hit = danger(CB, value)
    return (a_hit and not b_hit) if target == "a" else (b_hit and not a_hit)


def target_and_unit_word(value, target):
    """The scheduled tooth word, with guard and source-safe bits omitted."""
    require(target in ("a", "b"), "unknown target letter")
    if any(danger(speed, value) for speed in UNITS):
        return False
    a_hit = danger(CA, value)
    b_hit = danger(CB, value)
    return (a_hit and not b_hit) if target == "a" else (b_hit and not a_hit)


def B(value):
    return (P * value) % 1


def verify_factor_sparse_proof():
    # Qa at state j contains D_(13^3), while every literal terminal word at
    # state j+3 contains guard safety.  The former phase is literally the
    # latter base coordinate and 1/14 < 1/7.
    require(CA == P**3, "target-a/three-step identity changed")
    require(Fraction(1, 14) < Fraction(1, 7),
            "narrow target tooth stopped lying in the guard tooth")

    # For the target-b implication, put z=B^5 y=c+e with
    # c in {0,1/2}.  The first b tooth gives |e|<1/28.  On this range the
    # second b tooth has no nonzero wrap, so |26e|<1/14 and
    # |e|<1/(28*13).  Then |14e|<1/26<1/14, making speed 14 dangerous at
    # state five.  These exact endpoint comparisons exclude both wrap sides.
    require(CB == 2 * P**5, "target-b/five-step identity changed")
    require(26 * Fraction(1, 28) == Fraction(13, 14),
            "target-b wrap threshold changed")
    require(Fraction(1, 26) < Fraction(1, 14),
            "target-b contraction stopped forcing speed-14 danger")

    # Every length-six word is certified without interval enumeration.
    certificates = {}
    for word_tuple in product("ab", repeat=6):
        word = "".join(word_tuple)
        early_a = next((j for j in range(3) if word[j] == "a"), None)
        if early_a is not None:
            certificates[word] = f"a@{early_a}->guard@{early_a + 3}"
        else:
            require(word[:3] == "bbb", "six-word dichotomy failed")
            certificates[word] = "b@0,1->speed14@5"
    require(len(certificates) == 64, "length-six word census changed")
    counts = {
        label: sum(value == label for value in certificates.values())
        for label in sorted(set(certificates.values()))
    }
    require(counts == {
        "a@0->guard@3": 32,
        "a@1->guard@4": 16,
        "a@2->guard@5": 8,
        "b@0,1->speed14@5": 8,
    }, "length-six obstruction partition changed")
    return counts


def verify_sharp_witness():
    y0 = Fraction(132_661, CB)
    states = []
    margins = []
    value = y0
    for step in range(5):
        require(literal_word(value, "b"),
                f"sharp bbbbb witness failed at state {step}")
        distances = []
        for speed in ROW:
            phase = abs(centered(speed * value))
            threshold = Fraction(1, 7) if speed == GUARD else Fraction(1, 14)
            if speed != CB:
                distances.append(phase - threshold)
        require(min(distances) > 0,
                f"sharp witness acquired a safe-factor boundary at {step}")
        require(centered(CB * value) == 0,
                f"sharp witness left the target-b tooth centre at {step}")
        states.append(value)
        margins.append(min(distances))
        value = B(value)

    require(value == Fraction(1, 2), "sharp witness endpoint changed")
    require(danger(14, value), "state-five speed 14 is no longer dangerous")
    require(not literal_word(value, "a") and not literal_word(value, "b"),
            "sharp witness unexpectedly extended to a sixth semantic word")
    expected_states = (
        Fraction(132_661, 742_586),
        Fraction(18_417, 57_122),
        Fraction(841, 4_394),
        Fraction(165, 338),
        Fraction(9, 26),
    )
    expected_margins = (
        Fraction(186_041, 5_198_102),
        Fraction(3_314, 199_927),
        Fraction(1_493, 30_758),
        Fraction(66, 1_183),
        Fraction(15, 182),
    )
    require(tuple(states) == expected_states, "sharp witness orbit changed")
    require(tuple(margins) == expected_margins,
            "sharp witness safe margins changed")
    return y0, tuple(states), tuple(margins), value


def integer_danger(speed, numerator, modulus, denominator=14):
    """Exact half-open danger test for numerator/modulus."""
    residue = speed * numerator % modulus
    return (denominator * residue < modulus
            or denominator * residue >= (denominator - 1) * modulus)


def integer_target_and_unit(numerator, modulus, target):
    if any(integer_danger(speed, numerator, modulus) for speed in UNITS):
        return False
    a_hit = integer_danger(CA, numerator, modulus)
    b_hit = integer_danger(CB, numerator, modulus)
    return (a_hit and not b_hit) if target == "a" else (b_hit and not a_hit)


def integer_literal_word(numerator, modulus, target):
    return (
        integer_target_and_unit(numerator, modulus, target)
        and not integer_danger(GUARD, numerator, modulus, 7)
        and not integer_danger(C1, numerator, modulus)
    )


def verify_period_four_boundary():
    modulus = P**4 - 1
    rows = {}
    for schedule in ("baba", "abab"):
        tooth_only = []
        semantic = []
        for numerator in range(modulus):
            orbit = tuple(
                pow(P, step, modulus) * numerator % modulus
                for step in range(4)
            )
            if all(integer_target_and_unit(point, modulus, target)
                   for point, target in zip(orbit, schedule)):
                tooth_only.append(numerator)
            if all(integer_literal_word(point, modulus, target)
                   for point, target in zip(orbit, schedule)):
                semantic.append(numerator)
        rows[schedule] = (tuple(tooth_only), tuple(semantic))

    require(rows == {
        "baba": ((1_176, 27_384), ()),
        "abab": ((13_272, 15_288), ()),
    }, "period-four target lattice changed")

    orbit = tuple(
        Fraction(pow(P, step, modulus) * 1_176 % modulus, modulus)
        for step in range(4)
    )
    require(orbit == (
        Fraction(7, 170), Fraction(91, 170),
        Fraction(163, 170), Fraction(79, 170),
    ), "BABA orbit changed")
    defects = []
    for value, target in zip(orbit, "baba"):
        if target == "b":
            require(danger(GUARD, value, 7),
                    "B vertex lost its guard defect")
            defects.append("guard-danger")
        else:
            require(danger(C1, value),
                    "A vertex lost its unshifted-source defect")
            defects.append("unshifted-c1-danger")
    require(tuple(defects) == (
        "guard-danger", "unshifted-c1-danger",
        "guard-danger", "unshifted-c1-danger",
    ), "BABA defect schedule changed")
    return modulus, rows, orbit, tuple(defects)


def main():
    certificates = verify_factor_sparse_proof()
    y0, states, margins, endpoint = verify_sharp_witness()
    modulus, period_rows, orbit, defects = verify_period_four_boundary()

    print("LRC14 literal singleton-word one-step dilation nilpotence")
    print("status=PROVED + VERIFIED-EXACT; canonical typed row")
    print(f"row={ROW}; B(y)={{13y}}; Qa=E_c2; Qb=E_c3")
    print("literal_word_bits=guard-safe + five unit-safes + "
          "unshifted-c1-safe + exactly one target danger")
    print(f"length6_certificate_partition={certificates}")
    print("factor_sparse_a=Qa@j forces guard-danger@j+3")
    print("factor_sparse_b=Qb@0,Qb@1 force speed14-danger@5")
    print(f"sharp_word=bbbbb:y0={y0}:states={states}:"
          f"min_safe_margins={margins}:B5_endpoint={endpoint}:"
          "speed14_distance=0")
    print(f"period4_modulus={modulus}:rows={period_rows}:"
          f"BABA_orbit={orbit}:defects={defects}")
    print("verdict=literal Qa/Qb language is positive at five states and "
          "empty at six; no one-step semantic recurrent SCC")
    print("exit=enlarge the state by guard/source-clock debt or change the "
          "chronology; phase-scheduled tooth switching is not a terminal-word "
          "transition")
    print("SCOPE: canonical typed row; no scalar-row exclusion and no LRC14 "
          "conclusion")


if __name__ == "__main__":
    main()
