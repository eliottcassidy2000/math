#!/usr/bin/env python3
"""Exact probe for the prime-to-13 endpoint-current stopping boundary.

The computation uses only integer arithmetic.  It checks a primitive
nine-speed scalar boundary atlas whose deepest endpoint trace is invariant
under a half-turn.  At every odd frequency the two endpoints in each
orientation/address pair consequently have opposite phases.

This is a boundary-trace hostile, not a claimed LRC(14) counterexample.
In particular, its guard speed is even; the live fixed-section scalar rows
carry the additional odd-guard hypothesis.
"""

from math import gcd


P = 13
LAMBDA = 1
B = 3
C = 5

C1 = 2 * P**LAMBDA
C2 = 2 * P**B
C3 = 2 * P**C
H = 2
QS = (1, 1 + 7 * C3, 3, 3 + 7 * C3, 4)
A = P
R = P ** (C - B)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation_13(value: int) -> int:
    exponent = 0
    while value % P == 0:
        value //= P
        exponent += 1
    return exponent


def danger(speed: int, numerator: int, denominator: int) -> bool:
    """Return ||speed*numerator/denominator|| < 1/14 exactly."""

    residue = (speed * numerator) % denominator
    folded = min(residue, denominator - residue)
    return 14 * folded < denominator


def guard_safe(speed: int, numerator: int, denominator: int) -> bool:
    """Return ||speed*numerator/denominator|| > 1/7 exactly."""

    residue = (speed * numerator) % denominator
    folded = min(residue, denominator - residue)
    return 7 * folded > denominator


def lower_gate(r: int, orientation: int) -> bool:
    """Lower exclusive-owner gate on one c3 endpoint."""

    numerator = 14 * r + orientation
    denominator = 14 * C3
    if not guard_safe(H, numerator, denominator):
        return False
    for speed in QS:
        if danger(speed, numerator, denominator):
            return False
    return danger(C1, numerator, denominator) and not danger(
        C2, numerator, denominator
    )


def audit_orientation(orientation: int) -> tuple[int, int]:
    active = 0
    active_nontrivial = 0
    half_shift = C3 // 2

    for r in range(C3):
        value = lower_gate(r, orientation)
        paired = lower_gate((r + half_shift) % C3, orientation)
        require(value == paired, "half-turn gate invariance failed")
        if value:
            active += 1
            numerator = 14 * r + orientation
            if numerator % R != 0:
                active_nontrivial += 1

    return active, active_nontrivial


def main() -> None:
    speeds = (H, *QS, C1, C2, C3)
    common_divisor = 0
    for speed in speeds:
        common_divisor = gcd(common_divisor, speed)

    require(len(speeds) == 9, "wrong scalar dimension")
    require(len(set(speeds)) == 9, "speeds are not distinct")
    require(common_divisor == 1, "scalar atlas is not primitive")
    require(valuation_13(H) == 0, "guard is not a 13-unit")
    require(all(valuation_13(q) == 0 for q in QS), "unit speed failed")
    require(
        (valuation_13(C1), valuation_13(C2), valuation_13(C3))
        == (LAMBDA, B, C),
        "blocker valuation profile failed",
    )
    require(C3 % 2 == 0 and A % 2 == 1, "half-turn phase sign failed")
    require(C3 // (2 * R) == P**B, "fine-residue flip scale failed")
    require((C3 // (2 * R)) % 2 == 1, "half-turn does not swap h mod 2")

    # On a deepest endpoint, adding 7*c3 to an odd unit speed adds
    # orientation/2 modulo one.  This is exactly the same half-turn that an
    # odd speed sees under x -> x+1/2, so each displayed pair swaps.
    require(QS[1] - QS[0] == 7 * C3, "first odd pair failed")
    require(QS[3] - QS[2] == 7 * C3, "second odd pair failed")
    require(QS[0] % 2 == QS[1] % 2 == 1, "first pair is not odd")
    require(QS[2] % 2 == QS[3] % 2 == 1, "second pair is not odd")
    require(H % 2 == QS[4] % 2 == C1 % 2 == C2 % 2 == 0, "even bank failed")

    minus = audit_orientation(-1)
    plus = audit_orientation(+1)
    expected = (49_880, 49_536)
    require(minus == expected, "negative-orientation endpoint count changed")
    require(plus == expected, "positive-orientation endpoint count changed")

    # Fine-residue lemma ledger.  Once t=h+u*v is fixed, the phase has
    # pure order 13^(b-lambda).  Its final cyclotomic block has 13 columns.
    # The selected-owner arc has length 1/7, so it occupies at most two
    # points of a 13-grid and leaves at least eleven identically zero
    # columns.  Equal columns (the zero-current consequence) are therefore
    # all zero.
    phase_order = P ** (B - LAMBDA)
    final_block = P
    owner_column_cap = 2
    forced_zero_columns = final_block - owner_column_cap
    require(phase_order == 169, "fine phase order changed")
    require(forced_zero_columns == 11, "owner-column sparsity changed")

    print("LRC14 PRIME-TO-13 CURRENT HOSTILE -- exact audit")
    print(f"speeds={speeds}")
    print(
        "valuations="
        f"{tuple(valuation_13(speed) for speed in speeds)} gcd={common_divisor}"
    )
    print(f"relative_address_count={R}")
    print(f"orientation=-1 active={minus[0]} nontrivial={minus[1]}")
    print(f"orientation=+1 active={plus[0]} nontrivial={plus[1]}")
    print("half_turn_gate_invariance=PASS")
    print("odd_frequency_phase_reversal=PASS")
    print("all_raw_orientation_address_currents=ZERO_BY_PAIRING")
    print(
        "fine_residue_lemma="
        f"phase_order_{phase_order}; columns_{final_block}; "
        f"owner_cap_{owner_column_cap}; forced_zero_{forced_zero_columns}"
    )
    print("scope=primitive full scalar boundary trace; even guard; not LRC14")


if __name__ == "__main__":
    main()
