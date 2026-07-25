#!/usr/bin/env python3
"""Exact arithmetic audit for THM-2282's repeated-owner BV handoff."""

from fractions import Fraction


DELTA_5 = Fraction(961, 6930)
REPEATED_EXCLUSIVE = Fraction(5229541, 197927730)
OWNER_FLOOR = REPEATED_EXCLUSIVE / 3
TARGET_FLOOR = Fraction(2593, 90090)
HORIZON_FACTOR = 4 / OWNER_FLOOR
HANDOFF_FLOOR = OWNER_FLOOR * TARGET_FLOOR / 2
EXPIRATION_IMAGE_FLOOR = Fraction(169, 20) * OWNER_FLOOR


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def guard_blocker_cap(depth):
    """THM-2273's parity-sharp cap for mu(C_H intersect D_q)."""
    power = 13**depth
    if depth % 2:
        return Fraction(5, 49) + Fraction(5, 49 * power)
    return Fraction(5, 49) + Fraction(5, 294 * power)


def main():
    profiles = [(1, 1, deepest) for deepest in range(5, 20)]
    require(len(profiles) == 15, "repeated-profile census drift")
    require(
        OWNER_FLOOR == Fraction(5229541, 593783190),
        "selected-owner floor drift",
    )

    # The selected exclusive owner can be either shallow label (depth one)
    # or the deep label (depth c).  Depth one is the uniform worst target
    # debit over every possible selected label in all repeated rows.
    selected_depths = {
        depth
        for _, _, deepest in profiles
        for depth in (1, deepest)
    }
    caps = {depth: guard_blocker_cap(depth) for depth in selected_depths}
    require(max(caps, key=caps.get) == 1, "target-cap maximizing depth drift")
    require(caps[1] == Fraction(10, 91), "depth-one target cap drift")
    require(
        DELTA_5 - caps[1] == TARGET_FLOOR,
        "blocker-only target floor drift",
    )
    require(
        all(caps[depth] <= caps[1] for depth in selected_depths),
        "deep selected-owner cap exceeded depth-one cap",
    )

    # Var(1_E)<=2S and 13^k >= (4/e0)S imply that the Perron error
    # 2S/13^k is at most e0/2.
    require(
        HORIZON_FACTOR == Fraction(2375132760, 5229541),
        "BV horizon factor drift",
    )
    require(
        HORIZON_FACTOR * OWNER_FLOOR == 4,
        "BV threshold normalization drift",
    )
    require(
        HANDOFF_FLOOR
        == Fraction(13560199813, 106987855174200),
        "positive handoff floor drift",
    )
    require(HANDOFF_FLOOR > 0, "handoff floor lost positivity")

    # This is the same selected-owner floor at its prescribed expiration.
    # It crosses one oriented half-comb, but not one full danger comb; the
    # BV theorem deliberately makes no fixed-expiration claim.
    require(
        EXPIRATION_IMAGE_FLOOR == Fraction(5229541, 70270200),
        "expiration image floor drift",
    )
    require(
        EXPIRATION_IMAGE_FLOOR - Fraction(1, 14)
        == Fraction(210241, 70270200),
        "half-comb margin drift",
    )
    require(
        EXPIRATION_IMAGE_FLOOR < Fraction(1, 7),
        "unexpected full-comb crossing",
    )

    print("THM-2282 REPEATED-OWNER BV DELAYED HANDOFF -- exact audit")
    print(f"profiles: {len(profiles)}; deepest depths: 5..19")
    print(f"selected-owner depth set: {sorted(selected_depths)}")
    print(f"exclusive total/owner floor: {REPEATED_EXCLUSIVE} {OWNER_FLOOR}")
    print(f"target floor: {TARGET_FLOOR}")
    print(f"horizon factor: {HORIZON_FACTOR}")
    print(f"handoff floor: {HANDOFF_FLOOR}")
    print(f"expiration image/half-comb margin: {EXPIRATION_IMAGE_FLOOR} "
          f"{EXPIRATION_IMAGE_FLOOR - Fraction(1, 14)}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
