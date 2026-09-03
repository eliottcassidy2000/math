#!/usr/bin/env python3
"""Exact audit for the sharpened THM-4349 adaptive owner clock.

This companion checks only arithmetic and literal finite controls.  The
general proof that the clock certifies universal tail safety is analytic and
is recorded in the matching reflection.
"""

from fractions import Fraction
from itertools import combinations
from math import lcm


H72 = (1, 5, 11, 13, 17, 19, 23, 37, 41, 70, 72)


def odd_tail_box(h: int) -> tuple[int, ...]:
    return tuple(range(1, 12 * h, 2))


def lcm_many(values) -> int:
    ans = 1
    for value in values:
        ans = lcm(ans, value)
    return ans


def body_modulus(H: tuple[int, ...]) -> int:
    h = max(H)
    return lcm_many(tuple(2 * c for c in H) + odd_tail_box(h))


def height_modulus(h: int) -> int:
    return lcm_many(tuple(2 * c for c in range(1, h + 1)) + odd_tail_box(h))


def abs_residue(value: int, modulus: int) -> int:
    r = value % modulus
    return min(r, modulus - r)


def danger(z: int, r: int, lift: int, clock: int) -> bool:
    return 14 * abs_residue(z * (r + lift * clock), 2 * clock) < 2 * clock


def eligible_at(z: int, r: int, clock: int) -> bool:
    return 7 * abs_residue(z * r, clock) < clock


def nearest_integer(numerator: int, denominator: int) -> int:
    lower, remainder = divmod(numerator, denominator)
    if 2 * remainder == denominator:
        raise RuntimeError(("nearest-integer tie", numerator, denominator))
    return lower if 2 * remainder < denominator else lower + 1


def v2(value: int) -> int:
    ans = 0
    while value % 2 == 0:
        value //= 2
        ans += 1
    return ans


def audit_height(h: int) -> tuple[int, int, int, int]:
    D = height_modulus(h)
    Q = 7 * D
    odd_lcm = lcm_many(odd_tail_box(h))
    expected_D = (2 ** (1 + (h.bit_length() - 1))) * odd_lcm
    if D != expected_D:
        raise RuntimeError(("height modulus normal form", h, D, expected_D))

    # Every possible doubled core wall and odd-tail wall is a Q-grid point.
    for c in range(1, h + 1):
        if Q % (14 * c) != 0:
            raise RuntimeError(("core wall denominator", h, c, Q))
    for a in odd_tail_box(h):
        if Q % (7 * a) != 0:
            raise RuntimeError(("tail wall denominator", h, a, Q))

    # Full closed safe arc gives the exact sufficient long-block threshold.
    threshold_numerator = 1008 * h * h + 84 * h
    if Q <= threshold_numerator:
        raise RuntimeError(("long-block threshold", h, Q, threshold_numerator))

    # The two largest odd tails are coprime, and D is even because it also
    # contains a doubled body speed.  This is the missing factor two that lets
    # the witness clock itself clear the long-block threshold.
    lower = 14 * (12 * h - 1) * (12 * h - 3)
    if D % (2 * (12 * h - 1) * (12 * h - 3)) != 0:
        raise RuntimeError(("even coprime odd divisibility", h, D))
    if Q < lower:
        raise RuntimeError(("coprime odd lower bound", h, Q, lower))
    if lower - threshold_numerator != 42 * (24 * h * h - 18 * h + 1):
        raise RuntimeError(("threshold identity", h))

    L = lcm_many(range(1, 12 * h))
    simple = 7 * L
    if simple % Q != 0:
        raise RuntimeError(("simple envelope divisibility", h, simple, Q))

    A = 42 * h + 1
    old = 28 * A * A * lcm(L, 12 * h)
    if old % Q != 0:
        raise RuntimeError(("old clock divisibility", h, old, Q))
    return D, Q, simple // Q, old // Q


def audit_small_bodies(max_h: int = 18) -> int:
    checks = 0
    for h in range(11, max_h + 1):
        D_height = height_modulus(h)
        for lower_body in combinations(range(1, h), 10):
            H = lower_body + (h,)
            D = body_modulus(H)
            Q = 7 * D
            odd_lcm = lcm_many(odd_tail_box(h))
            expected_D = (2 ** (1 + max(v2(c) for c in H))) * odd_lcm
            if D != expected_D:
                raise RuntimeError(("body modulus normal form", H, D, expected_D))
            if D_height % D != 0:
                raise RuntimeError(("body-height divisibility", H))
            if Q <= 1008 * h * h + 84 * h:
                raise RuntimeError(("body long-block threshold", H, Q))
            for c in H:
                if Q % (14 * c) != 0:
                    raise RuntimeError(("body core wall", H, c))
            for a in odd_tail_box(h):
                if Q % (7 * a) != 0:
                    raise RuntimeError(("body tail wall", H, a))
            checks += 1
    return checks


def audit_wall_doubling() -> int:
    checks = 0
    for m in range(1, 101):
        for k in range(-2 * m, 2 * m + 1):
            for sign in (-1, 1):
                x = Fraction(14 * k + sign, 14 * m)
                y = 2 * x
                if (7 * m) % y.denominator != 0:
                    raise RuntimeError(("wall doubling", m, k, sign, x, y))
                checks += 1
    return checks


def audit_long_block_lemma() -> int:
    checks = 0
    for modulus in range(2, 401):
        for length in range(1, 33):
            for residue in range(modulus):
                premise = all(
                    7 * abs_residue(j * residue, modulus) < 2 * modulus
                    for j in range(1, length + 1)
                )
                if premise:
                    centered = residue if residue <= modulus // 2 else residue - modulus
                    if 7 * length * abs(centered) >= 2 * modulus:
                        raise RuntimeError(
                            ("long-block lemma", modulus, length, residue, centered)
                        )
                checks += 1
    return checks


def audit_adjacent_owner_rigidity() -> int:
    """Exhaust the local fact that repairs centering's lost top bit."""
    checks = 0
    for clock in range(2, 401, 2):
        for residue in range(clock):
            centered = (
                residue if residue <= clock // 2 else residue - clock
            )
            for r in range(clock - 1):
                if eligible_at(centered, r, clock) and eligible_at(
                    centered, r + 1, clock
                ):
                    n0 = nearest_integer(centered * r, clock)
                    n1 = nearest_integer(centered * (r + 1), clock)
                    if n0 != n1:
                        raise RuntimeError(
                            ("adjacent centered owner changed", clock, centered, r)
                        )
                    checks += 1
    return checks


def top_bit_control() -> tuple[
    int,
    int,
    int,
    tuple[bool, bool],
    tuple[bool, bool],
    tuple[bool, bool],
    tuple[bool, bool],
]:
    # At h=11 the height-uniform clock has a maximal 2-adic wall c=8.
    # A single point really does lose the top bit, but its adjacent point
    # changes that bit and prevents false complementarity on a long block.
    h = 11
    D = height_modulus(h)
    Q = 7 * D
    c = 8
    if v2(14 * c) != v2(Q):
        raise RuntimeError(("top-bit valuation", v2(14 * c), v2(Q)))
    r = Q // (14 * c)
    if r % 2 != 1:
        raise RuntimeError(("top-bit wall label parity", r))
    z0, z1 = 1, 1 + Q
    masks0_r = tuple(danger(z0, r, j, Q) for j in (0, 1))
    masks1_r = tuple(danger(z1, r, j, Q) for j in (0, 1))
    masks0_next = tuple(danger(z0, r + 1, j, Q) for j in (0, 1))
    masks1_next = tuple(danger(z1, r + 1, j, Q) for j in (0, 1))
    if masks0_r != (True, False) or masks1_r != (False, True):
        raise RuntimeError(("isolated top-bit literal masks", masks0_r, masks1_r))
    if masks0_next != (True, False) or masks1_next != (True, False):
        raise RuntimeError(
            ("adjacent top-bit repair masks", masks0_next, masks1_next)
        )
    return Q, c, r, masks0_r, masks1_r, masks0_next, masks1_next


def main() -> None:
    print("adaptive-clock LCM sharpening exact audit")
    print("height D_digits Q_digits simple_factor old_factor")
    for h in (1, 2, 11, 12, 40, 72, 128, 256):
        D, Q, simple_factor, old_factor = audit_height(h)
        print(h, len(str(D)), len(str(Q)), simple_factor, old_factor)

    # Exhaust all small heights as a boundary/parity control.
    for h in range(1, 257):
        audit_height(h)
    print("height_audit_1_256=PASS")

    body_checks = audit_small_bodies()
    print(f"eleven_body_checks_h11_18={body_checks}")

    wall_checks = audit_wall_doubling()
    print(f"wall_doubling_checks={wall_checks}")

    lemma_checks = audit_long_block_lemma()
    print(f"long_block_checks={lemma_checks}")

    adjacent_checks = audit_adjacent_owner_rigidity()
    print(f"adjacent_owner_rigidity_checks={adjacent_checks}")

    Q, c, r, masks0_r, masks1_r, masks0_next, masks1_next = top_bit_control()
    print(
        "same_clock_top_bit_repair="
        f"h=11 Q_digits={len(str(Q))} c={c} r_parity={r % 2} "
        f"at_r={masks0_r}/{masks1_r} "
        f"at_r_plus_1={masks0_next}/{masks1_next}"
    )

    D_body = body_modulus(H72)
    D_height = height_modulus(max(H72))
    if D_height % D_body != 0:
        raise RuntimeError("H72 body modulus does not divide height modulus")
    _, _, _, old_over_height = audit_height(max(H72))
    print(
        "H72="
        f"body_clock_digits={len(str(7 * D_body))} "
        f"height_clock_digits={len(str(7 * D_height))} "
        f"height_over_body={D_height // D_body} "
        f"old_over_body={old_over_height * (D_height // D_body)}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
