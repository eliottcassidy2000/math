#!/usr/bin/env python3
"""Exact deterministic audit for THM-1098.

No floating point is used in a proof-bearing assertion.  The script checks:
  * the strict 17/41 margin of the fixed 12-speed core;
  * arbitrary finite-atlas killing by one divisor-loaded coordinate;
  * the explicit odd-half witness near 17/41;
  * prefix blindness, q0, and exact q_min values on representative lcm rows;
  * fixed q=41 address with component-width bound tending to zero.
"""

from fractions import Fraction
from functools import reduce
from math import gcd, lcm


CORE = tuple(range(1, 12)) + (13,)
LEVEL_DEN = 14


def lonely_at(speeds: tuple[int, ...], numerator: int, denominator: int) -> bool:
    """Exact level-1/14 test; the fraction need not be reduced."""
    return all(
        LEVEL_DEN * min((v * numerator) % denominator,
                        denominator - ((v * numerator) % denominator))
        >= denominator
        for v in speeds
    )


def least_denominator(speeds: tuple[int, ...], cap: int) -> tuple[int, int] | None:
    for q in range(2, cap + 1):
        for p in range(1, q):
            if gcd(p, q) == 1 and lonely_at(speeds, p, q):
                return q, p
    return None


def nearest_odd_to(frac: Fraction) -> int:
    floor = frac.numerator // frac.denominator
    candidates = [z for z in range(floor - 2, floor + 4) if z % 2]
    return min(candidates, key=lambda z: abs(Fraction(z) - frac))


def explicit_half_witness(M: int) -> tuple[int, int]:
    """Return a/(2M), with a odd and within 1/(2M) of 17/41."""
    a = nearest_odd_to(Fraction(34 * M, 41))
    assert abs(Fraction(a, 2 * M) - Fraction(17, 41)) <= Fraction(1, 2 * M)
    assert a % 2 == 1
    return a, 2 * M


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def first_unblocked_modulus(speeds: tuple[int, ...], cap: int) -> int:
    """Least q for which no speed is divisible by q."""
    for q in range(2, cap + 1):
        if all(v % q != 0 for v in speeds):
            return q
    raise RuntimeError("first unblocked modulus exceeds cap")


print("THM-1098 EXACT DETERMINISTIC AUDIT")
print()

print("[1] fixed-core margin at t0=17/41")
core_residues = []
for v in CORE:
    r = (17 * v) % 41
    core_residues.append(min(r, 41 - r))
print("    speeds:   ", CORE)
print("    41*dist:  ", tuple(core_residues))
assert min(core_residues) == 3
margin = Fraction(3, 41) - Fraction(1, 14)
assert margin == Fraction(1, 574)
print("    min distance = 3/41 = 1/14 + 1/574")
print()

print("[2] arbitrary finite atlas is killed, but the row is explicitly lonely")
atlas_denominators = (21, 41, 53, 83, 89)
D = reduce(lcm, atlas_denominators, 1)
M_atlas = lcm(84, D)
if M_atlas < 3731:
    M_atlas *= (3731 + M_atlas - 1) // M_atlas
S_atlas = CORE + (M_atlas,)
assert reduce(gcd, S_atlas) == 1
assert is_covering(S_atlas)
for q in atlas_denominators:
    assert M_atlas % q == 0
    for p in range(1, q):
        if gcd(p, q) == 1:
            assert not lonely_at(S_atlas, p, q)
a, den = explicit_half_witness(M_atlas)
assert lonely_at(S_atlas, a, den)
red_den = den // gcd(a, den)
print("    atlas denominators:", atlas_denominators)
print("    M =", M_atlas)
print("    every reduced numerator in the atlas is killed by the M-runner")
print(f"    explicit lonely point = {a}/{den}; reduced denominator = {red_den}")
print()

print("[3] prefix rows S_B={1,...,11,13,lcm(1,...,B)}")
print("    B   digits(M)   q0   exact q_min   excess   explicit-denominator<=2M")
for B in (14, 23, 25, 30, 39, 41, 50, 60, 80, 100, 150, 200):
    M = reduce(lcm, range(1, B + 1), 1)
    speeds = CORE + (M,)
    assert M >= 3731 and M % 84 == 0
    assert reduce(gcd, speeds) == 1 and is_covering(speeds)
    for q in range(2, B + 1):
        assert M % q == 0
        for p in range(1, q):
            if gcd(p, q) == 1:
                assert not lonely_at(speeds, p, q)
    a, den = explicit_half_witness(M)
    assert lonely_at(speeds, a, den)
    red_den = den // gcd(a, den)
    assert B < red_den <= 2 * M
    exact = least_denominator(speeds, B + 250)
    assert exact is not None and exact[0] > B
    q0 = first_unblocked_modulus(speeds, B + 250)
    assert B < q0 <= exact[0]
    if B == 23:
        assert (q0, exact) == (25, (53, 22))
    print(f"    {B:3d} {len(str(M)):11d}  {q0:3d}   "
          f"{exact[0]:3d} at {exact[1]}/{exact[0]:<3d}  {exact[0] - q0:6d}   "
          f"yes ({red_den} <= 2M)")
print()

print("[4] fixed rational address, vanishing geometric thickness")
print("    k      M       17/41 lonely   component upper bound 6/(7M)")
for k in (0, 1, 2, 10, 100, 1000):
    M = 84 * (41 * k + 1)
    speeds = CORE + (M,)
    assert reduce(gcd, speeds) == 1 and is_covering(speeds)
    assert lonely_at(speeds, 17, 41)
    width = Fraction(6, 7 * M)
    print(f"    {k:4d} {M:8d}       yes          {width}")
print()
print("ALL EXACT ASSERTIONS PASSED")
