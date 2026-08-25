#!/usr/bin/env python3
"""Exact clock-sidecar audit for the Sun 2-4-6-8 counterexample.

This is a FINITE-EXACT reflection artifact, not a new coverage theorem.  It
uses THM-4027's certified binomial periods and exact local profiles, combines
prime-power factors by CRT, and contrasts them with the nonzero C_60 input
atlas over F_61.  Every audited local factor at the counterexample is positive.
"""

from __future__ import annotations

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path


N = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)
MODULI = (60, 420, 27_720)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def load_parent():
    root = Path(__file__).resolve().parents[1]
    path = root / "04-computation/sun_2468_local_density_thm4027.py"
    spec = spec_from_file_location("sun_local_thm4027", path)
    require(spec is not None and spec.loader is not None, "parent import spec")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prime_power_divisors(parent, q: int) -> tuple[int, ...]:
    return tuple(p**a for p, a in parent.factor_prime_powers(q))


def composite_target_profile(parent, q: int) -> tuple[tuple[int, ...], Fraction, int, int]:
    factors = Fraction(1)
    periods = [1] * len(DEGREES)
    for modulus in prime_power_divisors(parent, q):
        profile, local_periods = parent.local_profile(modulus)
        factors *= profile[N % modulus]
        periods = [a * b for a, b in zip(periods, local_periods)]
    expected_periods = tuple(parent.binomial_period(k, q) for k in DEGREES)
    require(tuple(periods) == expected_periods, (q, "CRT role periods"))
    denominator = 1
    for period in periods:
        denominator *= period
    accepted = factors * denominator / q
    require(accepted.denominator == 1, (q, "integral accepted phase count"))
    return tuple(periods), factors, denominator, accepted.numerator


def cyclic_convolve(a: list[int], b: list[int]) -> list[int]:
    q = len(a)
    out = [0] * q
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[(i + j) % q] += ai * bj
    return out


def nonzero_f61_atlas() -> tuple[tuple[int, ...], int, int, bool]:
    p = 61
    distributions = []
    supports = []
    for degree in DEGREES:
        counts = [0] * p
        for t in range(1, p):
            counts[comb(t, degree) % p] += 1
        require(sum(counts) == 60, (degree, "nonzero phase mass"))
        distributions.append(counts)
        supports.append({i for i, count in enumerate(counts) if count})
    total = [1] + [0] * (p - 1)
    for counts in distributions:
        total = cyclic_convolve(total, counts)
    target = N % p
    accepted = total[target]
    require(sum(total) == 60**4, "four-role nonzero atlas mass")
    two_role_full = {(a + b) % p for a in supports[0] for b in supports[1]} == set(range(p))
    require((target, accepted) == (21, 212_846), "F_61 target control")
    require(two_role_full, "S_2+S_4 coverage")
    return tuple(len(s) for s in supports), accepted, 60**4, two_role_full


def main() -> None:
    parent = load_parent()
    rows = tuple((q,) + composite_target_profile(parent, q) for q in MODULI)
    expected = (
        (60, (120, 720, 3600, 7200), Fraction(2486, 3375), 2_239_488_000_000, 27_493_171_200),
        (420, (840, 5040, 25200, 352800), Fraction(773146, 1157625), 37_639_074_816_000_000, 59_852_633_702_400),
        (27720, (55440, 332640, 1663200, 23284800), Fraction(19117792, 46690875), 714_191_507_917_848_576_000_000, 10_549_385_805_850_214_400),
    )
    require(rows == expected, "frozen composite controls")
    # Direct composite-modulus convolution is still cheap at 60 and 420 and
    # independently controls the CRT assembly used at 27720.
    for q, periods, factor, _, _ in rows[:2]:
        direct_profile, direct_periods = parent.local_profile(q)
        require(tuple(direct_periods) == periods, (q, "direct role periods"))
        require(direct_profile[N % q] == factor, (q, "direct target factor"))
    atlas = nonzero_f61_atlas()

    print("SUN_2468_CLOCK_HEIGHT_FIREWALL_20260824")
    print(f"target={N}")
    for q, periods, factor, denominator, accepted in rows:
        print(
            f"q={q};Nmod={N % q};role_periods={periods};"
            f"local_factor={factor};phase_denominator={denominator};accepted={accepted}"
        )
    print(
        "F61_nonzero_atlas="
        f"(target=21,image_sizes={atlas[0]},accepted={atlas[1]}/{atlas[2]},S2_plus_S4_full={atlas[3]})"
    )
    print("preserved=role-labelled congruence phase and exact periodic multiplicity")
    print("lost=integral lift quotient;anisotropic height coupling;exact triangular-square test")
    print("next_test=p61 phase tuple plus lift quotient and nearest-square gap on THM4026 survivors")
    print("scope=FINITE-EXACT local hostile;Sun conjecture already REFUTED;global leastness OPEN")


if __name__ == "__main__":
    main()
