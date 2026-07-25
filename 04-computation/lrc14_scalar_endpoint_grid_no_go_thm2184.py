#!/usr/bin/env python3
"""Exact audit for THM-2184's scalar endpoint-grid no-go addendum."""

from fractions import Fraction


INVOICE = Fraction(12493, 35640)
HARMONIC_SIX = sum(Fraction(1, j) for j in range(1, 7))
SC7_CONSTANT = Fraction(961, 28665)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def main() -> None:
    # Six distinct divisors of Q have six distinct complementary divisors.
    # The equality row also respects the odd-guard and thirteen-unit clauses.
    q_grid = 60
    guard = 15
    units = [60, 30, 20, 12, 10]
    coefficients = [guard, *units]
    require(guard % 2 == 1, "equality guard is not odd")
    require(len(set(coefficients)) == 6, "equality coefficients repeat")
    require(
        all(q_grid % a == 0 and a % 13 != 0 for a in coefficients),
        "equality coefficient is not a thirteen-unit divisor",
    )
    complementary = sorted(q_grid // a for a in coefficients)
    require(complementary == [1, 2, 3, 4, 5, 6], "wrong complements")
    require(
        Fraction(sum(coefficients), q_grid) == HARMONIC_SIX,
        "equality control misses the envelope",
    )
    require(HARMONIC_SIX == Fraction(49, 20), "wrong sixth harmonic number")

    # Independent finite hostile check of the divisor envelope.
    checked_q = 360
    for q in range(1, checked_q + 1):
        ds = divisors(q)
        if len(ds) < 6:
            continue
        six_largest = sorted(ds, reverse=True)[:6]
        require(
            Fraction(sum(six_largest), q) <= HARMONIC_SIX,
            f"divisor envelope failed at Q={q}",
        )

    grid_ratio = INVOICE / Fraction(7, 40)
    rho = 1 / grid_ratio
    require(grid_ratio == Fraction(12493, 6237), "wrong grid ratio")
    require(
        grid_ratio - 2 == Fraction(19, 6237), "wrong grid-ratio margin"
    )
    require(rho == Fraction(6237, 12493), "wrong reciprocal grid ratio")
    require(rho < Fraction(1, 2), "deep scale is not below half the grid")

    # THM-2192's unique-deepest alternative.  The bounded searches include
    # the symbolic minima 1+2+13 and 1+13+169 as positive controls.
    units_13 = [x for x in range(1, 200) if x % 13]
    multiples_13 = list(range(13, 200, 13))
    repeated_minimum = min(
        a + b + c
        for a in units_13
        for b in units_13
        if a != b
        for c in multiples_13
    )
    strict_minimum = min(
        a + b + c
        for a in units_13
        for b in multiples_13
        for c in range(169, 400, 169)
    )
    require(repeated_minimum == 16, "wrong repeated-minimum-depth floor")
    require(strict_minimum == 183, "wrong strictly-increasing-depth floor")

    # After division by J, failure of SC.7 is monotone in n.  The n=1 row
    # is therefore the hostile endpoint; n=2 is printed for transparency.
    gap_n1 = 16 * (1 - rho) - SC7_CONSTANT
    gap_n2 = 16 * (2 - rho) - 2 * SC7_CONSTANT
    require(
        gap_n1 == Fraction(219788159, 27547065), "wrong n=1 normalized gap"
    )
    require(gap_n1 > 0, "n=1 does not exclude SC.7")
    require(
        gap_n2 == Fraction(659617678, 27547065), "wrong n=2 normalized gap"
    )
    require(gap_n2 > gap_n1, "n=2 is not farther from SC.7")
    require(16 - SC7_CONSTANT > 0, "SC.7 gap is not monotone in n")

    # J is a multiple of fourteen, so k=13^d-nJ has this fixed signed
    # congruence.  Since J>2M, k is negative for every n>=1.
    residues = []
    for d in range(1, 13):
        signed_residue = pow(13, d, 14)
        absolute_residue = (-signed_residue) % 14
        expected_signed = 13 if d % 2 else 1
        expected_absolute = 1 if d % 2 else 13
        require(signed_residue == expected_signed, f"wrong k residue at d={d}")
        require(
            absolute_residue == expected_absolute,
            f"wrong absolute k residue at d={d}",
        )
        residues.append((d, signed_residue, absolute_residue))

    print("THM-2184 scalar endpoint-grid no-go exact audit")
    print(f"harmonic_six_envelope={HARMONIC_SIX}")
    print(
        "equality_control="
        f"Q:{q_grid};H:{guard};q:{','.join(map(str, units))};"
        f"complements:{','.join(map(str, complementary))}"
    )
    print(f"finite_divisor_hostile_check=Q<= {checked_q}:PASS")
    print(f"invoice_grid_ratio=J/M>={grid_ratio}=2+19/6237")
    print(f"deep_to_grid_ratio=M/J<={rho}<1/2")
    print(
        "unique_deepest_W_minima="
        f"repeated_min_depth:{repeated_minimum};strict_depths:{strict_minimum}"
    )
    print(f"SC7_normalized_gap_n1={gap_n1}")
    print(f"SC7_normalized_gap_n2={gap_n2}")
    print("k_mod_14=d_even:1,d_odd:13;abs_k_mod_14=d_even:13,d_odd:1")
    print("SC7_trigger_on_canonical_full_six_coefficient_grid=EMPTY")


if __name__ == "__main__":
    main()
