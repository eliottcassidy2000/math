#!/usr/bin/env python3
"""Clean-room arithmetic audit for the THM-3818/THM-4450 composition.

This script deliberately checks only the finite/arithmetic interfaces.  The
Haar localization and the actual-entry theorem remain inherited inputs.
"""

from fractions import Fraction
from itertools import combinations, product
from math import gcd


FIVE_RAYS = ((1, 11), (1, 23), (5, 11), (1, 37), (1, 25))
GENERAL_ODD_HOSTILE = (1, 9)
EXPECTED_FACTORS = {
    (1, 11): ((2, 2), (3, 1)),
    (1, 23): ((2, 3), (3, 1)),
    (5, 11): ((2, 4),),
    (1, 37): ((2, 1), (19, 1)),
    (1, 25): ((2, 1), (13, 1)),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factorization(n: int) -> tuple[tuple[int, int], ...]:
    answer: list[tuple[int, int]] = []
    prime = 2
    while prime * prime <= n:
        if n % prime == 0:
            exponent = 0
            while n % prime == 0:
                n //= prime
                exponent += 1
            answer.append((prime, exponent))
        prime += 1
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def inert_sum(sum_value: int) -> bool:
    """THM-3818 condition: every prime is 2 mod 3, exponent at most 2."""
    return all(prime % 3 == 2 and exponent <= 2
               for prime, exponent in factorization(sum_value))


def in_thm3818_atlas(p: int, q: int) -> bool:
    return (
        1 <= p < q
        and gcd(p, q) == 1
        and p + q <= 356
        and inert_sum(p + q)
    )


def gcd_many(values: tuple[int, ...]) -> int:
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def audit_atlas() -> list[tuple[int, int]]:
    atlas = [
        (p, q)
        for p in range(1, 178)
        for q in range(p + 1, 356)
        if in_thm3818_atlas(p, q)
    ]
    require(len(atlas) == 5855, "THM-3818 atlas count changed")
    require(max(p for p, _ in atlas) == 177, "THM-3818 max p changed")
    require(max(q for _, q in atlas) == 355, "THM-3818 max q changed")
    require(not set(FIVE_RAYS).intersection(atlas),
            "five-ray bank intersects THM-3818 atlas")
    for ray in FIVE_RAYS:
        require(factorization(sum(ray)) == EXPECTED_FACTORS[ray],
                f"factorization mismatch for {ray}")
        require(not inert_sum(sum(ray)), f"unexpected inert five-ray {ray}")
    require(in_thm3818_atlas(*GENERAL_ODD_HOSTILE),
            "general-odd hostile (1,9) should survive the inert atlas")
    return atlas


def audit_mass_gates() -> None:
    localization = Fraction(4, 91)

    # q=4, one v_2=1 tail: H=2C union {r}; THM-4450 gives
    # mu(G_H) >= mu(G_C)/2.  The localization boundary is inclusive.
    q4_gate = 2 * localization
    require(q4_gate == Fraction(8, 91), "q=4 gate arithmetic")
    require(q4_gate / 2 == localization, "q=4 transfer arithmetic")

    # q=2, one even tail: H=C union {r}; for the live odd 3-unit r,
    # mu(G_H) >= mu(G_C)-8/63.
    q2_loss = Fraction(8, 63)
    q2_gate = localization + q2_loss
    require(q2_gate == Fraction(20, 117), "q=2 one-even gate arithmetic")
    require(q2_gate - q2_loss == localization, "q=2 transfer arithmetic")

    # These strictly improve the old use of the full 3-unit pair cap 4/77.
    pair_cap = Fraction(4, 77)
    require(2 * pair_cap == Fraction(8, 77), "old q=4 gate arithmetic")
    require(pair_cap + q2_loss == Fraction(124, 693),
            "old q=2 one-even gate arithmetic")
    require(q4_gate < Fraction(8, 77), "q=4 gate is not an improvement")
    require(q2_gate < Fraction(124, 693),
            "q=2 one-even gate is not an improvement")


def audit_q4_unit_decoder() -> int:
    """Check the exact q=4 literal-unit simplification on a hostile bank.

    For H=2C union {r}, r odd, d=gcd(H) is odd.  Thus the normalized
    eleven-shape H/d contains 1 iff r=d, equivalently iff r divides every C.
    The proof is elementary; the finite bank is a regression control.
    """
    checked = 0
    for body in combinations(range(1, 16), 10):
        for r in range(1, 30, 2):
            h = tuple(2 * c for c in body) + (r,)
            d = gcd_many(h)
            has_literal_unit = d in h
            divisibility_test = all(c % r == 0 for c in body)
            require(has_literal_unit == divisibility_test,
                    f"q=4 unit equivalence failed at C={body}, r={r}")

            # The physical decoder core is 2H.  If a,b are odd and
            # g=gcd(a,b), its scale gcd with the pair equals the whole-row gcd.
            a, b = 21, 35
            pair_scale = gcd(a, b)
            physical_row = tuple(2 * value for value in h) + (a, b)
            core_scale = 2 * d
            require(gcd_many(physical_row) == gcd(core_scale, pair_scale),
                    "physical/shape scale gcd mismatch")
            checked += 1
    return checked


def audit_even_tail_obstruction() -> None:
    # Boolean distributivity proves, pointwise, for E=D_(2r), B=D_a union D_b,
    # tau(E)=E:
    #   (E union B) intersect (E union tau(B)) = E union (B intersect tau(B)).
    for e, b, tau_b in product((False, True), repeat=3):
        lhs = (e or b) and (e or tau_b)
        rhs = e or (b and tau_b)
        require(lhs == rhs, "even-tail Boolean identity failed")

    # At r=1 a D_2 tooth already has this physical width, larger than both
    # THM-4451 odd-tail strict-component caps.  Hence odd-tail transfer is
    # genuinely ill-typed for q=4's (2r,a,b) tail triple.
    tooth = Fraction(1, 14)
    require(tooth > Fraction(17, 693), "even tooth did not beat all-odd cap")
    require(tooth > Fraction(19, 1001),
            "even tooth did not beat odd-3-unit cap")


def main() -> None:
    atlas = audit_atlas()
    audit_mass_gates()
    unit_controls = audit_q4_unit_decoder()
    audit_even_tail_obstruction()

    print("THM3818_THM4450_ENTRY_COMPOSITION_AUDIT")
    print(f"atlas_count={len(atlas)} max_p={max(p for p, _ in atlas)} "
          f"max_q={max(q for _, q in atlas)}")
    for ray in FIVE_RAYS:
        print(f"ray={ray} sum={sum(ray)} factors={factorization(sum(ray))} "
              f"inert={inert_sum(sum(ray))} in_atlas={ray in atlas}")
    print("five_ray_atlas_intersection=empty")
    print("scope_hostile=(1, 9) sum=10 factors=((2, 1), (5, 1)) "
          "inert=True in_atlas=True not_3_unit=True")
    print("localization=4/91 q4_body_gate=8/91 "
          "q2_one_even_body_gate=20/117")
    print("old_full_cap_gates: q4=8/77 q2_one_even=124/693")
    print(f"q4_literal_unit_controls={unit_controls}")
    print("q4_unit_iff_r_divides_every_body_speed=PASS")
    print("even_tail_identity=D_2r union two_tail_cross_comb")
    print("r=1_tooth=1/14 > all_odd_cap=17/693 "
          "> odd_3_unit_cap=19/1001")
    print("PASS")


if __name__ == "__main__":
    main()
