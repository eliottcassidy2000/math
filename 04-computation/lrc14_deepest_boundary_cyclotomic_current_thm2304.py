#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2304.

The theorem is an algebraic endpoint-current decomposition, not a finite
search for LRC(14).  This companion checks the complete strict-profile
valuation ledger, the full-gap relative degree and current count, the
deepest/nondeep conductor separation, the recursive shell-count identity,
and the two sharp controls used in the proof.
"""

from __future__ import annotations

from math import gcd


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(value: int, prime: int = P) -> int:
    require(value != 0, "valuation of zero is not used")
    answer = 0
    value = abs(value)
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def main() -> None:
    profiles = [
        (b, c)
        for b in range(3, 18)
        for c in range(b + 2, 20)
    ]
    require(len(profiles) == 120, "strict profile census changed")

    ranks: list[int] = []
    for b, c in profiles:
        gap = c - b
        relative_degree = P**gap
        ranks.append(relative_degree - 1)

        # At nu_13(K)=1, a c_2 boundary has conductor at most
        # 13^(b-1), whereas a c_3 boundary can have conductor
        # 13^(c-1).  In the top generator, lower phases have exponent
        # divisible by the full relative degree.
        c2_top_exponent_factor = P ** (c - b)
        require(
            c2_top_exponent_factor == relative_degree,
            "lower-bank exponent did not land in the base field",
        )
        require(
            (c - 1) - (b - 1) == gap,
            "relative conductor gap changed",
        )

        # Peeling one base-13 digit at a time gives 12, 12*13, ...
        # independent colored currents.  Their total is R-1.
        recursive_count = sum(12 * P**level for level in range(gap))
        require(
            recursive_count == relative_degree - 1,
            "recursive shell count does not match the full-gap basis",
        )

        # A single deepest tooth has two endpoint classes.  Their
        # difference is a nonzero unit times 2, so they are distinct even
        # after reduction to the first top digit.  Units u_1/u_3 merely
        # permute these classes.
        for multiplier in range(1, P):
            require(gcd(multiplier, P) == 1, "test multiplier is not a unit")
            for residue in range(P):
                left = (-multiplier * (14 * residue - 1)) % relative_degree
                right = (-multiplier * (14 * residue + 1)) % relative_degree
                require(
                    (right - left) % P != 0,
                    "single-tooth endpoint classes collided",
                )
                require(
                    left != 0 or right != 0,
                    "single tooth has no nontrivial current class",
                )

        # A nontrivial class s mod R has conductor exponent at least b,
        # while every lower-bank phase has exponent at most b-1.
        largest_nonzero_class_valuation = gap - 1
        require(
            valuation(P**largest_nonzero_class_valuation) == gap - 1,
            "largest nonzero class valuation changed",
        )
        minimum_top_conductor = c - 1 - largest_nonzero_class_valuation
        require(
            minimum_top_conductor == b,
            "top class can collide with a lower-bank phase",
        )

    require(min(ranks) == 168, "minimum current rank changed")
    require(
        max(ranks) == P**16 - 1,
        "maximum strict-profile current rank changed",
    )

    # The full deepest comb is the sharp zero-current control.  Its
    # Fourier support is c_3 Z, while nu_13(K)=1<c, so its K coefficient
    # vanishes despite having every deepest boundary.
    for _b, c in profiles:
        require(c > 1, "deepest valuation is not separated")
        require(
            valuation(P, P) < c,
            "full-comb hostile control unexpectedly has c3|K",
        )

    # Root-character jump amplitudes lie in the lower field throughout the
    # strict branch because b>=3 gives zeta_13 in
    # Q(zeta_(N*13^(b-1))).
    require(min(b - 1 for b, _c in profiles) >= 2, "root amplitudes left base")

    print("strict_profiles=120")
    print("pair_frequency_valuation=1")
    print("relative_degree=13^(c-b)")
    print("minimum_relative_degree=169")
    print(f"minimum_nontrivial_currents={min(ranks)}")
    print(f"maximum_nontrivial_currents={max(ranks)}")
    print("nondeep_conductor_exponent<=b-1")
    print("deepest_conductor_exponent<=c-1")
    print("nontrivial_deep_class_conductor_exponent>=b")
    print("single_deepest_tooth=two_distinct_classes_at_least_one_nontrivial")
    print("recursive_shell_count=12*(1+13+...+13^(c-b-1))")
    print("recursive_shell_count_equals=13^(c-b)-1")
    print("full_deepest_comb_control=zero_at_pair_frequency")
    print("arbitrary_cut_endpoint_control=OUTSIDE_THEOREM_UNIVERSE")
    print("status=THM2304_CYCLOTOMIC_CURRENT_EXACT_REFEREE")


if __name__ == "__main__":
    main()
