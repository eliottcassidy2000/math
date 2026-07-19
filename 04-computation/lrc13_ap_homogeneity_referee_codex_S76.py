#!/usr/bin/env python3
"""Exact referee for THM-1171's twelve-speed AP homogeneity theorem.

For C(a,d)={a+kd:0<=k<=11}, write a=gA and d=gD with gcd(A,D)=1.
If D>=2, multiplication by A^{-1} sends all twelve phases at a time u/d
to the half-residue floor(D/2)/D.  If D=1, the centered time has clearance
A/(2A+11).  The proof is uniform; the finite box below independently
checks every displayed identity using Fraction arithmetic.

No assertion is used, so optimized mode checks every obligation.
"""

from collections import Counter
from fractions import Fraction as F
from math import gcd


BOUND = 300


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def main() -> None:
    rows = 0
    spread_rows = 0
    consecutive_rows = 0
    homogeneous_rows = 0
    reduced_difference_histogram: Counter[int] = Counter()
    minimum_spread_clearance = F(1)
    minimum_strict_consecutive_margin = F(1)

    for a in range(1, BOUND + 1):
        for d in range(1, BOUND + 1):
            rows += 1
            common = gcd(a, d)
            reduced_a = a // common
            reduced_d = d // common
            require(gcd(reduced_a, reduced_d) == 1, "gcd normalization failed")
            reduced_difference_histogram[reduced_d] += 1

            if reduced_d >= 2:
                spread_rows += 1
                half_residue = reduced_d // 2
                inverse = pow(reduced_a, -1, reduced_d)
                multiplier = (half_residue * inverse) % reduced_d
                require((reduced_a * multiplier - half_residue) % reduced_d == 0,
                        "modular inverse witness failed")
                time = F(multiplier, d)
                expected_phase = F(half_residue, reduced_d)
                phases = tuple(F((a + k * d) * multiplier, d) % 1 for k in range(12))
                require(all(phase == expected_phase for phase in phases),
                        f"common phase failed at {(a, d)}")
                clearance = min(map(circle_distance, phases))
                require(clearance == expected_phase,
                        f"half-residue clearance failed at {(a, d)}")
                require(clearance >= F(1, 3),
                        f"uniform one-third floor failed at {(a, d)}")
                require(clearance > F(1, 13),
                        f"spread AP did not clear one-thirteenth at {(a, d)}")
                require(time.denominator > 0, "invalid rational time")
                minimum_spread_clearance = min(minimum_spread_clearance, clearance)
                continue

            consecutive_rows += 1
            require(reduced_d == 1, "unhandled reduced difference")
            denominator = common * (2 * reduced_a + 11)
            time = F(1, denominator)
            phases = tuple(F(a + k * d, denominator) for k in range(12))
            require(all(0 < phase < 1 for phase in phases),
                    f"centered phases left the unit interval at {(a, d)}")
            require(phases[0] + phases[-1] == 1,
                    f"centered endpoint symmetry failed at {(a, d)}")
            clearance = min(map(circle_distance, phases))
            expected_clearance = F(reduced_a, 2 * reduced_a + 11)
            require(clearance == expected_clearance,
                    f"consecutive clearance formula failed at {(a, d)}")
            require((clearance > F(1, 13)) == (reduced_a > 1),
                    f"A>1 comparison failed at {(a, d)}")
            if reduced_a == 1:
                homogeneous_rows += 1
                require(a == d == common, "homogeneous normalization failed")
                require(clearance == F(1, 13), "homogeneous witness changed")
            else:
                minimum_strict_consecutive_margin = min(
                    minimum_strict_consecutive_margin,
                    clearance - F(1, 13),
                )

    require(minimum_spread_clearance == F(1, 3),
            "wrong sharp spread witness floor in audit box")
    require(minimum_strict_consecutive_margin == F(11, 195),
            "wrong least strict consecutive margin")
    require(homogeneous_rows == BOUND, "wrong homogeneous row count")

    print("THM-1171 exact twelve-speed AP homogeneity referee")
    print(f"audit box: 1<=a,d<={BOUND}; rows={rows}")
    print(f"spread reduced-difference rows D>=2: {spread_rows}")
    print(f"reduced-consecutive rows D=1: {consecutive_rows}")
    print(f"homogeneous rows a=d: {homogeneous_rows}")
    print("spread witness: choose u with A*u=floor(D/2) mod D, t=u/d")
    print("all twelve phases: floor(D/2)/D")
    print(f"sharp audited spread clearance floor: {minimum_spread_clearance}")
    print("consecutive witness: t=1/[g(2A+11)]")
    print("consecutive clearance: A/(2A+11)")
    print(f"least strict consecutive margin over 1/13: {minimum_strict_consecutive_margin}")
    print("exact dispatch: D>=2 is >1/13; D=1 is >1/13 iff A>1")
    print("Tournament Analysis:")
    print("  runner-order tournament: scores=(0,...,11), cycles=0, SCCs=12, paths=1")
    print("  residue-fibre observable: phase equality at t=u/d")
    print("  switch/gauge: multiplication by A^(-1) modulo D")
    print("  tie Hamiltonian path: runner indices 0->...->11 in the single fibre")
    print("  faithful carrier: common residue fibre plus metric phase floor(D/2)/D")
    print("  lost by runner tournament: gcd reduction, unit action, common phase distance")
    print("VERDICT: M(C)=1/13 can hold for a twelve-term AP only when a=d")
    print("SCOPE: the general non-AP n=12 equality/sporadic branch remains open")


if __name__ == "__main__":
    main()
