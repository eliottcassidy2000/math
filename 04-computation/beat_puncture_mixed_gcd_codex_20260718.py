#!/usr/bin/env python3
"""Exact beat-puncture audit for the six-on-seven slow-gap cone.

Mathematical predicate
----------------------
Let G_k be a closed safe gap of the slow carrier a, and suppose the six
faster combs b_1,...,b_6 cover it.  For a pair u,v among the seven speeds,
put q=u+v or q=|u-v|.  At every t=p/q, u and v have equal circle distance.
Consequently, whenever that common distance is at least 1/14 and t lies in
G_k, a complementary faster comb must kill t.  A failure is an exact lonely
witness, not a floating-point heuristic.

The script checks this all-pairs predicate, the weaker carrier-only predicate,
and the phase-free mixed-gcd count bound from THM-1192.  Every comparison is
integer-exact.

Tournament-analysis audit
-------------------------
The natural pair observable is the number of missed safe beat points for
{u,v}.  It is symmetric in u,v, so orienting it requires an arbitrary gauge;
the increasing-speed tie path would carry no mathematical information.  The
faithful object is instead the incidence hypergraph

    defining pair -> safe beat point -> complementary covering combs.

Thus tournament score vectors, cycles, SCCs, edge flips, and Hamiltonian-path
counts are deliberately not reported.  Alternate vertex sets considered were
speeds, beat denominators, beat points, and proof obligations.  Denominators
or obligations are better vertices than runners here, but the ternary
incidence is still lost by every binary tournament quotient.  This explicitly
challenges the default assumption that LRC tournament vertices should be
runners.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
from math import gcd
from random import Random

RADIUS_DENOMINATOR = 14


def ceil_div(n: int, d: int) -> int:
    return -((-n) // d)


def dangerous(speed: int, p: int, q: int) -> bool:
    """Return ||speed*p/q|| < 1/14, using only integers."""
    residue = (speed * p) % q
    return RADIUS_DENOMINATOR * min(residue, q - residue) < q


def gap_block(a: int, k: int, q: int) -> tuple[int, int]:
    """Integer p with p/q in G_k, including the safe endpoints."""
    denominator = RADIUS_DENOMINATOR * a
    lo = ceil_div(q * (RADIUS_DENOMINATOR * k + 1), denominator)
    hi = q * (RADIUS_DENOMINATOR * k + 13) // denominator
    return lo, hi


def window_count(q_reduced: int) -> int:
    """# {r mod Q : ||r/Q|| < 1/14} = 2 ceil(Q/14)-1."""
    return 2 * ceil_div(q_reduced, RADIUS_DENOMINATOR) - 1


def phase_free_cap(block_size: int, q_reduced: int) -> int:
    """Maximum supplied by full periods plus an arbitrary residual subset."""
    full, rem = divmod(block_size, q_reduced)
    count = window_count(q_reduced)
    return full * count + min(rem, count)


def pair_violation(
    a: int,
    faster: tuple[int, ...],
    k: int,
) -> tuple[int, int, str, int, int] | None:
    """First all-pairs beat puncture missed by the complementary combs."""
    speeds = (a,) + faster
    for left in range(7):
        for right in range(left + 1, 7):
            u, v = speeds[left], speeds[right]
            for sign, q in (("sum", u + v), ("difference", abs(u - v))):
                if q == 0:
                    continue
                lo, hi = gap_block(a, k, q)
                for p in range(lo, hi + 1):
                    if dangerous(u, p, q):
                        continue
                    killed_by_complement = any(
                        dangerous(c, p, q)
                        for index, c in enumerate(faster, 1)
                        if index not in (left, right)
                    )
                    if not killed_by_complement:
                        return u, v, sign, q, p
    return None


def packet_passes_all_pairs(a: int, faster: tuple[int, ...]) -> bool:
    """Whether at least one carrier gap survives every beat-puncture test."""
    return any(pair_violation(a, faster, k) is None for k in range(a))


def carrier_modes(
    a: int,
    faster: tuple[int, ...],
    k: int,
) -> tuple[bool, bool, bool]:
    """Phase-free count, exact sum count, and exact union for q=a+b_i."""
    phase_free_ok = True
    exact_sum_ok = True
    exact_union_ok = True
    for omitted, b in enumerate(faster):
        q = a + b
        lo, hi = gap_block(a, k, q)
        block_size = hi - lo + 1
        caps: list[int] = []
        counts: list[int] = []
        union_covers = True
        for index, c in enumerate(faster):
            if index == omitted:
                continue
            q_reduced = q // gcd(c, q)
            caps.append(phase_free_cap(block_size, q_reduced))
            counts.append(
                sum(dangerous(c, p, q) for p in range(lo, hi + 1))
            )
        for p in range(lo, hi + 1):
            if not any(
                dangerous(c, p, q)
                for index, c in enumerate(faster)
                if index != omitted
            ):
                union_covers = False
                break
        phase_free_ok &= block_size <= sum(caps)
        exact_sum_ok &= block_size <= sum(counts)
        exact_union_ok &= union_covers
    return phase_free_ok, exact_sum_ok, exact_union_ok


def harmonic_crowded(a: int, faster: tuple[int, ...]) -> bool:
    return a * sum((Fraction(1, b) for b in faster), Fraction()) > 1


def mixed_gcd(a: int, faster: tuple[int, ...]) -> bool:
    return gcd(a, *faster) == 1


def exhaustive_carrier_audit(max_a: int = 7, ratio: int = 3) -> tuple[int, int, int, int]:
    packets = phase_free = exact_sum = exact_union = 0
    for a in range(2, max_a + 1):
        for faster in combinations(range(a + 1, ratio * a + 1), 6):
            if not harmonic_crowded(a, faster) or not mixed_gcd(a, faster):
                continue
            packets += 1
            modes = [carrier_modes(a, faster, k) for k in range(a)]
            phase_free += any(mode[0] for mode in modes)
            exact_sum += any(mode[1] for mode in modes)
            exact_union += any(mode[2] for mode in modes)
    return packets, phase_free, exact_sum, exact_union


def exhaustive_all_pairs(max_a: int, ratio: int) -> tuple[int, int, list[tuple[int, int, int]]]:
    packets = survivors = 0
    per_a: list[tuple[int, int, int]] = []
    for a in range(2, max_a + 1):
        at_a = survive_at_a = 0
        for faster in combinations(range(a + 1, ratio * a + 1), 6):
            if not harmonic_crowded(a, faster) or not mixed_gcd(a, faster):
                continue
            packets += 1
            at_a += 1
            if packet_passes_all_pairs(a, faster):
                survivors += 1
                survive_at_a += 1
        per_a.append((a, at_a, survive_at_a))
    return packets, survivors, per_a


def random_all_pairs(
    samples: int,
    max_a: int,
    ratio: int,
    seed: int,
) -> tuple[int, int, int]:
    rng = Random(seed)
    accepted = survivors = attempts = 0
    while accepted < samples:
        attempts += 1
        if attempts > max(1000, 50 * samples):
            raise RuntimeError("random sampler could not fill the crowded cone")
        a = rng.randrange(3, max_a + 1)
        upper = ratio * a - 3
        span = upper - a
        values: set[int] = set()
        while len(values) < 6:
            # Square bias samples the crowded small-speed end while retaining
            # support all the way through the six-on-seven cutoff 6a-3.
            offset = min(span - 1, int(rng.random() ** 2 * span))
            values.add(a + 1 + offset)
        faster = tuple(sorted(values))
        if not harmonic_crowded(a, faster) or not mixed_gcd(a, faster):
            continue
        accepted += 1
        survivors += packet_passes_all_pairs(a, faster)
    return accepted, attempts, survivors


def check_window_formula(limit: int = 500) -> None:
    for q in range(1, limit + 1):
        brute = sum(dangerous(1, p, q) for p in range(q))
        assert brute == window_count(q), (q, brute, window_count(q))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=10)
    parser.add_argument("--ratio", type=int, default=3)
    parser.add_argument("--random", type=int, default=50_000)
    parser.add_argument("--random-max-a", type=int, default=60)
    parser.add_argument("--random-ratio", type=int, default=6)
    parser.add_argument("--seed", type=int, default=1192)
    parser.add_argument("--skip-carrier", action="store_true")
    parser.add_argument("--skip-exhaustive", action="store_true")
    args = parser.parse_args()

    check_window_formula()
    print("THM-1192 EXACT BEAT-PUNCTURE AUDIT")
    print("A(Q)=2*ceil(Q/14)-1 brute check, Q=1..500: PASS")

    if not args.skip_carrier:
        carrier = exhaustive_carrier_audit()
        print("carrier-only sum-beat audit (a<=7, b6<=3a):")
        print(f"  harmonic mixed-gcd packets: {carrier[0]}")
        print(f"  some gap passes phase-free count cap: {carrier[1]}")
        print(f"  some gap passes exact sum count: {carrier[2]}")
        print(f"  some gap passes exact carrier-beat union: {carrier[3]}")

        example = (5, (6, 7, 9, 11, 12, 15), 2)
        assert carrier_modes(*example)[2]
        violation = pair_violation(*example)
        assert violation == (6, 7, "sum", 13, 6), violation
        print("carrier-only survivor refuted by a fast-fast beat:")
        print(f"  (a,b,k)={example}; witness (u,v,kind,q,p)={violation}; t=6/13")

    if not args.skip_exhaustive:
        packets, survivors, per_a = exhaustive_all_pairs(args.max_a, args.ratio)
        print(f"all-pairs exhaustive audit (a<={args.max_a}, b6<={args.ratio}a):")
        for a, count, survive in per_a:
            print(f"  a={a}: packets={count}, survivors={survive}")
        print(f"  TOTAL packets={packets}, survivors={survivors}")

    if args.random:
        accepted, attempts, random_survivors = random_all_pairs(
            args.random,
            args.random_max_a,
            args.random_ratio,
            args.seed,
        )
        print(
            "seeded all-pairs audit "
            f"(a<={args.random_max_a}, b6<={args.random_ratio}a-3, seed={args.seed}):"
        )
        print(
            f"  accepted crowded mixed-gcd draws={accepted}, "
            f"attempts={attempts}, survivors={random_survivors}"
        )

    print("tournament audit: pair observable is symmetric; no canonical orientation")
    print("faithful object: pair -> safe beat point -> complementary-comb incidence")
    print("PROVED content: all-range inclusion and A(Q) count only")
    print("TELEMETRY content: every finite/random elimination count above")


if __name__ == "__main__":
    main()
