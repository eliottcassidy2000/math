#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2093.

The audit checks the five-level valuation flag, the private-coefficient
divisibility staircase on a finite exhaustive bank, sharp root-path valuation
budgets for every prime through 57, denominator clearing on synthetic stars,
and the stated uniform constant.  It deliberately uses explicit runtime
checks so that ``python -O`` runs the same audit.
"""

from __future__ import annotations

from math import gcd, isqrt, lcm
from random import Random


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(n: int, p: int) -> int:
    require(n != 0, "valuation of zero requested")
    n = abs(n)
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def primes_through(n: int) -> list[int]:
    out: list[int] = []
    for k in range(2, n + 1):
        if all(k % p for p in out if p * p <= k):
            out.append(k)
    return out


def audit_valuation_flags() -> int:
    rng = Random(2093)
    labels = ["x", "y", "2h0", "4h1", "8h2", "H"] + [f"V{i}" for i in range(7)]
    expected = [
        {"x", "y"},
        {"x", "y", "2h0"},
        {"x", "y", "2h0", "4h1"},
        {"x", "y", "2h0", "4h1", "8h2"},
        {"x", "y", "2h0", "4h1", "8h2", "H"},
    ]
    for trial in range(1000):
        odd = [2 * rng.randrange(1, 500) + 1 for _ in range(6)]
        h0, h1, h2, h, x, y = odd
        while True:
            q = [rng.randrange(1, 1000) for _ in range(7)]
            if gcd(*q) == 1:
                break
        values = [x, y, 2 * h0, 4 * h1, 8 * h2, 16 * h] + [32 * z for z in q]
        for k in range(1, 6):
            modulus = 1 << k
            support = {label for label, value in zip(labels, values) if value % modulus}
            require(support == expected[k - 1], f"support mismatch trial={trial} k={k}")
            retained = [value for label, value in zip(labels, values) if label not in support]
            retained_gcd = gcd(*retained)
            require(valuation(retained_gcd, 2) == k, f"deletion gcd mismatch trial={trial} k={k}")
    return 1000


def audit_private_staircase() -> int:
    checked = 0
    odds = range(1, 16, 2)
    for t in range(4):
        for h in odds:
            for q in odds:
                g = gcd(h, q)
                for u in odds:
                    forced = (1 << (4 - t)) * g // gcd(g, u)
                    for a in range(-10, 11):
                        for b in range(-10, 11):
                            numerator = -(16 * a * h + 32 * b * q)
                            denominator = (1 << t) * u
                            if numerator == 0 or numerator % denominator:
                                continue
                            c_big = numerator // denominator
                            require(c_big % forced == 0, "private divisor failed")
                            require(c_big % (1 << (4 - t)) == 0, "dyadic divisor failed")
                            c = c_big // (1 << (4 - t))
                            require(a * h + 2 * b * q + c * u == 0, "normalized relation failed")
                            require((a - c) % 2 == 0, "parity coupling failed")
                            checked += 1
    require(checked > 100_000, "staircase bank unexpectedly small")
    return checked


def audit_sharp_path_budgets() -> int:
    primes = primes_through(57)
    lambda_57 = lcm(*range(1, 58))
    reconstructed = 1
    for p in primes:
        f = 0
        power = 1
        while power * p <= 57:
            power *= p
            f += 1
        reconstructed *= p**f

        # The six relations -q_parent + p^f q_child=0 attain the full
        # six-edge loss budget.  The final child is 1, so the terminal gcd is
        # one while gcd(h,q_root) has valuation 6f.
        root = p ** (6 * f)
        h = root
        path = [root]
        for _ in range(6):
            child = path[-1] // power
            require(-path[-1] + power * child == 0, "path relation failed")
            path.append(child)
        require(path[-1] == 1, "sharp path did not reach a unit")
        require(gcd(*path) == 1, "sharp path is not primitive")
        require(valuation(gcd(h, root), p) == 6 * f, "sharp exponent mismatch")
    require(reconstructed == lambda_57, "LCM factor reconstruction failed")
    return len(primes)


def audit_denominator_clearing() -> int:
    rng = Random(2053)
    for trial in range(500):
        d0 = rng.randrange(1, 50)
        while True:
            a = rng.randrange(1, 100)
            b = rng.randrange(1, 100)
            if gcd(a, b) == 1:
                break
        h_anchor = d0 * a
        v_anchor = d0 * b
        triples: list[tuple[int, int, int, int]] = []
        speeds = [h_anchor, v_anchor]
        for _ in range(11):
            aa = -rng.randrange(1, 30)
            bb = -rng.randrange(1, 30)
            cc = d0
            speed = -(aa * h_anchor + bb * v_anchor) // cc
            require(aa * h_anchor + bb * v_anchor + cc * speed == 0, "star relation failed")
            triples.append((aa, bb, cc, speed))
            speeds.append(speed)

        # Force a primitive row without changing the relation pattern.
        row_gcd = gcd(*speeds)
        h_anchor //= row_gcd
        v_anchor //= row_gcd
        triples = [(aa, bb, cc, speed // row_gcd) for aa, bb, cc, speed in triples]
        speeds = [h_anchor, v_anchor] + [entry[3] for entry in triples]
        require(gcd(*speeds) == 1, "synthetic row not primitive")

        actual_d0 = gcd(h_anchor, v_anchor)
        aa0 = h_anchor // actual_d0
        bb0 = v_anchor // actual_d0
        coeffs = [abs(entry[2]) for entry in triples]
        ell = lcm(*coeffs)
        require(actual_d0 > 0, "zero anchor gcd")

        # If the synthetic coefficient choices acquired a common factor after
        # primitive normalization, replace them by the primitive triple
        # relations before clearing.
        primitive: list[tuple[int, int, int, int]] = []
        for aa, bb, cc, speed in triples:
            rel_gcd = gcd(abs(aa), abs(bb), abs(cc))
            aa1, bb1, cc1 = aa // rel_gcd, bb // rel_gcd, cc // rel_gcd
            require(aa1 * h_anchor + bb1 * v_anchor + cc1 * speed == 0, "primitive relation failed")
            primitive.append((aa1, bb1, cc1, speed))
        ell = lcm(*(abs(entry[2]) for entry in primitive))
        require(ell % actual_d0 == 0, f"D0 does not divide L on trial {trial}")

        scale = ell // actual_d0
        cleared = [scale * speed for speed in speeds]
        presented = [ell * aa0, ell * bb0]
        for aa, bb, cc, _speed in primitive:
            presented.append((-aa * ell // cc) * aa0 + (-bb * ell // cc) * bb0)
        require(cleared == presented, f"cleared presentation mismatch trial={trial}")
    return 500


def main() -> None:
    lambda_57 = lcm(*range(1, 58))
    q0 = 91**6
    finite_bound = 2912 * lambda_57**6 * q0**13

    flag_trials = audit_valuation_flags()
    staircase_relations = audit_private_staircase()
    sharp_primes = audit_sharp_path_budgets()
    clearing_trials = audit_denominator_clearing()

    require(lambda_57 == 164249358725037825439200, "Lambda_57 constant mismatch")
    require(q0 == 567869252041, "Q0 constant mismatch")
    require(len(str(finite_bound)) == 296, "finite-bound digit count mismatch")
    require(isqrt(2 * q0**24) < 2 * q0**12, "round-envelope arithmetic failed")

    print("LRC14 DYADIC GLOBAL-STAR AUDIT")
    print(f"lambda_57={lambda_57}")
    print(f"q0={q0}")
    print(f"valuation_flag_trials={flag_trials}")
    print(f"staircase_exact_relations={staircase_relations}")
    print(f"sharp_path_primes={sharp_primes}")
    print(f"denominator_clearing_trials={clearing_trials}")
    print(f"finite_bound_digits={len(str(finite_bound))}")
    print(f"finite_bound={finite_bound}")
    print("PASS")


if __name__ == "__main__":
    main()
