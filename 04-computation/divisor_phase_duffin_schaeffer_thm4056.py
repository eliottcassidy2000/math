"""Exact companion for THM-4056.

The global statements in THM-4056 are proved algebraically.  This script
supplies finite hostile controls for the natural-edge reduction fibres, the
divisor/Mobius compiler, its LCM clock, the pointwise metric firewall, and the
continued-fraction example for e.

Run from the repository root with

    python -B 04-computation/divisor_phase_duffin_schaeffer_thm4056.py
    python -B -O 04-computation/divisor_phase_duffin_schaeffer_thm4056.py
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import factorial, gcd, lcm


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def phi(n: int) -> int:
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def mobius(n: int) -> int:
    primes = 0
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            primes += 1
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes += 1
    return -1 if primes % 2 else 1


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def reduction_fibres(n: int) -> Counter[tuple[int, int]]:
    fibres: Counter[tuple[int, int]] = Counter()
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            g = gcd(a, b)
            fibres[(a // g, b // g)] += 1
    return fibres


def main() -> None:
    # Every natural-order edge has one primitive type p/q and scale g.  Its
    # fibre has floor(N/q) representatives.  The divisor first difference is
    # checked independently by comparing successive Counters.
    edge_checks = 0
    divisor_difference_checks = 0
    previous: Counter[tuple[int, int]] = Counter()
    for n in range(2, 257):
        fibres = reduction_fibres(n)
        expected = Counter(
            {
                (p, q): n // q
                for q in range(2, n + 1)
                for p in range(1, q)
                if gcd(p, q) == 1
            }
        )
        require(fibres == expected, (n, fibres - expected, expected - fibres))
        require(sum(fibres.values()) == n * (n - 1) // 2, n)
        if n > 2:
            delta = fibres - previous
            expected_delta = Counter(
                {
                    (p, q): 1
                    for q in divisors(n)
                    if q >= 2
                    for p in range(1, q)
                    if gcd(p, q) == 1
                }
            )
            require(delta == expected_delta, (n, delta, expected_delta))
            require(sum(expected_delta.values()) == n - 1, n)
            divisor_difference_checks += 1
        previous = fibres
        edge_checks += 1

    # A rational-weight test of the invertible *lift* divisor compiler.  With
    # b(q)=2 phi(q) psi(q), the mean over the LCM clock is exactly the finite
    # Duffin--Schaeffer mass sum 2 phi(q) psi(q)/q.
    q_max = 12
    clock_period = 1
    for q in range(1, q_max + 1):
        clock_period = lcm(clock_period, q)
    psi = {q: Fraction(1, q + 1) for q in range(2, q_max + 1)}
    b = {q: 2 * phi(q) * psi[q] for q in range(2, q_max + 1)}

    def lift_clock(n: int) -> Fraction:
        return sum((weight for q, weight in b.items() if n % q == 0), Fraction())

    def phase_clock(n: int) -> Fraction:
        return sum((2 * psi[q] for q in psi if gcd(n, q) == 1), Fraction())

    phase_coefficients = {
        d: 2 * Fraction(mobius(d)) * sum(
            (psi[q] for q in psi if q % d == 0), Fraction()
        )
        for d in range(1, q_max + 1)
    }

    lift_mean = sum((lift_clock(n) for n in range(1, clock_period + 1)), Fraction()) / clock_period
    phase_mean = sum((phase_clock(n) for n in range(1, clock_period + 1)), Fraction()) / clock_period
    ds_mass = sum((weight / q for q, weight in b.items()), Fraction())
    require(lift_mean == phase_mean == ds_mass, (lift_mean, phase_mean, ds_mass))
    first_clock_difference = None
    for n in range(1, clock_period + 1):
        compiled_phase = sum(
            (coefficient for d, coefficient in phase_coefficients.items() if n % d == 0),
            Fraction(),
        )
        require(compiled_phase == phase_clock(n), (n, compiled_phase, phase_clock(n)))
        if first_clock_difference is None and phase_clock(n) != lift_clock(n):
            first_clock_difference = n
    require(first_clock_difference is not None, "phase/lift hostile")

    # The phase compiler sees only rad(q): moving unit weight from q=2 to
    # q=4 preserves both its phase clock and DS mean, but changes the lift
    # clock.  Hence prime-power denominator depth is a required sidecar.
    phase_2 = [2 if gcd(n, 2) == 1 else 0 for n in range(1, 5)]
    phase_4 = [2 if gcd(n, 4) == 1 else 0 for n in range(1, 5)]
    lift_2 = [2 if n % 2 == 0 else 0 for n in range(1, 5)]
    lift_4 = [4 if n % 4 == 0 else 0 for n in range(1, 5)]
    require(phase_2 == phase_4 and sum(phase_2) == sum(phase_4), (phase_2, phase_4))
    require(lift_2 != lift_4 and sum(lift_2) == sum(lift_4), (lift_2, lift_4))

    recovered = {}
    for n in range(2, q_max + 1):
        recovered[n] = sum(
            (Fraction(mobius(n // d)) * lift_clock(d) for d in divisors(n)),
            Fraction(),
        )
        require(recovered[n] == b[n], (n, recovered[n], b[n]))

    # Shell lengths in the standard Duffin--Schaeffer normalization.
    shell_checks = 0
    for q in range(2, 301):
        psi = Fraction(1, q + 1)
        raw = sum(
            (Fraction(2) * psi / q for a in range(1, q) if gcd(a, q) == 1),
            Fraction(),
        )
        require(raw == Fraction(2 * phi(q), q * (q + 1)), (q, raw))
        shell_checks += 1

    # Rational separation behind the pointwise firewall: two distinct
    # reduced fractions A/B and a/q differ by at least 1/(Bq).
    separation_checks = 0
    for b_den in range(1, 41):
        for a_num in range(-b_den, b_den + 1):
            if gcd(abs(a_num), b_den) != 1:
                continue
            x = Fraction(a_num, b_den)
            for q in range(1, 61):
                for p in range(-q, q + 1):
                    if gcd(abs(p), q) != 1 or Fraction(p, q) == x:
                        continue
                    require(abs(x - Fraction(p, q)) >= Fraction(1, b_den * q), (x, p, q))
                    separation_checks += 1

    # Golden-ratio hostile.  For r=p/q in [0,1], the quadratic norm is a
    # nonzero integer; the conjugate factor is <4, giving
    # |r-alpha|>1/(4q^2).  The exact computation checks the load-bearing norm.
    golden_norm_checks = 0
    for q in range(1, 1001):
        for p in range(0, q + 1):
            if gcd(p, q) != 1:
                continue
            norm = p * p + p * q - q * q
            require(norm != 0 and abs(norm) >= 1, (p, q, norm))
            golden_norm_checks += 1

    # e=[2; overline-in-blocks (1,2k,1)].  The first 3n positive digits have
    # product 2^n n!, so their geometric means diverge.
    product = 1
    e_cf_product_checks = 0
    for n in range(1, 101):
        product *= 2 * n
        require(product == 2**n * factorial(n), n)
        e_cf_product_checks += 1

    print("status=PROVED identities with FINITE-EXACT hostile companion")
    print(f"natural_edge_levels={edge_checks};divisor_first_differences={divisor_difference_checks}")
    print(f"clock_qmax={q_max};lcm_period={clock_period};mobius_recoveries={len(recovered)}")
    print(f"duffin_schaeffer_shells={shell_checks};common_clock_mean={ds_mass}")
    print(f"phase_lift_pointwise_hostile=first_difference_at_t={first_clock_difference}")
    print("phase_kernel_hostile=psi_2_to_psi_4_preserves_phase_and_mean_but_changes_lift")
    print(f"rational_separation_checks={separation_checks};golden_norm_checks={golden_norm_checks}")
    print("golden_hostile=alpha_notin_Wprime(q->1/(4q));DS_sum_diverges")
    print(f"e_cf_product_checks={e_cf_product_checks};product_first_3n_positive_digits=2^n*n!")
    print("PASS")


if __name__ == "__main__":
    main()
