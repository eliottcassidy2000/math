#!/usr/bin/env python3
"""Complete the Berggren positive-cube slope atlas through denominator 101.

The inherited prime screen leaves eleven slopes.  This companion implements
the Lagrange--Mollin--Matthews continued-fraction algorithm for

    W^2 - K*h^2 = C,

then follows every solution class under the fundamental norm-one unit modulo
2*n.  The cube compiler requires h=n*u with u odd and W odd, equivalently

    h == n (mod 2*n),  W == 1 (mod 2).

No floating point arithmetic is used.  A SymPy implementation of the same
LMM algorithm is used only as an optional independent control.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt
from pathlib import Path


SURVIVORS = (
    (14, 23),
    (26, 29),
    (26, 47),
    (38, 47),
    (38, 53),
    (50, 53),
    (50, 71),
    (62, 95),
    (74, 95),
    (74, 101),
    (98, 101),
)

EXPECTED_ADMISSIBLE = ((14, 23), (26, 29), (26, 47), (98, 101))
EXPECTED_CLASS_COUNTS = {
    (14, 23): 3,
    (26, 29): 6,
    (26, 47): 9,
    (38, 47): 0,
    (38, 53): 9,
    (50, 53): 0,
    (50, 71): 0,
    (62, 95): 1,
    (74, 95): 1,
    (74, 101): 6,
    (98, 101): 6,
}
EXPECTED_BLOCKED_PERIODS = {
    (38, 53): (2,) * 9,
    (62, 95): (60,),
    (74, 95): (10,),
    (74, 101): (34,) * 6,
}

# If p occurs exactly once in C and K is a nonsquare modulo p, then
# W^2-K*h^2=C forces p|W,h and hence p^2|C, a contradiction.
PADIC_DESCENT_PRIMES = {(38, 47): 13, (50, 53): 11, (50, 71): 13}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factorint(n: int) -> dict[int, int]:
    """Trial-division factorization; all inputs here are below 40000."""
    n = abs(n)
    factors: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def divisors(n: int) -> tuple[int, ...]:
    values = [1]
    for p, exponent in factorint(n).items():
        values = [a * p**j for a in values for j in range(exponent + 1)]
    return tuple(sorted(values))


def square_factor(n: int) -> int:
    value = 1
    for p, exponent in factorint(n).items():
        value *= p ** (exponent // 2)
    return value


def square_roots_mod(d: int, modulus: int) -> tuple[int, ...]:
    return tuple(z for z in range(modulus) if (z * z - d) % modulus == 0)


def symmetric_residue(a: int, modulus: int) -> int:
    return a if a <= modulus // 2 else a - modulus


def pqa(p0: int, q0: int, d: int):
    """Exact P,Q,a,A,B,G continued-fraction generator."""
    root = isqrt(d)
    a2 = b1 = 0
    a1 = b2 = 1
    g1, g2 = q0, -p0
    p, q = p0, q0
    while True:
        digit = (p + root) // q
        a1, a2 = digit * a1 + a2, a1
        b1, b2 = digit * b1 + b2, b1
        g1, g2 = digit * g1 + g2, g1
        yield p, q, digit, a1, b1, g1
        p = digit * q - p
        q = (d - p * p) // q


def continued_fraction_length(p0: int, q0: int, d: int) -> int:
    """Aperiodic plus periodic length of (p0+sqrt(d))/q0."""
    seen: dict[tuple[int, int], int] = {}
    stream = pqa(p0, q0, d)
    index = 0
    while True:
        p, q, *_ = next(stream)
        state = (p, q)
        if state in seen:
            return index
        seen[state] = index
        index += 1


def pell_pm_one(d: int, sign: int) -> tuple[tuple[int, int], ...]:
    """Fundamental positive solution of x^2-d*y^2=sign, if it exists."""
    require(sign in (-1, 1), "Pell sign")
    root = isqrt(d)
    require(root * root != d, "Pell discriminant nonsquare")
    stream = pqa(0, 1, d)
    *_, previous_b, previous_g = next(stream)
    period_index = None
    for index, row in enumerate(stream):
        *_, digit, _a, b, g = row
        if digit == 2 * root:
            period_index = index
            break
        previous_b, previous_g = b, g
    require(period_index is not None, "Pell period closes")
    if period_index % 2:
        return ((previous_g, previous_b),) if sign == 1 else ()
    if sign == -1:
        return ((previous_g, previous_b),)
    for _ in range(period_index):
        row = next(stream)
    *_, b, g = row
    return ((g, b),)


def lmm_classes(d: int, norm: int) -> tuple[tuple[tuple[int, int], ...], dict[str, int]]:
    """One representative of every generalized-Pell solution class.

    This is the Lagrange--Mollin--Matthews algorithm.  Multiplication by
    powers of the fundamental norm-one unit, together with independent sign
    changes/conjugation, generates every integer solution from this list.
    """
    require(d > 0 and isqrt(d) ** 2 != d and norm != 0, "LMM input")
    solutions: list[tuple[int, int]] = []
    root_branches = 0
    negative_pell_fallbacks = 0
    for scale in divisors(square_factor(norm)):
        reduced_norm = norm // (scale * scale)
        modulus = abs(reduced_norm)
        for square_root in square_roots_mod(d, modulus):
            root_branches += 1
            z = symmetric_residue(square_root, modulus)
            stream = pqa(z, modulus, d)
            *_, previous_b, previous_g = next(stream)
            for _ in range(continued_fraction_length(z, modulus, d) - 1):
                _p, q, _digit, _a, b, g = next(stream)
                if abs(q) == 1:
                    if previous_g**2 - d * previous_b**2 == reduced_norm:
                        solutions.append((scale * previous_g, scale * previous_b))
                    else:
                        negative = pell_pm_one(d, -1)
                        if negative:
                            negative_pell_fallbacks += 1
                            a, b_neg = negative[0]
                            solutions.append(
                                (
                                    scale * (previous_g * a + previous_b * d * b_neg),
                                    scale * (previous_g * b_neg + previous_b * a),
                                )
                            )
                    break
                previous_b, previous_g = b, g
    for x, y in solutions:
        require(x * x - d * y * y == norm, "LMM representative norm")
    return tuple(solutions), {
        "root_branches": root_branches,
        "negative_pell_fallbacks": negative_pell_fallbacks,
    }


def conic_data(m: int, n: int) -> tuple[int, int]:
    numerator_k = 4 * m * m - n * n
    numerator_c = 4 * (2 * m * m + 2 * m * n - n * n)
    require(numerator_k % 3 == 0 and numerator_c % 3 == 0, "integral conic")
    return numerator_k // 3, numerator_c // 3


def unit_orbit(
    representative: tuple[int, int], d: int, unit: tuple[int, int], modulus: int
) -> tuple[int, tuple[int, ...]]:
    """Return orbit period and phases satisfying the odd compiler residue."""
    p, q = unit
    w, h = representative[0] % modulus, representative[1] % modulus
    start = (w, h)
    seen: dict[tuple[int, int], int] = {}
    hits: list[int] = []
    while (w, h) not in seen:
        seen[(w, h)] = len(seen)
        if w % 2 == 1 and h == modulus // 2:
            hits.append(len(seen) - 1)
        w, h = (p * w + d * q * h) % modulus, (q * w + p * h) % modulus
    require((w, h) == start, "unit orbit is purely periodic")
    return len(seen), tuple(hits)


def legendre_symbol(a: int, p: int) -> int:
    value = pow(a % p, (p - 1) // 2, p)
    return -1 if value == p - 1 else value


def main() -> None:
    ledger = []
    admissible = []
    sympy_matches = 0

    try:
        from sympy.solvers.diophantine.diophantine import diop_DN
    except ImportError:
        diop_DN = None

    for slope in SURVIVORS:
        m, n = slope
        k, constant = conic_data(m, n)
        require(k > 0 and constant > 0 and k % 2 == 1, f"{slope}: positive odd K")
        representatives, audit = lmm_classes(k, constant)
        require(len(representatives) == EXPECTED_CLASS_COUNTS[slope], f"{slope}: class count")

        unit_rows = pell_pm_one(k, 1)
        require(len(unit_rows) == 1, f"{slope}: positive fundamental unit")
        unit = unit_rows[0]
        require(unit[0] ** 2 - k * unit[1] ** 2 == 1, f"{slope}: unit norm")

        if diop_DN is not None:
            require(tuple(diop_DN(k, constant)) == representatives, f"{slope}: SymPy LMM control")
            require(tuple(diop_DN(k, 1)) == unit_rows, f"{slope}: SymPy Pell control")
            sympy_matches += 1

        orbit_rows = tuple(unit_orbit(rep, k, unit, 2 * n) for rep in representatives)
        good_classes = tuple(index for index, (_period, hits) in enumerate(orbit_rows) if hits)
        if good_classes:
            admissible.append(slope)

        if slope in EXPECTED_BLOCKED_PERIODS:
            require(
                tuple(period for period, _hits in orbit_rows) == EXPECTED_BLOCKED_PERIODS[slope],
                f"{slope}: blocked periods",
            )
            require(not good_classes, f"{slope}: blocked orbit has no target phase")

        if slope in PADIC_DESCENT_PRIMES:
            prime = PADIC_DESCENT_PRIMES[slope]
            require(constant % prime == 0 and constant % (prime * prime) != 0, f"{slope}: exact p order")
            require(legendre_symbol(k, prime) == -1, f"{slope}: K nonsquare")
            require(not representatives, f"{slope}: p-adic descent agrees with LMM")
        else:
            prime = 0

        representative_text = ",".join(f"{x}:{y}" for x, y in representatives)
        orbit_text = ",".join(
            f"{period}:" + ".".join(map(str, hits)) for period, hits in orbit_rows
        )
        ledger.append(
            f"{m}/{n}|K={k}|C={constant}|classes={len(representatives)}|"
            f"unit={unit[0]}:{unit[1]}|reps={representative_text}|"
            f"orbits={orbit_text}|good={'.'.join(map(str, good_classes))}|p={prime}|"
            f"branches={audit['root_branches']}|fallback={audit['negative_pell_fallbacks']}"
        )

    require(tuple(admissible) == EXPECTED_ADMISSIBLE, "complete admissible slope list")
    require(sympy_matches in (0, len(SURVIVORS)), "optional SymPy control is all-or-none")

    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    ledger_text = ";".join(ledger)
    print("BERGGREN POSITIVE-CUBE SLOPE ATLAS COMPLETION THROUGH 101")
    print("status=VERIFIED-EXACT+FINITE-EXACT")
    print("algorithm=Lagrange-Mollin-Matthews_PQa_plus_complete_unit_orbits_mod_2n")
    print(f"survivor_count={len(SURVIVORS)};sympy_independent_matches={sympy_matches}")
    print("admissible_slopes=" + repr(list(admissible)).replace(" ", ""))
    print("blocked_padic_descent=(38,47):13,(50,53):11,(50,71):13")
    print("blocked_unit_orbits=(38,53):9x2,(62,95):1x60,(74,95):1x10,(74,101):6x34")
    for row in ledger:
        print("row=" + row)
    print(f"ledger_sha256={sha256(ledger_text.encode('ascii')).hexdigest()}")
    print(f"source_lf_sha256={sha256(source).hexdigest()}")
    print(
        "scope=the_517_prime-screen exclusions are inherited from THM-3547;"
        "the four surviving slopes inherit explicit positive infinite rays;"
        "no asymptotic count or classification beyond n<=101"
    )


if __name__ == "__main__":
    main()
