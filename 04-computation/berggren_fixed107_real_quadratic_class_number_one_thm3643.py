#!/usr/bin/env python3
"""Exact/rigorous-ball companion for THM-3643.

The algebraic half is integer-only: it factors the fundamental discriminant,
constructs the complete continued fraction of sqrt(D), and reconstructs the
minimal Pell unit.  The analytic half evaluates the finite logarithmic
Dirichlet formula in Arb balls.  The final interval contains the sole integer
1, which proves the ordinary class number.
"""

from hashlib import sha256
import json
from math import isqrt

from flint import arb, ctx


D = 1_225_041
Q = 408_347
PRECISION_BITS = 160


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def sqrt_continued_fraction(n: int) -> tuple[int, tuple[int, ...]]:
    """Return a0 and the full period of sqrt(n), using exact recurrences."""
    a0 = isqrt(n)
    require(a0 * a0 != n, "D is nonsquare")
    m, den, a = 0, 1, a0
    period: list[int] = []
    while True:
        m = den * a - m
        num = n - m * m
        require(num > 0 and num % den == 0, "continued-fraction division")
        den = num // den
        a = (a0 + m) // den
        period.append(a)
        if a == 2 * a0:
            return a0, tuple(period)


def convergent(coefficients: tuple[int, ...]) -> tuple[int, int]:
    p_old2, p_old1 = 0, 1
    q_old2, q_old1 = 1, 0
    for a in coefficients:
        p = a * p_old1 + p_old2
        q = a * q_old1 + q_old2
        p_old2, p_old1 = p_old1, p
        q_old2, q_old1 = q_old1, q
    return p_old1, q_old1


def digest_ints(values: tuple[int, ...]) -> str:
    return sha256(",".join(map(str, values)).encode("ascii")).hexdigest()


def main() -> None:
    require(D == 3 * Q, "factorization")
    require(is_prime(Q), "prime cofactor")
    require(D % 8 == 1, "fundamental-discriminant congruence")

    a0, period = sqrt_continued_fraction(D)
    require(a0 == 1106, "continued-fraction head")
    require(len(period) == 1372, "continued-fraction period length")
    require(period[-1] == 2 * a0, "continued-fraction terminal")
    require(period[:-1] == tuple(reversed(period[:-1])), "period palindrome")

    # The period length is even, so the negative Pell equation is insoluble
    # and the +1 fundamental solution is convergent number L-1.
    x, y = convergent((a0,) + period[:-1])
    require(x * x - D * y * y == 1, "fundamental Pell identity")
    require(x.bit_length() == 2361 and y.bit_length() == 2351, "unit size")

    # O_K=Z[(1+sqrt(D))/2].  A half-integral unit would give odd x,y with
    # x^2-Dy^2=+/-4.  Since D=1 mod 8, its left side is 0 mod 8, impossible.
    # Thus the Pell unit above generates the full unit group modulo signs.

    # Quadratic residues modulo the prime Q.  Since D=1 mod 4,
    # chi_D(a)=(a/D)=(a/3)(a/Q), also for even a because D=1 mod 8.
    qr = bytearray(Q)
    for z in range(1, (Q + 1) // 2):
        qr[(z * z) % Q] = 1

    ctx.prec = PRECISION_BITS
    pi = arb.pi()
    half_log_sine_sum = arb(0)
    counts = {-1: 0, 0: 0, 1: 0}
    character_hash = sha256()
    for a in range(1, (D + 1) // 2):
        if a % 3 == 0 or a % Q == 0:
            chi = 0
        else:
            chi3 = 1 if a % 3 == 1 else -1
            chiq = 1 if qr[a % Q] else -1
            chi = chi3 * chiq
        counts[chi] += 1
        character_hash.update(bytes((chi + 1,)))
        if chi:
            half_log_sine_sum += chi * (2 * (pi * a / D).sin()).log()

    regulator = (arb(x) + arb(y) * arb(D).sqrt()).log()
    class_number_ball = -half_log_sine_sum / regulator
    lower_gate = arb("0.999999999999")
    upper_gate = arb("1.000000000001")
    require(class_number_ball.contains(1), "class-number ball misses 1")
    require(class_number_ball.lower() > lower_gate, "class-number lower gate")
    require(class_number_ball.upper() < upper_gate, "class-number upper gate")

    semantic = {
        "D": D,
        "q": Q,
        "period_length": len(period),
        "period_sha256": digest_ints(period),
        "pell_x_sha256": sha256(str(x).encode("ascii")).hexdigest(),
        "pell_y_sha256": sha256(str(y).encode("ascii")).hexdigest(),
        "character_counts": [counts[-1], counts[0], counts[1]],
        "character_sha256": character_hash.hexdigest(),
        "class_number_integer": 1,
    }
    semantic_hash = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("== THM-3643 fixed-107 real quadratic class-number certificate ==")
    print("D=1225041=3*408347;fundamental_discriminant=1;cofactor_prime=1")
    print(
        f"continued_fraction=a0:{a0};period:{len(period)};"
        f"period_sha256:{semantic['period_sha256']}"
    )
    print("negative_pell=insoluble_by_even_period;half_integral_units=excluded_mod8")
    print(
        f"pell_unit=x_bits:{x.bit_length()};y_bits:{y.bit_length()};"
        f"x_sha256:{semantic['pell_x_sha256']};y_sha256:{semantic['pell_y_sha256']}"
    )
    print(f"character_counts=(-1:{counts[-1]},0:{counts[0]},1:{counts[1]})")
    print(f"character_sha256={semantic['character_sha256']}")
    print(f"half_log_sine_sum={half_log_sine_sum}")
    print(f"regulator={regulator}")
    print(f"class_number_ball={class_number_ball}")
    print(f"semantic_sha256={semantic_hash}")
    print("status=PROVED FINITE-EXACT RIGOROUS-ARB;ordinary_class_number=1")
    print("scope=quadratic sidecar only;no Mordell rank equality or Selmer bound")


if __name__ == "__main__":
    main()
