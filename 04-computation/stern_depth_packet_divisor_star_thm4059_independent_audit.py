"""Independent hostile audit for THM-4059.

Rederives the inverse formula, packet symmetries, group-character and
arithmetic-multiplicativity failures, cyclic compiler, Möbius inverse, and
tournament lower-star degree counts.

Reproduce:
    python -B 04-computation/stern_depth_packet_divisor_star_thm4059_independent_audit.py
    python -B -O 04-computation/stern_depth_packet_divisor_star_thm4059_independent_audit.py
"""

from __future__ import annotations

from math import gcd


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def stern_depth(a: int, q: int) -> int:
    """Sum of canonical simple-CF coefficients minus one."""
    require(a > 0 and q > 0 and gcd(a, q) == 1, (a, q))
    coefficient_sum = 0
    while q:
        coefficient, remainder = divmod(a, q)
        coefficient_sum += coefficient
        a, q = q, remainder
    return coefficient_sum - 1


def depth_sign_direct(a: int, q: int) -> int:
    return -1 if stern_depth(a, q) % 2 else 1


def inverse_data(a: int, q: int) -> tuple[int, int]:
    inverse = pow(a, -1, q)
    carry, remainder = divmod(a * inverse - 1, q)
    require(remainder == 0, (a, q, inverse, carry, remainder))
    return inverse, carry


def depth_sign_inverse(a: int, q: int) -> int:
    inverse, carry = inverse_data(a, q)
    if q % 2:
        exponent = a + inverse
    else:
        exponent = carry + 1
    return -1 if exponent % 2 else 1


def packet(q: int) -> dict[int, int]:
    if q == 1:
        return {0: 1}
    return {
        a: depth_sign_direct(a, q)
        for a in range(1, q)
        if gcd(a, q) == 1
    }


def packet_sum(q: int) -> int:
    return sum(packet(q).values())


def moebius(n: int) -> int:
    require(n >= 1, n)
    prime_count = 0
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            prime_count += 1
            if n % p == 0:
                return 0
        while n % p == 0:
            n //= p
        p += 1
    if n > 1:
        prime_count += 1
    return -1 if prime_count % 2 else 1


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def audit_inverse_formula(q_max: int) -> int:
    checks = 0
    for q in range(2, q_max + 1):
        values = packet(q)
        for a, direct in values.items():
            inverse, carry = inverse_data(a, q)
            require(direct == depth_sign_inverse(a, q),
                    (a, q, direct, inverse, carry))
            require(values[inverse] == direct,
                    ("inverse symmetry", a, q, inverse, direct))
            require(values[q - a] == direct,
                    ("negation symmetry", a, q, q - a, direct))
            checks += 1

        total = sum(values.values())
        if q > 2:
            require((total - len(values)) % 4 == 0,
                    ("mod-four packet law", q, total, len(values)))
        if q % 2 == 1:
            odd_sum = sum(sign for a, sign in values.items() if a % 2 == 1)
            even_sum = sum(sign for a, sign in values.items() if a % 2 == 0)
            require(odd_sum == even_sum and odd_sum + even_sum == total,
                    ("odd/even half-shell", q, odd_sum, even_sum, total))
    return checks


def first_internal_character_failure(q_max: int) -> tuple[int, tuple[int, int, int]]:
    for q in range(2, q_max + 1):
        values = packet(q)
        for a in values:
            for b in values:
                product = a * b % q
                if values[product] != values[a] * values[b]:
                    return q, (a, b, product)
    raise AssertionError("no internal character hostile")


def first_odd_internal_character_failure(q_max: int) -> tuple[int, tuple[int, int, int]]:
    for q in range(3, q_max + 1, 2):
        values = packet(q)
        for a in values:
            for b in values:
                product = a * b % q
                if values[product] != values[a] * values[b]:
                    return q, (a, b, product)
    raise AssertionError("no odd internal character hostile")


def first_arithmetic_multiplicativity_failure(q_max: int) -> tuple[int, int]:
    sums = {q: packet_sum(q) for q in range(1, q_max * q_max + 1)}
    candidates: list[tuple[int, int, int]] = []
    for m in range(2, q_max + 1):
        for n in range(m + 1, q_max + 1):
            if gcd(m, n) == 1 and sums[m * n] != sums[m] * sums[n]:
                candidates.append((m * n, m, n))
    require(candidates, q_max)
    _, m, n = min(candidates)
    return m, n


def audit_compiler(n_max: int) -> tuple[int, list[tuple[int, int, int]]]:
    sums = {q: packet_sum(q) for q in range(1, n_max + 1)}
    star_rows: list[tuple[int, int, int]] = []
    checks = 0
    for n in range(1, n_max + 1):
        compiled_total = 0
        direct_star = 0
        for x in range(n):
            if x == 0:
                compiled_total += 1
                continue
            common = gcd(x, n)
            q = n // common
            a = x // common
            sign = depth_sign_direct(a, q)
            compiled_total += sign
            direct_star += sign

        divisor_total = sum(sums[q] for q in divisors(n))
        require(compiled_total == divisor_total,
                ("compiler", n, compiled_total, divisor_total))
        require(direct_star == divisor_total - 1,
                ("star", n, direct_star, divisor_total - 1))
        require((n - 1 + direct_star) % 2 == 0,
                ("degree parity", n, direct_star))
        lower_in = (n - 1 + direct_star) // 2
        lower_out = (n - 1 - direct_star) // 2
        require(lower_in >= 0 and lower_out >= 0 and lower_in + lower_out == n - 1,
                (n, lower_in, lower_out))

        recovered = sum(
            moebius(n // d) * (
                sum(packet_sum(q) for q in divisors(d))
            )
            for d in divisors(n)
        )
        require(recovered == sums[n], ("Moebius inverse", n, recovered, sums[n]))
        checks += 1
        if n <= 15:
            star_rows.append((n, sums[n], direct_star))
    return checks, star_rows


def audit_exact_hostiles() -> None:
    # The first even packet fails to be a group character at its identity.
    require(packet(2) == {1: -1}, packet(2))

    # The first odd failure is q=9: epsilon(2)^2 != epsilon(4).
    values9 = packet(9)
    require(values9[2] == -1 and values9[4] == -1,
            values9)
    require(values9[2] * values9[2] != values9[4], values9)

    # Arithmetic multiplicativity first fails at coprime 2 and 5.
    require(packet_sum(2) == -1, packet_sum(2))
    require(packet_sum(5) == 0, packet_sum(5))
    require(packet_sum(10) == -4, packet_sum(10))
    require(packet_sum(10) != packet_sum(2) * packet_sum(5),
            (packet_sum(2), packet_sum(5), packet_sum(10)))

    # Odd-only CRT hostile: q=15 is all positive although q=5 is balanced.
    require(packet_sum(3) == 2 and packet_sum(15) == 8,
            (packet_sum(3), packet_sum(15)))
    require(packet_sum(15) != packet_sum(3) * packet_sum(5),
            (packet_sum(3), packet_sum(5), packet_sum(15)))
    require(packet(15)[2] == 1, packet(15))
    require(packet(3)[2] * packet(5)[2] == -1,
            (packet(3), packet(5)))

    # The compiled color on C_4 is (+,-,-,-), not an additive character.
    compiled4 = []
    for x in range(4):
        if x == 0:
            compiled4.append(1)
        else:
            common = gcd(x, 4)
            compiled4.append(depth_sign_direct(x // common, 4 // common))
    require(tuple(compiled4) == (1, -1, -1, -1), compiled4)
    require(compiled4[2] != compiled4[1] * compiled4[1], compiled4)


def main() -> None:
    inverse_checks = audit_inverse_formula(2000)
    first_group_failure = first_internal_character_failure(100)
    first_odd_group_failure = first_odd_internal_character_failure(100)
    first_multiplicative_failure = first_arithmetic_multiplicativity_failure(30)
    compiler_checks, star_rows = audit_compiler(1000)
    audit_exact_hostiles()

    print("status=INDEPENDENT FINITE-EXACT hostile audit")
    print("definition=epsilon_q(a)=(-1)^D(a/q);epsilon_1(0)=1")
    print(f"modular_inverse_checks={inverse_checks};q_max=2000")
    print("inverse_formula=q_odd:(-1)^(a+abar);q_even:-(-1)^((a*abar-1)/q)")
    print("symmetries=epsilon(a)=epsilon(abar)=epsilon(q-a)")
    print("packet_congruence=S(q)=phi(q)_mod_4_for_q>2")
    print("odd_q_half_shells=odd_numerator_sum=even_numerator_sum=S(q)/2")
    print(f"first_group_character_failure={first_group_failure}")
    print(f"first_odd_group_character_failure={first_odd_group_failure}")
    print(f"first_arithmetic_multiplicativity_failure={first_multiplicative_failure}")
    print(f"compiler_levels={compiler_checks};identity=A=1*S;star=B=A-1")
    print(f"small_rows_(q,S,B)={star_rows}")
    print("compiled_C4_hostile=(+,-,-,-)_is_not_additive_character")
    print("PASS")


if __name__ == "__main__":
    main()

