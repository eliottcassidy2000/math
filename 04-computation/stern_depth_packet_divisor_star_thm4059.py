"""Exact companion for THM-4059.

Checks the modular-inverse Stern-depth parity formula, exact denominator
packets, inverse/complement symmetries, the mod-four packet law, the lawful
divisor-star convolution, and minimal noncharacter/nonmultiplicativity
hostiles.

Reproduce:
    python -B 04-computation/stern_depth_packet_divisor_star_thm4059.py
    python -B -O 04-computation/stern_depth_packet_divisor_star_thm4059.py
"""

from __future__ import annotations

from math import gcd


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def stern_depth(p: int, q: int) -> int:
    """Sum of Euclidean quotients minus one."""
    require(p > 0 and q > 0 and gcd(p, q) == 1, (p, q))
    total = 0
    while q:
        total += p // q
        p, q = q, p % q
    return total - 1


def inverse_packet_data(p: int, q: int) -> tuple[int, int]:
    require(q > 1 and gcd(p, q) == 1, (p, q))
    r = pow(p, -1, q)
    h = (p * r - 1) // q
    return r, h


def parity_formula(p: int, q: int) -> int:
    if q == 1:
        return (p - 1) & 1
    r, h = inverse_packet_data(p, q)
    return (p * q + p * r + r * h) & 1


def packet_direct(q: int) -> int:
    return sum((-1) ** stern_depth(a, q)
               for a in range(1, q) if gcd(a, q) == 1)


def packet_formula(q: int) -> int:
    total = 0
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        r, h = inverse_packet_data(a, q)
        if q & 1:
            exponent = a + r
        else:
            exponent = h + 1
        total += (-1) ** exponent
    return total


def star_direct(n: int) -> int:
    total = 0
    for a in range(1, n):
        g = gcd(a, n)
        total += (-1) ** stern_depth(a // g, n // g)
    return total


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def main() -> None:
    pair_bound = 700
    pair_checks = 0
    for p in range(1, pair_bound + 1):
        for q in range(1, pair_bound + 1):
            if gcd(p, q) != 1:
                continue
            depth = stern_depth(p, q)
            require((depth & 1) == parity_formula(p, q),
                    (p, q, depth, parity_formula(p, q)))
            if q > 1:
                r, h = inverse_packet_data(p, q)
                # The two columns are the unique Farey parents.
                require(0 <= h < p and 1 <= r < q, (p, q, r, h))
                require(h * (q - r) - r * (p - h) == -1,
                        (p, q, r, h))
            pair_checks += 1

    q_bound = 3000
    packets = {}
    for q in range(2, q_bound + 1):
        direct = packet_direct(q)
        formula = packet_formula(q)
        require(direct == formula, (q, direct, formula))
        packets[q] = direct
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            r = pow(a, -1, q)
            require(stern_depth(a, q) == stern_depth(r, q),
                    ("inverse exact depth", q, a, r))
            require(stern_depth(a, q) == stern_depth(q - a, q),
                    ("subtree mirror", q, a))
        if q > 2:
            require(direct % 2 == 0, ("complement pairs", q, direct))
            phi_q = sum(1 for a in range(1, q) if gcd(a, q) == 1)
            require((direct - phi_q) % 4 == 0,
                    ("packet mod four", q, direct, phi_q))

    star_bound = 3000
    for n in range(2, star_bound + 1):
        direct = star_direct(n)
        convolution = sum(packets[d] for d in divisors(n) if d > 1)
        require(direct == convolution, (n, direct, convolution))

    # Minimal hostile controls.
    require(stern_depth(2, 3) % 2 != stern_depth(2, 5) % 2,
            "endpoint parity alone")
    require(stern_depth(1, 8) % 2 != stern_depth(3, 8) % 2,
            "even-q inverse parity alone")
    chi9 = {a: (-1) ** stern_depth(a, 9)
            for a in range(1, 9) if gcd(a, 9) == 1}
    require(chi9[2] * chi9[2] != chi9[4], "not a U_9 character")
    require(packets[10] != packets[2] * packets[5],
            "packet arithmetic function not multiplicative")

    print("status=PROVED identities with FINITE-EXACT audit")
    print(f"coprime_pairs_checked={pair_checks};box={pair_bound}x{pair_bound}")
    print(f"packet_denominators_checked=2..{q_bound}")
    print(f"divisor_stars_checked=2..{star_bound}")
    print("parity=D(p/q)=p*q+p*r+r*h mod 2; p*r-q*h=1")
    print("odd_q_packet=sum_a (-1)^(a+a_inverse)")
    print("even_q_packet=-sum_a (-1)^((a*a_inverse-1)/q)")
    print("star(n)=sum_(d|n,d>1) packet(d)")
    print("packet_values_q2_to_q40="
          + repr(tuple(packets[q] for q in range(2, 41))))
    print("hostiles=endpoint_parity:(2/3,2/5);"
          "even_inverse_parity:(1/8,3/8);"
          "U9_noncharacter:2*2=4;"
          "nonmultiplicative:S10!=S2*S5")


if __name__ == "__main__":
    main()

