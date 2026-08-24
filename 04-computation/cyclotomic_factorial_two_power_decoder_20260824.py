#!/usr/bin/env python3
"""Exact controls for the cyclotomic factorial two-power decoder.

For F=a*x+b*y, let Pi_d average y -> zeta_d^r*y over r mod d and let
L(x^i*y^j)=i!*j!.  The normalized projected moments are

    M_m = L(Pi_d(F^(d*m))) / (d*m)!
        = sum_{j=0}^m (a^d)^(m-j) (b^d)^j.

The first two moments recover the unordered pair (a,b).  This companion
checks the coefficient identity, recurrence, decoder, the nonmultiplicative
projector boundary, and complete bounded two-cube collision atlases.  It
uses no Python assertions, so every gate remains active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, gcd, isqrt
import json
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_positive_nth_root(value: int, exponent: int) -> int | None:
    require(value >= 0 and exponent >= 1, "nth-root domain")
    if value in (0, 1):
        return value
    low, high = 1, 2
    while high ** exponent < value:
        high *= 2
    while low <= high:
        middle = (low + high) // 2
        power = middle ** exponent
        if power == value:
            return middle
        if power < value:
            low = middle + 1
        else:
            high = middle - 1
    return None


def projected_factorial_moment_raw(power: int, a: int, b: int, moment: int) -> int:
    """Apply the degree-divisible projector and the bivariate factorial functional."""
    degree = power * moment
    total = 0
    for y_degree in range(0, degree + 1, power):
        x_degree = degree - y_degree
        coefficient = comb(degree, y_degree) * a ** x_degree * b ** y_degree
        total += coefficient * factorial(x_degree) * factorial(y_degree)
    return total


def complete_homogeneous_two(u: int, v: int, moment: int) -> int:
    return sum(u ** (moment - index) * v ** index for index in range(moment + 1))


def normalized_projected_moment(power: int, a: int, b: int, moment: int) -> int:
    degree = power * moment
    raw = projected_factorial_moment_raw(power, a, b, moment)
    scale = factorial(degree)
    require(raw % scale == 0, "factorial normalization integrality")
    return raw // scale


def decode_pair(power: int, first: int, second: int) -> tuple[int, int] | None:
    """Recover positive distinct a<b from M_1,M_2, or reject the packet."""
    product = first * first - second
    if product <= 0:
        return None
    discriminant = first * first - 4 * product
    if discriminant <= 0:
        return None
    delta = isqrt(discriminant)
    if delta * delta != discriminant or (first - delta) % 2:
        return None
    u = (first - delta) // 2
    v = (first + delta) // 2
    if not 0 < u < v or u * v != product:
        return None
    a = exact_positive_nth_root(u, power)
    b = exact_positive_nth_root(v, power)
    if a is None or b is None or not a < b:
        return None
    return a, b


def factor(n: int) -> dict[int, int]:
    require(n >= 1, "factor domain")
    result: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            result[divisor] = result.get(divisor, 0) + 1
            n //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def admissible_cube_shell(shell: int) -> bool:
    factors = factor(shell)
    return shell >= 3 and bool(factors) and all(
        prime % 3 == 2 and exponent <= 2
        for prime, exponent in factors.items()
    )


def main() -> None:
    active_checks = 0

    # General all-power controls.  The proof is coefficientwise; this bounded
    # grid is a hostile implementation audit, not the source of the quantifier.
    power_range = range(2, 10)
    pair_cap = 80
    moment_cap = 8
    for power in power_range:
        for b in range(2, pair_cap + 1):
            for a in range(1, b):
                u, v = a ** power, b ** power
                moments = [1]
                for moment in range(1, moment_cap + 1):
                    measured = normalized_projected_moment(power, a, b, moment)
                    closed = complete_homogeneous_two(u, v, moment)
                    require(measured == closed, f"projected factorial identity {(power, a, b, moment)}")
                    moments.append(measured)
                    active_checks += 2
                require(decode_pair(power, moments[1], moments[2]) == (a, b),
                        f"two-moment decoder {(power, a, b)}")
                active_checks += 1
                product = moments[1] * moments[1] - moments[2]
                require(product == u * v, "second elementary symmetric coordinate")
                active_checks += 1
                for moment in range(2, moment_cap + 1):
                    require(
                        moments[moment]
                        == moments[1] * moments[moment - 1] - product * moments[moment - 2],
                        f"rank-two recurrence {(power, a, b, moment)}",
                    )
                    active_checks += 1

    # Pi_d is a lawful orbit average but is not an algebra homomorphism.
    # The x^d y^d coefficient changes from 2 to binomial(2d,d).
    projector_hostiles: list[tuple[int, int, int]] = []
    for power in power_range:
        product_cross = 2
        power_cross = comb(2 * power, power)
        require(power_cross != product_cross, "projector multiplicativity hostile disappeared")
        projector_hostiles.append((power, product_cross, power_cross))
        active_checks += 1

    # Direct functional hostile at the distinct pair (1,2), d=3.
    # The projected response, powers of the invariantized cubic, and ordinary
    # powers of the original cubic give three different normalized moments.
    boundary_a, boundary_b, boundary_power = 1, 2, 3
    boundary_projected = Fraction(
        projected_factorial_moment_raw(boundary_power, boundary_a, boundary_b, 2),
        factorial(2 * boundary_power),
    )
    boundary_fixed_invariant = Fraction(
        factorial(6) * boundary_a ** 6
        + 2 * boundary_a ** 3 * boundary_b ** 3 * factorial(3) ** 2
        + factorial(6) * boundary_b ** 6,
        factorial(6),
    )
    boundary_original = sum(Fraction(boundary_a ** (6 - k) * boundary_b ** k)
                            for k in range(7))
    require(boundary_projected == 73, "projected functional hostile value")
    require(boundary_fixed_invariant == Fraction(329, 5), "fixed invariant hostile value")
    require(boundary_original == 127, "ordinary factorial moment hostile value")
    require(len({boundary_projected, boundary_fixed_invariant, boundary_original}) == 3,
            "functional scope hostiles collided")
    active_checks += 4

    # Complete support-two atlas from THM-3743/3825.  One cube-sum coordinate
    # collides on the full atlas; the moment pair is injective everywhere.
    shell_cap = 356
    full_nodes = [
        (a, shell - a)
        for shell in range(3, shell_cap + 1)
        for a in range(1, (shell + 1) // 2)
        if a < shell - a and gcd(a, shell) == 1
    ]
    require(len(full_nodes) == 19_314, "support-two atlas count")
    cube_fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    moment_packets: set[tuple[int, int]] = set()
    selected_cube_values: set[int] = set()
    selected_nodes = 0
    for a, b in full_nodes:
        first = a ** 3 + b ** 3
        second = a ** 6 + a ** 3 * b ** 3 + b ** 6
        cube_fibres[first].append((a, b))
        require((first, second) not in moment_packets, "two-moment packet collision")
        moment_packets.add((first, second))
        require(decode_pair(3, first, second) == (a, b), "atlas decoder")
        if admissible_cube_shell(a + b):
            require(first not in selected_cube_values, "THM-3825 selected scalar collision")
            selected_cube_values.add(first)
            selected_nodes += 1
        active_checks += 4

    collisions = {
        value: sorted(pairs)
        for value, pairs in cube_fibres.items()
        if len(pairs) > 1
    }
    collision_histogram = Counter(len(pairs) for pairs in cube_fibres.values())
    require(len(collisions) == 28, "full-atlas cube collision count")
    require(collision_histogram == Counter({1: 19_258, 2: 28}), "full-atlas fibre histogram")
    require(len(moment_packets) == len(full_nodes), "moment-packet injection count")
    require(selected_nodes == len(selected_cube_values) == 5_855, "selected THM-3825 count")
    active_checks += 5

    # Complete positive-pair collision audit through coordinate 300.
    coordinate_cap = 300
    bounded_fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for b in range(2, coordinate_cap + 1):
        for a in range(1, b):
            bounded_fibres[a ** 3 + b ** 3].append((a, b))
    bounded_collisions = {
        value: sorted(pairs)
        for value, pairs in bounded_fibres.items()
        if len(pairs) > 1
    }
    first_collision_value = min(bounded_collisions)
    first_collision_pairs = bounded_collisions[first_collision_value]
    require(first_collision_value == 1729, "first positive distinct cube collision")
    require(first_collision_pairs == [(1, 12), (9, 10)], "taxicab hostile fibre")
    first_second_moments = {
        pair: pair[0] ** 6 + pair[0] ** 3 * pair[1] ** 3 + pair[1] ** 6
        for pair in first_collision_pairs
    }
    require(len(set(first_second_moments.values())) == 2, "second moment did not split 1729")
    for value, pairs in bounded_collisions.items():
        seconds = [a ** 6 + a ** 3 * b ** 3 + b ** 6 for a, b in pairs]
        require(len(seconds) == len(set(seconds)), f"second-moment collision at {value}")
        for pair, second in zip(pairs, seconds):
            require(decode_pair(3, value, second) == pair, f"bounded collision decode at {value}")
            active_checks += 2

    # Rejection controls: not every integer pair is a two-power packet.
    rejection_packets = [(3, 10, 20), (3, 1729, 1729 ** 2), (4, 17, 200)]
    for power, first, second in rejection_packets:
        require(decode_pair(power, first, second) is None, "invalid packet accepted")
        active_checks += 1

    semantic_rows = [
        {
            "value": value,
            "pairs": pairs,
            "second_moments": [
                a ** 6 + a ** 3 * b ** 3 + b ** 6 for a, b in pairs
            ],
        }
        for value, pairs in sorted(collisions.items())
    ]
    semantic_digest = sha256(
        "\n".join(json.dumps(row, sort_keys=True) for row in semantic_rows).encode("ascii")
    ).hexdigest()

    print("CYCLOTOMIC FACTORIAL TWO-POWER DECODER")
    print("THEOREM M_m=L(Pi_d((a*x+b*y)^(d*m)))/(d*m)! = h_m(a^d,b^d)")
    print("DECODER u+v=M_1; u*v=M_1^2-M_2; {u,v}={a^d,b^d}")
    print("ALL_POWER_CONTROL", "d=2..9,a<b<=80,m<=8")
    print("PROJECTOR_NONMULTIPLICATIVE", "d=3:x^3*y^3 coefficient 2 versus 20")
    print("PROJECTOR_HOSTILES", projector_hostiles)
    print("FUNCTIONAL_SCOPE_HOSTILE", "d=3,(a,b)=(1,2): projected=73,fixed_invariant=329/5,ordinary=127")
    print("SUPPORT_TWO_UNIVERSE", len(full_nodes))
    print("SUPPORT_TWO_CUBE_FIBRE_HISTOGRAM", dict(sorted(collision_histogram.items())))
    print("SUPPORT_TWO_COLLISION_VALUES", len(collisions))
    print("SUPPORT_TWO_MOMENT_PACKETS", len(moment_packets))
    print("THM3825_SELECTED_ONE_MOMENT_PACKETS", selected_nodes)
    print("FIRST_CUBE_COLLISION", first_collision_value, first_collision_pairs)
    print("FIRST_COLLISION_SECOND_MOMENTS", sorted(first_second_moments.items()))
    print("COORDINATE_CAP", coordinate_cap)
    print("BOUNDED_CUBE_COLLISION_VALUES", len(bounded_collisions))
    print("BOUNDED_MAX_FIBRE", max(map(len, bounded_collisions.values())))
    print("SEMANTIC_SHA256", semantic_digest)
    print("ACTIVE_CHECKS", active_checks)
    print("SCOPE uses (L o Pi_d)((F^d)^m), not original factorial moments L(f^m); no FC/HFC/JC/LRC consequence")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
