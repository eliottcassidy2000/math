#!/usr/bin/env python3
"""Exact companion for THM-3166.

The direct engine counts Hamiltonian paths on every vertex subset and then
unordered spanning path covers by their least-vertex component.  It audits the
order-join matching convolution, its falling-factorial conjugacy, the SCC
product, the Laguerre basis kernel, and fixed-depth exponential closed forms.
"""

from __future__ import annotations

import argparse
import hashlib
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
OUTPUT = ROOT / "05-knowledge/results/tournament_order_join_falling_factorial_transform_thm3166.out"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def edge_index(n: int, i: int, j: int) -> int:
    require(0 <= i < j < n, (n, i, j))
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def beats(n: int, code: int, i: int, j: int) -> bool:
    if i < j:
        return bool((code >> edge_index(n, i, j)) & 1)
    return not bool((code >> edge_index(n, j, i)) & 1)


def adjacency(n: int, code: int):
    return tuple(
        tuple(i != j and beats(n, code, i, j) for j in range(n))
        for i in range(n)
    )


@lru_cache(maxsize=None)
def path_cover_profile(n: int, code: int):
    """Return pc[d], using an endpoint DP and a least-block partition DP."""
    require(n >= 1 and 0 <= code < (1 << comb(n, 2)), (n, code))
    adj = adjacency(n, code)
    size = 1 << n
    endpoint = [[0] * n for _ in range(size)]
    for v in range(n):
        endpoint[1 << v][v] = 1
    hamilton = [0] * size
    for mask in range(1, size):
        if mask & (mask - 1):
            for last in range(n):
                bit = 1 << last
                if not mask & bit:
                    continue
                previous = mask ^ bit
                endpoint[mask][last] = sum(
                    endpoint[previous][v]
                    for v in range(n)
                    if previous & (1 << v) and adj[v][last]
                )
        hamilton[mask] = sum(endpoint[mask])

    covers = [[0] * (n + 1) for _ in range(size)]
    covers[0][0] = 1
    for mask in range(1, size):
        anchor = mask & -mask
        sub = mask
        while sub:
            if sub & anchor:
                weight = hamilton[sub]
                if weight:
                    rest = mask ^ sub
                    for d in range(n):
                        if covers[rest][d]:
                            covers[mask][d + 1] += weight * covers[rest][d]
            sub = (sub - 1) & mask
    profile = tuple(covers[-1])
    require(profile[0] == 0 and profile[n] == 1, (n, code, profile))
    return profile


def induced_code(n: int, code: int, vertices):
    vertices = tuple(vertices)
    out = 0
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            if beats(n, code, vertices[i], vertices[j]):
                out |= 1 << edge_index(len(vertices), i, j)
    return out


def join_code(n: int, first: int, m: int, second: int) -> int:
    total = n + m
    out = 0
    for i in range(total):
        for j in range(i + 1, total):
            value = (
                beats(n, first, i, j)
                if j < n
                else beats(m, second, i - n, j - n)
                if i >= n
                else True
            )
            if value:
                out |= 1 << edge_index(total, i, j)
    return out


def matching_convolution(first, second):
    out = [0] * (len(first) + len(second) - 1)
    for a in range(1, len(first)):
        for b in range(1, len(second)):
            for k in range(min(a, b) + 1):
                out[a + b - k] += (
                    first[a] * second[b] * comb(a, k) * comb(b, k) * factorial(k)
                )
    return tuple(out)


def falling_value(t: int, d: int) -> int:
    value = 1
    for j in range(d):
        value *= t - j
    return value


def path_colour_value(profile, t: int) -> int:
    return sum(profile[d] * falling_value(t, d) for d in range(1, len(profile)))


def backward_distribution(n: int, code: int):
    """Count permutations by backward consecutive adjacencies."""
    counts = [0] * n
    for order in permutations(range(n)):
        backward = sum(
            not beats(n, code, order[i], order[i + 1])
            for i in range(n - 1)
        )
        counts[backward] += 1
    return tuple(counts)


def negative_colour_value(backward, m: int) -> int:
    """The positive reciprocity coordinate (-1)^n Q_T(-m)."""
    n = len(backward)
    return sum(backward[b] * comb(m + b, n) for b in range(n))


def recover_backward_from_negative(profile):
    """Triangular inversion from Q_T(-1),...,Q_T(-n)."""
    n = len(profile) - 1
    recovered = [0] * n
    for m in range(1, n + 1):
        b = n - m
        value = (-1) ** n * path_colour_value(profile, -m)
        recovered[b] = value - sum(
            recovered[c] * comb(m + c, n)
            for c in range(b + 1, n)
        )
    return tuple(recovered)


def falling_power_rows(n: int):
    rows = [(1,)]
    for d in range(1, n + 1):
        old = rows[-1]
        row = [0] * (d + 1)
        for j, value in enumerate(old):
            row[j + 1] += value
            row[j] -= (d - 1) * value
        rows.append(tuple(row))
    return tuple(rows)


def power_polynomial(profile):
    rows = falling_power_rows(len(profile) - 1)
    out = [0] * len(profile)
    for d in range(1, len(profile)):
        for j, value in enumerate(rows[d]):
            out[j] += profile[d] * value
    return tuple(out)


def polynomial_product(first, second):
    out = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            out[i + j] += a * b
    return tuple(out)


def polynomial_value(coefficients, t: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = value * t + coefficient
    return value


def profile_from_power(coefficients):
    n = len(coefficients) - 1
    values = [polynomial_value(coefficients, j) for j in range(n + 1)]
    profile = [0] * (n + 1)
    for d in range(n + 1):
        difference = sum(
            (-1) ** (d - j) * comb(d, j) * values[j]
            for j in range(d + 1)
        )
        require(difference % factorial(d) == 0, (d, difference))
        profile[d] = difference // factorial(d)
    return tuple(profile)


def scc_components(n: int, code: int):
    reach = [[i == j or beats(n, code, i, j) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    unseen = set(range(n))
    components = []
    while unseen:
        first = min(unseen)
        component = tuple(v for v in sorted(unseen) if reach[first][v] and reach[v][first])
        unseen.difference_update(component)
        components.append(component)
    components.sort(key=lambda block: -sum(beats(n, code, block[0], v)
                                           for other in components if other != block
                                           for v in other))
    return tuple(components)


def repeated_closed_form(profile, repetitions: int, d: int) -> int:
    numerator = sum(
        (-1) ** (d - j) * comb(d, j) * path_colour_value(profile, j) ** repetitions
        for j in range(d + 1)
    )
    require(numerator % factorial(d) == 0, (repetitions, d, numerator))
    return numerator // factorial(d)


def laguerre_kernel(a: int, b: int):
    require(1 <= a <= b, (a, b))
    out = [0] * (a + b + 1)
    for j in range(a + 1):
        require(factorial(a) % factorial(j) == 0, (a, j))
        out[b + j] = factorial(a) // factorial(j) * comb(b, a - j)
    return tuple(out)


def basis_matching_kernel(a: int, b: int):
    out = [0] * (a + b + 1)
    for k in range(min(a, b) + 1):
        out[a + b - k] = comb(a, k) * comb(b, k) * factorial(k)
    return tuple(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    # All labelled tournaments through order five.
    tournament_checks = 0
    negative_reciprocity_checks = 0
    negative_inverse_checks = 0
    scc_checks = 0
    for n in range(1, 6):
        for code in range(1 << comb(n, 2)):
            profile = path_cover_profile(n, code)
            power = power_polynomial(profile)
            require(profile_from_power(power) == profile, ("transform inverse", n, code))
            require(path_colour_value(profile, 1) == profile[1], ("endpoint", n, code))
            backward = backward_distribution(n, code)
            for m in range(1, 8):
                require(
                    (-1) ** n * path_colour_value(profile, -m)
                    == negative_colour_value(backward, m),
                    ("negative-colour reciprocity", n, code, m),
                )
                negative_reciprocity_checks += 1
            require(recover_backward_from_negative(profile) == backward,
                    ("negative-colour inverse", n, code))
            negative_inverse_checks += 1
            tournament_checks += 1

            product = (1,)
            for component in scc_components(n, code):
                component_profile = path_cover_profile(
                    len(component), induced_code(n, code, component)
                )
                product = polynomial_product(product, power_polynomial(component_profile))
            require(product == power, ("SCC product", n, code, scc_components(n, code)))
            scc_checks += 1
    require(tournament_checks == scc_checks == 1099, (tournament_checks, scc_checks))
    require(negative_reciprocity_checks == 7693, negative_reciprocity_checks)
    require(negative_inverse_checks == 1099, negative_inverse_checks)

    # Every ordered pair of labelled tournaments through order four.
    small = tuple(
        (n, code, path_cover_profile(n, code))
        for n in range(1, 5)
        for code in range(1 << comb(n, 2))
    )
    require(len(small) == 75, len(small))
    join_checks = 0
    for n, first_code, first_profile in small:
        for m, second_code, second_profile in small:
            joined = path_cover_profile(n + m, join_code(n, first_code, m, second_code))
            convolved = matching_convolution(first_profile, second_profile)
            require(joined == convolved,
                    ("join convolution", n, first_code, m, second_code))
            product = polynomial_product(
                power_polynomial(first_profile), power_polynomial(second_profile)
            )
            require(power_polynomial(joined) == product, ("join product", n, m))
            require(profile_from_power(product) == joined, ("product inverse", n, m))
            for t in range(0, n + m + 2):
                require(path_colour_value(joined, t) ==
                        path_colour_value(first_profile, t) *
                        path_colour_value(second_profile, t),
                        ("integer colour product", n, m, t))
            join_checks += 1
    require(join_checks == 5625, join_checks)

    laguerre_checks = 0
    for a in range(1, 21):
        for b in range(a, 21):
            require(laguerre_kernel(a, b) == basis_matching_kernel(a, b), (a, b))
            laguerre_checks += 1
    require(laguerre_checks == 210, laguerre_checks)

    repeated_checks = 0
    for _, _, profile in small:
        running = (1,)
        for repetitions in range(1, 6):
            running = matching_convolution(running, profile) if len(running) > 1 else profile
            for d in range(1, len(running)):
                require(repeated_closed_form(profile, repetitions, d) == running[d],
                        ("repeated closed form", profile, repetitions, d))
                repeated_checks += 1

    k1 = path_cover_profile(1, 0)
    c3 = path_cover_profile(3, 5)
    transitive3 = path_cover_profile(3, 7)
    source_c3 = path_cover_profile(4, join_code(1, 0, 3, 5))
    c3_sink = path_cover_profile(4, join_code(3, 5, 1, 0))
    require(k1 == (0, 1) and c3 == (0, 3, 3, 1), (k1, c3))
    require(source_c3 == c3_sink == matching_convolution(k1, c3),
            (source_c3, c3_sink))
    require(power_polynomial(transitive3) == (0, 0, 0, 1), transitive3)
    require(power_polynomial(c3) == (0, 2, 0, 1), c3)

    lines = [
        "TOURNAMENT ORDER-JOIN FALLING-FACTORIAL TRANSFORM -- THM-3166 EXACT COMPANION",
        "status=PROVED+VERIFIED-EXACT;direct subset Hamilton-path/path-cover engine",
        f"labelled_tournament_controls={tournament_checks};orders=1..5;transform_inverse=PASS;endpoint=PASS",
        f"negative_colour_reciprocity_controls={negative_reciprocity_checks};m=1..7;positive_binomial_sum=PASS",
        f"negative_colour_inverse_controls={negative_inverse_checks};Q(-1..-n)_to_backward_distribution=PASS",
        f"SCC_factorization_controls={scc_checks};full_profile_transform_product=PASS",
        f"ordered_join_controls={join_checks};all ordered pairs of 75 labelled tournaments through order 4",
        "join_kernel=sum_k binom(a,k)binom(b,k)k!*x^(a+b-k);direct_profiles=PASS",
        f"Laguerre_basis_controls={laguerre_checks};a<=b<=20;kernel=a!*x^b*L_a^(b-a)(-x)",
        f"repeated_join_closed_form_controls={repeated_checks};75 seeds;r=1..5;all depths=PASS",
        "transform=Q_T(t)=sum_d pc_T(d)*(t)_d;Q_(A join B)=Q_A*Q_B",
        "negative_reciprocity=(-1)^n*Q_T(-m)=sum_pi binom(m+b_T(pi),n);m>=1",
        "negative_endpoint=Q_T(-1)=(-1)^n*H(T);negative_values_recover_full_backward_jet",
        "inverse=pc_T(d)=Delta^d Q_T(0)/d!;H(T)=Q_T(1)",
        "repeated=pc_(A^join r)(d)=1/d!*sum_j(-1)^(d-j)binom(d,j)Q_A(j)^r",
        f"hostile_order_loss=K1_join_C3 and C3_join_K1 share profile {source_c3}",
        f"cyclic_boundary=Q_transitive3={power_polynomial(transitive3)};Q_C3={power_polynomial(c3)}",
        "scope=order-join/SCC chains;not cyclic substitution;not arbitrary tournament complexity collapse",
        "reproduction=python3 04-computation/tournament_order_join_falling_factorial_transform_thm3166.py",
        f"source_sha256={sha256(HERE)}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
