#!/usr/bin/env python3
"""Exact companion for THM-3130's divisor-antichain response laws."""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc_grid_witness_antichain_response_thm3130.py"
OUTPUT = ROOT / "05-knowledge/results/lrc_grid_witness_antichain_response_thm3130.out"
MAXIMUM_Q = 84


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def prime_powers(n):
    answer = []
    p = 2
    while p * p <= n:
        if n % p:
            p += 1
            continue
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        answer.append((p, exponent))
        p += 1
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def phi(n):
    answer = n
    for p, _ in prime_powers(n):
        answer = answer // p * (p - 1)
    return answer


def lcm(left, right):
    return left // gcd(left, right) * right


def principal_tail_direct(q, m):
    require(q % m == 0, ("nondivisor tail", q, m))
    return sum(phi(d) for d in divisors(q) if d % m == 0)


def principal_tail_product(q, m):
    require(q % m == 0, ("nondivisor product tail", q, m))
    answer = 1
    residual = m
    for p, exponent in prime_powers(q):
        b = 0
        while residual % p == 0:
            residual //= p
            b += 1
        require(b <= exponent, (q, m, p, b, exponent))
        answer *= p**exponent if b == 0 else p**exponent - p ** (b - 1)
    require(residual == 1, ("tail residual", q, m, residual))
    return answer


def enumerate_upsets(ds):
    answer = []
    for mask in range(1 << len(ds)):
        if all(
            not (mask >> i) & 1
            or all((mask >> j) & 1 for j, e in enumerate(ds) if e % d == 0)
            for i, d in enumerate(ds)
        ):
            answer.append(mask)
    return tuple(answer)


def minimal_elements(mask, ds):
    return tuple(
        d
        for i, d in enumerate(ds)
        if (mask >> i) & 1
        and not any(
            e != d and d % e == 0 and (mask >> j) & 1
            for j, e in enumerate(ds)
        )
    )


def subset_lcms(generators):
    rows = []
    for mask in range(1, 1 << len(generators)):
        value = 1
        size = 0
        for index, generator in enumerate(generators):
            if (mask >> index) & 1:
                value = lcm(value, generator)
                size += 1
        rows.append((value, 1 if size % 2 else -1))
    return tuple(rows)


def response_coefficients(ds, generators, ceiling):
    coefficients = [0] * len(ds)
    for multiple, sign in subset_lcms(generators):
        require(ceiling % multiple == 0, (ceiling, generators, multiple))
        for index, d in enumerate(ds):
            if ceiling % d == 0 and d % multiple == 0:
                coefficients[index] += sign
    return tuple(coefficients)


def totient_mass(mask, ds):
    return sum(phi(d) for index, d in enumerate(ds) if (mask >> index) & 1)


def lost_mask(mask, ds, g):
    return sum(
        1 << index
        for index, d in enumerate(ds)
        if (mask >> index) & 1 and g % d == 0
    )


def run_audit():
    upset_count = static_vector_checks = tail_checks = 0
    update_checks = supermodular_checks = 0
    maximum_antichain = (0, None)
    semantic_rows = []

    for q in range(2, MAXIMUM_Q + 1):
        ds = divisors(q)
        upsets = enumerate_upsets(ds)
        weights = tuple(phi(d) for d in ds)
        require(sum(weights) == q, ("totient divisor identity", q))
        for m in ds:
            direct = principal_tail_direct(q, m)
            product = principal_tail_product(q, m)
            require(direct == product, ("principal tail", q, m, direct, product))
            tail_checks += 1

        q_rows = []
        for mask in upsets:
            # THM-3055 states have 1 hit and hence exclude the bottom divisor.
            if mask == 0 or mask & 1:
                continue
            generators = minimal_elements(mask, ds)
            if len(generators) > maximum_antichain[0]:
                maximum_antichain = (len(generators), (q, generators))
            coefficients = response_coefficients(ds, generators, q)
            indicator = tuple((mask >> index) & 1 for index in range(len(ds)))
            require(coefficients == indicator,
                    ("static response vector", q, generators, coefficients, indicator))
            static_vector_checks += 1
            inclusion_mass = sum(
                sign * principal_tail_product(q, multiple)
                for multiple, sign in subset_lcms(generators)
            )
            require(inclusion_mass == totient_mass(mask, ds),
                    ("static totient response", q, generators))
            upset_count += 1

            update_rows = []
            for g in ds:
                lost = lost_mask(mask, ds, g)
                selected = tuple(m for m in generators if g % m == 0)
                loss_coefficients = response_coefficients(ds, selected, g)
                expected = tuple(
                    1 if ((lost >> index) & 1) else 0 for index in range(len(ds))
                )
                # response_coefficients is supported only on divisors of g.
                require(loss_coefficients == expected,
                        ("insertion loss vector", q, generators, g))
                loss_formula = sum(
                    sign * principal_tail_product(g, multiple)
                    for multiple, sign in subset_lcms(selected)
                )
                loss_direct = totient_mass(lost, ds)
                require(loss_formula == loss_direct,
                        ("insertion loss mass", q, generators, g))
                update_checks += 1
                update_rows.append((g, loss_direct))
            q_rows.append((generators, totient_mass(mask, ds), tuple(update_rows)))

        # Exact diminishing-loss form of supermodularity on every nested pair.
        live = tuple(mask for mask in upsets if mask and not mask & 1)
        for large in live:
            for small in live:
                if small & ~large:
                    continue
                for g in ds:
                    require(
                        totient_mass(lost_mask(small, ds, g), ds)
                        <= totient_mass(lost_mask(large, ds, g), ds),
                        ("supermodular loss", q, small, large, g),
                    )
                    supermodular_checks += 1
        semantic_rows.append((q, ds, tuple(q_rows)))

    # Positive antichain/update control: q=12, U generated by {3,4}.
    require(principal_tail_product(12, 3) == 8, "q12 T(3)")
    require(principal_tail_product(12, 4) == 6, "q12 T(4)")
    require(principal_tail_product(12, 12) == 4, "q12 T(12)")
    require(8 + 6 - 4 == 10, "q12 static inclusion-exclusion")
    require(principal_tail_product(6, 3) == 4, "q12 insertion loss at gcd 6")

    # Sharp hostile: normalized witness distributions need not move upward.
    mu_before = Fraction(phi(2) + phi(6), phi(2) + phi(3) + phi(6))
    mu_after = Fraction(phi(6), phi(3) + phi(6))
    require((mu_before, mu_after, mu_after - mu_before)
            == (Fraction(3, 5), Fraction(1, 2), Fraction(-1, 10)),
            "q6 normalized-transport hostile")

    semantic = (
        tuple(semantic_rows), upset_count, static_vector_checks, tail_checks,
        update_checks, supermodular_checks, maximum_antichain,
        (mu_before, mu_after),
    )
    semantic_sha = hashlib.sha256(repr(semantic).encode()).hexdigest()
    return {
        "upset_count": upset_count,
        "static_vector_checks": static_vector_checks,
        "tail_checks": tail_checks,
        "update_checks": update_checks,
        "supermodular_checks": supermodular_checks,
        "maximum_antichain": maximum_antichain,
        "semantic_sha": semantic_sha,
    }


def render(result):
    return "\n".join((
        "LRC grid-witness divisor-antichain response and supermodular loss",
        f"source_sha256={lf_sha256(SOURCE)}",
        f"exact_universe=q2..{MAXIMUM_Q};nonempty_upsets_excluding_bottom={result['upset_count']}",
        f"principal_tail_product_checks={result['tail_checks']}",
        f"static_antichain_response_vectors={result['static_vector_checks']}",
        f"speed_insertion_loss_checks={result['update_checks']}",
        f"nested_supermodular_loss_checks={result['supermodular_checks']}",
        f"maximum_antichain_control={result['maximum_antichain']}",
        "q12_control=M:(3,4);T12:8+6-4=10;g6_loss:T6(3)=4",
        "q6_hostile=U:(2,3,6)->(3,6);upset:(2,6);mass:3/5->1/2;delta:-1/10",
        "scope=exact_grid witness response;not continuous safe mass or normalized Hasse transport",
        f"semantic_sha256={result['semantic_sha']}",
        "all_exact_checks=PASS",
    )) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_audit()
    payload = render(result)
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload)


if __name__ == "__main__":
    main()
