#!/usr/bin/env python3
"""Independent exact certificates for the unit-power-sum detector boundary.

Pure Python, standard library only. Universe for the complete census: all
2^((p-1)/2) labelled tournament connection sets, for every odd prime p <= 31.
Construction checks: every bad-regime prime <= --bound (default 1000).
No numerical Fourier transforms, random searches, or external state writes.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path


def factor(n):
    result = {}
    q = 2
    while q * q <= n:
        while n % q == 0:
            result[q] = result.get(q, 0) + 1
            n //= q
        q += 1
    if n > 1:
        result[n] = 1
    return result


def primes(bound):
    return [p for p in range(3, bound + 1, 2)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def primitive_root(p):
    qs = factor(p - 1)
    return next(g for g in range(2, p)
                if all(pow(g, (p - 1) // q, p) != 1 for q in qs))


def good_regime(p):
    factors = factor(p - 1)
    e = factors.pop(2)
    return (len(factors) <= 1 or
            (e == 1 and len(factors) == 2 and
             all(a == 1 for a in factors.values())))


def point_sign(xs):
    return 1 if all(x == 0 for x in xs) else -1


def construct(p):
    """Return (primitive root, signs in exponent order, mechanism)."""
    n = p - 1
    factors = factor(n)
    e = factors.pop(2)
    odd = sorted(factors.items())
    moduli = [q ** a for q, a in odd]
    if good_regime(p):
        raise ValueError("No rigid traceless construction in a good regime")
    g = primitive_root(p)
    signs = []
    for i in range(n):
        t = i % (2 ** e)
        x = [i % d for d in moduli]
        if e >= 2:
            half = 2 ** (e - 1)
            base = (point_sign(x[1:]) if t % half == 0
                    else point_sign([x[0]] + x[2:]))
            f = base if t < half else -base
            mechanism = "independent antipodal two-adic slices"
        elif len(odd) >= 3:
            base = (point_sign([x[1]]) if x[2] == 0
                    else point_sign([x[0]]))
            base *= point_sign(x[3:])
            f = base if t == 0 else -base
            mechanism = "third full primary-coordinate selector"
        else:
            repeated = next(j for j, (_, a) in enumerate(odd) if a >= 2)
            other = 1 - repeated
            r, b = odd[repeated]
            step = r ** (b - 1)
            y = x[repeated]
            base = (point_sign([y // step]) if y % step == 0
                    else point_sign([x[other]]))
            f = base if t == 0 else -base
            mechanism = "low-digit selector with full top-digit branch"
        signs.append(f)
    return g, signs, mechanism


def exact_certificate(p, g, signs, include_set=False):
    n = p - 1
    selected = sorted(pow(g, i, p) for i, s in enumerate(signs) if s == 1)
    selected_set = set(selected)
    assert len(selected) == n // 2
    assert selected_set.isdisjoint({p - x for x in selected})
    unit_indices = [k for k in range(1, n) if gcd(k, n) == 1]
    # Direct power sums over the actual field elements, not the structural
    # decomposition used to produce the construction.
    sums = [[k, sum(pow(s, k, p) for s in selected) % p] for k in unit_indices]
    # Independently inspect every multiplier, not only prime-order subgroups.
    stabilizer = [a for a in range(1, p)
                  if {(a * x) % p for x in selected} == selected_set]
    # A separate exponent-shift inspection gives the same answer.
    periods = [d for d in range(n)
               if all(signs[(i + d) % n] == signs[i] for i in range(n))]
    assert sorted(pow(g, d, p) for d in periods) == stabilizer
    cert = {"p": p, "factorization_p_minus_1": factor(n),
            "primitive_root": g, "size": len(selected), "tournament": True,
            "unit_power_sums": sums, "stabilizer": stabilizer,
            "exponent_periods": periods}
    if include_set:
        cert["set"] = selected
    return cert


def exhaustive_census(p):
    """Gray-code census in the fixed gauge {±1},...,{±(p-1)/2}."""
    n, h = p - 1, (p - 1) // 2
    indices = [k for k in range(1, n) if gcd(k, n) == 1]
    cols = [[pow(a, k, p) for k in indices] for a in range(1, h + 1)]
    # Mask bit 0 means choose +a; bit 1 means choose -a.
    sums = [sum(col[j] for col in cols) % p for j in range(len(indices))]
    traceless = rigid = 0
    distribution = Counter()
    mask = 0
    for ordinal in range(1 << h):
        if ordinal:
            bit = (ordinal & -ordinal).bit_length() - 1
            mask ^= 1 << bit
            coeff = -2 if mask & (1 << bit) else 2
            sums = [(s + coeff * col) % p
                    for s, col in zip(sums, cols[bit])]
        if any(sums):
            continue
        traceless += 1
        selected = {p - a if mask & (1 << (a - 1)) else a
                    for a in range(1, h + 1)}
        stab = [a for a in range(1, p)
                if {(a * x) % p for x in selected} == selected]
        distribution[len(stab)] += 1
        rigid += len(stab) == 1
    return {"p": p, "universe": 1 << h, "traceless": traceless,
            "rigid_traceless": rigid,
            "traceless_stabilizer_order_distribution": dict(distribution)}


def meet_in_middle_61():
    """Complete p=61 census, with 2^15 vectors on each side."""
    p, n, h = 61, 60, 30
    indices = [k for k in range(1, n) if gcd(k, n) == 1]
    cols = [[pow(a, k, p) for k in indices] for a in range(1, h + 1)]

    def vectors(part):
        width = len(part)
        sums = [sum(col[j] for col in part) % p for j in range(len(indices))]
        mask = 0
        for ordinal in range(1 << width):
            if ordinal:
                bit = (ordinal & -ordinal).bit_length() - 1
                mask ^= 1 << bit
                coeff = -2 if mask & (1 << bit) else 2
                sums = [(s + coeff * col) % p
                        for s, col in zip(sums, part[bit])]
            yield tuple(sums), mask

    table = defaultdict(list)
    for v, mask in vectors(cols[:15]):
        table[v].append(mask)
    count = rigid = 0
    distribution = Counter()
    for v, high in vectors(cols[15:]):
        for low in table.get(tuple((-x) % p for x in v), []):
            mask = low | high << 15
            selected = {p - a if mask & (1 << (a - 1)) else a
                        for a in range(1, h + 1)}
            # Fresh direct sums on each matched connection set.
            assert all(sum(pow(a, k, p) for a in selected) % p == 0
                       for k in indices)
            stab = [a for a in range(1, p)
                    if {(a * x) % p for x in selected} == selected]
            count += 1
            rigid += len(stab) == 1
            distribution[len(stab)] += 1
    return {"p": p, "universe": 1 << h, "method": "exact 15+15 meet-in-middle",
            "traceless": count, "rigid_traceless": rigid,
            "traceless_stabilizer_order_distribution": dict(distribution)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("verification.json"))
    args = parser.parse_args()
    small = [exhaustive_census(p) for p in primes(31)]
    assert all(row["rigid_traceless"] == 0 for row in small)
    construction_rows = []
    examples = {}
    selected_examples = {61, 127, 151, 181, 211, 241, 271, 331, 421, 541}
    for p in primes(args.bound):
        if good_regime(p):
            continue
        g, signs, mechanism = construct(p)
        cert = exact_certificate(p, g, signs, p in selected_examples)
        assert all(value == 0 for _, value in cert["unit_power_sums"])
        assert cert["stabilizer"] == [1]
        construction_rows.append({"p": p, "p_minus_1": factor(p - 1),
                                  "primitive_root": g, "mechanism": mechanism,
                                  "unit_indices_checked": len(cert["unit_power_sums"]),
                                  "stabilizer": cert["stabilizer"]})
        if p in selected_examples:
            examples[str(p)] = cert | {"mechanism": mechanism}
    mitm = meet_in_middle_61()
    assert mitm["traceless"] == 1444 and mitm["rigid_traceless"] == 360
    output = {
        "complete_small_prime_census": small,
        "complete_p61_meet_in_middle_census": mitm,
        "construction_universe": f"Every bad-regime prime p <= {args.bound}",
        "construction_count": len(construction_rows),
        "construction_checks": construction_rows,
        "named_examples": examples,
        "first_bad_prime": next(p for p in primes(args.bound) if not good_regime(p)),
        "first_bad_prime_3_mod_4": next(p for p in primes(args.bound)
                                     if p % 4 == 3 and not good_regime(p)),
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: v for k, v in output.items()
                      if k not in {"named_examples", "construction_checks"}}, indent=2))


if __name__ == "__main__":
    main()
