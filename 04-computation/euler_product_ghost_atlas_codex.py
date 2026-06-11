#!/usr/bin/env python3
"""
Euler-product ghost atlas: eta, primes, Mobius, Liouville, and random exponents.

codex-2026-06-11-P3

For integer exponents b_m, study

    P_b(q) = prod_{m>=1} (1-q^m)^{b_m}.

The logarithm has ghost coordinates

    -q d/dq log P_b(q) = sum_{n>=1} G_b(n) q^n,
    G_b(n) = sum_{d|n} d b_d.

This is the ordinary-power-series cousin of Euler products for zeta/L-functions
and the Witt-vector view of partition functions.  The goal is not to prove a new
number-theory theorem here; it is to locate more cancellation gates and produce
proof-facing observables beyond raw coefficients.
"""

from __future__ import annotations

from collections import Counter
from math import comb, log
import random


N = 180
SEED = 20260611


def primes_upto(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            for k in range(p * p, n + 1, p):
                sieve[k] = False
    return [i for i in range(n + 1) if sieve[i]]


PRIMES = set(primes_upto(N))


def mobius_values(n: int) -> list[int]:
    mu = [1] * (n + 1)
    is_prime = [True] * (n + 1)
    for p in range(2, n + 1):
        if is_prime[p]:
            for k in range(p, n + 1, p):
                is_prime[k] = False
                mu[k] *= -1
            p2 = p * p
            for k in range(p2, n + 1, p2):
                mu[k] = 0
    mu[0] = 0
    return mu


def liouville_values(n: int) -> list[int]:
    lam = [1] * (n + 1)
    for p in primes_upto(n):
        pp = p
        while pp <= n:
            for k in range(pp, n + 1, pp):
                lam[k] *= -1
            pp *= p
    lam[0] = 0
    return lam


MU = mobius_values(N)
LAMBDA = liouville_values(N)


def exponent_schedule(name: str) -> list[int]:
    rng = random.Random(SEED)
    b = [0] * (N + 1)
    if name == "eta_all":
        for m in range(1, N + 1):
            b[m] = 1
    elif name == "prime_only":
        for m in PRIMES:
            b[m] = 1
    elif name == "mobius":
        b = MU[:]
    elif name == "liouville":
        b = LAMBDA[:]
    elif name == "squarefree_only":
        for m in range(1, N + 1):
            b[m] = 1 if MU[m] else 0
    elif name == "random_pm":
        for m in range(1, N + 1):
            b[m] = 1 if rng.randrange(2) else -1
    elif name == "random_prime_pm":
        for m in PRIMES:
            b[m] = 1 if rng.randrange(2) else -1
    else:
        raise ValueError(name)
    return b


def multiply_factor(coeffs: list[int], m: int, exp: int) -> list[int]:
    if exp == 0:
        return coeffs
    factor = [0] * (N + 1)
    if exp > 0:
        for j in range(exp + 1):
            if j * m <= N:
                factor[j * m] = (-1) ** j * comb(exp, j)
    else:
        # (1-x)^(-r) = sum_j C(r+j-1,j) x^j
        r = -exp
        j = 0
        while j * m <= N:
            factor[j * m] = comb(r + j - 1, j)
            j += 1
    out = [0] * (N + 1)
    for i, ai in enumerate(coeffs):
        if ai == 0:
            continue
        for j, fj in enumerate(factor):
            if i + j > N:
                break
            if fj:
                out[i + j] += ai * fj
    return out


def product_coeffs(b: list[int]) -> list[int]:
    coeffs = [1] + [0] * N
    for m in range(1, N + 1):
        coeffs = multiply_factor(coeffs, m, b[m])
    return coeffs


def ghost(b: list[int]) -> list[int]:
    g = [0] * (N + 1)
    for d in range(1, N + 1):
        if b[d] == 0:
            continue
        db = d * b[d]
        for k in range(d, N + 1, d):
            g[k] += db
    return g


def slope(coeffs: list[int]) -> float:
    start = N // 2
    xs = list(range(start, N + 1))
    ys = [log(abs(coeffs[i]) or 1) for i in xs]
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def summarize(name: str) -> dict[str, object]:
    b = exponent_schedule(name)
    coeffs = product_coeffs(b)
    g = ghost(b)
    support = sum(1 for c in coeffs if c)
    g_support = sum(1 for x in g[1:] if x)
    b_hist = Counter(b[1:])
    return {
        "name": name,
        "b_hist": dict(sorted(b_hist.items())),
        "support": support,
        "max_abs": max(abs(c) for c in coeffs),
        "slope": slope(coeffs),
        "ghost_support": g_support,
        "ghost_max": max(abs(x) for x in g),
        "coeffs": coeffs,
        "ghost": g,
    }


def tournament(rows: list[dict[str, object]]) -> dict[str, object]:
    # Lower slope wins; tie by lower ghost support, then listed order.
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    keys = [(float(r["slope"]), int(r["ghost_support"]), i) for i, r in enumerate(rows)]
    for i in range(n):
        for j in range(i + 1, n):
            if keys[i] <= keys[j]:
                adj[i][j] = True
            else:
                adj[j][i] = True
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    cycles = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                if (adj[a][b] and adj[b][c] and adj[c][a]) or (
                    adj[a][c] and adj[c][b] and adj[b][a]
                ):
                    cycles += 1
    order = [rows[i]["name"] for _key0, _key1, i in sorted(keys)]
    return {"score_hist": dict(sorted(Counter(scores).items())), "cycles": cycles, "order": order}


def main() -> None:
    names = [
        "eta_all",
        "prime_only",
        "squarefree_only",
        "mobius",
        "liouville",
        "random_prime_pm",
        "random_pm",
    ]
    rows = [summarize(name) for name in names]
    print("Euler-product ghost atlas")
    print(f"N={N}, SEED={SEED}")
    print()
    print("[1] Product coefficient and ghost fingerprints")
    for r in rows:
        print(
            f"  {r['name']:16s} b_hist={r['b_hist']} "
            f"coeff_support={r['support']:3d}/{N+1} slope={r['slope']:.5f} "
            f"max|a|={r['max_abs']} ghost_support={r['ghost_support']:3d}/{N} "
            f"max|G|={r['ghost_max']}"
        )
        nz = [(i, r["coeffs"][i]) for i in range(N + 1) if r["coeffs"][i]][:10]
        gz = [(i, r["ghost"][i]) for i in range(1, N + 1) if r["ghost"][i]][:10]
        print(f"    coeff first nonzero: {nz}")
        print(f"    ghost first nonzero: {gz}")
    print()
    print("[2] Tournament Analysis over exponent schedules")
    ta = tournament(rows)
    print("  observable: lower coefficient slope wins; tie by lower ghost support")
    print("  score histogram:", ta["score_hist"])
    print("  directed 3-cycles:", ta["cycles"])
    print("  order:", " -> ".join(ta["order"]))
    print()
    print("[3] Reading")
    print("  b_m are Witt/product coordinates; G(n)=sum_{d|n} d b_d are ghost")
    print("  coordinates.  Eta_all has dense ghost sigma_1 but sparse coefficients.")
    print("  Prime_only has sparse exponents but dense ghost on multiples of primes.")
    print("  Mobius/Liouville convert multiplicative signs into ordinary product")
    print("  gates; coefficient growth and ghost support separate cleanly.")


if __name__ == "__main__":
    main()
