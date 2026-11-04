#!/usr/bin/env python3
"""Exact controls for valuation interfaces and finite fixed-support compilation.

Standard library. Writes no source/receipt files. Verification remains active
under -O. The JSON printed to stdout is the reproducible receipt.
"""
from itertools import combinations, product
from math import isqrt, lcm, prod
from fractions import Fraction
import hashlib
import json
from pathlib import Path


def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def valuation(n, r):
    require(n > 0, ("valuation-domain", n, r))
    v = 0
    while n % r == 0:
        n //= r
        v += 1
    return v


def order(p, r):
    x = p % r
    d = 1
    while x != 1:
        x = x * p % r
        d += 1
    return d


def data(p, r):
    d = order(p, r)
    return d, 0 if d == 1 else valuation(p ** d - 1, r)


def allowed(d, role):
    return valuation(d, 2) <= role


def sigma_mod(p, n, modulus):
    """1+p+...+p^(n-1), via independent binary concatenation."""
    total, power = 0, 1
    block_sum, block_power = 1, p % modulus
    while n:
        if n & 1:
            total = (total + power * block_sum) % modulus
            power = power * block_power % modulus
        block_sum = block_sum * (1 + block_power) % modulus
        block_power = block_power * block_power % modulus
        n //= 2
    return total


def value(p, r, n):
    d, c = data(p, r)
    if d == 1:
        return valuation(n, r)
    return c + valuation(n, r) if n % d == 0 else 0


def exact_row(p, targets, role, wanted):
    """Return least divisibility witness n, or a checkable rejection reason."""
    L = 1 + role
    positive = {}
    for r, t in zip(targets, wanted):
        d, c = data(p, r)
        if t == 0:
            continue
        if not allowed(d, role):
            return None, ("role", r, d)
        if d == 1:
            s, k = t, r ** t
        elif t < c:
            return None, ("baseline", r, c, t)
        else:
            s, k = t - c, d * r ** (t - c)
        positive[r] = s
        L = lcm(L, k)
    for r, t in zip(targets, wanted):
        d, _ = data(p, r)
        if t:
            if valuation(L, r) != positive[r]:
                return None, ("overflow", r, positive[r], valuation(L, r), L)
        elif allowed(d, role):
            k = r if d == 1 else d
            if L % k == 0:
                return None, ("forced-edge", r, k, L)
    if L == 1:
        # Avoid all finitely named requirements while making a=n-1 positive.
        L = max([3] + list(targets) + [order(p, r) for r in targets]) + 1
        while not prime(L):
            L += 1
    return L, None


def box_bound(S, q):
    h = len(S) - 1
    Cs = {}
    for r in S:
        Cs[r] = sum(c for p in S if p != r
                    for d, c in [data(p, r)]
                    if allowed(d, int(p == q)))
    C = max(Cs.values(), default=0)
    B = max(2, 2 * h, C + 1)
    while 3 ** (B - C) <= (B + 1) ** h:
        B += 1
    return B, C, Cs


def graph_strong(S, q):
    edges = {(p, r) for p in S for r in S if p != r
             and allowed(order(p, r), int(p == q))}
    for p in S:
        reached = {p}
        while True:
            more = {v for u, v in edges if u in reached} - reached
            if not more:
                break
            reached |= more
        if reached != set(S):
            return False
    return True


def finite_support(S, q):
    B, C, Cs = box_bound(S, q)
    domains = [list(range(1, B, 4)) if p == q else list(range(2, B, 2))
               for p in S]
    table = {(p, a): tuple(0 if p == r else value(p, r, a + 1) for r in S)
             for p, domain in zip(S, domains) for a in domain}
    for (p, a), row in table.items():
        direct_sigma = sum(p ** j for j in range(a + 1))
        require(row == tuple(valuation(direct_sigma, r) for r in S),
                ("direct-component-row", S, q, p, a, row))
    count = balanced = perfect = 0
    obstruction_histogram = {str(r): 0 for r in S}
    for exponents in product(*domains):
        count += 1
        incoming = tuple(sum(table[p, a][j] for p, a in zip(S, exponents))
                         for j in range(len(S)))
        failure = next((j for j, a in enumerate(exponents) if incoming[j] != a), None)
        if failure is not None:
            obstruction_histogram[str(S[failure])] += 1
            continue
        balanced += 1
        N = prod(p ** a for p, a in zip(S, exponents))
        sigma = prod((p ** (a + 1) - 1) // (p - 1) for p, a in zip(S, exponents))
        perfect += int(sigma == 2 * N)
    return dict(support=S, euler=q, strong=graph_strong(S, q), bound_exclusive=B,
                baseline_C=C, baselines=Cs, candidates=count,
                valuation_balanced=balanced, perfect=perfect,
                first_failed_budget=obstruction_histogram,
                abundancy_ceiling=str(prod(Fraction(p, p - 1) for p in S)))


def main():
    source_bytes = Path(__file__).read_bytes()
    receipt = {}
    # Exhaustive requested valuation vectors 0..3 for every 3-target choice
    # from odd primes below 20, every other such source, and both roles.
    P = [p for p in range(3, 20, 2) if prime(p)]
    rows = realized = rejected = 0
    reasons = {}
    for p in P:
        for R in combinations([r for r in P if r != p], 3):
            for role in (0, 1):
                for t in product(range(4), repeat=3):
                    rows += 1
                    n, why = exact_row(p, R, role, t)
                    if n is None:
                        rejected += 1
                        reasons[why[0]] = reasons.get(why[0], 0) + 1
                        continue
                    realized += 1
                    require(n >= 2 and (n % 2 == 1 if role == 0 else n % 4 == 2),
                            ("witness-role", p, R, role, t, n))
                    for r, v in zip(R, t):
                        residue = sigma_mod(p, n, r ** (v + 1))
                        require(residue % (r ** v) == 0 and residue != 0,
                                ("witness-valuation", p, R, role, t, n, r, residue))
    receipt["requested_rows"] = dict(total=rows, realized=realized,
                                      rejected=rejected, rejection_reasons=reasons)

    # Independent complete-period enumeration: all exact 0/1 valuation triples
    # at targets 3,5,7, sources 11,13,17,19, both roles. Sum update uses no LTE.
    periodic = 0
    period_rows = []
    R = (3, 5, 7)
    for p in (11, 13, 17, 19):
        Q = lcm(4, *(order(p, r) * r ** 2 for r in R))
        observed = {0: set(), 1: set()}
        sums = [0] * 3
        powers = [1] * 3
        for n in range(1, Q + 1):
            for j, r in enumerate(R):
                sums[j] = (sums[j] + powers[j]) % (r * r)
                powers[j] = powers[j] * p % (r * r)
            if all(sums):
                t = tuple(int(s % r == 0) for s, r in zip(sums, R))
                if n % 2:
                    observed[0].add(t)
                elif n % 4 == 2:
                    observed[1].add(t)
            periodic += 1
        require(sums == [0, 0, 0] and powers == [1, 1, 1], ("period", p, Q))
        for role in (0, 1):
            predicted = {t for t in product(range(2), repeat=3)
                         if exact_row(p, R, role, t)[0] is not None}
            require(predicted == observed[role], ("independent-period", p, role,
                                                 predicted, observed[role]))
        period_rows.append(dict(source=p, period=Q,
                                feasible=[len(observed[i]) for i in (0, 1)]))
    receipt["independent_periods"] = dict(iterations=periodic, cases=period_rows)

    # Sharp interaction: positive support {11,23} is realizable, but (2,1) is not.
    require(data(3, 11) == (5, 2) and data(3, 23) == (11, 1), "3-order-baselines")
    require(exact_row(3, (11, 23), 0, (2, 1))[0] is None, "interaction-not-rejected")
    n, why = exact_row(3, (11, 23), 0, (3, 1))
    require(n == 55 and why is None, ("sharp-interaction", n, why))
    sigma = (3 ** 55 - 1) // 2
    require(valuation(sigma, 11) == 3 and valuation(sigma, 23) == 1,
            "direct-sharp-valuations")
    require(value(3, 11, 5) == 2 and value(3, 23, 5) == 0, "single-edge-11")
    require(value(3, 11, 11) == 0 and value(3, 23, 11) == 1, "single-edge-23")
    receipt["sharp_interaction"] = dict(source=3, targets=[11, 23],
                                         forbidden=[2, 1], realized=[3, 1],
                                         exponent=54, sigma=str(sigma))

    # Denominator exception and a disallowed even order are mandatory hostiles.
    require(exact_row(23, (11,), 0, (0,))[0] is not None, "order-one-empty")
    require(value(23, 11, 3) == 0 and value(23, 11, 11) == 1, "order-one-gate")
    require(exact_row(3, (5,), 1, (1,))[0] is None, "order-four-Euler-gate")
    receipt["hostile_controls"] = ["order-one denominator", "Euler order-four exclusion",
                                    "joint order valuation overflow", "individual antecedents insufficient"]

    # Every role in this fixed five-prime support passes the maximal graph.
    # Exact budgets then exclude its entire mathematically bounded exponent box.
    supports = [finite_support((3, 5, 7, 11, 13), q) for q in (5, 13)]
    require(all(x["strong"] and x["valuation_balanced"] == 0 for x in supports), supports)
    receipt["fixed_supports"] = supports
    require(Path(__file__).read_bytes() == source_bytes, "source changed during verification")
    receipt["source_sha256"] = hashlib.sha256(source_bytes).hexdigest()
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
