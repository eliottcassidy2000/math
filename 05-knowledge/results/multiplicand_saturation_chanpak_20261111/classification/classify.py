"""Exact finite classification plus deliberately small synchronous controls.

All finite-time sets use integer arithmetic.  PROOF.md supplies the
all-seed theorem; these controls do not infer asymptotics from a census.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import hashlib
import heapq
import itertools
import json
from math import gcd, prod
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]


def need(condition, label):
    if not condition:
        raise AssertionError(label)


def factor(n):
    if not isinstance(n, int) or n < 1:
        raise ValueError("Seeds must be positive integers")
    out = {}
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            out[divisor] = out.get(divisor, 0) + 1
            n //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    r = 0
    for col in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        divisor = a[r][col]
        a[r] = [x / divisor for x in a[r]]
        for i in range(r + 1, len(a)):
            multiplier = a[i][col]
            if multiplier:
                a[i] = [x - multiplier * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def classify(seeds):
    seeds = sorted(set(seeds))
    if not seeds:
        raise ValueError("An empty seed set has no maximum")
    factors = [factor(n) for n in seeds]
    primes = sorted(set().union(*(set(f) for f in factors)))
    vectors = [[f.get(p, 0) for p in primes] for f in factors]
    dimension = rank(vectors)
    result = dict(seeds=seeds, primes=primes, valuation_vectors=vectors,
                  rational_rank=dimension,
                  eventually_saturated=(dimension <= 1))
    if dimension == 1:
        first = next(v for v in vectors if any(v))
        divisor = reduce(gcd, first)
        primitive = [x // divisor for x in first]
        coordinate = next(i for i, x in enumerate(primitive) if x)
        exponents = []
        for v in vectors:
            quotient, remainder = divmod(v[coordinate], primitive[coordinate])
            need(remainder == 0, "primitive ray scalar is integral")
            need(v == [quotient * x for x in primitive], "common ray")
            exponents.append(quotient)
        result.update(common_integer_base=prod(p ** v for p, v in zip(primes, primitive)),
                      seed_exponents=exponents)
    return result


def next_state(live):
    return live | {x * y for x, y in itertools.combinations(sorted(live), 2)}


def bounded_monoid(seeds, limit):
    """Untimed multiplicative monoid, including its identity."""
    nonunits = sorted(s for s in set(seeds) if s > 1)
    seen, todo = {1}, [1]
    while todo:
        x = heapq.heappop(todo)
        for seed in nonunits:
            y = x * seed
            if y <= limit and y not in seen:
                seen.add(y)
                heapq.heappush(todo, y)
    return seen


def bounded_closure(seeds, limit):
    """Independent mixed-seed-ideal route; no generation-time recurrence."""
    nonunits = sorted(s for s in set(seeds) if s > 1)
    out = {s for s in seeds if s <= limit}
    monoid = bounded_monoid(nonunits, limit)
    for a, b in itertools.combinations(nonunits, 2):
        out.update(a * b * x for x in monoid if a * b * x <= limit)
    return out


def valuation(n, prime):
    out = 0
    while n % prime == 0:
        out += 1
        n //= prime
    return out


def pure_power_exponent(n, base):
    out = 0
    while n > 1 and n % base == 0:
        out += 1
        n //= base
    return out if n == 1 and out else None


def run_case(seeds, rounds, kind):
    live = set(seeds)
    rows = []
    for t in range(rounds + 1):
        maximum = max(live)
        full = bounded_closure(seeds, maximum)
        need(live <= full, "every birth belongs to independent closure")
        holes = sorted(full - live)
        row = dict(t=t, maximum=maximum, cardinality=len(live),
                   saturated_to_max=(not holes),
                   first_hole_below_max=(holes[0] if holes else None),
                   missing_count_below_max=len(holes), live=sorted(live))
        if kind == "aliased_singleton_layer":
            layer = {x for x in live if valuation(x, 3) == 1}
            expected = {b * 2 ** k for b in (3, 12) for k in range(t + 1)}
            need(layer == expected, "exact aliased minimum positive layer")
            target = 12 * 2 ** (t + 1)
            need(target not in live, "next layer maximum remains missing")
            need(target in next_state(live), "next layer maximum arrives next round")
            need(max(layer) == 12 * 2 ** t, "exact layer clock")
            row.update(first_layer=sorted(layer), next_layer_target=target,
                       next_layer_target_below_max=(target < maximum))
            if t >= 2:
                need(target < maximum, "persistent singleton-layer failure")
        if kind == "two_populated_rays":
            expected_max_exponent = 2 ** t + 1
            ray_info = {}
            for base in (2, 3):
                exponents = sorted(e for x in live
                                   if (e := pure_power_exponent(x, base)) is not None)
                need(exponents == list(range(1, expected_max_exponent + 1)),
                     "exact saturated ray restriction")
                ray_info[str(base)] = dict(exponents=exponents, maximum=base ** expected_max_exponent)
            target = 2 ** (expected_max_exponent + 1)
            need(target not in live, "slow-ray target is missing")
            need(target < 3 ** expected_max_exponent <= maximum,
                 "slow-ray hole beneath fast-ray maximum")
            row.update(ray_states=ray_info, next_slow_ray_target=target)
        rows.append(row)
        if t < rounds:
            live = next_state(live)
    if kind == "rank_one_transient":
        need(rows[3]["first_hole_below_max"] == 2 ** 14,
             "rank one can have a temporary nonsaturated frontier")
        need(rows[-1]["saturated_to_max"], "rank-one positive final control")
    return dict(classifier=classify(seeds), control=kind, rounds=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, nargs="+")
    args = parser.parse_args()
    if args.seeds is not None:
        print(json.dumps(classify(args.seeds), indent=2))
        return

    cases = [
        run_case([2, 3, 12], 4, "aliased_singleton_layer"),
        run_case([2, 3, 4, 9], 4, "two_populated_rays"),
        run_case([2, 4, 6, 9], 4, "mixed_ray_population"),
        run_case([4, 8], 6, "rank_one_transient"),
        run_case([6, 36, 216], 4, "composite_common_base"),
        run_case([1, 2, 4], 4, "unit_is_inert"),
        run_case([1], 2, "unit_only"),
        run_case([16], 2, "singleton_nonunit"),
    ]
    expected_ranks = [2, 2, 2, 1, 1, 1, 0, 1]
    need([case["classifier"]["rational_rank"] for case in cases] == expected_ranks,
         "finite exact classifier controls")
    provenance = [
        HERE / "PROOF.md", HERE / "classify.py",
        REPO / "05-knowledge/results/summand_multiplicand_growth_20261111/summand/TOP_SUFFIX.md",
        REPO / "05-knowledge/results/summand_multiplicand_growth_20261111/summand/CLOSURE.md",
    ]
    receipt = dict(status="PASS", scope="eight fixed finite controls; no asymptotic inference",
                   universe=dict(cases=len(cases), states=sum(len(c["rounds"]) for c in cases),
                                 maximum_round=6, filters="none; all strict births retained"),
                   hashes={str(path.relative_to(REPO)): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in provenance}, cases=cases)
    destination = HERE / "receipt.json"
    destination.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: {receipt['universe']['states']} complete finite states; eight classifier controls")
    print(f"Receipt: {destination}")


if __name__ == "__main__":
    main()
