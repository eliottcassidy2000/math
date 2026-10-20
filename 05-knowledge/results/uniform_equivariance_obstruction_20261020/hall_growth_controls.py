#!/usr/bin/env python3
"""Independent exact controls for the two-base Hall-growth addendum."""
import argparse
import itertools
import json
import math
from pathlib import Path


def valuation(a, p):
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def factorial_valuation(n, p):
    v = 0
    while n:
        n //= p
        v += n
    return v


def binomial_valuation(n, k, p):
    return factorial_valuation(n, p) - factorial_valuation(k, p) - factorial_valuation(n-k, p)


def digits_fit(n, k, p):
    while n or k:
        if k % p > n % p:
            return False
        n //= p
        k //= p
    return True


def graph_controls(n, k):
    edges = [(a, b) for a in range(n) for b in range(a+1, n) if (a < k) != (b < k)]
    edge_set = set(edges)
    all_pairs = list(itertools.combinations(range(n), 2))
    orders = 0
    odd_auts = 0
    for p in itertools.permutations(range(n)):
        if all(((a, b) in edge_set) == (tuple(sorted((p[a], p[b]))) in edge_set) for a, b in all_pairs):
            orders += 1
            odd_auts += sum(p[a] > p[b] for a, b in edges) % 2
    expected = math.factorial(k)*math.factorial(n-k)*(2 if k > 0 and 2*k == n else 1)
    assert orders == expected
    even = odd_auts == 0
    assert even == (2*k != n or k % 2 == 0 or k == 0)
    degrees = [sum(v in edge for edge in edges) for v in range(n)]
    euler = all(d % 2 == 0 for d in degrees)
    assert euler == (k == 0 or (k % 2 == 0 and (n-k) % 2 == 0))
    full = all(valuation(orders, p) == factorial_valuation(n, p) for p in (3, 5))
    assert full == (math.gcd(math.comb(n, k), 15) == 1)
    return {"n": n, "k": k, "automorphism_order": orders, "sign_odd_automorphisms": odd_auts,
            "sign_even": even, "eulerian": euler, "both_full_ambient_sylows": full}


def bounds(n):
    admissible = [k for k in range(n//2+1) if math.gcd(math.comb(n, k), 15) == 1]
    b = [k for k in admissible if 2*k < n or k % 2 == 0]
    e = [k for k in admissible if k == 0 or (k % 2 == 0 and (n-k) % 2 == 0)]
    return {"n": n, "tournament_lower_bound": len(b), "switching_lower_bound": len(e),
            "tournament_k": b, "switching_k": e}


def family_case(a, b, balanced=False):
    three, five = 3**a, 5**b
    n, limit = three*five-1, min(three, five)
    assert 2*(limit-1) < n
    for k in range(limit):
        for p in (3, 5):
            assert binomial_valuation(n, k, p) == 0
            assert digits_fit(n, k, p)
        if n <= 5000:
            assert math.gcd(math.comb(n, k), 15) == 1
    count_even = sum(k % 2 == 0 for k in range(limit))
    assert count_even == (limit+1)//2
    if balanced:
        assert five <= three < 3*five
        assert 3*limit*limit > n+1
    return {"a": a, "b": b, "n": n, "tournament_lower_bound": limit,
            "switching_lower_bound": count_even, "tested_k_count": limit,
            "exact_binomial_control": n <= 5000, "square_root_case": balanced}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="hall_growth_controls.json")
    args = parser.parse_args()
    report = {"scope": "All permutations of complete bipartite graphs through n=8; all binomials through n=250; full declared prime-power prefixes"}
    report["graph_controls"] = [graph_controls(n, k) for n in range(1, 9) for k in range(n//2+1)]
    report["graph_permutation_checks"] = sum(math.factorial(n)*(n//2+1) for n in range(1, 9))
    digit_checks = 0
    for n in range(251):
        for k in range(n+1):
            coefficient = math.comb(n, k)
            for p in (3, 5):
                assert (coefficient % p != 0) == digits_fit(n, k, p)
                assert valuation(coefficient, p) == binomial_valuation(n, k, p)
                digit_checks += 1
    report["digit_and_exact_binomial_checks"] = digit_checks
    report["general_bounds"] = [bounds(n) for n in range(5, 251)]
    report["extra_bounds"] = [bounds(674)]
    expected = {5:(1,1),6:(1,1),7:(2,1),8:(4,2),14:(6,3),44:(18,9),224:(36,18),674:(27,14)}
    for n, pair in expected.items():
        row = bounds(n)
        assert (row["tournament_lower_bound"], row["switching_lower_bound"]) == pair
    report["prime_power_grid"] = [family_case(a,b) for a in range(1,6) for b in range(1,6)]
    report["diagonal_cases"] = [family_case(r,r) for r in range(1,9)]
    report["square_root_cases"] = []
    for b in range(1,9):
        a = 1
        while 3**a < 5**b:
            a += 1
        row = family_case(a,b,balanced=True)
        report["square_root_cases"].append(row)
        print(json.dumps(row, sort_keys=True))
    report["all_checks_passed"] = True
    Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print("ALL CHECKS PASSED")

if __name__ == "__main__":
    main()
