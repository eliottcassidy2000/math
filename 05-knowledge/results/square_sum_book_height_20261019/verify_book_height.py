#!/usr/bin/env python3
"""Exact finite checks for BOOK_HEIGHT_ASYMPTOTIC.md; standard library only.

No `assert` statements: all verification stays active under Python -O.
The asymptotic theorem is proved in the accompanying document, not inferred
from these finite tests. JSON output is deterministic in normal and -O modes.
"""
from collections import Counter
from itertools import combinations
from math import isqrt, prod
import json
import random

COUNTS = Counter()


def check(condition, category, detail=None):
    COUNTS[category] += 1
    if not condition:
        raise RuntimeError(f"{category}: {detail!r}")


def square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def primes_at_least_five(k):
    out = []
    p = 5
    while len(out) < k:
        if prime(p):
            out.append(p)
        p += 2
    return out


def all_divisors(n):
    small = [a for a in range(1, isqrt(n) + 1) if n % a == 0]
    return small + [n // a for a in reversed(small) if a * a != n]


def squarefree_divisors(ps):
    ds = [1]
    for p in ps:
        ds += [d * p for d in ds]
    return sorted(ds)


def tau(n):
    return len(all_divisors(n))


def capacity(d):
    if d % 4 == 2:
        return 0
    return tau(d if d % 2 else d // 4) // 2


def factor_neighbors(A, B, N):
    """Exact bijection, independently compared with direct adjacency sets."""
    if not 0 < A < B:
        raise ValueError("positive increasing anchors required")
    d = B - A
    answer = []
    for u in range(1, isqrt(d) + 1):
        if d % u:
            continue
        v = d // u
        if u >= v or (u - v) % 2:
            continue
        r, s = (v - u) // 2, (v + u) // 2
        X = r * r - A
        if 1 <= X <= N and X not in (A, B):
            answer.append((X, r, s, u, v))
    return sorted(answer)


def test_codegrees():
    receipts = []
    for N in (8, 16, 32, 64, 128, 256, 512):
        Amax = min(N, 96)
        direct = {}
        for A in range(1, Amax + 1):
            direct[A] = {X for X in range(1, N + 1)
                         if X != A and square(X + A)}
        maximum = 0
        maximum_over_adjacent = 0
        cases = 0
        for A, B in combinations(range(1, Amax + 1), 2):
            records = factor_neighbors(A, B, N)
            actual = direct[A] & direct[B]
            encoded = {rec[0] for rec in records}
            check(encoded == actual, "factorization_vs_direct", (N, A, B))
            check(len(records) == len(encoded), "factorization_injectivity")
            check(len(actual) <= capacity(B - A), "parity_capacity")
            if (B - A) % 4 == 2:
                check(not actual, "difference_two_mod_four")
            for X, r, s, u, v in records:
                check(X + A == r * r and X + B == s * s,
                      "factorization_root_identities")
                check((s - r, s + r) == (u, v), "factorization_inverse")
            maximum = max(maximum, len(actual))
            if square(A + B):
                maximum_over_adjacent = max(maximum_over_adjacent, len(actual))
            cases += 1
        receipts.append({"N": N, "anchor_maximum": Amax,
                         "anchor_pairs": cases, "max_codegree": maximum,
                         "max_codegree_adjacent_anchors": maximum_over_adjacent})
    # Test the parity capacity itself before positivity/range restrictions.
    for d in range(1, 5001):
        count = sum(1 for u in range(1, isqrt(d) + 1)
                    if d % u == 0 and u < d // u and (u - d // u) % 2 == 0)
        check(count == capacity(d), "exact_unrestricted_parity_capacity", d)
    return receipts


def exact_h(k):
    # Binary search avoids any floating-point fractional power.
    lo, hi = 1, k
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** 4 >= k ** 3:
            hi = mid
        else:
            lo = mid + 1
    return lo


def parameters(k):
    ps = primes_at_least_five(k)
    Q = prod(ps)
    P = ps[-1]
    h = exact_h(k)
    S = isqrt(9 * Q) + 1
    while S % 2 != 1 or S % 3 == 0:
        S += 1
    U = (S * S - 9 * Q) // 2
    V = U + 9 * Q
    cutoff = P ** (2 * h)
    ceiling = 10 * Q * cutoff
    return {"k": k, "primes": ps, "Q": Q, "P": P, "h": h,
            "S": S, "U": U, "V": V, "cutoff": cutoff,
            "ceiling": ceiling}


def retained(a, par):
    return 4 * a * a < par["Q"] and a * a * par["cutoff"] >= par["Q"]


def page(a, par):
    Q, U = par["Q"], par["U"]
    c = Q // a
    r, s = 3 * (c - a) // 2, 3 * (c + a) // 2
    return r * r - U, r, s


def check_parameters(par):
    k, Q, P, h = (par[name] for name in ("k", "Q", "P", "h"))
    S, U, V = (par[name] for name in ("S", "U", "V"))
    check(Q >= 5005 and Q % 6 in (1, 5) and Q % 5 == 0,
          "parameter_moduli")
    check(h ** 4 >= k ** 3 and (h - 1) ** 4 < k ** 3, "exact_h")
    check(h * h >= 2 * k, "finite_chebyshev_constant")
    check(S * S > 9 * Q and S % 2 == 1 and S % 3 != 0,
          "anchor_root_choice")
    for T in range(max(1, S - 6), S):
        check(T % 2 != 1 or T % 3 == 0 or T * T <= 9 * Q,
              "anchor_root_minimality")
    check(S - 6 <= 0 or (S - 6) ** 2 <= 9 * Q, "anchor_root_gap_six")
    check(0 < U < V <= 10 * Q, "anchor_positivity_height")
    check(U + V == S * S and V - U == 9 * Q, "anchor_identities")
    check(U % 3 == V % 3 == 2, "anchor_residues")
    check(par["ceiling"] == 10 * Q * P ** (2 * h), "height_ceiling")


def check_page(a, par):
    Q, U, V, ceiling = (par[name] for name in ("Q", "U", "V", "ceiling"))
    X, r, s = page(a, par)
    check(Q % a == 0 and retained(a, par), "page_domain")
    check(0 < X < ceiling, "page_positive_height")
    check(X + U == r * r and X + V == s * s, "page_root_identities")
    check(r > 0 and s > r and r % 3 == s % 3 == 0, "page_root_domain")
    check(X % 3 == 1 and X not in (U, V), "page_residue_separation")
    check(16 * r * r > 81 * Q, "page_positivity_gap")
    check(4 * X < 9 * Q * par["cutoff"], "page_refined_height")
    return X


def test_all_divisors():
    receipts = []
    for k in range(4, 19):
        par = parameters(k)
        check_parameters(par)
        Q = par["Q"]
        ds = squarefree_divisors(par["primes"])
        check(len(ds) == len(set(ds)) == 2 ** k, "divisor_enumeration")
        lower = [a for a in ds if 4 * a * a < Q]
        upper = [a for a in ds if a * a > 4 * Q]
        central = [a for a in ds if 4 * a * a >= Q and a * a <= 4 * Q]
        check(len(lower) == len(upper), "reciprocal_outside_counts")
        check(len(central) <= 2 ** (k - 1), "central_count")
        check(len(lower) >= 2 ** (k - 2), "lower_quarter_count")
        for a in ds:
            check(Q % a == 0, "enumerated_divisor")
            check((4 * a * a < Q) == ((Q // a) ** 2 > 4 * Q),
                  "reciprocal_outside_bijection")
            if a % 5 != 0:
                b = 5 * a
                a_central = 4 * a * a >= Q and a * a <= 4 * Q
                b_central = 4 * b * b >= Q and b * b <= 4 * Q
                check(not (a_central and b_central), "toggle_five_central_exclusion")
        discarded = [a for a in ds if a * a * par["cutoff"] < Q]
        # Exact rational form of the asserted one-sided Chebyshev estimate.
        check(4 * par["h"] ** 2 * len(discarded) <= k * 2 ** k,
              "finite_chebyshev_tail_bound")
        aa = [a for a in lower if retained(a, par)]
        check(len(aa) >= 2 ** (k - 3), "retained_eighth_count")
        previous_X = None
        for a in aa:
            X = check_page(a, par)
            if previous_X is not None:
                check(X < previous_X, "page_strict_monotonicity")
            previous_X = X
        receipts.append({"k": k, "Q": str(Q), "P": par["P"],
                         "h": par["h"], "divisors": len(ds),
                         "central": len(central), "lower": len(lower),
                         "tail_discarded": len(discarded),
                         "retained": len(aa), "guaranteed": 2 ** (k - 3),
                         "ceiling_digits": len(str(par["ceiling"]))})
    return receipts


def test_large_parameter_samples():
    # Full divisor enumeration is exponential. These sampled checks exercise
    # the genuinely active lower cutoff at much larger k, without claiming
    # to replace the mathematical concentration proof.
    rng = random.Random(20261019)
    receipts = []
    for k in (32, 64, 128, 256):
        par = parameters(k)
        check_parameters(par)
        ps, Q = par["primes"], par["Q"]
        ds = {1, Q, prod(ps[::2]), prod(ps[1::2])}
        for _ in range(600):
            a = prod(p for p in ps if rng.randrange(2))
            ds.add(a)
            ds.add(Q // a)
        leaves = []
        discarded = 0
        for a in sorted(ds):
            if a * a * par["cutoff"] < Q:
                discarded += 1
            if retained(a, par):
                leaves.append(check_page(a, par))
        check(len(leaves) > 0, "large_sample_has_pages")
        check(len(leaves) == len(set(leaves)), "large_sample_distinct_pages")
        if k >= 64:
            check(discarded > 0, "large_sample_exercises_tail")
        for X, Y in combinations(leaves[:100], 2):
            check(not square(X + Y), "large_sample_page_nonedge")
        receipts.append({"k": k, "samples": len(ds),
                         "retained_samples": len(leaves),
                         "tail_discarded_samples": discarded,
                         "ceiling_digits": len(str(par["ceiling"]))})
    return receipts


def build_book(m):
    k = max(4, (m - 1).bit_length() + 3)
    par = parameters(k)
    ds = squarefree_divisors(par["primes"])
    aa = [a for a in ds if retained(a, par)][:m]
    check(len(aa) == m, "all_m_supply")
    leaves = [page(a, par)[0] for a in aa]
    return par, aa, [par["U"], par["V"]] + leaves


def test_induced_books():
    sizes = list(range(1, 33)) + [48, 64, 96, 128, 192, 256, 512, 1024]
    receipts = []
    certificate = None
    for m in sizes:
        par, aa, labels = build_book(m)
        check(len(set(labels)) == m + 2, "book_distinct_labels")
        check(min(labels) >= 1 and max(labels) <= par["ceiling"],
              "book_label_domain")
        edge_count = 0
        for i, j in combinations(range(m + 2), 2):
            expected = i < 2
            observed = square(labels[i] + labels[j])
            check(observed == expected, "full_induced_adjacency", (m, i, j))
            edge_count += observed
        check(edge_count == 2 * m + 1, "book_edge_count")
        receipts.append({"m": m, "k": par["k"],
                         "max_label": str(max(labels)),
                         "ceiling": str(par["ceiling"]),
                         "edges": edge_count})
        if m == 8:
            certificate = {"m": m, "k": par["k"], "Q": str(par["Q"]),
                           "P": par["P"], "h": par["h"],
                           "S": str(par["S"]),
                           "labels": [str(x) for x in labels],
                           "divisors": [str(a) for a in aa],
                           "roots": [[str(page(a, par)[1]), str(page(a, par)[2])]
                                     for a in aa]}
    return receipts, certificate


def main():
    codegrees = test_codegrees()
    divisor_cases = test_all_divisors()
    large_samples = test_large_parameter_samples()
    books, certificate = test_induced_books()
    output = {"status": "PASS", "arithmetic": "exact integer; standard library",
              "checks_active_under_python_O": True,
              "total_checks": sum(COUNTS.values()),
              "checks_by_category": dict(sorted(COUNTS.items())),
              "codegree_cases": codegrees,
              "complete_divisor_enumerations": divisor_cases,
              "large_parameter_samples": large_samples,
              "fully_checked_books": books,
              "example_certificate": certificate,
              "scope": "Finite tests only; unbounded theorem is proved in the companion text."}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
