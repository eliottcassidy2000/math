"""Exact finite Selmer interfaces; no oracle for arbitrary Selmer groups.

Run: python 04-computation/experiments/bsd_selmer_growth_20261013a.py
The arithmetic input for E_n is the cited Monsky matrix theorem. The BSD
consequence additionally uses the owner's supplied low-corank BSD preprint.
All checks use explicit exceptions and remain active under python -O.
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

CHECKS = 0


def check(predicate, message):
    global CHECKS
    CHECKS += 1
    if not predicate:
        raise RuntimeError(message)


def rank(rows, width):
    """Bit-row Gaussian elimination over F_2."""
    if type(width) is not int or width < 0:
        raise ValueError("width must be an exact nonnegative integer")
    basis = {}
    for value in rows:
        if type(value) is not int or value < 0 or value >= 1 << width:
            raise ValueError("row outside declared width")
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def kernel_count(rows, width):
    """Independent exhaustive vectors, not another elimination."""
    return sum(all((row & v).bit_count() % 2 == 0 for row in rows)
               for v in range(1 << width))


def factors(n):
    if type(n) is not int or n < 1 or n % 2 == 0:
        raise ValueError("expected a positive odd squarefree integer")
    remaining = n
    out = []
    p = 3
    while p * p <= remaining:
        if remaining % p == 0:
            remaining //= p
            out.append(p)
            if remaining % p == 0:
                raise ValueError("n is not squarefree")
        p += 2
    if remaining > 1:
        out.append(remaining)
    return out


def legendre_bit(a, p):
    r = pow(a % p, (p - 1) // 2, p)
    if r == 0:
        raise ValueError("Legendre symbol must be nonzero")
    return 0 if r == 1 else 1


def monsky(n):
    ps = factors(n)
    k = len(ps)
    a = [[0] * k for _ in range(k)]
    for i, p in enumerate(ps):
        for j, q in enumerate(ps):
            if i != j:
                a[i][j] = legendre_bit(q, p)
        a[i][i] = sum(a[i]) % 2
    d2 = [legendre_bit(2, p) for p in ps]
    dm2 = [legendre_bit(-2, p) for p in ps]
    rows = []
    for i in range(k):
        rows.append(sum((a[i][j] ^ (d2[i] if i == j else 0)) << j
                        for j in range(k)) | (d2[i] << (k + i)))
    for i in range(k):
        rows.append((d2[i] << i) | sum(
            (a[i][j] ^ (dm2[i] if i == j else 0)) << (k + j)
            for j in range(k)))
    return ps, a, rows


def phi3_rows(a):
    k = len(a)
    return [sum(((sum(a[i][t] * a[t][j] for t in range(k)) % 2)
                 ^ a[i][j] ^ int(i == j)) << j for j in range(k))
            for i in range(k)]


def tournament_a(k, code):
    a = [[0] * k for _ in range(k)]
    for b, (i, j) in enumerate(combinations(range(k), 2)):
        a[i][j] = (code >> b) & 1
        a[j][i] = 1 ^ a[i][j]
    for i in range(k):
        a[i][i] = sum(a[i]) % 2
    return a


def switched(a, mask):
    k = len(a)
    b = [[0 if i == j else a[i][j] ^ ((mask >> i) & 1) ^ ((mask >> j) & 1)
          for j in range(k)] for i in range(k)]
    for i in range(k):
        b[i][i] = sum(b[i]) % 2
    return b


def delete_sink(a):
    k = len(a)
    check(all(a[k - 1][j] == 0 for j in range(k - 1)), "last vertex not sink")
    b = [row[:-1] for row in a[:-1]]
    for i in range(k - 1):
        b[i][i] ^= 1
    return b


def all3_matrix(a):
    k = len(a)
    packed = [sum(value << j for j, value in enumerate(row)) for row in a]
    return [(row ^ (1 << i)) | (1 << (k + i))
            for i, row in enumerate(packed)] + [
                (1 << i) | (row << k) for i, row in enumerate(packed)]


def abstract_u(s, finite_exponents, m):
    if s < 0 or m < 0 or any(e < 1 for e in finite_exponents):
        raise ValueError("invalid cofinitely generated group data")
    return m * s + sum(min(m, e) for e in finite_exponents)


def low_corank_from_exact_logs(m, selmer_log, torsion_log,
                              next_selmer_log=None, next_torsion_log=None):
    """Only for authenticated EXACT orders of full Selmer and rational torsion.

    This function checks the numerical implication, not authentication.
    Isogeny Selmer orders, rank bounds, or two upper bounds are not inputs.
    """
    if any(type(v) is not int for v in (m, selmer_log, torsion_log)):
        raise ValueError("levels and logarithms must be exact integers")
    if (next_selmer_log is None) != (next_torsion_log is None):
        raise ValueError("both next-level orders required together")
    if m < 1 or min(selmer_log, torsion_log) < 0:
        raise ValueError("invalid receipt")
    u = selmer_log - torsion_log
    if u < 0:
        raise ValueError("negative normalized Selmer size")
    if next_selmer_log is None:
        return u < 2 * m
    if any(type(v) is not int for v in (next_selmer_log, next_torsion_log)):
        raise ValueError("next-level logarithms must be exact integers")
    if min(next_selmer_log, next_torsion_log) < 0:
        raise ValueError("both exact next-level orders required")
    v = next_selmer_log - next_torsion_log
    if v < u:
        raise ValueError("normalized sizes must increase")
    return v - u <= 1


def main():
    out = {"status": "FINITE-EXACT; BSD applications CONDITIONAL on supplied Theorem 1.1",
           "universe": "all positive odd squarefree n <= 10000; E_n: y^2=x^3-n^2*x",
           "input": "Monsky full 2-Selmer matrix, Watkins et al. section 4.1",
           "normalization": "u_m=log_q#Sel_(q^m)-log_q#E(Q)[q^m]"}
    # Abstract groups verify the whole growth interface and delayed detection.
    models = 0
    for s in range(5):
        for length in range(5):
            for es in product(range(1, 5), repeat=length):
                us = [abstract_u(s, es, m) for m in range(22)]
                differences = [us[m + 1] - us[m] for m in range(21)]
                check(all(d >= s for d in differences), "increment below corank")
                check(all(a >= b for a, b in zip(differences, differences[1:])),
                      "increments not decreasing")
                for m in range(1, 21):
                    cert = low_corank_from_exact_logs(m, us[m] + 2, 2)
                    cert2 = low_corank_from_exact_logs(m, us[m] + 2, 2,
                                                     us[m + 1] + 2, 2)
                    check(not cert or s <= 1, "false one-level certificate")
                    check(not cert2 or s <= 1, "false two-level certificate")
                check(any(u < 2 * m for m, u in enumerate(us) if m) == (s <= 1),
                      "eventual one-level criterion")
                check(any(d <= 1 for d in differences) == (s <= 1),
                      "eventual two-level criterion")
                models += 1
    out["abstract_group_models"] = models
    k = 10
    us = [abstract_u(1, (k, k), m) for m in range(1, 23)]
    out["delayed_detection_control"] = {
        "s": 1, "finite_exponents": [k, k],
        "first_average_certificate_level": next(m for m, u in enumerate(us, 1) if u < 2 * m),
        "first_increment_certificate_pair": [k, k + 1],
        "first_ten_levels_equal_corank3_without_finite_part": True}
    check(all(abstract_u(1, (k, k), m) == 3 * m for m in range(k + 1)),
          "truncation hostile")

    histogram = Counter()
    certified = []
    examples = {}
    arithmetic_resonances = 0
    for n in range(1, 10001, 2):
        try:
            ps, a, rows = monsky(n)
        except ValueError:
            continue
        width = 2 * len(ps)
        nullity = width - rank(rows, width)
        check(kernel_count(rows, width) == 1 << nullity, "Monsky independent kernel check")
        histogram[nullity] += 1
        examples.setdefault(str(nullity), {"n": n, "primes": ps, "matrix_rows": rows})
        if nullity <= 1:
            certified.append(n)
        if ps and all(p % 8 == 3 for p in ps):
            check(all(a[i][j] ^ a[j][i] == 1
                      for i in range(len(ps)) for j in range(i)), "reciprocity tournament")
            spectral_nullity = len(ps) - rank(phi3_rows(a), len(ps))
            check(nullity == spectral_nullity, "arithmetic Phi3 correspondence")
            check(nullity % 2 == 0, "F4 kernel dimension")
            arithmetic_resonances += 1
    out["monsky_nullity_histogram"] = dict(sorted(histogram.items()))
    out["curves_examined"] = sum(histogram.values())
    out["low_corank_certified_count"] = len(certified)
    out["low_corank_certified_n"] = certified
    out["smallest_n_by_nullity"] = examples
    out["all_primes_3mod8_checked"] = arithmetic_resonances
    # Independent primary positive and hostile controls.
    for n, expected in [(1, 0), (3, 0), (5, 1), (7, 1), (17, 2), (627, 2)]:
        ps, a, rows = monsky(n)
        check(2 * len(ps) - rank(rows, 2 * len(ps)) == expected, "named curve control")
    for n in (0, 2, 9, 45, 75, -1, True, 1.0):
        try:
            monsky(n)
        except ValueError:
            check(True, "invalid arithmetic input rejected")
        else:
            check(False, "invalid arithmetic input accepted")
    for args in [(True, 0, 0), (1.0, 0, 0), (1, 0.5, 0), (1, 0, False),
                 (1, 0, 0, None, 0), (1, 0, 0, 0, None), (1, 0, 0, 0.0, 0)]:
        try:
            low_corank_from_exact_logs(*args)
        except ValueError:
            check(True, "invalid Selmer receipt rejected")
        else:
            check(False, "invalid Selmer receipt accepted")

    atlas = {}
    triangle_examples = {}
    for k in range(1, 7):
        hist = Counter()
        for code in range(1 << (k * (k - 1) // 2)):
            a = tournament_a(k, code)
            r = k - rank(phi3_rows(a), k)
            matrix_nullity = 2 * k - rank(all3_matrix(a), 2 * k)
            check(r == matrix_nullity, "Schur kernel formula")
            check(r % 2 == 0, "irreducible quadratic requires even nullity")
            hist[r] += 1
            if k % 2 == 0:
                # The marked last vertex becomes a sink under this UNIQUE
                # switch with last bit zero; retain mask for exact decoding.
                mask = sum(a[k - 1][j] << j for j in range(k - 1))
                normalized = switched(a, mask)
                smaller = delete_sink(normalized)
                check(switched(normalized, mask) == a, "lossless switch receipt")
                check(r == k - 1 - rank(phi3_rows(smaller), k - 1),
                      "even-order deletion preserves Phi3 kernel")
                if k <= 4:
                    for switch_mask in range(1 << (k - 1)):
                        b = switched(a, switch_mask)
                        check(k - rank(phi3_rows(b), k) == r, "all small even switches")
            if k <= 4:
                check(kernel_count(phi3_rows(a), k) == 1 << r, "small atlas exhaustive kernel")
            if k == 3:
                outdegrees = sorted(sum(a[i][j] for j in range(k) if j != i)
                                    for i in range(k))
                cyclic = outdegrees == [1, 1, 1]
                check(r == (2 if cyclic else 0), "triangle resonance classification")
        atlas[str(k)] = dict(sorted(hist.items()))
        if k % 2 == 0:
            check(dict(hist) == {int(r): (1 << (k - 1)) * count
                                 for r, count in atlas[str(k - 1)].items()},
                  "exact even-order distribution recursion")
    for ps in combinations([p for p in range(3, 100, 8)
                            if all(p % d for d in range(2, int(p ** 0.5) + 1))], 3):
        n = ps[0] * ps[1] * ps[2]
        _, a, rows = monsky(n)
        r = 6 - rank(rows, 6)
        name = "cyclic" if r == 2 else "transitive"
        candidate = {"n": n, "primes": list(ps), "A": a, "pure_2_selmer_dimension": r}
        if name not in triangle_examples or n < triangle_examples[name]["n"]:
            triangle_examples[name] = candidate
    out["labelled_tournament_atlas_k1_to6"] = atlas
    # At odd order switching loses the arithmetic obstruction.
    cyclic_a = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    check(3 - rank(phi3_rows(cyclic_a), 3) == 2, "odd switch original")
    check(3 - rank(phi3_rows(switched(cyclic_a, 1)), 3) == 0, "odd switch hostile")
    out["switching_scope"] = "Phi3 nullity invariant at even order; false already at order3"
    out["least_triangle_examples_with_primes_below100"] = triangle_examples
    out["checks"] = CHECKS
    payload = json.dumps(out, indent=2, sort_keys=True) + "\n"
    target = Path(__file__).resolve().parents[2] / "05-knowledge/results/bsd_selmer_growth_20261013a.json"
    target.write_text(payload, encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in out.items() if k != "low_corank_certified_n"},
                     indent=2, sort_keys=True))
    print("output_sha256=" + hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
