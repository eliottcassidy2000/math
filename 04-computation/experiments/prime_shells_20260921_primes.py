"""Exact, optimization-safe harmonic/exceptional-prime audit; standard library only.

Run with --output PATH; the default output is the sibling .json. No network,
randomness, assert statements, conjectural primality filters, or orbit claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path

CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(message)


def sieve(n: int) -> bytearray:
    a = bytearray(b"\1") * (n + 1)
    a[:2] = b"\0\0"
    for p in range(2, isqrt(n) + 1):
        if a[p]:
            a[p * p::p] = b"\0" * ((n - p * p) // p + 1)
    return a


def segmented_count(lo: int, hi: int) -> int:
    """Independent block sieve, inclusive bounds, trial-generated base primes."""
    bases = [p for p in range(2, isqrt(hi) + 1)
             if all(p % d for d in range(2, isqrt(p) + 1))]
    answer = 0
    for left in range(lo, hi + 1, 32768):
        right = min(hi, left + 32767)
        composite = bytearray(right - left + 1)
        for p in bases:
            start = max(p * p, ((left + p - 1) // p) * p)
            for value in range(start, right + 1, p):
                composite[value - left] = 1
        answer += sum(not composite[v - left] for v in range(max(2, left), right + 1))
    return answer


def harmonic(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction())


def rational_mod(x: Fraction, modulus: int) -> int:
    return x.numerator * pow(x.denominator, -1, modulus) % modulus


def eta(p: int) -> int:
    """H_(p-1)/p^2 modulo p, by direct reciprocal sum modulo p^3."""
    h = sum(pow(k, -1, p**3) for k in range(1, p)) % p**3
    check(h % (p * p) == 0, f"Wolstenholme theorem at {p}")
    return h // (p * p)


def short_eta(p: int) -> int:
    """Independent cited McIntosh-1995 Thm 2 formula, p>=11."""
    total = sum(pow(pow(k, -1, p), 3, p)
                for k in range(p // 6 + 1, p // 4 + 1)) % p
    return -total * pow(63, -1, p) % p


def q2(p: int) -> int:
    value = pow(2, p - 1, p * p) - 1
    check(value % p == 0, f"Fermat at {p}")
    return value // p


def wilson(p: int) -> int:
    value = 1
    for k in range(1, p):
        value = value * k % (p * p)
    check((value + 1) % p == 0, f"Wilson at {p}")
    return (value + 1) // p % p


def fib_pair(n: int, modulus: int) -> tuple[int, int]:
    a, b = 0, 1
    for bit in bin(n)[2:]:
        c = a * (2 * b - a) % modulus
        d = (a * a + b * b) % modulus
        a, b = (d, (c + d) % modulus) if bit == "1" else (c, d)
    return a, b


def fibonacci_quotient(p: int) -> int:
    character = 1 if pow(5, (p - 1) // 2, p) == 1 else -1
    value = fib_pair(p - character, p * p)[0]
    check(value % p == 0, f"Fibonacci prime congruence at {p}")
    return value // p


def bernoulli_rational(bound: int) -> list[Fraction]:
    values = [Fraction(1)]
    for n in range(1, bound + 1):
        values.append(-sum((comb(n + 1, k) * values[k]
                            for k in range(n)), Fraction()) / (n + 1))
    return values


def bernoulli_mod(p: int) -> list[int]:
    """Akiyama-Tanigawa triangular recurrence; B1=+1/2 convention only here."""
    row, result = [], []
    for m in range(p - 2):
        row.append(pow(m + 1, -1, p))
        for j in range(m, 0, -1):
            row[j - 1] = j * (row[j - 1] - row[j]) % p
        result.append(row[0])
    return result


def bernoulli_witness(p: int, r: int) -> int:
    """Hart-Harvey-Ong (2016), Eq.(1), direct moment, no FFT needed here."""
    c = next(c for c in range(2, p) if pow(c, r, p) != 1)
    invc = pow(c, -1, p)
    # Double f_c to avoid rational arithmetic inside the sum.
    moment = sum(pow(x, r - 1, p) *
                 (2 * (c * (x * invc % p) // p) - (c - 1))
                 for x in range(1, p)) % p
    return moment * r * pow(2 * (pow(c, r, p) - 1), -1, p) % p


def valuation(n: int, p: int) -> int:
    if n == 0:
        raise ValueError("valuation requires nonzero input")
    answer = 0
    while n % p == 0:
        n //= p
        answer += 1
    return answer


def affine_binomial_controls() -> list[dict]:
    rows = []
    for p in (5, 7, 11, 13, 17, 16843):
        q = comb(2 * p - 1, p - 1)
        r = valuation(q - 1, p)
        check(r == (4 if p == 16843 else 3), "exact binomial valuation")
        check(valuation(q, 2) == p.bit_count() - 1 >= 1, "even multiplier")
        c = (q - 1) // p
        check(c % 2 == 1, "odd affine offset")
        levels = []
        for s in range(1, 6):
            period = p ** max(0, s - r + 1)
            modulus = 2 * p**s
            # Exact closed iterate computed modulo p*modulus before division.
            for n in (1, 3, modulus - 1):
                a = pow(q, period, p * modulus)
                check(((p * n + 1) * (a - 1) // p) % modulus == 0,
                      "claimed affine period returns")
                if period > 1:
                    a = pow(q, period // p, p * modulus)
                    check(((p * n + 1) * (a - 1) // p) % modulus != 0,
                          "minimal affine period")
            if p in (5, 7, 11) and s <= 4:
                # Every odd residue, direct iteration, no closed-form reuse.
                unseen = set(range(1, modulus, 2))
                while unseen:
                    start = min(unseen)
                    n, length = start, 0
                    while n in unseen:
                        unseen.remove(n)
                        n = (q * n + c) % modulus
                        length += 1
                    check(n == start and length == period, "all direct affine cycles")
            levels.append({"s": s, "period": period})
        rows.append({"p": p, "v_p_Q_minus_1": r, "v2_Q": valuation(q, 2), "levels": levels})
    check(comb(9, 4) == 126 and 5 * (126 * 1 + 25) + 1 == 126 * (5 * 1 + 1),
          "p5 affine intertwining")
    check((5 * (126 + 25) + 1) // 4 == 189 != (5 + 1) // 2,
          "odd core not preserved by binomial multiplier")
    return rows


def main() -> dict:
    flags = sieve(2124679)
    primes = [p for p, yes in enumerate(flags) if yes]
    for n in range(10001):
        check(bool(flags[n]) == (n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))),
              f"sieve/trial disagreement {n}")
    ranks = {p: sum(flags[:p + 1]) for p in (1093, 2851, 3511, 12101, 16843, 2124679)}
    between = sum(flags[16844:2124679])
    check(between == segmented_count(16844, 2124678) == 155559, "strict prime interval")
    check(ranks[16843] == 1944 and ranks[2124679] == 157504, "prime ranks")
    check(comb(9, 2) == 36 != 5**3 + 1, "literal binomial hostile")
    check(comb(9, 4) == comb(9, 5) == 5**3 + 1 == 126, "binomial repair")
    for n in range(1, 101):
        check((comb(2 * n - 1, n - 1) == n**3 + 1) == (n == 5), "unique binomial cube-plus-one control")
    check(harmonic(4) == Fraction(25, 12), "H4")
    check(harmonic(6) == Fraction(49, 20), "H6")
    harmonic_totals = {}
    for n in (4, 6):
        total = sum((harmonic(j) for j in range(1, n + 1)), Fraction())
        check(total == (n + 1) * harmonic(n) - n, "sum of harmonic numbers")
        harmonic_totals[n] = str(total)
    bernoulli = bernoulli_rational(297)
    small_rows, irregular = [], {}
    for p in primes:
        if p < 5:
            continue
        if p > 300:
            break
        h, hhalf = harmonic(p - 1), harmonic((p - 1) // 2)
        e = eta(p)
        beta_exact = (comb(2 * p - 1, p - 1) - 1) // p**3
        check((comb(2 * p - 1, p - 1) - 1) % p**3 == 0, "binomial divisibility")
        check(e == rational_mod(h / (p * p), p), "rational/direct harmonic")
        check(beta_exact % p == 2 * e % p, "beta=2eta")
        check(rational_mod(hhalf, p) == -2 * q2(p) % p, "Eisenstein")
        check(e == rational_mod(-bernoulli[p - 3] / 3, p), "Bernoulli eta")
        if p >= 11:
            check(e == short_eta(p), "short cubic sum")
        bmod = bernoulli_mod(p)
        for k in range(2, p - 1, 2):
            check(bmod[k] == rational_mod(bernoulli[k], p), "Bernoulli independent recurrence")
        indices = [k for k in range(2, p - 1, 2) if bmod[k] == 0]
        if indices:
            irregular[p] = indices
        small_rows.append({"p": p, "eta": e, "beta": beta_exact % p,
                           "half_harmonic": rational_mod(hhalf, p), "q2": q2(p)})
    check([r["eta"] for r in small_rows[:7]] == [3, 6, 6, 7, 10, 14, 18], "first user sequence")
    check([r["beta"] for r in small_rows[:7]] == [1, 5, 1, 1, 3, 9, 13], "second user sequence")
    check(irregular[37] == [32] and 24 == next(r["eta"] for r in small_rows if r["p"] == 37),
          "irregular does not imply Wolstenholme")
    check(harmonic(10).numerator == 7381 == 11 * 11 * 61 and isqrt(7381)**2 != 7381,
          "harmonic numerator need not be a square")
    for p in (1093, 2851):
        b = bernoulli_mod(p)
        check(all(b[k] for k in range(2, p - 1, 2)), f"complete regularity at {p}")
    # The larger lists are author-dataset rows; independently certify every positive index.
    cited_irregular = {1093: [], 2851: [], 3511: [1416, 1724], 12101: [7718],
                      16843: [16840], 2124679: [701898, 2124676]}
    for p, indices in cited_irregular.items():
        for r in indices:
            check(bernoulli_witness(p, r) == 0, f"irregular witness {p},{r}")
    selected = []
    for p in cited_irregular:
        e = eta(p)
        check(e == short_eta(p), f"large independent short formula {p}")
        half = sum(pow(k, -1, p) for k in range(1, (p + 1) // 2)) % p
        q = q2(p)
        check(half == -2 * q % p, f"large half harmonic {p}")
        selected.append({"p": p, "prime_rank": ranks[p], "eta": e, "beta": 2 * e % p,
                         "q2": q, "half_harmonic": half, "wilson_quotient": wilson(p),
                         "fibonacci_quotient": fibonacci_quotient(p),
                         "cited_complete_irregular_indices": cited_irregular[p]})
    check({r["p"] for r in selected if r["eta"] == 0} == {16843, 2124679}, "selected W positives")
    check({r["p"] for r in selected if r["q2"] == 0} == {1093, 3511}, "selected Wieferich positives")
    check(all(r["wilson_quotient"] and r["fibonacci_quotient"] for r in selected),
          "selected Wilson/WSS negatives")
    # Bounded searches, not new records or all-prime classifications.
    bounded = [p for p in primes if 5 <= p <= 20000]
    near_eta = {j: [] for j in range(4)}
    for p in bounded:
        e = eta(p) if p < 11 else short_eta(p)
        if e < 4:
            near_eta[e].append(p)
    check(near_eta[0] == [16843], "W census to 20000")
    wp = [p for p in primes if 3 <= p <= 10000 and q2(p) == 0]
    wilp = [p for p in primes if 3 <= p <= 10000 and wilson(p) == 0]
    fwp = [p for p in primes if 3 <= p <= 10000 and p != 5 and fibonacci_quotient(p) == 0]
    check(wp == [1093, 3511], "Wieferich census to 10000")
    check(wilp == [5, 13, 563], "Wilson positive controls/census")
    check(fwp == [], "WSS census to 10000")
    # Fast-doubling independent recurrence controls, including n=0.
    for modulus in (5, 25, 49, 121):
        a, b = 0, 1
        for n in range(101):
            check(fib_pair(n, modulus) == (a, b), "Fibonacci independent recurrence")
            a, b = b, (a + b) % modulus
    affine_rows = affine_binomial_controls()
    return {"status": "FINITE-EXACT", "checks": CHECKS,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "universes": {"prime_interval": [16843, 2124679], "interval_endpoints": "excluded",
                          "harmonic_independent_checks_prime_bound": 300,
                          "near_eta_prime_bound": 20000, "other_prime_search_bound": 10000,
                          "irregular_census_prime_bound": 300,
                          "regularity_complete_additional": [1093, 2851]},
            "strict_prime_count": between, "prime_ranks": ranks,
            "H4": str(harmonic(4)), "H6": str(harmonic(6)),
            "sums_of_harmonic_numbers": harmonic_totals,
            "first_prime_rows": small_rows, "irregular_census": irregular,
            "selected_primes": selected, "near_eta_census": near_eta,
            "affine_binomial_controls": affine_rows,
            "bounded_wieferich": wp, "bounded_wilson": wilp, "bounded_fibonacci_wieferich": fwp,
            "scope": "All positive irregular witnesses checked. Larger complete index lists are cited, not independently exhaustively recomputed. No new-prime record or Collatz implication."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    result = main()
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"PASS: {result['checks']} exact checks; wrote {args.output}")
