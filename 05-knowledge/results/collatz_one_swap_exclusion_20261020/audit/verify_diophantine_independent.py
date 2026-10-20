#!/usr/bin/env python3
"""Independent closure audit using a different logarithm series.

Reads, but does not modify, the closure report. Does not import its verifier.
"""

from fractions import Fraction as F
from math import factorial, gcd
from pathlib import Path
import argparse
import json


def require(value, label):
    if not value:
        raise RuntimeError(label)


def log_reciprocal_one_minus(t, terms=240):
    """Exact enclosure for -log(1-t), using the ordinary power series."""
    lower = sum((t**k / k for k in range(1, terms + 1)), F(0))
    upper = lower + t ** (terms + 1) / ((terms + 1) * (1 - t))
    return lower, upper


root = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser()
parser.add_argument(
    "--closure-report",
    type=Path,
    default=root / "session13-collatz-diophantine/closure.json",
    help="Path to the production closure.json, useful after publication or relocation.",
)
report_path = parser.parse_args().closure_report
reported = json.loads(report_path.read_text())
require(reported["global_p_cutoff_exclusive"] == 10**12, "current cutoff")
require(reported["legendre_threshold"] == 100, "current Legendre threshold")

# log 2=-log(1-1/2); log(3/2)=-log(1-1/3); log 3 is their sum.
l2, u2 = log_reciprocal_one_minus(F(1, 2))
l32, u32 = log_reciprocal_one_minus(F(1, 3))
l3, u3 = l2 + l32, u2 + u32
require(u2 < F(7, 10) and u3 < F(11, 10), "Matveev heights")
require(512 < 23**2, "2^(9/2)<23")
constant = F(7, 5) * 30**5 * 23 * F(7, 10) * F(11, 10)
require(constant == 602494200 < 10**9, "Matveev constant")

# log(5/4)<1/4 gives log(2*10^12)<0.7+12*(3*0.7+0.25)=28.9<30.
require(F(7, 10) + 12 * (3 * F(7, 10) + F(1, 4)) < 30, "cutoff log")
require(9 * 10**9 * 31 < 10**12, "Matveev cutoff boundary")
exp_lower = sum((F(100, 9) ** j / factorial(j) for j in range(9)), F(0))
require(F(200) / exp_lower < l2, "Legendre boundary")

alpha_lo, alpha_hi = l3 / u2, u3 / l2
lo, hi = alpha_lo, alpha_hi
hm2, hm1, km2, km1 = 0, 1, 1, 0
digits, pairs, signs = [], [], []
for _ in range(100):
    digit = lo.numerator // lo.denominator
    require(digit == hi.numerator // hi.denominator, "CF digit unique")
    digits.append(digit)
    h, k = digit * hm1 + hm2, digit * km1 + km2
    require(gcd(h, k) == 1, "convergent reduced")
    hm2, hm1, km2, km1 = hm1, h, km1, k
    if 100 <= k < 10**12:
        low, high = h * l2 - k * u3, h * u2 - k * l3
        if high < 0:
            signs.append("negative")
        else:
            require(low > F(1, 2 ** min(k // 9, 256)), "upper convergent exclusion")
            signs.append("positive excluded")
        pairs.append((h, k))
    if k >= 10**12:
        break
    lo, hi = 1 / (hi - digit), 1 / (lo - digit)
else:
    raise RuntimeError("CF not exhausted")
require(digits == reported["continued_fraction_prefix"], "independent CF prefix")
require(
    pairs == [(row["q"], row["p"]) for row in reported["excluded_large_convergents"]],
    "all reported convergents",
)
require(reported["terminal_convergent"] == {"q": h, "p": k}, "terminal convergent")

# Re-enumerate the residual words, with q determined by exact integer powers
# instead of logarithm intervals, and carry accumulated right-to-left.
small_rows = []
swaps = shapes = 0
for p in range(2, 100):
    q = (3**p).bit_length()
    if gcd(p, q) != 1:
        small_rows.append((p, q, 0, None))
        continue
    shapes += 1
    d = 2**q - 3**p
    word = [((i + 1) * p) // q - (i * p) // q for i in range(q)]
    count, max_g = 0, 1
    for i in range(q):
        j = (i + 1) % q
        if word[i] == word[j]:
            continue
        edited = word.copy()
        edited[i], edited[j] = edited[j], edited[i]
        b, threes = 0, 1
        for pos in range(q - 1, -1, -1):
            if edited[pos]:
                b += (1 << pos) * threes
                threes *= 3
        content = gcd(d, b)
        require(content < d, ("small word", p, q, i))
        max_g = max(max_g, content)
        swaps += 1
        count += 1
    small_rows.append((p, q, count, max_g))
require(
    small_rows
    == [
        (row["p"], row["q"], row["swaps"], row.get("maximum_content"))
        for row in reported["small_shapes"]
    ],
    "all small shapes and counts",
)
require(swaps == reported["small_swap_occurrences"] == 3548, "small swap count")
require(shapes == reported["small_positive_coprime_shapes"] == 59, "small shape count")
print(
    json.dumps(
        {
            "status": "PASS",
            "independent_logarithm_series": "-log(1-t), 240 terms",
            "matveev_constant_upper": str(constant),
            "global_p_cutoff_exclusive": 10**12,
            "legendre_threshold": 100,
            "continued_fraction_prefix": digits,
            "checked_large_convergents": len(pairs),
            "positive_large_convergents": signs.count("positive excluded"),
            "negative_large_convergents": signs.count("negative"),
            "terminal_convergent": {"q": h, "p": k},
            "small_shapes": shapes,
            "small_swap_occurrences": swaps,
        },
        indent=2,
        sort_keys=True,
    )
)
