"""Exact controls for the discrepancy-capacity proof, not Collatz convergence.

Universe: mechanical prefixes of lengths 1..128; all odd starts through
1999 for at most 100 odd steps; signed controls stated explicitly below.
All assertions use require so optimized Python retains the checks.
"""
from fractions import Fraction
import json
from math import comb
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def v2(n):
    require(n > 0, "valuation requires a positive integer")
    return (n & -n).bit_length() - 1


def word_source(word):
    carry, total, power = 0, 0, 1
    for k in word:
        carry = 3 * carry + 2 ** total
        total += k
        power *= 3
    modulus = 2 ** (total + 1)
    source = ((2 ** total - carry) * pow(power, -1, modulus)) % modulus
    return source, modulus, carry, total


def check_orbit(n, a=3, b=1, steps=100):
    initial, q, qsum, seen, nodes = n, Fraction(1), Fraction(0), set(), []
    for j in range(steps + 1):
        require(n * q == initial + Fraction(b, a) * qsum,
                "exact discrepancy/carry identity")
        nodes.append(n)
        if n in seen or j == steps:
            break
        seen.add(n)
        z = a * n + b
        require(z > 0, "signed control left positive domain")
        k = v2(z)
        require(k >= 1, "odd source guard lost")
        qsum += q
        q *= Fraction(2 ** k, a)
        n = z // 2 ** k
    return nodes


def unit_odd_count(m):
    return (m + 5) // 6 + (m + 1) // 6


def contraction_density_controls():
    rows = []
    for length in (4, 8, 16, 64, 128, 256):
        threshold = 7 * length // 4
        bad = sum((Fraction(comb(s-1, length-1), 2 ** s)
                   for s in range(length, threshold)), Fraction(0))
        require(bad <= Fraction(32, length), "Chebyshev source-density bound")
        if length <= 16:
            # Independent convolution of geometric weights, truncated exactly
            # where larger cumulative exponents cannot reenter the bad set.
            weights = {0: Fraction(1)}
            for _ in range(length):
                nxt = {}
                for s, weight in weights.items():
                    for k in range(1, threshold-s):
                        nxt[s+k] = nxt.get(s+k, Fraction(0)) + weight / 2 ** k
                weights = nxt
            require(sum(weights.values(), Fraction(0)) == bad, "geometric convolution")
        if length <= 8:
            hits = 0
            for source in range(1, 2 ** threshold, 2):
                n, total = source, 0
                for _ in range(length):
                    k = v2(3*n+1)
                    total += k
                    n = (3*n+1) // 2 ** k
                hits += total < threshold
            require(Fraction(hits, 2 ** (threshold-1)) == bad,
                    "exhaustive cylinder-period density")
        rows.append({"L": length, "bad_event": f"K_L<{threshold}",
                     "density": str(bad), "Chebyshev_bound": str(Fraction(32, length)),
                     "good_leading_factor_at_most": str(Fraction(81, 128) ** (length//4))})
    return rows


def main():
    # floor(j log_2(3)) needs no floating-point logarithm.
    cumulative = [(3 ** j).bit_length() - 1 for j in range(129)]
    word = [cumulative[j] - cumulative[j-1] for j in range(1, 129)]
    require(set(word) == {1, 2}, "mechanical alphabet")
    qvalues = [Fraction(2 ** cumulative[j], 3 ** j) for j in range(129)]
    require(all(Fraction(1, 2) < q <= 1 for q in qvalues), "mechanical strip")
    for length in range(1, 33):
        weights = [sum(word[i:i+length]) for i in range(129-length)]
        require(max(weights)-min(weights) <= 1, "finite mechanical balance")

    prefix_rows, last_source, last_modulus = [], 0, 1
    for length in range(1, 129):
        source, modulus, carry, total = word_source(word[:length])
        require(source > 0 and source % 2 == 1, "positive odd word source")
        require(source % last_modulus == last_source, "nested cylinders")
        require(source >= last_source, "least representatives are monotone")
        n = source
        for k in word[:length]:
            require(v2(3*n+1) == k, "exact finite-word valuations")
            n = (3*n+1) // 2 ** k
        require(2 ** total * n == 3 ** length * source + carry, "word carry")
        if length <= 12 or length in (16, 32, 64, 128):
            prefix_rows.append({"length": length, "source": source,
                                "modulus": modulus, "K": total,
                                "carry": carry, "endpoint": n})
        last_source, last_modulus = source, modulus

    # Independent iteration checks the additive identity, including repetitions.
    for n in range(1, 2000, 2):
        check_orbit(n)
    for m in range(3000):
        require(unit_odd_count(m) == sum(x % 6 in (1, 5) for x in range(1, m+1)),
                "exact unit-odd counting function")
        require(3 * unit_odd_count(m) <= m + 3, "density bound with endpoint error")

    signed = {}
    for a, b, starts in ((3, -1, (1, 5, 17, 95)), (5, 1, (1, 3, 7)),
                         (5, -1, (1, 3, 9)), (3, -5, (5, 15, 25)),
                         (3, 5, (1, 5, 7))):
        signed[f"a={a},b={b}"] = [check_orbit(n, a, b, 50) for n in starts]
    try:
        check_orbit(3, 3, -5, 5)
    except RuntimeError as error:
        require(str(error) == "signed control left positive domain", "unexpected guard error")
    else:
        raise RuntimeError("3n-5 hostile at3 should leave the positive domain via1")

    # Hostiles: a positive fixed point has UNBOUNDED discrepancy, so it
    # does not contradict the theorem; the negative variant has q -> 0.
    require(Fraction(4, 3) ** 20 > 9, "positive fixed-point strip escape")
    require(sum((Fraction(2, 3) ** j for j in range(80)), Fraction(0)) < 3,
            "negative fixed-point summability")
    result = {
        "status": "PASS", "collatz_convergence_proved": False,
        "scope": "Exact finite controls; infinite claims have proofs in the companion note.",
        "mechanical_prefix_lengths": [1, 128], "balance_lengths_checked": [1, 32],
        "mechanical_word_first_64": word[:64], "prefix_controls": prefix_rows,
        "ordinary_orbit_controls": {"odd_sources": [1, 1999], "max_steps": 100},
        "signed_controls": signed,
        "fixed_factor_descent_density_controls": contraction_density_controls(),
        "proved_capacity_ratio_threshold": 9,
        "stronger_proved_statement": "No positive Collatz orbit has q bounded above and bounded away from zero",
        "general_capacity_threshold": "2*a*a/(b*phi(a)), b>0, gcd(a,b)=1",
        "negative_carry_consequence": "sum(q_j)<=a*n0/abs(b); q_j tends to zero",
    }
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print("PASS: exact discrepancy identities, nested mechanical guards, unit counts and signed controls")


if __name__ == "__main__":
    main()
