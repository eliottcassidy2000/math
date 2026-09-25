"""Exact decoder/carry controls; no assertion of universal Collatz convergence.

Run from the repository root with Python 3, normally or with -O.
Universe: all binary words of lengths 1..12; all starts 3..4096 for
the root-certificate control (10000 shortcut steps permitted per start).
Only the standard library, integers, and Fraction are used.
"""

from fractions import Fraction
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def collatz(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def ceiling(n):
    return (3 * n + n % 2) // 2


def carries(word):
    r = c = a = 0
    for j, e in enumerate(word):
        r = (3 ** e) * r + e * (2 ** j)
        c = 3 * c + e * (2 ** j)
        a += e
    return a, r, c


def direct_carries(word):
    r = sum(e * 2 ** j * 3 ** sum(word[j + 1:])
            for j, e in enumerate(word))
    c = sum(e * 2 ** j * 3 ** (len(word) - 1 - j)
            for j, e in enumerate(word))
    return r, c


def parity_trace(n, length, step):
    result = []
    for _ in range(length):
        result.append(n % 2)
        n = step(n)
    return tuple(result), n


def rational_ceiling_trace(n, length):
    result = []
    for _ in range(length):
        require(n.denominator % 2 == 1, "non-2-adic-integral fraction")
        e = n.numerator % 2
        result.append(e)
        n = (3 * n + e) / 2
    return tuple(result), n


def main():
    word_count = swap_count = lift_count = 0
    for length in range(1, 13):
        q = 2 ** length
        collatz_residues = set()
        ceiling_residues = set()
        by_weight = {}
        for word in product((0, 1), repeat=length):
            word_count += 1
            a, r, c = carries(word)
            require((r, c) == direct_carries(word), "independent carry sum")
            p = 3 ** a
            rc = (-r * pow(p, -1, q)) % q
            rm = (-c * pow(3 ** length, -1, q)) % q
            collatz_residues.add(rc)
            ceiling_residues.add(rm)
            by_weight.setdefault(a, []).append(r)
            for lift in (0, 1, 3):
                nc, nm = rc + lift * q, rm + lift * q
                wc, outc = parity_trace(nc, length, collatz)
                wm, outm = parity_trace(nm, length, ceiling)
                require(wc == word and wm == word, "cylinder decoding")
                require(q * outc == p * nc + r, "Collatz affine formula")
                require(q * outm == 3 ** length * nm + c,
                        "ceiling affine formula")
                wnext, outnext = parity_trace(nc + q, length, collatz)
                require(wnext == word and outnext - outc == p,
                        "dyadic translation law")
                lift_count += 1
            for j in range(length - 1):
                if word[j:j + 2] != (0, 1):
                    continue
                swap = word[:j] + (1, 0) + word[j + 2:]
                _, rs, cs = carries(swap)
                suffix = word[j + 2:]
                require(r - rs == 2 ** j * 3 ** sum(suffix),
                        "Collatz adjacent exchange")
                require(c - cs == -(2 ** j * 3 ** len(suffix)),
                        "ceiling adjacent exchange")
                swap_count += 1
        require(len(collatz_residues) == len(ceiling_residues) == q,
                "full cylinder bijection")
        for a, values in by_weight.items():
            require(len(values) == len(set(values)), "carry collision")
            lo = 3 ** a - 2 ** a
            hi = 2 ** (length - a) * lo
            require(min(values) == lo and max(values) == hi,
                    "fixed-clock carry extrema")

    require(carries((0, 1)) == (1, 2, 2), "01 hostile")
    require(carries((1, 0)) == (1, 1, 3), "10 positive control")
    require(Fraction(4 * 4 - 1, 3) == 5, "10 root decode")
    require(Fraction(4 * 4 - 2, 3) == Fraction(14, 3), "01 noninteger decode")

    target_image = Fraction(-4, 15)
    require(rational_ceiling_trace(target_image, 18)[0]
            == parity_trace(4, 18, collatz)[0], "root4 Mahler image")
    bad = (1, 0, 1, 0, 1)
    _, _, badcarry = carries(bad)
    require(Fraction(2 * badcarry, 3 ** len(bad)) == Fraction(266, 243),
            "inherited 10101 unsafe tail")
    samples = max_depth = 0
    for start in range(3, 4097):
        n, bits = start, []
        for _ in range(10000):
            if n == 4:
                break
            bits.append(n % 2)
            n = collatz(n)
        require(n == 4, "finite root control did not terminate")
        word = tuple(bits)
        length = len(word)
        a, r, c = carries(word)
        require(3 ** a * start + r == 4 * 2 ** length,
                "root certificate equality")
        image = (2 ** length * target_image - c) / (3 ** length)
        require(image < 0 and image.denominator % 2 == 1,
                "root-basin image loses positive ordinary domain")
        expected, _ = parity_trace(start, length + 18, collatz)
        observed, _ = rational_ceiling_trace(image, length + 18)
        require(expected == observed, "same-word conjugacy on root controls")
        require(expected[length + 2:length + 7] == bad,
                "root certificate supplies Mahler hostile suffix")
        samples += 1
        max_depth = max(max_depth, length)

    # A unit gap is neither the root-decoder test nor the general cycle gate.
    require(2 ** 11 - 3 ** 7 == -139, "signed-cycle hostile gap")
    minus_cycle = (17, 25, 37, 55, 41, 61, 91)
    valuation_word = []
    for j, n in enumerate(minus_cycle):
        m, exponent = 3 * n - 1, 0
        while m % 2 == 0:
            m //= 2
            exponent += 1
        require(m == minus_cycle[(j + 1) % len(minus_cycle)],
                "known minus-cycle control")
        valuation_word.append(exponent)
    k = 0
    b = 0
    for j, exponent in enumerate(valuation_word):
        b += 3 ** (len(valuation_word) - 1 - j) * 2 ** k
        k += exponent
    require(k == 11 and b == 2363 == 17 * 139, "carry cancels nonunit gap")

    print("STATUS FINITE-EXACT; general statements proved in companion note")
    print("binary_words_lengths_1_to_12", word_count)
    print("cylinder_lift_controls", lift_count)
    print("adjacent_exchange_controls", swap_count)
    print("all_cylinder_bijections_and_fixed_weight_carry_extrema PASS")
    print("clock_2_1_collatz_carries 01:2 10:1; ceiling_carries 01:2 10:3")
    print("clock_2_1_root_candidates 01:14/3 10:5")
    print("root4_same_word_ceiling_image", target_image)
    print("root_controls_start_3_to_4096", samples, "maximum_depth", max_depth)
    print("all_root_images_negative_and_same_word PASS")
    print("root_suffix_10101_mahler_partial_tail", Fraction(266, 243))
    print("signed_minus_cycle_exponents", valuation_word,
          "power_gap", -139, "carry", b)
    print("SCOPE no universal root certificate or Mahler Z-number conclusion")


if __name__ == "__main__":
    main()
