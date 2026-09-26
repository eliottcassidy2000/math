"""Exact residue projection controls; analytic displays are diagnostics only."""
from fractions import Fraction
from itertools import product
from math import cos, floor, log2, pi, sin
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def data(word):
    odd, carry, good = 0, 0, True
    for j, bit in enumerate(word):
        if bit:
            odd += 1
            carry = 3 * carry + (1 << j)
        good = good and 3 ** odd > 2 ** (j + 1)
    modulus = 1 << len(word)
    residue = (-carry * pow(3 ** odd, -1, modulus)) % modulus
    return odd, carry, residue, good


def actual_word(n, length):
    out = []
    for _ in range(length):
        bit = n & 1
        out.append(bit)
        n = (3 * n + 1) // 2 if bit else n // 2
    return tuple(out)


def first_before(word, first_index, second_index):
    ones = [j for j, bit in enumerate(word) if bit]
    zeros = [j for j, bit in enumerate(word) if not bit]
    return ones[first_index - 1] < zeros[second_index - 1]


def residue_controls():
    totals = []
    counts = 0
    for length in range(1, 15):
        modulus = 1 << length
        residues = set()
        positive = 0
        by_count = {}
        for word in product((0, 1), repeat=length):
            odd, carry, residue, good = data(word)
            require(actual_word(residue, length) == word, "fixed source lift")
            require(residue not in residues, "parity cylinder bijection")
            residues.add(residue)
            if good:
                positive += 1
                cyclic = by_count.setdefault(odd, {})
                exponent = carry % modulus
                require(exponent not in cyclic, "one atom per residue at fixed odd count")
                cyclic[exponent] = word
                require((carry + 3 ** odd * residue) % modulus == 0, "Fourier projector congruence")
            counts += 1
        require(len(residues) == modulus, "complete universe")
        totals.append((length, positive))
    words = [w for w in product((0, 1), repeat=6) if sum(w) == 5 and data(w)[3]]
    rows = [("".join(map(str, w)), data(w)[1], data(w)[2], first_before(w, 4, 1)) for w in words]
    require([r[2] for r in rows] == [27, 39, 47, 31], "27/31 source control")
    require(sum(r[3] for r in rows) == 2, "uniform poset balanced pair")
    for source in (27, 31):
        selected = [w for w in words if data(w)[2] == source]
        require(len(selected) == 1, "fixed integer has a singleton extension law")
    return {"all_words_checked": counts, "positive_prefix_counts": totals, "six_event_rows": rows,
            "uniform_A4_before_B1": "1/2", "source27_probability": 0, "source31_probability": 1}


def analytic_diagnostics():
    def vf(x):
        t = log2(x)
        u = t - floor(t)
        return floor(t) + cos(pi * u / 2) ** 2
    def vg(x):
        t = log2(x)
        return t + sin(2 * pi * t) / (6 * pi)
    rows = []
    for n in (1, 3, 19, 27, 31, 41, 703, 10087):
        target = (3 * n + 1) // 2
        rows.append({"n": n, "Tn": target, "delta_log2F": round(vf(target) - vf(n), 12),
                     "delta_VG": round(vg(target) - vg(n), 12)})
    # Exact all-odd paths; no numerical logarithm is used in their proof.
    for k in range(2, 301):
        n = 2 ** k - 1
        for j in range(k):
            require(n % 2 == 1, "all-odd growing source")
            n = (3 * n + 1) // 2
            require(n == 3 ** (j + 1) * 2 ** (k - j - 1) - 1, "all-odd affine identity")
        require(Fraction(n, 2 ** k - 1) > Fraction(3, 2) ** k, "bounded correction obstruction")
    # At u=1/4 the sine series is exactly 1-1/3+1/5-... = pi/4.
    # Hence the user's printed expression is 0, while {u}=1/4.
    return {"float_diagnostics_not_proofs": rows, "Fourier_quarter_exact": {"printed": "0", "fractional_part": "1/4"},
            "all_odd_lengths_checked": [2, 300], "dyadic_jump_ratio_exact": 4,
            "VG_derivative_in_log2_coordinate": "1+(1/3)cos(2*pi*t), in [2/3,4/3]"}


def main():
    print(json.dumps({"scope": "FINITE-EXACT finite parity universes; analytic floats diagnostic only",
                      "residue_projection": residue_controls(), "analytic": analytic_diagnostics()}, indent=2))


if __name__ == "__main__":
    main()
