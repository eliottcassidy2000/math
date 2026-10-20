#!/usr/bin/env python3
"""Independent exact audit of the negative one-adjacent-swap mechanism.

No repository mathematical modules or production certificates are imported.
All pass/fail decisions use integers or fractions; checks survive python -O.
"""

from fractions import Fraction
from hashlib import sha256
from math import factorial, gcd
from pathlib import Path
import json


F = Fraction
checks = 0


def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def log_interval(t, terms=300):
    """Enclose -log(1-t) by its positive ordinary power series."""
    total = sum((t**j / j for j in range(1, terms + 1)), F(0))
    tail = t ** (terms + 1) / ((terms + 1) * (1 - t))
    return total, total + tail


def carry(word):
    """Accumulate odd-position contributions from right to left."""
    total, suffix = 0, 1
    for pos in range(len(word) - 1, -1, -1):
        if word[pos]:
            total += (1 << pos) * suffix
            suffix *= 3
    return total


def affine_carry(word):
    total, denominator = 0, 1
    for bit in word:
        total = (3 if bit else 1) * total + bit * denominator
        denominator *= 2
    return total


def mechanical(p, q):
    return tuple(((i + 1) * p) // q - (i * p) // q for i in range(q))


def swaps(word):
    for i, bit in enumerate(word):
        j = (i + 1) % len(word)
        if bit != word[j]:
            edited = list(word)
            edited[i], edited[j] = edited[j], edited[i]
            yield i, tuple(edited)


def shortcut(x):
    return (3 * x + 1) // 2 if x % 2 else x // 2


def integer_orbit(word, clock):
    b = carry(word)
    require(b % clock == 0, ("integer start", word))
    start = b // clock
    x, entries = start, []
    for bit in word:
        require(x % 2 == bit, ("parity", word, x))
        entries.append(x)
        x = shortcut(x)
    require(x == start, ("cycle closure", word))
    return entries


def determinant_bareiss(matrix):
    """Fraction-free determinant with row pivoting, independent of SymPy."""
    a = [list(row) for row in matrix]
    n, sign, previous = len(a), 1, 1
    if n == 0:
        return 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((j for j in range(k + 1, n) if a[j][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * a[i][j] - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss exact division")
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def sparse_norm(p, a, epsilon):
    """Determinant of multiplication by Q in Z[Z]/(Z^p-2)."""
    terms = ((0, 1), (a + 1, 1), (a, -1)) if epsilon == 1 else (
        (0, 1), (p - a, 1), (1, -1)
    )
    matrix = [[0] * p for _ in range(p)]
    for col in range(p):
        for degree, coefficient in terms:
            quotient, remainder = divmod(degree + col, p)
            matrix[remainder][col] += coefficient * (2**quotient)
    value = abs(determinant_bareiss(matrix))
    require(value > 0, ("nonzero norm", p, a, epsilon))
    for prime in (2, 3):
        while value % prime == 0:
            value //= prime
    return value


# Different logarithm expansion from the production atanh certificate.
L2, U2 = log_interval(F(1, 2))
L32, U32 = log_interval(F(1, 3))
L3, U3 = L2 + L32, U2 + U32
alpha_lower, alpha_upper = L3 / U2, U3 / L2
require(U2 < F(7, 10) and U3 < F(11, 10), "Matveev height choices")
require(2**9 < 23**2, "2^(9/2) < 23")
constant_upper = F(7, 5) * 30**5 * 23 * F(7, 10) * F(11, 10)
require(constant_upper == 602494200 < 10**9, "Matveev constant")
require(3**5 < 2**8, "negative residue ratio below 3/5")
require(F(4**3) < F(5, 2) ** 5, "4^(3/5) < 5/2")
require(L32 > F(1, 3), "log(3/2) > 1/3")
require(F(4, 9) < F(1, 2), "p>=4 lower critical clock")
require(F(7, 10) + 12 * (3 * F(7, 10) + F(1, 4)) < 30,
        "log(2*10^12) < 30")
require(6 * 10**9 * 31 < 10**12, "global cutoff")
exp_lower = sum((F(100, 6) ** j / factorial(j) for j in range(20)), F(0))
require(F(400) / exp_lower < L2, "Legendre at p=100")

# Certified CF interval transforms and independent positive absolute gaps.
lo, hi = alpha_lower, alpha_upper
hm2, hm1, km2, km1 = 0, 1, 1, 0
digits, convergents = [], []
for index in range(100):
    digit = lo.numerator // lo.denominator
    require(digit == hi.numerator // hi.denominator, ("unique CF digit", index))
    digits.append(digit)
    q, p = digit * hm1 + hm2, digit * km1 + km2
    hm2, hm1, km2, km1 = hm1, q, km1, p
    require(gcd(p, q) == 1, "reduced convergent")
    if 100 <= p < 10**12:
        signed_lo, signed_hi = q * L2 - p * U3, q * U2 - p * L3
        if signed_hi < 0:
            side, abs_lower = "below", -signed_hi
        else:
            require(signed_lo > 0, ("separated gap", q, p))
            side, abs_lower = "above", signed_lo
        exponent = min(p // 6, 256)
        rational_barrier = F(2, 2**exponent)
        require(abs_lower > rational_barrier, ("absolute CF exclusion", q, p))
        convergents.append({
            "q": q, "p": p, "side": side, "cap_exponent": exponent,
            "absolute_gap_lower": str(abs_lower),
            "rational_exclusion_barrier": str(rational_barrier),
        })
    if p >= 10**12:
        terminal = {"q": q, "p": p}
        break
    lo, hi = 1 / (hi - digit), 1 / (lo - digit)
else:
    raise RuntimeError("No terminal convergent")

# Exact residual universe p=2,...,99, using bit_length rather than logarithms.
small_rows, small_hits = [], []
small_swaps = small_shapes = 0
for p in range(2, 100):
    q = (3**p).bit_length() - 1
    require(2**q < 3**p < 2 ** (q + 1), ("floor critical clock", p))
    require(int(p * alpha_lower) == int(p * alpha_upper) == q,
            ("independent critical clock agreement", p))
    row = {"p": p, "q": q, "swaps": 0}
    if gcd(p, q) != 1:
        row["status"] = "nonprimitive; repetition obstruction applies"
        small_rows.append(row)
        continue
    small_shapes += 1
    word, clock, maximum = mechanical(p, q), 2**q - 3**p, 1
    for index, edited in swaps(word):
        b = carry(edited)
        require(b == affine_carry(edited), ("two carries", p, q, index))
        g = gcd(-clock, b)
        require(g * g < 6**p, ("small sqrt6 gate", p, q, index))
        maximum = max(maximum, g)
        if g == -clock:
            small_hits.append({
                "p": p, "q": q, "swap_index": index,
                "word": "".join(map(str, edited)), "clock": clock,
                "carry": b, "orbit": integer_orbit(edited, clock),
            })
        else:
            require(g < -clock, ("proper content", p, q, index))
        row["swaps"] += 1
        small_swaps += 1
    row["status"] = "checked all unequal cyclic swaps"
    row["maximum_content"] = maximum
    small_rows.append(row)
require([(x["p"], x["q"], x["word"]) for x in small_hits] ==
        [(2, 3, "101"), (2, 3, "110")], "complete small exception set")

# A hostile census is separate from the proof's finite residual universe.
# Enumerate ALL negative primitive clocks through q=160, not only critical ones.
norm_cache = {}
census_shapes = census_swaps = census_crossings = 0
census120_swaps = census120_shapes = 0
census_hits = []
norm_checks = 0
for q in range(3, 161):
    for p in range(2, q):
        if gcd(p, q) != 1 or 2**q >= 3**p:
            continue
        census_shapes += 1
        if q <= 120:
            census120_shapes += 1
        r, clock = q - p, 2**q - 3**p
        require(5 * r < 3 * p, ("negative residue range", p, q))
        word = mechanical(p, q)
        labels = {q - 1: 0}
        for j in range(1, p):
            labels[j + (j * r + p - 1) // p - 1] = j
        require(set(labels) == {i for i, bit in enumerate(word) if bit},
                ("mechanical one labels", p, q))
        for position, j in labels.items():
            a = (-j * r) % p
            plus_legal = word[(position + 1) % q] == 0
            minus_legal = word[(position - 1) % q] == 0
            require(plus_legal == (a < r), ("plus legality", p, q, j))
            require(minus_legal == (a >= p - r), ("minus legality", p, q, j))
        for index, edited in swaps(word):
            source = index if word[index] else (index + 1) % q
            epsilon = 1 if word[index] else -1
            a = (-labels[source] * r) % p
            require((a, epsilon) not in ((0, -1), (p - 1, 1)),
                    ("large endpoints forbidden", p, q, index))
            b = carry(edited)
            require(b == affine_carry(edited), ("census independent carries", p, q, index))
            g = gcd(-clock, b)
            require(g * g < 6**p, ("all-negative sqrt6 gate", p, q, index))
            if (a, epsilon) in ((0, 1), (p - 1, -1)):
                require(g == 1, ("free endpoints", p, q, index))
            if p <= 18:
                key = (p, a, epsilon)
                if key not in norm_cache:
                    norm_cache[key] = sparse_norm(*key)
                require(norm_cache[key] % g == 0, ("full norm divisibility", key, q, g))
                norm_checks += 1
            if g == -clock:
                census_hits.append({
                    "p": p, "q": q, "swap_index": index,
                    "word": "".join(map(str, edited)),
                })
            census_swaps += 1
            census_crossings += index == q - 1
            census120_swaps += q <= 120
require(census120_swaps == 47892, "independent q<=120 count")
require([(x["p"], x["q"], x["word"]) for x in census_hits] ==
        [(2, 3, "101"), (2, 3, "110")], "complete hostile census exception set")

# Independent negative-cycle boundary control: derive its word from the map.
x, minus17_orbit, minus17_bits = -17, [], []
for _ in range(100):
    minus17_orbit.append(x)
    minus17_bits.append(x % 2)
    x = shortcut(x)
    if x == -17:
        break
else:
    raise RuntimeError("-17 control did not close")
require("".join(map(str, minus17_bits)) == "11110111000", "-17 parity control")
control_p, control_q = sum(minus17_bits), len(minus17_bits)
control_clock = 2**control_q - 3**control_p
require(control_clock == -139 and carry(minus17_bits) == -17 * control_clock,
        "-17 affine control")
control_base = mechanical(control_p, control_q)
distances = [
    sum(minus17_bits[(i + shift) % control_q] != control_base[i]
        for i in range(control_q))
    for shift in range(control_q)
]
require(min(distances) == 4, "-17 remains outside the one-swap family")

# Repetition-factor hostile checks over arbitrary nonconstant bases.
repeated_swaps = repeated_crossings = negative_repeated_swaps = 0
for q0 in range(2, 9):
    for number in range(1, 2**q0 - 1):
        base = tuple((number >> i) & 1 for i in range(q0))
        p0 = sum(base)
        clock0 = 2**q0 - 3**p0
        for repetitions in range(2, 6):
            word = base * repetitions
            clock = 2 ** (q0 * repetitions) - 3 ** (p0 * repetitions)
            H = sum(2 ** (q0 * (repetitions - 1 - j)) * 3 ** (p0 * j)
                    for j in range(repetitions))
            require(H > 1 and gcd(H, 6) == 1 and clock == H * clock0,
                    "repetition multiplier")
            require(carry(word) == H * carry(base), "repeated carry")
            for index, edited in swaps(word):
                b = carry(edited)
                require(gcd(H, b) == 1, "entire repetition factor destroyed")
                require(gcd(abs(clock), b) == gcd(abs(clock0), b) < abs(clock),
                        "exact content reduction for either clock sign")
                repeated_swaps += 1
                repeated_crossings += index == len(word) - 1
                negative_repeated_swaps += clock < 0
require(repeated_swaps == 25088 and repeated_crossings == 1016,
        "independent repetition census counts")

report = {
    "status": "PASS: independent exact negative one-swap audit",
    "scope": "One cyclic adjacent unequal-bit swap of a mechanical word; not arbitrary negative cycles",
    "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    "explicit_checks": checks,
    "log_series": "-log(1-t), 300 terms; log3=log2+log(3/2)",
    "matveev_source": "https://arxiv.org/pdf/math/0403046",
    "matveev_locator": "Theorem 9.4, printed page 16",
    "matveev_prefactor_upper": str(constant_upper),
    "global_weight_cutoff_exclusive": 10**12,
    "legendre_threshold": 100,
    "continued_fraction_digits": digits,
    "large_convergents": convergents,
    "large_below_count": sum(row["side"] == "below" for row in convergents),
    "large_above_count": sum(row["side"] == "above" for row in convergents),
    "terminal_convergent": terminal,
    "small_primitive_shapes": small_shapes,
    "small_swap_occurrences": small_swaps,
    "small_integer_hits": small_hits,
    "small_rows": small_rows,
    "hostile_census": {
        "maximum_q": 160,
        "primitive_shapes": census_shapes,
        "swap_occurrences": census_swaps,
        "crossing_cut_swaps": census_crossings,
        "through_q120_shapes": census120_shapes,
        "through_q120_swaps": census120_swaps,
        "integer_hits": census_hits,
        "norm_determinants": len(norm_cache),
        "norm_gate_checks_through_p18": norm_checks,
    },
    "minus17_boundary": {
        "orbit": minus17_orbit,
        "parity_word": "".join(map(str, minus17_bits)),
        "mechanical_word": "".join(map(str, control_base)),
        "p": control_p, "q": control_q, "clock": control_clock,
        "hamming_distances_over_rotations": distances,
        "minimum_hamming_distance": min(distances),
    },
    "arbitrary_repeated_bases": {
        "base_lengths": [2, 8], "repetitions": [2, 5],
        "swap_occurrences": repeated_swaps,
        "crossing_cut_swaps": repeated_crossings,
        "negative_clock_swap_occurrences": negative_repeated_swaps,
    },
    "rational_log_intervals": {
        "log2_lower": str(L2), "log2_upper": str(U2),
        "log3_lower": str(L3), "log3_upper": str(U3),
    },
}
print(json.dumps(report, indent=2, sort_keys=True))
