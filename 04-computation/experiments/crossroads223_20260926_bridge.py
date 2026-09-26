"""Exact finite controls for the rational bridge theorem; stdlib only.

Run normally and with -O. Integer acceptance and all matrix minors are
FINITE-EXACT; sine/eigenvalue/analytic-bound probes are VERIFIED numerical.
No assertion statement carries a validation gate.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import ceil, exp, log, pi, sin
from pathlib import Path
import hashlib

C = log(2) / log(3)
D = 1 - C
LOG_SCALE = -log(2) - C * log(C) - D * log(D)
A = pi * pi * C * D / 2
KAPPA = 1.5 * (pi * pi * C * D) ** (1 / 3) * log(3) ** (2 / 3)
THETA = (pi * pi * C * D / log(3)) ** (1 / 3)


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def endpoint(length):
    """ceil(length log_3 2), using integers exclusively."""
    e, power = 0, 1
    while power <= 1 << length:
        power *= 3
        e += 1
    return e


def matrix(length, width):
    e = endpoint(length)
    answer = []
    for start in range(1, width):
        states = {start: 1}
        for time in range(1, length + 1):
            nxt = defaultdict(int)
            for a, count in states.items():
                for target in (a, a + 1):
                    if 0 < length * target - e * time < width * length:
                        nxt[target] += count
            states = nxt
        answer.append([states.get(e + j, 0) for j in range(1, width)])
    return answer


def mul(left, right):
    size = len(left)
    return [[sum(left[i][k] * right[k][j] for k in range(size))
             for j in range(size)] for i in range(size)]


def trace(mat):
    return sum(row[i] for i, row in enumerate(mat))


def log_mu(p, width):
    b, t = 1 - p, pi / width
    beta = log(b * sin(p * t) / (p * sin(b * t)))
    return log(b) - p * beta + log(sin(t)) - log(sin(b * t)), beta


def direct_bridges(length, width):
    e = endpoint(length)
    counts = 0
    words = set()
    outputs = defaultdict(int)
    for chosen in combinations(range(length), e):
        bits = tuple(int(j in chosen) for j in range(length))
        prefixes, u = [0], 0
        for j, bit in enumerate(bits, 1):
            u += bit
            prefixes.append(length * u - e * j)
        starts = sum(all(0 < length * s + y < width * length
                         for y in prefixes)
                     for s in range(1, width))
        if not starts:
            continue
        counts += starts
        words.add(bits)
        rotation = min(range(length), key=prefixes.__getitem__)
        rotated = bits[rotation:] + bits[:rotation]
        outputs[rotated] += 1
        u = 0
        for j, bit in enumerate(rotated, 1):
            u += bit
            check(3 ** u > 1 << j, 'rotation must be strictly Collatz bad')
            check(3 ** u < 3 ** (width + 1) * (1 << j), 'rotated hard cap')
        check(starts <= width - 1, 'start multiplicity')
    check(all(n <= length for n in outputs.values()), 'rotation multiplicity')
    check(len(outputs) * (width - 1) * length >= counts, 'total multiplicity')
    return counts, len(words), len(outputs)


def actual_peak(length):
    """Independent exact tree, retaining a rational peak at every node."""
    states = [(0, Fraction(1))]
    for j in range(1, length + 1):
        nxt = []
        for u, peak in states:
            for target in (u, u + 1):
                if 3 ** target > 1 << j:
                    newpeak = max(peak, Fraction(3 ** target, 1 << j)) if j < length else peak
                    nxt.append((target, newpeak))
        states = nxt
    return sum((1 / peak for _, peak in states), Fraction()) / (1 << length)


def main():
    matrices = minors = trace_powers = direct_words = rotations = 0
    numeric_trace_checks = numeric_vector_checks = 0
    largest_lower_to_trace = 0.0
    for length in range(3, 41):
        e = endpoint(length)
        check(3 ** (e - 1) < 1 << length < 3 ** e, 'exact endpoint rounding')
        check(5 * e <= 4 * length, 'endpoint probability <=4/5')
        for width in range(2, 10):
            mat = matrix(length, width)
            matrices += 1
            for i, k in combinations(range(width - 1), 2):
                for j, ell in combinations(range(width - 1), 2):
                    check(mat[i][j] * mat[k][ell] >= mat[i][ell] * mat[k][j], 'TP2 minor')
                    minors += 1
            power = mat
            for n in range(2, 5):
                power = mul(power, mat)
                check(trace(power) <= trace(mat) ** n, 'trace power bound')
                trace_powers += 1
            if length <= 14:
                total, count, output = direct_bridges(length, width)
                check(total == trace(mat), 'matrix trace vs direct bridges')
                direct_words += count
                rotations += output
            p = e / length
            logw = e * log(p) + (length - e) * log(1 - p)
            lm, beta = log_mu(p, width)
            ratio = exp(length * lm - log(trace(mat)) - logw)
            check(ratio <= 1 + 2e-12, 'numerical spectral trace lower')
            largest_lower_to_trace = max(largest_lower_to_trace, ratio)
            numeric_trace_checks += 1
            f = [exp(beta * x) * sin(pi * x / width) for x in range(1, width)]
            for i, row in enumerate(mat):
                actual = sum(exp(log(count) + logw + (j - i) * log(p / (1-p))) * f[j]
                             for j, count in enumerate(row) if count)
                lower = exp(length * lm) * f[i]
                check(actual + 2e-12 * lower >= lower, 'numerical Pf lower')
                numeric_vector_checks += 1

    larger = []
    for length in (64, 128, 223, 256, 512, 1024):
        width = ceil(THETA * length ** (1/3))
        mat = matrix(length, width)
        e = endpoint(length)
        p = e / length
        logw = e * log(p) + (length-e) * log(1-p)
        lm, _ = log_mu(p, width)
        logtrace = log(trace(mat)) + logw
        check(length * lm <= logtrace + 1e-10, 'large trace spectral lower')
        entropy_loss = -length * log(2) - logw - length * LOG_SCALE
        check(entropy_loss >= -log(4)-1e-10, 'constant entropy rounding loss')
        lm_c, _ = log_mu(C, width)
        check(lm >= lm_c-1e-13, 'mu monotonicity in p')
        lower = (-(width+1)*log(3) + length*lm - log((width-1)*length)
                 - length*log(2) - logw)
        explicit = (length*LOG_SCALE - KAPPA*length**(1/3)
                    - log(108)-1-(4/3)*log(length))
        check(lower >= explicit-1e-9, 'explicit logarithmic error constant')
        larger.append((length, width, round(logtrace-length*lm, 9)))

    peak_checks = []
    for length in range(8, 17):
        value = actual_peak(length)
        exactlog = log(value.numerator)-log(value.denominator)
        lower = (length*LOG_SCALE-KAPPA*length**(1/3)-log(108)-1-(4/3)*log(length))
        check(lower <= exactlog, 'finite actual peak lower')
        peak_checks.append(length)

    hostile = [[0, 1], [1, 0]]
    check(trace(hostile) == 0 and trace(mul(hostile, hostile)) == 2,
          'non-TP2 hostile must fail trace-power property')
    print('source_sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('FINITE-EXACT universe: L=3..40, M=2..9; explicit endpoint rounding')
    print('matrices', matrices, 'all_2x2_minors', minors, 'trace_power_checks', trace_powers)
    print('direct bridge controls L=3..14; surviving words summed over widths', direct_words,
          'rotated outputs', rotations)
    print('VERIFIED numerical trace/vector probes', numeric_trace_checks, numeric_vector_checks)
    print('largest spectral_lower / trace_probability', format(largest_lower_to_trace, '.12g'))
    print('larger (L,M,log_trace_minus_L_log_mu)', larger)
    print('exact peak trees with VERIFIED analytic-lower comparison', peak_checks)
    print('non-TP2 hostile detected; no inferred theorem from numerical probes')
    print('PASS')


if __name__ == '__main__':
    main()
