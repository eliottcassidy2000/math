"""Independent exact counts and numerical audits for the sine barrier proof.

Run with ordinary Python or -O. All gates are explicit exceptions.
Integer DP counts and word realizations are FINITE-EXACT. Transcendental
comparisons are VERIFIED numerical probes; the accompanying note proves
the inequalities analytically for every parameter in its stated domain.
"""
from collections import defaultdict
from itertools import product
from math import exp, log, log2, pi, sin
import hashlib
from pathlib import Path

C = log(2) / log(3)
D = 1 - C
ENTROPY = -C * log2(C) - D * log2(D)
LOG_SCALE = -(1 - ENTROPY) * log(2)
A = pi * pi * C * D / 2
KAPPA = 1.5 * (pi * pi * C * D) ** (1 / 3) * log(3) ** (2 / 3)


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def floors(length):
    """floor(j log_3 2), with exact integer acceptance decisions."""
    answer = [0]
    power, exponent = 1, 0
    for j in range(1, length + 1):
        while power * 3 <= 1 << j:
            power *= 3
            exponent += 1
        answer.append(exponent)
    return answer


def counts_all(length, m):
    fl = floors(length)
    states = {0: 1}
    counts = [1]
    for j in range(length):
        nxt = defaultdict(int)
        zone = m - 1 if j == 0 else m + fl[j]
        for u, mass in states.items():
            if u >= zone:
                if u > fl[j + 1]:
                    nxt[u] += 2 * mass
            else:
                if u > fl[j + 1]:
                    nxt[u] += mass
                if u + 1 > fl[j + 1]:
                    nxt[u + 1] += mass
        states = dict(nxt)
        counts.append(sum(states.values()))
    return counts


def word_survives(word, m):
    u = 0
    for j, bit in enumerate(word):
        zone = 3 ** u >= (3 ** (m - 1)) * (1 << j)
        if bit and not zone:
            u += 1
        if 3 ** u <= 1 << (j + 1):
            return False
    return True


def integer_realization(word, m):
    """Reconstruct the source residue from independent affine equations."""
    residue, u, carry = 0, 0, 0
    decisions = []
    for j, bit in enumerate(word):
        image = (3 ** u * residue + carry) // (1 << j)
        if image % 2 != bit:
            residue += 1 << j
        check(((3 ** u * residue + carry) // (1 << j)) % 2 == bit,
              'source residue parity')
        flip = bool(bit and 3 ** u >= 3 ** (m - 1) * (1 << j))
        decisions.append(flip)
        if flip:
            carry -= 1 << j
        elif bit:
            carry = 3 * carry + (1 << j)
            u += 1
    start = residue + (1 << len(word)) * 10 ** 8
    value, ups, survived = start, 0, True
    for j, (bit, flip) in enumerate(zip(word, decisions)):
        check(value % 2 == bit, 'ordinary integer word')
        if flip:
            value = (value - 1) // 2
        elif bit:
            value = (3 * value + 1) // 2
            ups += 1
        else:
            value //= 2
        check(value > 0, 'positive realization')
        expected = 3 ** ups > 1 << (j + 1)
        check((value > start) == expected, 'affine descent sign')
        survived &= expected
    check(survived == word_survives(word, m), 'realization survival')


def spectral_parameters(m, c=C):
    d = 1 - c
    theta = pi / (m + 4)
    beta = log(d * sin(c * theta) / (c * sin(d * theta)))
    mu = d * exp(-c * beta) * sin(theta) / sin(d * theta)
    return theta, beta, mu


def main():
    print('ROBIN SINE CERTIFICATE: exact finite counts plus labelled numerical audits')
    print('source_sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    check(3 ** 5 < 2 ** 8 and 3 ** 2 > 2 ** 3, '5/8<c<2/3')
    exact = 0
    for length in range(1, 13):
        for m in range(2, 9):
            brute = sum(word_survives(w, m) for w in product((0, 1), repeat=length))
            check(brute == counts_all(length, m)[length], (length, m, brute))
            exact += 2 ** length
    print('FINITE-EXACT direct word checks', exact)

    realized = 0
    for length in (8, 12, 16):
        for m in (2, 3, 5, 8):
            for seed in (0, 1, 3, 17, 223, 65535):
                word = tuple((seed >> (j % 16)) & 1 for j in range(length))
                integer_realization(word, m)
                realized += 1
    print('FINITE-EXACT positive integer realizations', realized)

    numeric = 0
    worst_kernel, worst_eigen = 0.0, 0.0
    for c in (5 / 8, C, 2 / 3):
        d = 1 - c
        for m in tuple(range(2, 129)) + (223, 256, 512, 1024):
            theta, beta, mu = spectral_parameters(m, c)
            def f(x):
                return exp(beta * x) * sin(theta * (x + 1))
            check(beta < 0 and m * abs(beta) < 1, 'beta bound')
            eigen_gap = log(mu) + c * d * theta * theta / 2
            worst_eigen = max(worst_eigen, eigen_gap)
            check(eigen_gap <= 2e-15, ('spectral coefficient', c, m, eigen_gap))
            # Include both sides of the zone edge and the limiting top state.
            grid = [j * (m - c) / 400 for j in range(401)]
            grid += [m - 1 - 1e-9, m - 1, m - 1 + 1e-9]
            for x in grid:
                if x >= m - 1:
                    bf = 2 * d * f(x - c) if x > c else 0.0
                else:
                    bf = c * f(x + d) + (d * f(x - c) if x > c else 0.0)
                ratio = bf / (mu * f(x))
                worst_kernel = max(worst_kernel, ratio)
                check(ratio <= 1 + 3e-12, ('kernel', c, m, x, ratio))
                check(f(x) >= exp(-1) * f(0) - 1e-14, 'endpoint comparison')
                numeric += 1
    print('VERIFIED real-state kernel probes', numeric,
          'max_ratio', format(worst_kernel, '.12f'),
          'max_eigen_excess', format(worst_eigen, '.3g'))

    strongest = (0.0, None)
    count_gates = 0
    print('VERIFIED count bound: L m log(N/2^L) log(bound) ratio')
    for m in tuple(range(2, 41)) + (64, 100, 223):
        counts = counts_all(4096, m)
        for length, n in enumerate(counts[1:], 1):
            if not n:
                continue
            actual = log(n) - length * log(2)
            bound = 1 + length * LOG_SCALE - A * length / (m + 4) ** 2
            ratio = exp(actual - bound)
            if ratio > strongest[0]:
                strongest = (ratio, (length, m))
            check(actual <= bound + 1e-10, ('finite count inequality', length, m))
            count_gates += 1
            if length in (64, 256, 1024, 4096) and m in (4, 16, 64, 223):
                print(length, m, format(actual, '.9f'), format(bound, '.9f'),
                      format(ratio, '.8g'))
    print('VERIFIED finite-count comparisons', count_gates,
          'max ratio', strongest)
    print('A', format(A, '.12f'), 'kappa', format(KAPPA, '.12f'))
    print('PASS: no uniform Robin-ratio or globally consistent pairing claim.')


if __name__ == '__main__':
    main()
