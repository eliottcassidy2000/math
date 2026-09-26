"""Exact critical-strip word counts and cyclic-rotation certificates.

Reproduce from repository root with:
python3 04-computation/experiments/crossroads_20260926_ballot.py
All acceptance decisions use integer powers, never floating logarithms.
Logarithms below describe integer counts, not validity gates.
"""
from itertools import product
from math import comb, isqrt, log2


def ceil_critical(q, b):
    k = 0
    while q ** k < 2 ** b:
        k += 1
    return k


def count_words(q, length, cap=None, endpoint=None):
    row = {0: 1}
    for j in range(1, length + 1):
        nxt = {}
        for o, count in row.items():
            for bit in (0, 1):
                u = o + bit
                if q ** u < 2 ** j:
                    continue
                if cap is not None and q ** u > cap * 2 ** j:
                    continue
                nxt[u] = nxt.get(u, 0) + count
        row = nxt
    return row.get(endpoint, 0) if endpoint is not None else sum(row.values())


def good(q, word, cap=None):
    a = 1
    for j, bit in enumerate(word, 1):
        a *= q ** bit
        if a < 2 ** j or (cap is not None and a > cap * 2 ** j):
            return False
    return True


def good_rotation(q, word):
    # Locate a minimum cumulative logarithm using exact rational comparison.
    a, minimum_a, minimum_j, cut = 1, 1, 0, 0
    for j, bit in enumerate(word[:-1], 1):
        a *= q ** bit
        if a * 2 ** minimum_j < minimum_a * 2 ** j:
            minimum_a, minimum_j, cut = a, j, j
    rotated = word[cut:] + word[:cut]
    assert good(q, rotated)
    return rotated


def source_residue(q, word, shift=1):
    residue, modulus = 0, 1
    for j, bit in enumerate(word):
        y = residue
        for _ in range(j):
            y = (q * y + shift) // 2 if y & 1 else y // 2
        if (y & 1) != bit:
            residue += modulus
        modulus *= 2
    return residue, modulus


def main():
    print('FINITE-EXACT: critical-strip counts, rotation and realization controls')
    checks = 0
    for q in (3, 5):
        for length in range(1, 13):
            words = list(product((0, 1), repeat=length))
            for cap in (None, 2, 4, 8, 16):
                actual = sum(good(q, w, cap) for w in words)
                assert actual == count_words(q, length, cap)
                checks += len(words)
            k = ceil_critical(q, length)
            inputs = [w for w in words if sum(w) == k]
            outputs = {good_rotation(q, w) for w in inputs}
            assert len(outputs) * length >= comb(length, k)
            assert len(outputs) <= count_words(q, length, endpoint=k)
            checks += len(inputs)
    print('independent enumeration / rotation word-checks:', checks)

    print('q L count_bad count_cap2 count_cap4 count_cap8 count_cap16')
    for q in (3, 5):
        for length in (8, 16, 32, 64, 128, 256):
            counts = [count_words(q, length, cap) for cap in (None, 2, 4, 8, 16)]
            print(q, length, *counts)

    print('block construction: q L b t r good_blocks log2(lower_count)/L log2(cap)')
    for q in (3, 5):
        for length in (16, 64, 256, 1024, 4096):
            b = max(1, isqrt(length))
            t, r = divmod(length, b)
            k = ceil_critical(q, b)
            block_count = count_words(q, b, endpoint=k)
            assert block_count * b >= comb(b, k)
            lower_count = block_count ** t
            cap = q ** (t + b + r)
            if length <= 256:
                assert count_words(q, length, cap) >= lower_count
            print(q, length, b, t, r, block_count,
                  format(log2(lower_count) / length, '.9f'),
                  format(log2(cap), '.9f'))

    realizations = 0
    for q in (3, 5):
        for shift in (-1, 1):
            for b in range(2, 10):
                k = ceil_critical(q, b)
                block = good_rotation(q, (0,) * (b-k) + (1,) * k)
                word = block * 7
                residue, modulus = source_residue(q, word, shift)
                n = residue + modulus
                y = n
                for bit in word:
                    assert (y & 1) == bit
                    y = (q * y + shift) // 2 if bit else y // 2
                    assert y > 0
                realizations += 1
    print('guarded positive integer realizations on both signs and multipliers:', realizations)
    print('PASS; finite critical strips do not establish an infinite integer survivor.')


if __name__ == '__main__':
    main()
