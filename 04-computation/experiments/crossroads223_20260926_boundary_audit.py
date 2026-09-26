"""Independent integer control for THM4487's repaired gamma=1 boundary."""
from itertools import combinations
from pathlib import Path
import hashlib


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def main():
    words_checked = realizations = 0
    for length in range(4, 17):
        tail_length = length - 1
        e = 0
        while 3 ** e <= 1 << tail_length:
            e += 1
        e += 1
        rotated_words = set()
        for chosen in combinations(range(tail_length), e):
            bits = tuple(int(i in chosen) for i in range(tail_length))
            # Maximize i-e_i log2(3), by exactly minimizing 3^e_i/2^i.
            from fractions import Fraction
            slopes, ups = [Fraction(1)], 0
            for i, bit in enumerate(bits, 1):
                ups += bit
                slopes.append(Fraction(3 ** ups, 1 << i))
            rot = min(range(tail_length), key=slopes.__getitem__)
            full = (1,) + bits[rot:] + bits[:rot]
            rotated_words.add(full)
            ups = 0
            for i, bit in enumerate(full, 1):
                ups += bit
                check(2 * 3 ** ups >= 3 * (1 << i), 'prefix margin >=3/2')
            words_checked += 1
        for word in rotated_words:
            for sign in (-1, 1):
                residue = carry = ups = 0
                for j, bit in enumerate(word):
                    value = (3 ** ups * residue + carry) // (1 << j)
                    if value % 2 != bit:
                        residue += 1 << j
                    if bit:
                        carry = 3 * carry + sign * (1 << j)
                        ups += 1
                start = residue + (1 << length)
                value = start
                for bit in word:
                    check(value % 2 == bit, 'Terras realization parity')
                    value = (3 * value + sign) // 2 if bit else value // 2
                    check(value >= start, 'actual endpoint-gamma no dip on both sheets')
                realizations += 1
    print('source_sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('FINITE-EXACT universe: lengths4..16, prescribed surplus-count tails')
    print('all rotated inputs', words_checked, 'distinct full-word sheet realizations', realizations)
    print('integer >=3/2 prefix margins and actual no-dip condition: PASS')
    print('This tests the repaired boundary only, not all of THM4487.')


if __name__ == '__main__':
    main()
