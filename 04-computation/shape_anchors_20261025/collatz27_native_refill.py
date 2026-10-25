#!/usr/bin/env python3
"""Native 27-word exhaustion/refill; standard library, no producer imports."""
from hashlib import sha256
from pathlib import Path
import argparse
import json

Q = 2 ** 59
A = 3 ** 37
B = 23 * Q - 27 * A
D = Q - A
C = 23 * D - 2 * B
WORD = '11011111010110111011110100111011011111100111100010101000100'


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(x):
    check(x != 0, 'zero deviation')
    x = abs(x)
    return (x & -x).bit_length() - 1


def bits(n):
    return valuation(D * n - B)


def shortcut(n):
    return (3 * n + 1) // 2 if n & 1 else n // 2


def ordinary(n):
    return 3 * n + 1 if n & 1 else n // 2


def guarded_replay(n, copies):
    source = n
    ordinary_clock = 0
    for _ in range(copies):
        block_start = n
        for i, bit in enumerate(WORD):
            check(n & 1 == int(bit), 'wrong native parity')
            ordinary_clock += 1 + (n & 1)
            n = shortcut(n)
            if i < 58:
                check(n > block_start, 'earlier-than-59 descent')
        check(n < block_start, 'missing block descent')
    check(copies == 0 or n < source, 'original source lost')
    return n, ordinary_clock


def t_progression(h):
    modulus = 1 << (59 * h + 2)
    residue = ((1 << (59 * h + 1)) - C) * pow(A * D, -1, modulus) % modulus
    check(residue > 0 and residue & 1, 'positive odd residue failed')
    return residue, modulus


def direct_phase_check(n):
    """Independent fixed-denominator rational orbit, not affine block powers."""
    ell = bits(n)
    native = n
    numerator = B
    clock = 0
    for j in range(ell):
        check((native & 1) == (numerator & 1), 'early phase departure')
        check((numerator & 1) == int(WORD[j % 59]), 'rational periodic word')
        clock += 1 + (native & 1)
        native = shortcut(native)
        numerator = (3 * numerator + D) // 2 if numerator & 1 else numerator // 2
    check((native & 1) != (numerator & 1), 'missing exact phase departure')
    check(clock == 96 * (ell // 59) + ell % 59 + WORD[:ell % 59].count('1'), 'ordinary phase clock')
    return {'n': n, 'ell': ell, 'departure_value': native, 'ordinary_clock': clock}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    check(len(WORD) == 59 and WORD.count('1') == 37, 'word shape')
    check(C == 31 * A - 23 * Q, 'cofactor identity')
    base, base_clock = guarded_replay(27, 1)
    check(base == 23 and base_clock == 96, 'base word')
    numerator = B
    for bit in WORD:
        check(numerator & 1 == int(bit), 'rational cycle parity')
        numerator = (3 * numerator + D) // 2 if numerator & 1 else numerator // 2
    check(numerator == B, 'rational cycle does not close')

    fixtures = []
    boundary_count = 0
    same_residue = []
    for h in range(1, 65):
        residue, modulus = t_progression(h)
        for lift in (0, 1, 2, (1 << 80) + 1):
            t = residue + modulus * lift
            source = 27 + Q * t
            middle, clock = guarded_replay(source, 1)
            check(middle == 23 + A * t and not middle & 1, 'actual even endpoint')
            after = shortcut(middle)
            check(after == middle // 2 and clock + 1 == 97, 'native even step/clock')
            check((bits(source), bits(middle), bits(after)) == (59, 0, 59 * h), 'exact fuel transition')
            check((source > middle > after > 0), 'positivity or source descent')
            terminal, suffix_clock = guarded_replay(after, h)
            check(bits(terminal) == 0, 'suffix not exactly exhausted')
            check(suffix_clock == 96 * h, 'suffix ordinary clock')
            if h <= 4 and lift == 0:
                # Independent ordinary map, with no shortcut clock arithmetic.
                ordinary_state = source
                for _ in range(96):
                    ordinary_state = ordinary(ordinary_state)
                check(ordinary_state == middle, 'ordinary 96 endpoint')
                check(ordinary(ordinary_state) == after, 'ordinary 97 endpoint')
                fixtures.append({'h': h, 't': t, 'source': source,
                                 'at_shortcut_59': middle, 'at_shortcut_60': after,
                                 'fuel': [1, 0, h], 'bit_depth': [59, 0, 59 * h]})
            boundary_count += 1
        # Force identical mod-seven state at all three boundaries across h.
        lift = (-residue * pow(modulus, -1, 7)) % 7
        t = residue + modulus * lift
        source = 27 + Q * t
        middle = 23 + A * t
        after = middle // 2
        check(t % 7 == 0, 'CRT residue')
        check((source % 7, middle % 7, after % 7) == (6, 2, 1), 'residue transport')
        check((bits(source), bits(middle), bits(after)) == (59, 0, 59 * h), 'same residue fuel')
        if h <= 3:
            same_residue.append({'h': h, 'source': source, 'middle': middle, 'after': after,
                                 'residues': [6, 2, 1], 'bit_depth': [59, 0, 59 * h]})
    phase_count = 0
    for n in range(1, 513):
        direct_phase_check(n)
        phase_count += 1
    for ell in range(129):
        modulus = 1 << (ell + 1)
        residue = (B + (1 << ell)) * pow(D, -1, modulus) % modulus
        n = residue or modulus
        check(bits(n) == ell, 'phase boundary source')
        direct_phase_check(n)
        phase_count += 1
    # Mutation controls: an unguarded affine word and an omitted even reset.
    check((A * 28 + B) % Q != 0, 'unguarded affine hostile')
    check(bits(fixtures[0]['at_shortcut_59']) // 59 != fixtures[0]['h'], 'omitted reset hostile')
    result = {'status': 'PROVED explicit native refill family; FINITE-EXACT replay',
              'constants': {'q': Q, 'a': A, 'b': B, 'd': D, 'c': C},
              'refill_checks': boundary_count, 'same_residue_checks': 64,
              'independent_phase_checks': phase_count, 'fixtures': fixtures,
              'same_residue_fixtures': same_residue,
              'phase_controls': [direct_phase_check(27), direct_phase_check(27 + Q)],
              'scope': 'No convergence, unrestricted coverage, or globally decreasing fuel rank'}
    payload = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({'refill_checks': boundary_count, 'independent_phase_checks': phase_count,
                      'sha256': sha256(payload.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
