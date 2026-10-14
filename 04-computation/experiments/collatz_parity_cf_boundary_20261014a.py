"""Exact finite controls for parity CFs and a stopping-gated CF exponent.

Run: python 04-computation/experiments/collatz_parity_cf_boundary_20261014a.py
No finite experiment decides an irrationality exponent or universal Collatz.
The corresponding note supplies the all-height proofs. All checks remain on
under python -O. Output is deterministic JSON on stdout.
"""
from fractions import Fraction
from itertools import product
import json
from math import isqrt

CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def step(n):
    if type(n) is not int:
        raise ValueError('integer state required')
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def bits(n, length):
    if type(n) is not int or type(length) is not int or length < 0:
        raise ValueError('integer state and nonnegative exact length required')
    out = []
    for _ in range(length):
        out.append(n % 2)
        n = step(n)
    return tuple(out)


def cylinder(word):
    """T^L(x)=(P*x+B)/Q on the word's unique residue c mod Q."""
    if type(word) is not tuple or any(type(b) is not int or b not in (0, 1) for b in word):
        raise ValueError('exact binary tuple required')
    P, B, Q = 1, 0, 1
    for bit in word:
        if bit:
            P, B = 3 * P, 3 * B + Q
        Q *= 2
    c = (-B * pow(P, -1, Q)) % Q
    return P, B, Q, c


def cf_matrix(digits):
    if type(digits) not in (tuple, list) or any(type(a) is not int or a < 1 for a in digits):
        raise ValueError('positive exact CF digits required')
    a, b, c, d = 1, 0, 0, 1
    for digit in digits:
        a, b, c, d = a * digit + b, a, c * digit + d, c
    return a, b, c, d


def convergent(digits):
    a, b, c, d = cf_matrix(digits)
    return Fraction(c, a)


def hitting_time(n, limit=4096):
    if type(n) is not int or n < 1 or type(limit) is not int or limit < 0:
        raise ValueError('positive source and nonnegative exact horizon required')
    for t in range(limit + 1):
        if n == 1:
            return t
        n = step(n)
    raise RuntimeError('finite control exceeded its horizon')


def gated_code(n, length=14, exponent=3):
    """a0=n, a1=2; by stage j inspect T through time 2^j.

    Emits q_j^(exponent-2) while unresolved and 1 after an actual visit to1.
    Stage j+1 is total computable; a capped run is only a prefix.
    """
    if (type(n) is not int or n < 1 or type(exponent) is not int or exponent < 3
            or type(length) is not int or length < 1):
        raise ValueError('positive exact source/length and integer exponent >=3 required')
    digits = [2]
    q_prev, q = 1, 2
    state, time, hit, switched_at = n, 0, (0 if n == 1 else None), None
    rows = []
    for j in range(1, length):
        while hit is None and time < 2 ** j:
            state = step(state)
            time += 1
            if state == 1:
                hit = time
        digit = 1 if hit is not None else q ** (exponent - 2)
        if hit is not None and switched_at is None:
            switched_at = j
        rows.append((j, q, digit, hit))
        digits.append(digit)
        q_prev, q = q, digit * q + q_prev
    return digits, rows, hit, switched_at


def main():
    hostiles = (
        lambda: step(True), lambda: bits(1, True), lambda: bits(1, -1),
        lambda: cylinder((2,)), lambda: cylinder((True,)), lambda: cylinder([1]),
        lambda: cf_matrix([0]), lambda: cf_matrix([True]), lambda: cf_matrix([1.0]),
        lambda: hitting_time(0), lambda: gated_code(True), lambda: gated_code(1, 0),
        lambda: gated_code(1, True), lambda: gated_code(1, exponent=2),
    )
    for hostile in hostiles:
        try:
            hostile()
        except ValueError:
            require(True, 'malformed interface rejected')
        else:
            require(False, 'malformed interface accepted')
    by_length = []
    for L in range(1, 13):
        seen = set()
        for word in product((0, 1), repeat=L):
            P, B, Q, c = cylinder(word)
            require(c not in seen, 'distinct parity cylinders')
            seen.add(c)
            require(P == 3 ** sum(word) and Q == 2 ** L, 'affine slope')
            for lift in (1, 2, 3):
                n = c + lift * Q
                require(bits(n, L) == word, 'cylinder positive lift')
                end = n
                for _ in range(L):
                    end = step(end)
                require(Q * end == P * n + B, 'affine endpoint')
            digits = tuple(1 + bit for bit in word)
            a, b, cc, d = cf_matrix(digits)
            require(a * d - b * cc == (-1) ** L, 'unimodular CF prefix')
            require(convergent(digits) == 1 / (digits[0] + convergent(digits[1:])),
                    'finite Gauss deletion')
        require(seen == set(range(2 ** L)), 'complete residue partition')
        by_length.append({'length': L, 'cylinders': len(seen)})
    # Shortcut cycles of length dividing two: parity algebra gives the proof;
    # this census is a bounded consistency check, not an all-height argument.
    period_two = [n for n in range(1, 10001) if step(step(n)) == n]
    require(period_two == [1, 2], 'two-cycle finite control')
    require(step(-1) == -1, 'signed hostile fixed point')
    for L in range(1, 65):
        n = 2 ** L - 1
        require(bits(n, L) == (1,) * L, 'positive shadow of signed fixed point')
        require(cylinder((1,) * L)[3] == n, 'Mersenne cylinder')
    gate_rows = []
    for n in range(1, 2049):
        digits, rows, hit, switch = gated_code(n)
        require(hit == hitting_time(n), 'authenticated gate time')
        require(switch is not None, 'bounded source control switched')
        for j, q, digit, found in rows:
            require(digit == (1 if found is not None else q), 'actual gate rule')
        if n in (1, 2, 3, 27, 97, 703):
            gate_rows.append({'n': n, 'hit_time': hit, 'switch_stage': switch,
                             'digits': digits})
    for exponent in (3, 4, 5):
        # Abstract never-hit control, not a claimed positive Collatz orbit.
        q0, q1 = 1, 2
        for j in range(1, 10):
            a = q1 ** (exponent - 2)
            q2 = a * q1 + q0
            require(q2 > q1 ** (exponent - 1), 'unresolved spike lower growth')
            require(q2 < 2 * q1 ** (exponent - 1), 'unresolved spike upper growth')
            q0, q1 = q1, q2
    # A field label might retain more than mu, but finite absence is not proof.
    field_census = []
    for length in range(1, 17):
        matches = 0
        nonalternating = 0
        for word in product((1, 2), repeat=length):
            a, b, c, d = cf_matrix(word)
            disc = (a + d) ** 2 - 4 * ((-1) ** length)
            if disc % 3 == 0 and isqrt(disc // 3) ** 2 == disc // 3:
                matches += 1
                nonalternating += int(any(word[j] == word[(j + 1) % length]
                                          for j in range(length)))
        field_census.append({'length': length, 'sqrt3_words': matches,
                             'nonalternating': nonalternating})
    # The apparent field pattern through16 fails at18. Retain the source guard.
    field_word = tuple(map(int, '111111121121211122'))
    field_bits = tuple(a - 1 for a in field_word)
    matrix = cf_matrix(field_word)
    require(matrix == (29872, 12553, 18471, 7762), 'field hostile matrix')
    a, b, c, d = matrix
    require(a * d - b * c == 1, 'field hostile determinant')
    require((a + d) ** 2 - 4 == 3 * 21728 ** 2, 'field hostile discriminant')
    require(any(field_word[j] == field_word[(j + 1) % 18] for j in range(18)),
            'field hostile not alternating')
    P, B, Q, residue = cylinder(field_bits)
    source = Fraction(B, Q - P)
    require((P, Q, B) == (243, 262144, 402560), 'field hostile parity carrier')
    require(source == Fraction(402560, 261901), 'field hostile rational source')
    require(source.denominator > 1, 'integer guard rejects field hostile')
    orbit = []
    state = source
    for bit in field_bits:
        orbit.append(state)
        require(state > 0 and state.denominator % 2 == 1, 'positive2-adic state')
        require(state.numerator % 2 == bit, 'rational parity guard')
        state = state / 2 if bit == 0 else (3 * state + 1) / 2
    require(state == source and len(set(orbit)) == 18, 'primitive rational cycle')
    return {'status': 'FINITE-EXACT; infinite statements proved in note',
            'checks': CHECKS, 'parity_partition': by_length,
            'two_cycle_census_1_to_10000': period_two,
            'stopping_gate_universe': 'all sources1..2048,14 CF digits, horizon4096',
            'gate_samples': gate_rows,
            'quadratic_field_exploration_only': field_census,
            'quadratic_field_hostile': {'word': ''.join(map(str, field_word)),
                'matrix': matrix, 'source': str(source),
                'rational_orbit': [str(x) for x in orbit]}}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
