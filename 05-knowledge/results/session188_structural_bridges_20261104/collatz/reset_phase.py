#!/usr/bin/env python3
"""Exact Collatz one-reset phase compiler. Standard library only.

Production certificates use modular arithmetic and explicit checks. Expanding
an actual source may cost exponentially many bits; the API keeps it separate.
"""
from itertools import product
from fractions import Fraction
import json


def require(ok, message):
    if not ok:
        raise ValueError(message)


def validate_word(word):
    require(len(word) > 0 and all(isinstance(a, int) and a >= 1 for a in word),
            'word must be nonempty with positive integer letters')


def word_carry(word, modulus=None):
    validate_word(word)
    A = B = 0
    for a in word:
        B = 3 * B + (2 ** A if modulus is None else pow(2, A, modulus))
        if modulus is not None:
            B %= modulus
        A += a
    return A, B


def log4_ternary(unit, depth):
    """Return j mod 3**depth with 4**j == unit mod 3**(depth+1).

    Each stage tests precisely three lifts, never an exponential interval.
    """
    require(isinstance(depth, int) and depth >= 0, 'negative logarithm depth')
    modulus = 3 ** (depth + 1)
    unit %= modulus
    require(unit % 3 == 1, 'target must be 1 modulo 3')
    j, period, mod = 0, 1, 3
    for _ in range(depth):
        mod *= 3
        candidates = [j + digit * period for digit in range(3)
                      if pow(4, j + digit * period, mod) == unit % mod]
        require(len(candidates) == 1, 'ternary lift is not unique')
        j = candidates[0]
        period *= 3
    return j


def reset_phase(word):
    """Least even beta in [0, 2*3**L) for word followed by one reset to 1."""
    validate_word(word)
    L = len(word)
    modulus = 3 ** (L + 1)
    A, B = word_carry(word, modulus)
    target = (1 + 3 * B * pow(pow(2, A, modulus), -1, modulus)) % modulus
    beta = 2 * log4_ternary(target, L)
    require(beta % 6 in (2, 4), 'invalid reset phase')
    return beta


def completion(word, expand=False):
    """Return a compact certificate for the least first-hit reset exponent.

    The source is exactly (2**(A+b) - 2**A - 3*B) / 3**(L+1).
    No source is expanded unless expand=True.
    """
    validate_word(word)
    L = len(word)
    modulus = 3 ** (L + 1)
    A, B_residue = word_carry(word, modulus)
    beta, period = reset_phase(word), 2 * 3 ** L
    b = beta if beta >= 4 else beta + period
    cert = dict(word=list(word), length=L, cost=A,
                carry_residue=B_residue, carry_modulus=modulus,
                source_formula='(2**(A+b) - 2**A - 3*B) / 3**(L+1)',
                carry_formula='sum(3**(L-1-i) * 2**sum(word[:i]) for i in range(L))',
                phase=beta, period=period, first_hit_reset=b)
    if expand:
        _, B = word_carry(word)
        cert['carry'] = B
        numerator = 2 ** (A + b) - 2 ** A - 3 * B
        denominator = 3 ** (L + 1)
        require(numerator % denominator == 0, 'nonintegral source')
        n = numerator // denominator
        require(n > 1 and n % 2 == 1, 'bad source')
        cert['source'] = n
    return cert


def prefix_for_phase(length, beta):
    """Compile any allowable beta to a word over {1,2,3,4}.

    Only residues of at most length+1 ternary digits are constructed.
    The reversed list is the actual forward prefix. No integer 2**beta
    is materialized.
    """
    require(isinstance(length, int) and length >= 1, 'length must be positive')
    period = 2 * 3 ** length
    require(isinstance(beta, int) and 0 <= beta < period and beta % 6 in (2, 4),
            'phase must lie in range and equal 2 or 4 modulo 6')
    modulus = 3 ** length
    residue = (pow(2, beta, 3 * modulus) - 1) // 3
    reverse_word = []
    for remaining in range(length, 0, -1):
        require(residue % 3 in (1, 2), 'dead inverse residue')
        low = 2 if residue % 3 == 1 else 1
        choices = (low, low + 2) if remaining > 1 else (low,)
        selected = None
        for a in choices:
            value = ((2 ** a) * residue - 1) // 3
            if remaining == 1 or value % 3 != 0:
                selected = (a, value)
                break
        require(selected is not None, 'four-letter controller got stuck')
        a, residue = selected
        reverse_word.append(a)
        modulus //= 3
        if modulus > 1:
            residue %= modulus
    return tuple(reversed(reverse_word))


def periodic_phase(word, repetitions):
    """Compile the phase for word**repetitions from one rational center."""
    validate_word(word)
    require(isinstance(repetitions, int) and repetitions >= 1, 'bad repetition count')
    L = len(word)
    modulus = 3 ** (repetitions * L + 1)
    A, B = word_carry(word, modulus)
    D = (pow(2, A, modulus) - pow(3, L, modulus)) % modulus
    target = ((D + 3 * B) * pow(D, -1, modulus)) % modulus
    return 2 * log4_ternary(target, repetitions * L)


def actual_step(n):
    require(n > 0 and n % 2 == 1, 'source must be a positive odd integer')
    value = 3 * n + 1
    a = (value & -value).bit_length() - 1
    return value >> a, a


def direct_replay(n, expected):
    states = [n]
    for a in expected:
        n, actual = actual_step(n)
        require(a == actual, 'actual valuation differs from certificate')
        states.append(n)
    require(n == 1 and all(x > 1 for x in states[:-1]), 'not a first-hit root route')
    return states


def run_tests():
    checks = 0
    words = sources = 0
    phase_counts = []
    # Exhaustive words, separate direct arithmetic oracle, and brute log tables.
    for L in range(1, 7):
        M, P = 3 ** (L + 1), 2 * 3 ** L
        brute_logs = {pow(2, b, M): b for b in range(P)}
        phases = set()
        for word in product(range(1, 5), repeat=L):
            cert = completion(word, expand=True)
            A, B, beta, b = (cert[k] for k in ('cost', 'carry', 'phase', 'first_hit_reset'))
            target = (1 + 3 * B * pow(pow(2, A, M), -1, M)) % M
            require(beta == brute_logs[target], 'brute logarithm mismatch')
            require(b >= 4 and b - P < 4, 'first-hit exponent is not minimal')
            direct_replay(cert['source'], list(word) + [b])
            # A second member of every exponent progression, independently replayed.
            b2 = b + P
            n2 = (2 ** (A + b2) - 2 ** A - 3 * B) // 3 ** (L + 1)
            direct_replay(n2, list(word) + [b2])
            D = 2 ** A - 3 ** L
            C = Fraction(D + 3 * B, D)
            numerator_abs = abs(C.numerator)
            is_signed_power_two = (C.denominator == 1 and numerator_abs > 0
                                   and numerator_abs & (numerator_abs - 1) == 0)
            require(is_signed_power_two == (all(a == 1 for a in word) or all(a == 2 for a in word)),
                    'rational-address classification failed')
            checks += 1
            phases.add(beta)
            words += 1
            sources += 2
            checks += 3 + 2 * (L + 1)
        expected = {b for b in range(2, P, 2) if b % 3 != 0}
        require(phases == expected, 'four-letter phase spectrum incomplete')
        phase_counts.append(dict(length=L, phases=len(phases), words=4 ** L))
        checks += 1
    # Full converse controller for every allowed phase through depth ten.
    compiled_phases = 0
    for L in range(1, 11):
        for beta in range(2, 2 * 3 ** L, 2):
            if beta % 3 == 0:
                continue
            word = prefix_for_phase(L, beta)
            require(len(word) == L and all(1 <= a <= 4 for a in word), 'bad bounded prefix')
            require(reset_phase(word) == beta, 'converse phase mismatch')
            compiled_phases += 1
            checks += 2
    # Alphabet three obstruction, exact exhaustive depth two.
    phases3 = {reset_phase(w) for w in product(range(1, 4), repeat=2)}
    require(phases3 == {2, 4, 8, 10, 16}, 'alphabet-three obstruction changed')
    m = (2 ** 14 - 1) // 3
    require(m % 9 == 7 and ((4 * m - 1) // 3) % 3 == 0, 'obstruction arithmetic failed')
    checks += 2
    # Analytic extremal formulas, without large output integers.
    for L in range(1, 61):
        P = 2 * 3 ** L
        require(reset_phase((1,) * (L - 1) + (3,)) == P - 2, 'extremal phase failed')
        require(reset_phase((1,) * L) == 3 ** L + 1, 'all-one formula failed')
        require(completion((2,) * L)['first_hit_reset'] == P + 2, 'root phase formula failed')
        checks += 3
    periodic_tests = 0
    for word in [(1,), (2,), (1, 2), (1, 1, 3), (2, 1, 1, 2), (3, 1, 2, 4)]:
        A, B = word_carry(word)
        L, D = len(word), 2 ** A - 3 ** len(word)
        previous = None
        for r in range(1, 26):
            Ar, Br = word_carry(word * r)
            beta = periodic_phase(word, r)
            require(beta == reset_phase(word * r), 'center compiler mismatch')
            require(Fraction(D + 3 * B, D) - Fraction(2 ** Ar + 3 * Br, 2 ** Ar)
                    == Fraction(3 ** (r * L + 1) * B, D * 2 ** Ar),
                    'center approximation identity failed')
            if previous is not None:
                old_beta, old_P = previous
                require(beta % old_P == old_beta and beta >= old_beta, 'phase coherence failed')
                require(beta == old_beta or beta >= old_P, 'new ternary block is too small')
                checks += 1
            previous = beta, 2 * 3 ** (r * L)
            periodic_tests += 1
            checks += 3
    # Exact S179 shadow precision for expanded periodic examples.
    shadow_tests = 0
    for word in [(1,), (2,), (1, 2), (1, 1, 3), (3, 1)]:
        A, B = word_carry(word)
        D = 2 ** A - 3 ** len(word)
        for r in range(1, 4):
            cert = completion(word * r)
            if cert['cost'] + cert['first_hit_reset'] > 20000:
                continue
            cert = completion(word * r, expand=True)
            n, b = cert['source'], cert['first_hit_reset']
            z = D * n - B
            v = (abs(z) & -abs(z)).bit_length() - 1
            require(b != word[0], 'test requires mismatch at reset')
            require(v == r * A + min(b, word[0]), 'template mismatch precision failed')
            direct_replay(n, list(word) * r + [b])
            shadow_tests += 1
            checks += 2
    # Huge binary-encoded letters must not expand the carry or source by default.
    compact = completion((10 ** 20, 1, 3))
    require('source' not in compact and 'carry' not in compact, 'compact path expanded large integers')
    require(compact['carry_modulus'] == 81 and 0 <= compact['carry_residue'] < 81,
            'compact carry is not bounded by ternary precision')
    require(periodic_phase((10 ** 20, 1, 3), 5) == reset_phase((10 ** 20, 1, 3) * 5),
            'large-letter periodic compiler mismatch')
    checks += 3
    # Rejection controls use checks rather than optimization-removable asserts.
    rejected = 0
    invalid = [lambda: reset_phase(()), lambda: reset_phase((0,)),
               lambda: prefix_for_phase(0, 2), lambda: prefix_for_phase(2, 3),
               lambda: prefix_for_phase(2, 6), lambda: prefix_for_phase(2, 18),
               lambda: log4_ternary(2, 3), lambda: periodic_phase((1,), 0)]
    for f in invalid:
        try:
            f()
        except ValueError:
            rejected += 1
        else:
            raise RuntimeError('invalid input accepted')
    checks += rejected
    examples = []
    for w in [(1,), (1, 1), (1, 3), (2, 2), (1, 1, 3)]:
        c = completion(w, expand=True)
        examples.append(dict(word=list(w), phase=c['phase'], reset=c['first_hit_reset'],
                             source=c['source'], states=direct_replay(c['source'], list(w) + [c['first_hit_reset']])))
    return dict(status='PASS', exhaustive_words=words, literal_root_routes=sources,
                inverse_phase_compilations=compiled_phases, periodic_center_tests=periodic_tests,
                template_precision_tests=shadow_tests, rejected_controls=rejected,
                explicit_check_count=checks, phase_spectra=phase_counts, examples=examples)


if __name__ == '__main__':
    print(json.dumps(run_tests(), indent=2, sort_keys=True))
