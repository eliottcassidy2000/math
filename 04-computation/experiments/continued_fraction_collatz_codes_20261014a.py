"""Exact CF/state codecs and local-growth controls; no convergence oracle.

The infinite irrationality-exponent statements are proved in the companion
note. These finite checks authenticate words, guards and integer inequalities.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from math import factorial, gcd
import json

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def integer(n, minimum=1):
    need(type(n) is int and n >= minimum, 'exact integer in domain')
    return n


def word(values, empty=True):
    need(type(values) is tuple and (empty or len(values) > 0), 'literal tuple')
    for value in values:
        integer(value)
    return values


def valuation(n):
    integer(n)
    return (n & -n).bit_length()-1


def step(n, clock='standard'):
    integer(n)
    need(type(clock) is str and clock in ('standard', 'shortcut', 'odd'), 'clock')
    if clock == 'odd':
        need(n % 2 == 1, 'odd source')
        if n == 1:
            return 1  # declared infinite-code convention only
        value = 3*n+1
        return value >> valuation(value)
    if n % 2 == 0:
        return n//2
    return 3*n+1 if clock == 'standard' else (3*n+1)//2


def orbit_digits(source, count, clock='standard', marked=False):
    integer(source), integer(count)
    need(type(marked) is bool, 'exact marking flag')
    current = source
    digits = []
    for index in range(1, count+1):
        # Validate the clock/domain even when only one digit was requested.
        successor = step(current, clock)
        digits.append(index*current+1 if marked else current)
        current = successor
    return tuple(digits)


def decode_marked(digits, clock='standard'):
    word(digits, empty=False)
    states = []
    for index, value in enumerate(digits, 1):
        need((value-1) % index == 0, 'retained index divisibility')
        state = (value-1)//index
        integer(state)
        successor = step(state, clock)
        if states:
            need(state == step(states[-1], clock), 'actual consecutive states')
        states.append(state)
    return tuple(states)


def continuant(digits):
    """(p_n,q_n,p_(n-1),q_(n-1)) for [0;digits], including empty prefix."""
    word(digits)
    p_old, p, q_old, q = 1, 0, 0, 1
    for a in digits:
        p_old, p = p, a*p+p_old
        q_old, q = q, a*q+q_old
    return p, q, p_old, q_old


def cylinder(digits):
    word(digits, empty=False)
    p, q, r, s = continuant(digits)
    return tuple(sorted((Fraction(p, q), Fraction(p+r, q+s))))


def growth_bill(digits, bound, width):
    """Validate a finite local-growth word and the last-window theorem."""
    word(digits, empty=False), integer(bound), integer(width)
    need(len(digits) >= width+1, 'enough past digits')
    need(all(b <= bound*a for a, b in zip(digits, digits[1:])), 'local ratio')
    q = continuant(digits[:-1])[1]
    left = digits[-1]**width
    right = bound**(width*(width+1)//2)*q
    need(left <= right, 'proved last-window inequality')
    return left, right


def exponent_gate(digits, bound, denominator):
    """If the exact gate holds, a_next^d <= q; false is inconclusive."""
    integer(denominator)
    width = 2*denominator
    growth_bill(digits, bound, width)
    q = continuant(digits[:-1])[1]
    if q < bound**(denominator*(2*denominator+1)):
        return False
    need(digits[-1]**denominator <= q, 'exact exponent gate consequence')
    return True


def native_cell(valuations):
    """Unique odd cell plus a strict no-earlier-ROOT lower cut."""
    word(valuations)
    p, q, b, lower = 1, 1, 0, 1
    for a in valuations:
        lower = max(lower, (q-b)//p+1)
        p, q, b = 3*p, q*(1 << a), 3*b+q
    modulus = 2*q
    residue = ((q-b)*pow(p, -1, modulus)) % modulus
    least = residue+max(0, (lower-residue+modulus-1)//modulus)*modulus
    return residue, modulus, least


@dataclass(frozen=True)
class OddCode:
    source: int
    valuations: tuple
    states: tuple
    digits: tuple


def odd_code(source, valuations):
    integer(source), word(valuations)
    need(source % 2 == 1, 'odd source')
    states, current = [source], source
    for a in valuations:
        need(current > 1, 'strict first-hit ROOT boundary')
        raw = 3*current+1
        need(valuation(raw) == a, 'exact native valuation')
        current = raw >> a
        states.append(current)
    values = tuple(states)
    digits = tuple(index*x+1 for index, x in enumerate(values, 1))
    return OddCode(source, valuations, values, digits)


def audit_code(receipt):
    need(type(receipt) is OddCode, 'literal OddCode')
    integer(receipt.source)
    word(receipt.valuations), word(receipt.states, empty=False), word(receipt.digits, empty=False)
    need(odd_code(receipt.source, receipt.valuations) == receipt, 'source-bound receipt')
    need(decode_marked(receipt.digits, 'odd') == receipt.states, 'lossless state decoder')
    return receipt


def discharge(receipt, supplied_child_word):
    audit_code(receipt), word(supplied_child_word)
    suffix = odd_code(receipt.states[-1], supplied_child_word)
    need(suffix.states[-1] == 1, 'supplied child ROOT required')
    return audit_code(odd_code(receipt.source, receipt.valuations+supplied_child_word))


def reject(call, label):
    try:
        call()
    except (ValueError, TypeError):
        check(True, label)
    else:
        check(False, label)


def main():
    model_controls = {}
    for clock, bound in (('standard', 4), ('shortcut', 2), ('odd', 2)):
        sources = range(1, 129) if clock != 'odd' else range(1, 256, 2)
        total = 0
        for source in sources:
            raw = orbit_digits(source, 65, clock)
            marked = orbit_digits(source, 65, clock, True)
            check(decode_marked(marked, clock) == raw, 'exact marked recovery')
            for digits, c in ((raw, bound), (marked, 2*bound)):
                check(all(b <= c*a for a, b in zip(digits, digits[1:])), 'all local ratios')
                for depth in (8, 16, 32, 64):
                    for width in (1, 2, 4, 8):
                        left, right = growth_bill(digits[:depth+1], c, width)
                        check(left <= right, 'window product versus continuant')
                for d in (1, 2, 3):
                    gate = exponent_gate(digits, c, d)
                    if gate:
                        check(digits[-1]**d <= continuant(digits[:-1])[1], 'quantified digit gate')
            total += 1
        model_controls[clock] = {'sources': total, 'digits_each': 65, 'raw_ratio_bound': bound}

    native_words = 0
    for length in range(4):
        for values in product(range(1, 5), repeat=length):
            residue, modulus, least = native_cell(values)
            check(least > 0 and least % 2 == 1 and least % modulus == residue, 'native cell source')
            for lift in range(8):
                code = audit_code(odd_code(least+modulus*lift, values))
                p, q, pp, qq = continuant(code.digits)
                check(abs(p*qq-pp*q) == 1 and gcd(p, q) == 1, 'unimodular continuant')
                lo, hi = cylinder(code.digits)
                longer = continuant(code.digits+(2, 3))
                check(lo < Fraction(longer[0], longer[1]) < hi, 'infinite-tail cylinder')
            native_words += 1
    control = odd_code(7, (1, 1))
    root = discharge(control, (2, 3, 4))
    check(root.states == (7, 11, 17, 13, 5, 1), 'supplied ROOT word transport')
    check(discharge(odd_code(1, ()), ()).states == (1,), 'empty actual ROOT seed')
    check(native_cell((2,))[2] == 9, 'formal ROOT padding excluded by lower cut')

    # Exact continuant bounds prove the exponents; finite ratios do not estimate them.
    engineering = []
    for power in (2, 3, 4):
        count = 7
        exponents = tuple(power**j for j in range(count))
        values = tuple(1 << e for e in exponents)
        q = continuant(values)[1]
        s = sum(exponents)
        check((1 << s) <= q <= (1 << (s+count)), 'geometric exponent product bounds')
        engineering.append({'digit_power': power, 'proved_mu': power+1,
                            'next_log_ratio_lower': str(Fraction(power**count, s+count)),
                            'next_log_ratio_upper': str(Fraction(power**count, s))})
    count = 7
    exponents = tuple(factorial(j) for j in range(1, count+1))
    q = continuant(tuple(1 << e for e in exponents))[1]
    check((1 << sum(exponents)) <= q <= (1 << (sum(exponents)+count)), 'factorial denominator bounds')
    check(sum(exponents)+count <= 3*factorial(count), 'factorial ratio diverges at least linearly')
    # Exponential envelope alone is weaker than a local ratio bound.
    sparse = tuple((1 << k) if k & (k-1) == 0 else 1 for k in range(1, 513))
    q = continuant(sparse[:-1])[1]
    check(all(a <= (1 << k) for k, a in enumerate(sparse, 1)), 'global exponential envelope')
    check(sparse[-1] > 4*sparse[-2], 'sparse spike violates local growth')
    check(sparse[-1]**2 > q, 'positive subsequence ratio despite exponential envelope')
    growing = tuple(1 << k for k in range(65))
    check(growth_bill(growing, 2, 8)[0] <= growth_bill(growing, 2, 8)[1], 'divergent synthetic orbit has mu2')

    # Literal tuples retain the clock; a terminating rational alone does not.
    check(Fraction(*continuant((2, 1))[:2]) == Fraction(*continuant((3,))[:2]), 'finite CF last-one ambiguity')
    check(orbit_digits(1, 6, 'standard') == (1, 4, 2, 1, 4, 2), 'standard ROOT cycle')
    check(orbit_digits(1, 6, 'shortcut') == (1, 2, 1, 2, 1, 2), 'shortcut ROOT cycle')
    check(orbit_digits(1, 6, 'odd') == (1,)*6, 'declared absorbed odd code')
    hostile = [lambda: orbit_digits(True, 3), lambda: orbit_digits(3, 3.0),
               lambda: orbit_digits(2, 3, 'odd'), lambda: orbit_digits(3, 3, 'wrong'),
               lambda: orbit_digits(3, 3, marked=1), lambda: decode_marked((4, 99)),
               lambda: continuant((2, True)), lambda: native_cell((1.0,)),
               lambda: odd_code(1, (2,)), lambda: odd_code(7, (2,)),
               lambda: audit_code(replace(control, source=7.0)),
               lambda: audit_code(replace(control, states=(7, 11, 17.0))),
               lambda: audit_code(replace(control, digits=(8, 23, 53))),
               lambda: discharge(control, (4,)), lambda: discharge(control, ()),
               lambda: growth_bill((1, 10), 4, 1),
               lambda: exponent_gate((1,)*9, 2, True)]
    for index, call in enumerate(hostile):
        reject(call, 'hostile '+str(index))

    report = {'status': 'PROVED asymptotic statements in companion note; finite exact controls here',
              'models': model_controls, 'native_words': native_words, 'native_lifts_each': 8,
              'engineering_controls': engineering, 'factorial_proved_mu': 'infinity',
              'source_7_receipt': {'states': root.states, 'valuations': root.valuations, 'digits': root.digits},
              'sparse_envelope': 'a_k<=2^k does not imply mu2; this example has mu>=5/2',
              'hostiles': len(hostile), 'checks': CHECKS,
              'scope': 'No finite data are labelled an irrationality-exponent estimate; no target ROOT oracle.'}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
