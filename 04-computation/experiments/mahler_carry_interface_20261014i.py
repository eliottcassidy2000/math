"""Source-bound finite separation of Mahler and Collatz carry codes.

No Z-number or Collatz convergence oracle. Stdout is deterministic JSON.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from math import gcd
import json


def integer(n, least=None):
    if type(n) is not int or (least is not None and n < least):
        raise ValueError('exact integer outside domain')
    return n


def rational(x):
    if type(x) not in (int, Fraction):
        raise ValueError('exact rational required')
    x = Fraction(x)
    if x.denominator % 2 == 0:
        raise ValueError('odd denominator required')
    return x


def word(w):
    if type(w) is not tuple:
        raise ValueError('binary tuple required')
    for e in w:
        if type(e) is not int or e not in (0, 1):
            raise ValueError('exact binary digits required')
    return w


def model(name):
    if type(name) is not str or name not in ('mahler', 'collatz'):
        raise ValueError('known map required')
    return name


def parity(x):
    return rational(x).numerator % 2


def step(x, name):
    x, name = rational(x), model(name)
    e = parity(x)
    return ((3 if name == 'mahler' or e else 1) * x + e) / 2


def prefix(x, length, name):
    x, length, name = rational(x), integer(length, 0), model(name)
    out = []
    for _ in range(length):
        out.append(parity(x))
        x = step(x, name)
    return tuple(out), x


def carrier(w, name):
    w, name = word(w), model(name)
    p, q, c = 1, 1, 0
    for e in w:
        mult = 3 if name == 'mahler' or e else 1
        p, c, q = mult * p, mult * c + e * q, 2 * q
    return p, q, c


def native_cell(w, name):
    p, q, c = carrier(w, name)
    return (-c * pow(p, -1, q)) % q, q


def periodic_source(w, name):
    w = word(w)
    if not w:
        raise ValueError('nonempty period required')
    p, q, c = carrier(w, name)
    return Fraction(c, q - p)


def periodic_mahler_tails(w):
    w = word(w)
    if not w:
        raise ValueError('nonempty period required')
    return tuple(-2 * periodic_source(w[i:] + w[:i], 'mahler')
                 for i in range(len(w)))


def v2(n):
    integer(n)
    if n == 0:
        raise ValueError('nonzero integer required')
    n = abs(n)
    return (n & -n).bit_length() - 1


@dataclass(frozen=True)
class Transport:
    source: int
    length: int
    code: tuple
    mahler_residue: int
    collatz_residue: int
    modulus: int


def transport(source, length):
    integer(source, 1)
    integer(length, 0)
    w, _ = prefix(source, length, 'mahler')
    a, q = native_cell(w, 'mahler')
    b, _ = native_cell(w, 'collatz')
    return Transport(source, length, w, a, b, q)


def audit_transport(packet, source):
    integer(source, 1)
    if type(packet) is not Transport:
        raise ValueError('transport packet required')
    integer(packet.source, 1)
    integer(packet.length, 0)
    word(packet.code)
    for n in (packet.mahler_residue, packet.collatz_residue, packet.modulus):
        integer(n, 0)
    if packet != transport(source, packet.length):
        raise ValueError('transport is not bound to supplied source')
    return packet


def lattice_bound(x):
    x = rational(x)
    u, v = x.numerator, x.denominator
    b = max(abs(u), 2 * v)
    k = (3 * b + v) // 2
    return b, k, 2 * k + 1


@dataclass(frozen=True)
class RationalObserver:
    source: Fraction
    kind: str
    index: int
    start: int
    code: tuple
    states: tuple
    state_bound: int


def observe_rational(x):
    """Stop at the second bit of 11 or a repeated rational Collatz state."""
    x = rational(x)
    _, k, bound = lattice_bound(x)
    current, seen, states, bits = x, {}, [], []
    for j in range(bound + 1):
        e = parity(current)
        if bits and bits[-1] == e == 1:
            return RationalObserver(x, '11', j, -1, tuple(bits + [e]),
                                    tuple(states + [current]), bound)
        if current in seen:
            return RationalObserver(x, 'cycle', j, seen[current], tuple(bits),
                                    tuple(states + [current]), bound)
        if abs(current * x.denominator) > k:
            raise ArithmeticError('finite no-11 lattice bound failed')
        seen[current] = j
        states.append(current)
        bits.append(e)
        current = step(current, 'collatz')
    raise ArithmeticError('finite lattice observer failed')


def audit_observer(packet, source):
    x = rational(source)
    if type(packet) is not RationalObserver or type(packet.source) is not Fraction:
        raise ValueError('exact observer required')
    rational(packet.source)
    for n in (packet.index, packet.start, packet.state_bound):
        integer(n)
    if type(packet.kind) is not str or type(packet.states) is not tuple:
        raise ValueError('exact observer fields required')
    for z in packet.states:
        if type(z) is not Fraction:
            raise ValueError('exact state required')
        rational(z)
    word(packet.code)
    if packet != observe_rational(x):
        raise ValueError('wrong observer or source')
    return packet


def separation_bound(source, x):
    integer(source, 1)
    return 2 * lattice_bound(x)[2] + (source + 1).bit_length()


def height_bound(source, height):
    integer(source, 1)
    integer(height, 1)
    return 14 * height + 2 + (source + 1).bit_length()


@dataclass(frozen=True)
class Separation:
    source: int
    rational_source: Fraction
    bound: int
    kind: str
    index: int
    mahler_word: tuple
    collatz_word: tuple


def compare(source, x):
    """Total finite test: first mismatch or first shared forbidden 11."""
    integer(source, 1)
    x = rational(x)
    bound = separation_bound(source, x)
    a, z, wm, wt = Fraction(source), x, [], []
    for i in range(bound):
        wm.append(parity(a))
        wt.append(parity(z))
        kind = 'mismatch' if wm[-1] != wt[-1] else (
            'shared11' if i and wm[-2:] == [1, 1] else None)
        if kind:
            return Separation(source, x, bound, kind, i, tuple(wm), tuple(wt))
        a, z = step(a, 'mahler'), step(z, 'collatz')
    raise ArithmeticError('proved separation bound failed')


def audit_separation(packet, source, x):
    """Check supplied literal words; no unbounded orbit search."""
    integer(source, 1)
    x = rational(x)
    if type(packet) is not Separation or type(packet.rational_source) is not Fraction:
        raise ValueError('exact separation packet required')
    for n in (packet.source, packet.bound, packet.index):
        integer(n, 0)
    rational(packet.rational_source)
    if type(packet.kind) is not str:
        raise ValueError('exact event kind required')
    word(packet.mahler_word)
    word(packet.collatz_word)
    if (packet.source != source or packet.rational_source != x
            or packet.bound != separation_bound(source, x)
            or not 0 <= packet.index < packet.bound
            or len(packet.mahler_word) != packet.index + 1
            or len(packet.collatz_word) != packet.index + 1):
        raise ValueError('source, horizon, or shape mismatch')
    wm, wt = packet.mahler_word, packet.collatz_word
    if wm != prefix(source, len(wm), 'mahler')[0] or wt != prefix(x, len(wt), 'collatz')[0]:
        raise ValueError('non-native supplied word')
    for i in range(packet.index):
        if wm[i] != wt[i] or (i and wm[i - 1] == wm[i] == 1):
            raise ValueError('not the first event')
    expected = 'mismatch' if wm[-1] != wt[-1] else (
        'shared11' if packet.index and wm[-2:] == (1, 1) else None)
    if packet.kind != expected or expected is None:
        raise ValueError('false separation event')
    return packet


def cycle_mismatch_bound(source, observer):
    """If matching reaches a rational cycle, its next mismatch is exact."""
    integer(source, 1)
    audit_observer(observer, observer.source)
    if observer.kind != 'cycle':
        raise ValueError('cycle observer required')
    i, j = observer.start, observer.index
    if prefix(source, i, 'mahler')[0] != observer.code[:i]:
        raise ValueError('preperiod already mismatches')
    period = observer.code[i:j]
    a = prefix(source, i, 'mahler')[1]
    p, q, c = carrier(period, 'mahler')
    d = p - q
    numerator = d * a + c
    if numerator.denominator != 1 or numerator <= 0:
        raise ArithmeticError('positive integer numerator required')
    h = v2(numerator.numerator)
    return dict(start=i, period=list(period), numerator=numerator.numerator,
                mismatch_index=i + h, bound=separation_bound(source, observer.source))


def delayed_family(index):
    integer(index, 0)
    length = 18 * index + 8
    source, rem = divmod((1 << length) - 9, 19)
    if rem or source <= 0:
        raise ArithmeticError('denominator-19 phase failed')
    return source, length


def main():
    checks = 0

    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise AssertionError(label)

    def rejects(f):
        try:
            f()
        except ValueError:
            check(True, 'rejected')
        else:
            check(False, 'hostile accepted')

    cells = 0
    for length in range(11):
        native = {name: set() for name in ('mahler', 'collatz')}
        for w in product((0, 1), repeat=length):
            for name in native:
                r, q = native_cell(w, name)
                native[name].add(r)
                p, _, c = carrier(w, name)
                for lift in (1, 2, 3):
                    source = r + q * lift
                    got, end = prefix(source, length, name)
                    check(got == w and end * q == p * source + c, 'native cylinder and full carrier')
                cells += 1
        check(all(len(z) == 1 << length for z in native.values()), 'both native cell bijections')

    periodic = 0
    for length in range(1, 9):
        for w in product((0, 1), repeat=length):
            for name in ('mahler', 'collatz'):
                z = periodic_source(w, name)
                got, end = prefix(z, length, name)
                check(got == w and end == z, 'periodic rational source')
            p, q, c = carrier(w, 'mahler')
            check(0 <= c <= p - q, 'Mahler fixed point interval')
            periodic += 1

    rationals = sorted({Fraction(u, v) for v in range(1, 16, 2)
                        for u in range(-15, 16) if gcd(u, v) == 1})
    observers = {x: observe_rational(x) for x in rationals}
    cycle_count = sum(z.kind == 'cycle' for z in observers.values())
    cycle_claims = 0
    events = {'mismatch': 0, 'shared11': 0}
    largest_index = -1
    largest_example = None
    for x, observation in observers.items():
        check(audit_observer(observation, x) == observation, 'observer audit')
        if observation.kind == 'cycle':
            period = observation.code[observation.start:observation.index]
            doubled = period + period
            check(all(doubled[j:j+2] != (1, 1) for j in range(len(doubled)-1)), 'including cyclic seam')
            check(prefix(x, 5 * observation.state_bound, 'collatz')[1].denominator <= x.denominator,
                  'fixed denominator lattice')
        for source in range(1, 129):
            packet = compare(source, x)
            check(audit_separation(packet, source, x) == packet, 'literal separation audit')
            check(packet.bound <= height_bound(source, max(abs(x.numerator), x.denominator)),
                  'uniform height bound')
            events[packet.kind] += 1
            if packet.index > largest_index:
                largest_index = packet.index
                largest_example = dict(A=source, x=str(x), kind=packet.kind, index=packet.index,
                                       bound=packet.bound)
            i = observation.start
            if observation.kind == 'cycle' and prefix(source, i, 'mahler')[0] == observation.code[:i]:
                claim = cycle_mismatch_bound(source, observation)
                h = claim['mismatch_index']
                wm, _ = prefix(source, h + 1, 'mahler')
                wt, _ = prefix(x, h + 1, 'collatz')
                check(wm[:h] == wt[:h] and wm[h] != wt[h] and h < claim['bound'],
                      'independent exact periodic mismatch index')
                cycle_claims += 1

    delayed = []
    for index in range(6):
        a, length = delayed_family(index)
        p = compare(a, Fraction(1, 5))
        check(p.kind == 'mismatch' and p.index == length, 'arbitrarily delayed rational comparison')
        check(v2(19 * a + 9) == length, 'independent exact denominator phase')
        check(p.mahler_word[:length] == ((1, 0, 0) * ((length + 2)//3))[:length], 'safe periodic prefix')
        for j in range(length):
            tail = sum((Fraction(2, 3) ** (k + 1) * p.mahler_word[j+k]
                        for k in range(length - j)), Fraction(0))
            check(tail < Fraction(18, 19), 'every finite suffix tail strictly safe')
        delayed.append(dict(index=index, A=a, matching_bits=length, rational_collatz_source='1/5'))
        tr = transport(a, length)
        check(audit_transport(tr, a) == tr, 'source-retaining transport')
        check((tr.collatz_residue * 5 - 1) % tr.modulus == 0, 'different native source same code')

    check(periodic_mahler_tails((1, 0, 0)) == (Fraction(18, 19), Fraction(8, 19), Fraction(12, 19)),
          'strictly safe denominator-19 rational period')
    check(prefix(Fraction(-9, 19), 99, 'mahler')[0] == prefix(Fraction(1, 5), 99, 'collatz')[0],
          'positive ordinary-integer source gate is essential')
    check(periodic_source((0, 1), 'mahler') == Fraction(-2, 5)
          and periodic_source((0, 1), 'collatz') == 2
          and max(periodic_mahler_tails((0, 1))) == Fraction(6, 5), 'no11 is insufficient for real safety')
    check(compare(1, 1).kind == 'mismatch', 'same initial integer is not same itinerary')

    good = compare(13, Fraction(1, 5))
    for bad in (replace(good, source=True), replace(good, bound=float(good.bound)),
                replace(good, index=True), replace(good, rational_source=0.2),
                replace(good, mahler_word=(True,) + good.mahler_word[1:]),
                replace(good, collatz_word=good.collatz_word[:-1]),
                replace(good, kind='shared11')):
        rejects(lambda bad=bad: audit_separation(bad, 13, Fraction(1, 5)))
    rejects(lambda: audit_separation(good, 32, Fraction(1, 5)))
    rejects(lambda: audit_separation(good, 13, Fraction(3, 5)))
    for bad in (True, 1.0, 0, -1, Fraction(13)):
        rejects(lambda bad=bad: compare(bad, Fraction(1, 5)))
    for bad in (True, 0.2, Fraction(1, 2)):
        rejects(lambda bad=bad: compare(13, bad))
    rejects(lambda: carrier((0, True), 'mahler'))
    rejects(lambda: native_cell([1, 0], 'collatz'))
    rejects(lambda: periodic_source((), 'mahler'))
    rejects(lambda: audit_transport(replace(transport(13, 8), length=8.0), 13))
    rejects(lambda: audit_observer(replace(observers[Fraction(1, 5)], index=3.0), Fraction(1, 5)))
    rejects(lambda: cycle_mismatch_bound(13, observe_rational(3)))
    rejects(lambda: delayed_family(True))

    print(json.dumps(dict(status='PROVED finite source-bound obstruction; FINITE-EXACT controls',
        checks=checks, native_cells=cells, native_lengths=list(range(11)), periodic_words=periodic,
        rational_universe=dict(numerator=[-15, 15], odd_denominator=[1, 15], reduced_count=len(rationals),
                               ordinary_Mahler_sources=[1, 128], no11_cycle_observers=cycle_count),
        event_counts=events, periodic_mismatch_claims=cycle_claims, largest_observed_event=largest_example,
        delayed_safe_family=delayed,
        theorem='Positive ordinary Mahler source plus all-depth no11 implies irrational transported Collatz 2-adic source',
        limits='No Z-number existence/nonexistence; no Collatz convergence; finite safe prefixes do not certify infinite safety'),
        indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
