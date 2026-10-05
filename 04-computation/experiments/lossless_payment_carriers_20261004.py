"""Exact source-preserving affine, quadratic and guarded H/G carriers.

Run with python -B, also python -B -O. No import-time computation or writes.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


def need(ok, message):
    if not ok:
        raise ValueError(message)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'exact positive odd source')


def vp(n, p):
    need(type(n) is int and n != 0, 'nonzero integer valuation')
    n, k = abs(n), 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def vf(x, p):
    x = F(x)
    need(x != 0, 'nonzero rational valuation')
    return vp(x.numerator, p)-vp(x.denominator, p)


@dataclass(frozen=True)
class Carrier:
    P: int
    Q: int
    B: int

    def __post_init__(self):
        need(all(type(x) is int for x in (self.P, self.Q, self.B))
             and self.P > 0 and self.Q > 0 and self.B >= 0,
             'positive integer slope data and nonnegative carry')

    def at(self, n):
        return F(self.P*n+self.B, self.Q)


I, H, G = Carrier(1, 1, 0), Carrier(729, 1024, 669), Carrier(9, 8, 5)
LETTERS = {'H': H, 'G': G}
R, S, D = F(729, 1024), F(9, 8), F(67, 32)
V = (1, 2, 1, 1, 1, 2)


def compose(first, second):
    need(type(first) is Carrier and type(second) is Carrier, 'typed carriers')
    return Carrier(second.P*first.P, second.Q*first.Q,
                   second.P*first.B+second.B*first.Q)


def word(w):
    need(type(w) is str and all(c in LETTERS for c in w), 'H/G string')


def compile_word(w):
    word(w)
    result = I
    for c in w:
        result = compose(result, LETTERS[c])
    return result


def decode(data):
    need(type(data) is Carrier, 'typed carrier')
    u, v = vp(data.P, 3), vp(data.Q, 2)
    need(data.P == 3**u and data.Q == 2**v, 'H/G prime support')
    need((2*v-3*u) % 2 == 0, 'integral H count')
    h, g = (2*v-3*u)//2, 5*u-3*v
    need(h >= 0 and g >= 0, 'nonnegative H/G counts')
    total = (F(data.B, data.Q)+5-5*F(data.P, data.Q))/D
    reverse = []
    while h:
        need(total > 0, 'positive H translation')
        value = vf(total, 3)
        need(value >= 0 and value % 2 == 0, 'trailing G valuation')
        k = value//2
        need(k <= g, 'trailing G count within budget')
        reverse.extend('G'*k+'H')
        total = (total/S**k-1)/R
        h, g = h-1, g-k
    need(total == 0, 'pure G terminal translation')
    reverse.extend('G'*g)
    result = ''.join(reversed(reverse))
    need(compile_word(result) == data, 'exact carrier reconstruction')
    return result


def cylinder(data):
    """Exact legal positive odd source class for a decoded H/G program."""
    decode(data)
    modulus = 2*data.Q
    return ((data.Q-data.B)*pow(data.P, -1, modulus)) % modulus, modulus


def decode_carry(data):
    """Independent unshifted decoder: terminal B has v3 zero(G) or one(H)."""
    need(type(data) is Carrier, 'typed carrier')
    result, current = [], data
    while current != I:
        need(current.B > 0, 'nonempty positive carry')
        value = vp(current.B, 3)
        need(value in (0, 1), 'terminal carry valuation')
        c = 'G' if value == 0 else 'H'
        last = LETTERS[c]
        need(current.P % last.P == 0 and current.Q % last.Q == 0,
             'terminal slope divisibility')
        p, q = current.P//last.P, current.Q//last.Q
        numerator = current.B-last.B*q
        need(numerator >= 0 and numerator % last.P == 0, 'terminal carry divisibility')
        current = Carrier(p, q, numerator//last.P)
        result.append(c)
    w = ''.join(reversed(result))
    need(compile_word(w) == data, 'independent carry reconstruction')
    return w


def credit(w):
    word(w)
    current = minimum = 0
    for c in w:
        current += 2 if c == 'H' else -1
        minimum = min(minimum, current)
    return minimum, current


def credit_join(left, right):
    for value in (left, right):
        need(type(value) is tuple and len(value) == 2
             and all(type(x) is int for x in value), 'integer credit summary')
        need(value[0] <= min(0, value[1]), 'prefix-minimum summary')
    return min(left[0], left[1]+right[0]), left[1]+right[1]


def literal(n, exponents):
    odd(n)
    for a in exponents:
        need(n != 1, 'strict first-hit prefix')
        z = 3*n+1
        need(vp(z, 2) == a, 'actual valuation guard')
        n = z >> a
    return n


def stage(n, c):
    odd(n)
    need(c in LETTERS, 'H/G stage')
    data = LETTERS[c]
    value = data.at(n)
    need(value.denominator == 1 and value.numerator % 2 == 1,
         'stage oddness guard')
    result = value.numerator
    if c == 'G':
        need(literal(n, (1, 2)) == result, 'G actual word')
    else:
        join = literal(n, V)
        need(join == 4*result+1, 'H common-future branch')
        a = vp(3*join+1, 2)
        need(a >= 3 and literal(join, (a,)) == literal(result, (a-2,)),
             'H common-future witness')
    return result


def monitor(n, w, require_credit=False):
    odd(n)
    word(w)
    if require_credit:
        need(credit(w)[0] >= 0, 'prefix credit guard')
    data, current, states = I, n, [n]
    for c in w:
        current = stage(current, c)
        data = compose(data, LETTERS[c])
        need(data.at(n) == current, 'immutable original source')
        states.append(current)
    debt = (data.P-data.Q)*n+data.B
    need(F(debt, data.Q) == current-n, 'original debt identity')
    return tuple(states), data, debt


def lift(a, b, dimension):
    a, b = F(a), F(b)
    if dimension == 3:
        return ((a, 0, b), (0, 1, 0), (0, 0, 1))
    if dimension == 5:
        return ((a*a, 0, 2*a*b, 0, b*b), (0, 1, 0, 0, 0),
                (0, 0, a, 0, b), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1))
    need(dimension == 6, 'observer dimension 3, 5 or 6')
    return ((a*a, 0, 0, 2*a*b, 0, b*b), (0, a, 0, 0, b, 0),
            (0, 0, 1, 0, 0, 0), (0, 0, 0, a, 0, b),
            (0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1))


def vector(x, n, dimension):
    if dimension == 3:
        return x, n, 1
    if dimension == 5:
        return x*x, n*n, x, n, 1
    need(dimension == 6, 'observer dimension')
    return x*x, x*n, n*n, x, n, 1


def mv(m, v):
    return tuple(sum(a*b for a, b in zip(row, v)) for row in m)


def mm(a, b):
    return tuple(tuple(sum(x*y for x, y in zip(row, col))
                       for col in zip(*b)) for row in a)


def energy(n):
    odd(n)
    if n == 1:
        return F(0), 0
    k = vp(n-1, 2)
    return F(3, 4)**k*(n-1)**2, k


def main():
    words = [''.join(w) for length in range(9) for w in product('HG', repeat=length)]
    source_checks = credit_checks = 0
    for w in words:
        data = compile_word(w)
        need(decode(data) == w and decode_carry(data) == w, 'two independent word decoders')
        residue, modulus = cylinder(data)
        for n in (residue+modulus, residue+2*modulus):
            states, actual, debt = monitor(n, w)
            need(actual == data, 'guard sufficiency')
            source_checks += 1
            if w and credit(w)[0] == 0:
                need(debt < 0, 'credit program pays immutable source')
                credit_checks += 1
        for cut in range(len(w)+1):
            need(compose(compile_word(w[:cut]), compile_word(w[cut:])) == data,
                 'affine cut composition')
            need(credit_join(credit(w[:cut]), credit(w[cut:])) == credit(w),
                 'credit cut composition')
        if w:
            bad = residue+modulus+modulus//2
            try:
                monitor(bad, w)
            except ValueError:
                pass
            else:
                raise ValueError('one-bit guard hostile accepted')
    print('H/G lengths 0..8:', len(words), 'round trips by two independent decoders;', source_checks,
          'independently replayed legal sources;', credit_checks, 'credit-paid controls')
    small = [w for w in words if len(w) <= 3]
    for a, b, c in product(small, repeat=3):
        need(credit_join(credit_join(credit(a), credit(b)), credit(c)) ==
             credit_join(credit(a), credit_join(credit(b), credit(c))),
             'credit monoid associativity')
    print('Credit associativity:', len(small)**3, 'triples; HG/GH carry and prefix hostile')
    need(compile_word('HG').P == compile_word('GH').P and
         compile_word('HG').B != compile_word('GH').B and
         credit('HG') == (0, 1) and credit('GH') == (-1, 1), 'order hostile')
    observer_checks = composition_checks = 0
    values = (F(-2), F(-1, 2), F(0), F(1), F(3, 2))
    for a, b, x, n in product(values, repeat=4):
        for dim in (3, 5, 6):
            need(mv(lift(a, b, dim), vector(x, n, dim)) == vector(a*x+b, n, dim),
                 'polynomial observer transport')
            observer_checks += 1
    for a, b, c, d in product(values, repeat=4):
        for dim in (3, 5, 6):
            need(mm(lift(c, d, dim), lift(a, b, dim)) == lift(c*a, c*b+d, dim),
                 'observer composition')
            composition_checks += 1
    print('Polynomial controls:', observer_checks, 'state transports;', composition_checks,
          'composition identities over five rational values')
    data = compile_word('GGGH')
    residue, modulus = cylinder(data)
    n = residue or modulus
    states, _, debt = monitor(n, 'GGGH')
    need(states[-1] < states[-2] and debt > 0, 'local descent versus original debt')
    print('Local/global hostile GGGH:', states, '; original debt numerator', debt)
    need(15 < 17 and energy(15) > energy(17), 'numeric versus proper rank hostile')
    need(literal(171, (1,)) == 257 and energy(257) < energy(171), 'rank pays growth')
    print('Rank controls: E(15)=147 > E(17)=81; 171->257 lowers E 21675->6561')
    invalid = 0
    for action in (lambda: Carrier(True, 1, 0), lambda: Carrier(1, 1.0, 0),
                   lambda: decode(Carrier(7, 8, 1)), lambda: decode(Carrier(729, 1024, 670)),
                   lambda: monitor(True, ''), lambda: monitor(1.0, ''),
                   lambda: monitor(155, 'GH', True), lambda: compile_word('X')):
        try:
            action()
        except ValueError:
            invalid += 1
        else:
            raise ValueError('invalid control accepted')
    print('Invalid type/carrier/credit controls rejected:', invalid)
    print('PASS: exact guarded program memory and source payment; no universal completion claim')


if __name__ == '__main__':
    main()
