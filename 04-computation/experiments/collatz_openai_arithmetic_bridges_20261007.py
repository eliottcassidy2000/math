"""Exact transfer controls; no audit of the accepted OpenAI manuscripts.

The finite computations check our affine maps and the supplied idoneal list.
They do not estimate an asymptotic correlation or certify new ROOT sources.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt, prod

from collatz_floor_transport_deadlines_20261005 import step, word_weight

IDONEAL = (1,2,3,4,5,6,7,8,9,10,12,13,15,16,18,21,22,24,25,28,30,
           33,37,40,42,45,48,57,58,60,70,72,78,85,88,93,102,105,112,
           120,130,133,165,168,177,190,210,232,240,253,273,280,312,
           330,345,357,385,408,462,520,760,840,1320,1365,1848)


def positive_integer(n):
    if type(n) is not int or n < 1:
        raise ValueError('exact positive integer required')


def reduced_forms(n):
    """Primitive positive reduced forms of discriminant -4n, one per class."""
    positive_integer(n)
    out = []
    for a in range(1, isqrt(4*n//3)+1):
        for b in range(-a, a+1):
            if (b*b+4*n) % (4*a):
                continue
            c = (b*b+4*n)//(4*a)
            if a > c or gcd(gcd(a,b),c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            out.append((a,b,c))
    return tuple(out)


def ambiguous(form):
    a,b,c = form
    return b == 0 or abs(b) == a or a == c


def carrier(word):
    if type(word) is not tuple or not word:
        raise ValueError('nonempty exact valuation tuple required')
    P,Q,B = 1,1,0
    for a in word:
        positive_integer(a)
        P,B,Q = 3*P,3*B+Q,Q*(1 << a)
    return P,Q,B


def cylinder(word):
    P,Q,B = carrier(word)
    r = (Q-B)*pow(P,-1,2*Q) % (2*Q)
    e = (P*r+B)//Q
    return r,2*Q,e,2*P,B


def replay(n, word):
    positive_integer(n)
    if n % 2 != 1:
        raise ValueError('odd source required')
    for a in word:
        n, actual = step(n)
        if actual != a:
            raise ValueError('wrong source valuation')
    return n


def root_word(n):
    word = []
    for _ in range(100):
        if n == 1:
            return tuple(word)
        n,a = step(n)
        word.append(a)
    raise ValueError('finite control cap exhausted')


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(label)

    found = []
    for n in range(1,1849):
        forms = reduced_forms(n)
        check(bool(forms), ('principal form exists',n))
        check(all(b*b-4*a*c == -4*n for a,b,c in forms), ('discriminant',n))
        if all(map(ambiguous,forms)):
            found.append(n)
    check(tuple(found) == IDONEAL, 'entire supplied list, bounded range')
    check(reduced_forms(11) == ((1,0,11),(3,-2,4),(3,2,4)), '11 is not idoneal')
    check(len(reduced_forms(105)) == 8, '105 has eight classes, one per genus')
    print('Primitive reduced forms: every n=1..1848; exactly the supplied 65 pass ambiguity.')
    print('  n11: three classes; n105: eight classes. Idoneal is not class number one.')

    words = 0
    for length in range(1,5):
        for w in product((1,2,3),repeat=length):
            r,a,e,c,B = cylinder(w)
            check(a*e-c*r == 2*B > 0, ('nonproportional forms',w))
            check(r % 2 == e % 2 == 1, ('odd endpoints',w))
            for t in (1,2,7,31,100):
                check(replay(r+a*t,w) == e+c*t, ('actual first-hit prefix',w,t))
            words += 1
    check(cylinder((1,2)) == (11,16,13,18,5), 'word12 interface')
    check(cylinder((1,1,2)) == (7,32,13,54,19), 'word112 interface')
    print('Legal affine-cylinder controls: %d words, five actual sources each; determinant=2B.' % words)
    print('  word12: n=11+16t, endpoint=13+18t, determinant10; t>=1 avoids ROOT padding.')

    weights = {n:word_weight(root_word(n)) for n in (3,5,15)}
    check(weights == {3:F(1,6),5:F(1,3),15:F(1,252)}, 'exact weight control')
    check(weights[15] != weights[3]*weights[5], 'Collatz weight is not multiplicative')
    check(step(5)[0] == 1, 'absorbing invariant measure cannot give 5 positive mass')
    print('Weight multiplicativity obstruction: W15=1/252 != W3*W5=1/18.')

    # A finite fixed-prime filter is an exact translation of coprimality.
    r,a,_,_,_ = cylinder((1,2))
    primes = (3,5,7,11)
    M = prod(primes)
    forbidden = {p:-r*pow(a,-1,p) % p for p in primes}
    shift = next(b for b in range(M) if all(b % p == forbidden[p] for p in primes))
    for t in range(2*M):
        check((gcd(r+a*t,M) == 1) == (gcd(t-shift,M) == 1), 'translated prime filter')
    check(shift == (-r*pow(a,-1,M)) % M, 'CRT shift')
    print('Jacobsthal interface: word12, prime filter{3,5,7,11}, all 2310 parameters checked.')

    # Same positive local intersection count, opposite real slope directions.
    for w, expected_anchor in (((1,),F(-1)),((1,2),F(-5)),((1,4),F(5,23))):
        P,Q,B = carrier(w)
        anchor = F(B,Q-P)
        check(anchor == expected_anchor, ('fixed point',w))
        check(P != Q and F(P,Q)*anchor+F(B,Q) == anchor, 'transverse graph intersection')
    print('Serre toy interface: every listed affine graph/diagonal intersection has length1;')
    print('  anchors -1,-5,5/23 retain neither positive-integer existence nor descent.')
    for bad in (True,0,-1,1.0):
        try:
            reduced_forms(bad)
        except ValueError:
            checks += 1
        else:
            raise RuntimeError('type hostile accepted')
    print('PASS: %d exact transfer controls; accepted external theorems were not audited.' % checks)


if __name__ == '__main__':
    main()
