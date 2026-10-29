#!/usr/bin/env python3
"""Exact finite certificate for INTEGER_RETURNS.md; no optional dependencies.

Run normally and with -O. Checks use explicit exceptions, never assert.
The exhaustive finite reduction is a proof ingredient; random large tests are
only hostile controls. Neither is a proof of convergence of iterated returns.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random

CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(message)


BAD = ((1,), (3, 1, 1, 1), (3, 3, 1, 1, 1, 1, 1),
       (3, 1, 3, 1, 1, 1, 1), (3, 1, 1, 3, 1, 1, 1))


def affine(word):
    c = d = 0
    for a in word:
        c = 3*c + (1 << d)
        d += a
    return len(word), d, c


def guard(word):
    j, d, c = affine(word)
    modulus = 1 << (d+1)
    return (((1 << d)-c)*pow(3**j, -1, modulus)) % modulus, modulus


def predicted_bad(n):
    return (n % 4 == 3 or n % 128 == 125 or
            n % 4096 in (1357, 2269, 2877))


def block(n, checked=False):
    start, q, word, peak = n, 1, [], 1
    k = max(1, (n-1).bit_length())
    while q:
        before, oldq = n, q
        t = 3*n+1
        a = (t & -t).bit_length()-1
        n = t >> a
        q += 2 if a % 3 == 0 else -1
        word.append(a)
        peak = max(peak, q)
        if len(word) > 69*k:
            raise ArithmeticError(('clock bound failed', start, word))
        if checked:
            f = lambda x: 32 if x == 3 else x**3
            # Keep the actual endpoint even at q=0; this is stronger than
            # assigning the absorbing state potential zero.
            require(128*f(n)*4**q <= 125*f(before)*4**oldq,
                    ('patched integer potential', start, before, a))
            if start >= 7:
                require(1372*n**3*4**q <= 1331*before**3*4**oldq,
                        ('unpatched potential', start, before, a))
    word = tuple(word)
    if checked:
        j, d, c = affine(word)
        require((1 << d)*n == 3**j*start+c, 'affine endpoint')
        require(j % 3 == 1, 'tree size')
        require(d <= 2*j+start.bit_length()-1, 'integer division clock')
        require(j+d <= 208*k, 'standard integer clock')
        r, modulus = guard(word)
        require(start % modulus == r, 'source guard')
        require((n > start) == predicted_bad(start), 'five-class prediction')
        require((n == start) == (start == 1), 'fixed endpoint')
    return n, word, peak


def second_block(n):
    """Independent queue execution: literal divisions and two counters."""
    start, branches, leaves, j = n, 0, 0, 0
    while leaves <= 2*branches:
        n = 3*n+1
        a = 0
        while n % 2 == 0:
            n //= 2
            a += 1
        if a % 3:
            leaves += 1
        else:
            branches += 1
        j += 1
        if j > 69*max(1, (start-1).bit_length()):
            raise ArithmeticError('independent clock bound')
    return n, branches


def max_carries():
    cs, witnesses = [1], [(1,)]
    for m in range(1, 15):
        options = []
        for i in range(m):
            for j in range(m-i):
                k = m-1-i-j
                c, d = 1, 3
                for a in (i, j, k):
                    c = 3**(3*a+1)*c + (1 << d)*cs[a]
                    d += 5*a+1
                options.append((c, (i, j, k)))
        c, parts = max(options)
        cs.append(c)
        witnesses.append((3,)+sum((witnesses[a] for a in parts), ()))
        require(affine(witnesses[-1]) == (3*m+1, 5*m+1, c),
                'maximum attained')
    return cs


@lru_cache(None)
def all_words(m):
    if m == 0:
        return ((1,),)
    result = []
    for i in range(m):
        for j in range(m-i):
            k = m-1-i-j
            result.extend((3,)+u+v+w for u in all_words(i)
                          for v in all_words(j) for w in all_words(k))
    return tuple(result)


def main():
    require(2*1331**23 < 1372**23, '23-step halving')
    require(11*1331**15 < 7*1372**15, '15-branch descent')
    cs = max_carries()
    thresholds = [F(cs[m], 2**(5*m+1)-3**(3*m+1)) for m in range(3, 15)]
    cutoff = max(thresholds).__floor__()
    require(cutoff == 178169, 'finite cutoff')
    for m in range(9):
        require(max(affine(w)[2] for w in all_words(m)) == cs[m],
                ('independent tree enumeration', m))

    # Exhaustive remaining universe, no inherited filters: every odd source
    # 1..178169. Both implementations must reject every long expanding block.
    histogram = Counter()
    for n in range(1, cutoff+1, 2):
        h, w, peak = block(n, True)
        h2, m2 = second_block(n)
        require((h2, m2) == (h, (len(w)-1)//3), 'independent finite replay')
        require(not (n > 1 and h >= n and len(w) >= 10),
                ('long nondecreasing exception', n))
        if h > n:
            require(w in BAD, 'exact failure word')
            histogram[str(w)] += 1

    # Larger finite and very-large-source controls do not enter the proof.
    max_length = (0, 0)
    for n in range(1, 1000002, 2):
        h, w, _ = block(n)
        require((h > n) == predicted_bad(n), ('million-source control', n))
        if len(w) > max_length[0]:
            max_length = (len(w), n)
    rng = random.Random(20261029)
    for bits in (32, 64, 128, 512, 2048, 4096):
        for _ in range(50):
            block(rng.getrandbits(bits) | 1, True)
    for k in range(2, 1025):
        block(2**k-1, True)
        block(2**k+1, True)

    rows = []
    p = s = b = F(0)
    for w in BAD:
        j, d, c = affine(w)
        r, modulus = guard(w)
        p += F(1, 2**d)
        s += F(3**j, 4**d)
        b += F(c, 4**d)
        require(3**j > 2**d, 'each failure branch expands')
        for t in (0, 1, 17, 10**100):
            h, actual, _ = block(r+t*modulus, True)
            require(actual == w, 'failure route realization')
        rows.append(dict(word=w, J=j, D=d, C=c, residue=r, modulus=modulus))
    require(p == F(1059, 2048), 'failure mass')
    require(s == F(3235233, 4194304), 'failure slope mass')
    require(b == F(1270035, 4194304), 'failure carry mass')
    require(2*s**3 < 1, 'source refinement payment')
    prefix_checks = []
    for k in range(1, 6):
        mass = slope = carry = F(0)
        for letters in product(BAD, repeat=k):
            w = sum(letters, ())
            j, d, c = affine(w)
            mass += F(1, 2**d)
            slope += F(3**j, 4**d)
            carry += F(c, 4**d)
            r, modulus = guard(w)
            x = r+modulus
            start = x
            for expected in letters:
                x, actual, _ = block(x)
                require(actual == expected, 'same-source concatenation')
            require((1 << d)*x == 3**j*start+c, 'total original-source endpoint')
        require(mass == p**k and slope == s**k, 'prefix mass/slope recurrence')
        require(carry == b*(s**k-p**k)/(s-p), 'prefix carry recurrence')
        prefix_checks.append(dict(k=k, words=5**k))

    result = dict(status='PASS; auxiliary entry proved, iterated descent OPEN',
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  explicit_checks=CHECKS, finite_cutoff=cutoff,
                  exhaustive_odd_sources=(cutoff+1)//2,
                  maximum_carries=cs, thresholds=list(map(str, thresholds)),
                  failure_routes=rows, finite_failure_histogram=dict(histogram),
                  control_maximum_J_source=max_length,
                  failure_mass=str(p), failure_slope_mass=str(s),
                  failure_carry_mass=str(b), source_refinement_factor=str(2*s**3),
                  prefix_checks=prefix_checks)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
