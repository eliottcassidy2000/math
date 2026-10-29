#!/usr/bin/env python3
"""Exact least source of k consecutive growing returns, 1<=k<=40.

Two branch-and-bound implementations: direct affine cylinders and canonical
source/endpoint lifts. Pruning follows inclusion of arithmetic progressions.
This is a finite horizon result, not a conjectural recurrence extrapolation.
"""
from hashlib import sha256
from pathlib import Path
import json

B = ((1, 1, 1, 3), (4, 6, 179, 125), (7, 11, 16177, 1357),
     (7, 11, 12289, 2269), (7, 11, 9697, 2877))
WORDS = ((1,), (3, 1, 1, 1), (3, 3, 1, 1, 1, 1, 1),
         (3, 1, 3, 1, 1, 1, 1), (3, 1, 1, 3, 1, 1, 1))


def require(ok, msg):
    if not ok:
        raise ArithmeticError(msg)


def direct(k):
    best, witness = 2**(k+1)-1, (0,)*k
    stack, visited = [(0, 0, 0, 0, ())], 0
    while stack:
        depth, J, D, C, word = stack.pop()
        visited += 1
        for i, (j, d, c, residue) in enumerate(B):
            JJ, DD, CC = J+j, D+d, 3**j*C+2**D*c
            modulus = 2**(DD+1)
            r = ((2**DD-CC)*pow(3**JJ, -1, modulus)) % modulus
            if r >= best:
                continue
            if depth+1 == k:
                best, witness = r, word+(i,)
            else:
                stack.append((depth+1, JJ, DD, CC, word+(i,)))
    return best, witness, visited


def lifted(k):
    best, witness = 2**(k+1)-1, (0,)*k
    # State: block depth, odd clock, valuation clock, least source, endpoint
    # at that source, labelled word. Empty word has source=endpoint=1.
    stack, visited = [(0, 0, 0, 1, 1, ())], 0
    inv = {j: pow(3**j, -1, 2048) for j in range(7*k+1)}
    while stack:
        depth, J, D, r, h, word = stack.pop()
        visited += 1
        for i, (j, d, c, target) in enumerate(B):
            t = (((target-h)//2)*inv[J]) % (2**d)
            rr = r+2**(D+1)*t
            if rr >= best:
                continue
            hh = (3**j*(h+2*3**J*t)+c)//2**d
            if depth+1 == k:
                best, witness = rr, word+(i,)
            else:
                stack.append((depth+1, J+j, D+d, rr, hh, word+(i,)))
    return best, witness, visited


def validate(n, word):
    for i in word:
        before = n
        q, actual = 1, []
        while q:
            n = 3*n+1
            a = 0
            while n % 2 == 0:
                n //= 2
                a += 1
            actual.append(a)
            q += 2 if a % 3 == 0 else -1
        require(tuple(actual) == WORDS[i], 'native witness word')
        require(n > before, 'native witness grows')
    return n


def main():
    result = []
    for k in range(1, 41):
        a = direct(k)
        b = lifted(k)
        require(a == b, ('independent least-source recurrence', k))
        n, word, visits = a
        endpoint = validate(n, word)
        validate(2**(k+1)-1, (0,)*k)  # authenticated initial upper bound
        result.append(dict(k=k, least_source=n, word=list(word), endpoint=endpoint,
                           visited_prefixes=visits))
    require(result[14]['least_source'] == 56575, 'Mersenne-only hostile')
    require(all(a['least_source'] <= b['least_source'] for a, b in zip(result, result[1:])),
            'nested failure cylinders')
    require(result[-1]['least_source'] == 406360948735, 'horizon forty')
    print(json.dumps(dict(status='FINITE-EXACT; horizons 1..40, all positive odd sources',
                          source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                          rows=result), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
