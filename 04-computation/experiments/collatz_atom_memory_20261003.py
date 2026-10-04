"""Exact audits of signed atoms, legal-word summaries, and certified splices.

Finite universes are printed. All arithmetic is exact; no orbit is presumed
to terminate. Certified seeds are checked by bounded explicit iteration.
"""
from itertools import product
import json


def v2(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def edge(n, sign=1):
    assert n > 0 and n % 2 == 1 and sign in (-1, 1)
    q = 3*n + sign
    k = v2(q)
    return q >> k, k


def summary(word):
    A, p, S = 0, 0, 0
    for a in word:
        assert a >= 1
        S = 3*S + 2**A
        A += a
        p += 1
    return A, p, S


def compose(u, v):
    A, p, S = u
    B, q, R = v
    return A+B, p+q, 3**q*S + 2**A*R


def splice(b, word):
    A, p, S = summary(word)
    modulus = 3**p
    phases = []
    for j in range(modulus):
        h = (4**j * (3*b+1) - 1)//3
        if (2**A*h-S) % modulus == 0:
            phases.append(j)
    assert len(phases) == 1
    j0 = phases[0]
    R = 4**modulus
    K, remainder = divmod((R-1)*(2**A+3*S), 3**(p+1))
    assert remainder == 0
    return j0, R, K


def main():
    for n in range(1, 10001, 2):
        for sign in (-1, 1):
            m, k = edge(n, sign)
            assert m % 3 != 0
            assert edge(2*n+sign, -sign) == (m, k+1)
            assert edge(4*n+sign, sign) == (m, k+2)
            assert edge(4*n-sign, sign) == (6*n-sign, 1)
    print('Signed atoms: 10000 edges; mirror/doubling lift, atom z lift, opposite-sheet control')
    inverse_count = 0
    for m in range(1, 501, 2):
        if m % 3 == 0:
            continue
        for k in range(1, 31):
            sign = 1 if 2**k*m % 3 == 1 else -1
            n = (2**k*m-sign)//3
            assert edge(n, sign) == (m, k)
            inverse_count += 1
    print(f'Inverse signed chart: {inverse_count} target/exponent pairs')
    assert edge(5, 1)[0] == 1
    assert edge(5, -1)[0] == 7 and edge(7, -1)[0] == 5
    print('Sheet-loss hostile: plus 5 -> 1, minus 5 -> 7 -> 5')
    minus_prefix = [9]
    for _ in range(5):
        minus_prefix.append(edge(minus_prefix[-1], -1)[0])
    assert minus_prefix == [9,13,19,7,5,7]
    previous = None
    for parameter in (0,1):
        numerator = 56*4**(9*parameter+1)+19
        assert numerator % 27 == 0
        n = numerator//27
        if previous is not None:
            assert n == 262144*previous-184471
        previous = n
        for _ in range(2):
            n, k = edge(n, -1)
            assert k == 1
        assert edge(n, -1)[0] == 7
    assert edge(1) == (1,2)
    print('Minus-splice control: parameters 0/1 join cycle 5 <-> 7; root-padding control 1 -> 1')
    for x in range(1, 10001):
        n, steps = x, 0
        while n > 1:
            n = (n+1)//2
            steps += 1
        assert steps == (x-1).bit_length()
    print('Binary parent: exact rank bit_length(x-1), 1 <= x <= 10000')
    for b in (1, 5, 27):
        n = b
        for _ in range(200):
            if n == 1:
                break
            n = edge(n)[0]
        assert n == 1
    checks = 0
    for p in range(1, 5):
        for word in product(range(1, 5), repeat=p):
            A, _, S = summary(word)
            for split in range(p+1):
                assert compose(summary(word[:split]), summary(word[split:])) == (A,p,S)
            exact_r = ((2**A-S)*pow(3**p, -1, 2**(A+1))) % 2**(A+1)
            for b in (1, 5, 27):
                # Independent full permutation test of sibling residues.
                h_residues = {((pow(4,j,3**(p+1))*(3*b+1)-1)//3) % 3**p
                              for j in range(3**p)}
                assert len(h_residues) == 3**p
                j0, R, K = splice(b, word)
                previous = None
                for parameter in (0, 1):
                    j = j0 + 3**p*parameter
                    h = (4**j*(3*b+1)-1)//3
                    numerator = 2**A*h-S
                    assert numerator % 3**p == 0
                    n = numerator//3**p
                    assert n > 0 and n % 2 == 1
                    assert n % 2**(A+1) == exact_r
                    if previous is not None:
                        assert n == R*previous+K
                    previous = n
                    current = n
                    for a in word:
                        current, actual = edge(current)
                        assert a == actual
                    assert current == h and edge(h)[0] == edge(b)[0]
                    checks += 1
    print(f'Splices: {checks} instances; word lengths 1..4, letters 1..4, seeds 1/5/27, parameters 0/1')
    print('Each splice checked by forward valuation replay, residue permutation, exact cylinder, and recurrence')
    assert summary((1,2)) == (3,2,5)
    assert summary((2,1)) == (3,2,7)
    print('Order-loss hostile: (1,2) and (2,1) have same clock (A,p) but carries 5 and 7')
    j0, R, K = splice(1, (1,1))
    assert (j0,R,K) == (4,262144,184471)
    record = {'seed':1, 'head':[1,1], 'summary':{'A':2,'p':2,'S':5},
              'sibling_phase':4, 'sibling_period':9,
              'source_initial':151, 'recurrence_multiplier':R,
              'recurrence_offset':K, 'first_route':[151,227,341,1],
              'scope':'all integer parameters >= 0; certified plus-sheet family'}
    print('EXAMPLE_JSON=' + json.dumps(record, separators=(',',':')))
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
