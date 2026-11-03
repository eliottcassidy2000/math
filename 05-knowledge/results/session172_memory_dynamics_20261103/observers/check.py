#!/usr/bin/env python3
"""Standalone exact checks. No assertions, external packages, or repository imports."""
import itertools
import json
from fractions import Fraction


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def step(x, b, coeff):
    a = x % b
    u, v = coeff[a]
    require((u*x+v) % b == 0, 'branch integrality')
    return (u*x+v)//b


def word(x, k, b, coeff):
    out = []
    for _ in range(k):
        out.append(x % b)
        x = step(x, b, coeff)
    return tuple(out)


def affine(w, b, coeff):
    U, C, power = 1, 0, 1
    for a in w:
        u, v = coeff[a]
        U, C = u*U, u*C+power*v
        power *= b
    return U, C


def valuation_mod(x, b, cap):
    k = 0
    while k < cap and x % b == 0:
        x //= b
        k += 1
    return k


def audit(name, b, coeff, depth):
    counters = {'cylinders': 0, 'precision_checks': 0, 'metric_pairs': 0}
    for k in range(1, depth+1):
        modulus = b**k
        bank = {}
        for x in range(modulus):
            w = word(x, k, b, coeff)
            require(w not in bank, 'distinct residues have distinct words')
            bank[w] = x
            U, C = affine(w, b, coeff)
            require((-C*pow(U, -1, modulus)) % modulus == x, 'word residue')
            y = x
            for _ in range(k):
                y = step(y, b, coeff)
            require(modulus*y == U*x+C, 'affine iteration')
            require(word(x+modulus, k, b, coeff) == w, 'cylinder preservation')
            delta = step(x+modulus, b, coeff)-step(x, b, coeff)
            require(delta % (b**(k-1)) == 0, 'lower precision update')
            require(delta % modulus != 0, 'same precision obstruction')
            counters['cylinders'] += 1
            counters['precision_checks'] += 1
        require(len(bank) == modulus, 'all words realized')
        require(all(w in bank for w in itertools.product(range(b), repeat=k)), 'surjectivity')
    # Exhaustive finite metric check on a lower-resolution set.
    k = min(depth, 5)
    words = [word(x, k, b, coeff) for x in range(b**k)]
    for x, wx in enumerate(words):
        for y, wy in enumerate(words):
            common = next((j for j in range(k) if wx[j] != wy[j]), k)
            require(common == valuation_mod(x-y, b, k), 'observer metric')
            counters['metric_pairs'] += 1
    return {'name': name, 'base': b, 'depth': depth, **counters}


def main():
    systems = [audit('ceil_half', 2, ((1,0),(1,1)), 10),
               audit('collatz_3', 2, ((1,0),(3,1)), 10),
               audit('collatz_5', 2, ((1,0),(5,1)), 10),
               audit('ternary_affine', 3, ((1,0),(2,1),(4,1)), 7),
               audit('composite_base', 4, ((1,0),(3,1),(1,2),(3,3)), 5)]
    cycle = [13,33,83,208,104,52,26]
    require([step(x,2,((1,0),(5,1))) for x in cycle] == cycle[1:]+cycle[:1], '5x+1 cycle')
    w = tuple(x % 2 for x in cycle)
    U3,C3 = affine(w,2,((1,0),(3,1)))
    U5,C5 = affine(w,2,((1,0),(5,1)))
    require((U3,C3) == (27,19), '3x+1 affine constants')
    require((U5,C5) == (125,39), '5x+1 affine constants')
    rational = Fraction(C3,2**len(w)-U3)
    require(rational == Fraction(19,101), 'rational image cycle')
    # Fraction with odd denominator has an unambiguous parity.
    x = rational
    observed = []
    for _ in range(7):
        a = (x.numerator * pow(x.denominator,-1,2)) % 2
        observed.append(a)
        x = ((3*x+1) if a else x)/2
    require(tuple(observed) == w and x == rational, 'rational orbit exact')
    for n in range(2,10001):
        require(step(n,2,((1,0),(1,1))) < n, 'easy positive descent')
    return {'status':'PASS', 'systems':systems, 'five_cycle':cycle,
            'shared_periodic_word':list(w), 'three_rational_cycle_start':str(rational),
            'easy_descent_checked_through':10000}

if __name__ == '__main__':
    print(json.dumps(main(), sort_keys=True, indent=2))
