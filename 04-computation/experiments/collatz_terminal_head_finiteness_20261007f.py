"""Exact algebra controls for the fixed-middle-head S-unit specialization.

Finiteness is a cited theorem, not inferred from this bounded enumeration.
An endpoint equation is necessary; strict first-hit legality is still checked.
"""
from itertools import product
import collatz_marked_completion_20261007c as marked
import collatz_uncovered_join_routes_20261007 as routes

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def head_constant(word):
    routes.letters(word)
    need(not word or word[0] >= 2, 'maximal initial run is outside the middle head')
    p, q, b = routes.carrier(word)
    invariant = 3*b-3*p+q
    need(invariant % 4 == 2, 'nonzero invariant has exact dyadic valuation one')
    return invariant//2


def equation(exponent, word, terminal):
    need(type(exponent) is int and exponent >= 2, 'positive nonroot Mersenne exponent')
    need(type(terminal) is int and terminal >= 1, 'positive exact terminal valuation')
    constant = head_constant(word)
    return (1 << (sum(word)+terminal-1))-3**(exponent+len(word)) == constant


def first_hit(exponent, word, terminal):
    if not equation(exponent, word, terminal):
        return False
    try:
        return marked.replay((1 << exponent)-1, (1,)*(exponent-1)+word+(terminal,)) == 1
    except ValueError:
        return False


def main():
    heads = [()]
    for length in range(1, 5):
        heads += [w for w in product(range(1, 5), repeat=length) if w[0] >= 2]
    solutions, strict = [], []
    for w in heads:
        C = head_constant(w)
        invariant, cost = -2, 0
        for a in w:
            cost += a
            invariant = 3*invariant+(1 << cost)
        need(invariant == 2*C, 'independent carry-invariant recurrence')
        p, q, b = routes.carrier(w)
        for E in range(2, 13):
            for c in range(1, 25):
                eq = equation(E, w, c)
                x = 2*3**(E-1)-1
                need(eq == (3*(p*x+b)+q == q*(1 << c)),
                     'independent unnormalized terminal equation')
                if eq:
                    solutions.append((E, w, c))
                    if first_hit(E, w, c):
                        strict.append((E, w, c))
    need(equation(3, (2, 3), 4) and first_hit(3, (2, 3), 4), 'source7 ROOT control')
    need(equation(2, (4, 2), 2) and not first_hit(2, (4, 2), 2),
         'arithmetic endpoint does not authorize ROOT padding')
    for args in ((True, (), 4), (2, (), True), (2, (1,), 4)):
        try:
            equation(*args)
        except ValueError:
            need(True, 'malformed exact type or head rejected')
        else:
            need(False, 'hostile accepted')
    print('CITED input: Beukers-Schlickewei Theorem1.1, at most2^(8r+8) solutions.')
    print('PROVED specialization: fixed middle head => at most2^32 Mersenne exponents.')
    print('Normalized group generators: (2,1), (1,3), (1/C,-1/C); rank<=3.')
    print('FINITE-EXACT universe:', len(heads), 'heads; E2..12; terminal1..24.')
    print('Endpoint-equation solutions:', len(solutions), '; strict first-hit solutions:', len(strict))
    print('Hostile: M2 head(4,2) terminal2 satisfies equation but pads ROOT.')
    print('No height bound, exhaustive solution algorithm, or universal obstruction to variable heads.')
    print('Checks:', CHECKS)


if __name__ == '__main__':
    main()
