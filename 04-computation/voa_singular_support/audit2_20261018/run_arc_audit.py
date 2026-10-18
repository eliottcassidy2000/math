#!/usr/bin/env python3
"""Driver for the HYP-9298 audit.  usage: python run_arc_audit.py case Wmax [P]
cases: q7, q7c (engine coefficient 7/39), q8, q8c, q10, q11, hostile, cusp, q10scan, q11scan"""
import sys, random
from fractions import Fraction
from arc_audit import Arc, a2_product, cusp_ci

case, Wmax = sys.argv[1], int(sys.argv[2])
P = int(sys.argv[3]) if len(sys.argv) > 3 else 2147483647


def report(name, rels, ref, W=Wmax):
    A = Arc(rels, W, P)
    hs = A.hs_list(W)
    diff = [w for w in range(W + 1) if hs[w] != ref[w]]
    tag = 'EQUAL through %d' % W if not diff else 'first differs at %d (HS %d vs ref %d)' % (diff[0], hs[diff[0]], ref[diff[0]])
    print(f'{name}: {tag}', flush=True)
    return hs


if case == 'q7':
    ref = a2_product(7, Wmax)
    hs = report('q=7 (xy, y^2-x^3)', [{(1, 1): 1}, {(0, 2): 1, (3, 0): -1}], ref)
    print(' HS ', hs); print(' ref', ref)
elif case == 'q7c':
    report('q=7 (xy, y^2-(7/39)x^3)', [{(1, 1): 1}, {(0, 2): 1, (3, 0): Fraction(-7, 39)}], a2_product(7, Wmax))
elif case == 'q8':
    ref = a2_product(8, Wmax)
    hs = report('q=8 (y^2-x^3, x^2y)', [{(0, 2): 1, (3, 0): -1}, {(2, 1): 1}], ref)
    print(' HS ', hs); print(' ref', ref)
elif case == 'q8c':
    report('q=8 (y^2-(16/279)x^3, x^2y)', [{(0, 2): 1, (3, 0): Fraction(-16, 279)}, {(2, 1): 1}], a2_product(8, Wmax))
elif case == 'q10':
    ref = a2_product(10, Wmax)
    report('q=10 engine (xy^2-(2/123)x^4, y^3-(8/41)x^3y)', [{(1, 2): 1, (4, 0): Fraction(-2, 123)}, {(0, 3): 1, (3, 1): Fraction(-8, 41)}], ref)
elif case == 'q11':
    ref = a2_product(11, Wmax)
    report('q=11 engine (y^3-(11/138)x^3y, x^2y^2-(11/1035)x^5)', [{(0, 3): 1, (3, 1): Fraction(-11, 138)}, {(2, 2): 1, (5, 0): Fraction(-11, 1035)}], ref)
elif case == 'q10scan':
    ref = a2_product(10, Wmax)
    print('q=10 product:', ref)
    random.seed(7)
    rlist = [Fraction(1, 12), Fraction(1), Fraction(1, 6), Fraction(-1, 12), Fraction(1, 3), Fraction(2), Fraction(-1),
             Fraction(1, 24), Fraction(5, 12), 'generic1', 'generic2']
    for r in rlist:
        if isinstance(r, str):
            rr = random.randrange(2, P - 1)
            rels = [{(1, 2): 1, (4, 0): -rr}, {(0, 3): 1, (3, 1): -1}]
        else:
            rels = [{(1, 2): 1, (4, 0): -r}, {(0, 3): 1, (3, 1): -1}]
        A = Arc(rels, Wmax, P)
        hs = A.hs_list(Wmax)
        diff = [w for w in range(Wmax + 1) if hs[w] != ref[w]]
        print(f' ratio beta/gamma={r}: ', 'EQUAL' if not diff else
              f'first differs at {diff[0]} (HS {hs[diff[0]]} vs product {ref[diff[0]]}); HS-product = {[hs[w]-ref[w] for w in range(Wmax+1)]}', flush=True)
elif case == 'q11scan':
    ref = a2_product(11, Wmax)
    print('q=11 product:', ref)
    random.seed(11)
    rlist = [Fraction(15, 2), Fraction(1), Fraction(7), Fraction(8), Fraction(15, 4), Fraction(15), Fraction(-15, 2),
             'generic1', 'generic2']
    for r in rlist:
        if isinstance(r, str):
            rr = random.randrange(2, P - 1)
            rels = [{(0, 3): 1, (3, 1): -rr}, {(2, 2): 1, (5, 0): -1}]
        else:
            rels = [{(0, 3): 1, (3, 1): -r}, {(2, 2): 1, (5, 0): -1}]
        A = Arc(rels, Wmax, P)
        hs = A.hs_list(Wmax)
        diff = [w for w in range(Wmax + 1) if hs[w] != ref[w]]
        print(f' ratio d/e={r}: ', 'EQUAL' if not diff else
              f'first differs at {diff[0]} (HS {hs[diff[0]]} vs product {ref[diff[0]]}); HS-product = {[hs[w]-ref[w] for w in range(Wmax+1)]}', flush=True)
elif case == 'hostile':
    report('q=7 hostile (xy, y^2)', [{(1, 1): 1}, {(0, 2): 1}], a2_product(7, Wmax))
    report('q=7 orbit (xy, y^2+x^3)', [{(1, 1): 1}, {(0, 2): 1, (3, 0): 1}], a2_product(7, Wmax))
    report('q=8 hostile (y^2, x^2y)', [{(0, 2): 1}, {(2, 1): 1}], a2_product(8, Wmax))
    report('q=7 other hostile (xy, x^3)', [{(1, 1): 1}, {(3, 0): 1}], a2_product(7, Wmax))
elif case == 'cusp':
    ref = cusp_ci(Wmax)
    hs = report('cusp (y^2-x^3) vs prod_{n>=2}(1-t^n)^-1/((1-t^3)(1-t^4)(1-t^5))', [{(0, 2): 1, (3, 0): -1}], ref)
    print(' HS ', hs); print(' ref', ref)
    print(' HS-ref', [a - b for a, b in zip(hs, ref)])

if case == 'q13':
    ref = a2_product(13, Wmax)
    report('q=13 TRUE coefficient (xy^3-(13/448)x^4y, y^4-(13/56)x^3y^2+(169/150528)x^6)',
           [{(1, 3): 1, (4, 1): Fraction(-13, 448)}, {(0, 4): 1, (3, 2): Fraction(-13, 56), (6, 0): Fraction(169, 150528)}], ref)
    report('q=13 AUTHOR coefficient (xy^3-(13/448)x^4y, y^4-(13/56)x^3y^2-(28529/3278)x^6)',
           [{(1, 3): 1, (4, 1): Fraction(-13, 448)}, {(0, 4): 1, (3, 2): Fraction(-13, 56), (6, 0): Fraction(-28529, 3278)}], ref)
