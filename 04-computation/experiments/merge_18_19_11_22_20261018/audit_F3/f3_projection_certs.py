#!/usr/bin/env python3
"""Audit F3: explicit line-projection colourings of G_q (chi of the quotient circulant only, no alpha), written as
certificates in the session's file format and re-verified with the independent checker of f3_kappa_check.py.
Usage: python3 f3_projection_certs.py q [q ...]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from f3_projection import qnr, circle, lines, functional_for_line, chromatic_circulant
from f3_kappa_check import check_certificate

out = []
for q in [int(a) for a in sys.argv[1:]]:
    r = qnr(q)
    mu = circle(q, r)
    best = None
    for v in lines(q):
        l = functional_for_line(q, v)
        S = sorted({l(a, b) for (a, b) in mu})
        if 0 in S:
            continue
        k, col = chromatic_circulant(q, S)
        if best is None or k < best[0]:
            best = (k, v, col)
        break   # all admissible lines are equivalent under the rotations z -> lambda z (lambda in mu_(q+1))
    k, v, col = best
    l = functional_for_line(q, v)
    fn = os.path.join(HERE, 'kappa_%d_colouring_k%d_projection.txt' % (q, k))
    with open(fn, 'w') as fo:
        fo.write('# q=%d r=%d (t^2=r); vertex a*q+b for a+bt; colours 0..%d; proper colouring of Cay(F_q^2, mu_(q+1)) '
                 '[audit F3: colour = phi(%d*a - %d*b mod q), phi a %d-colouring of the quotient circulant]\n'
                 % (q, r, k - 1, v[1], v[0], k))
        fo.write(' '.join(str(col[l(a, b)]) for a in range(q) for b in range(q)) + '\n')
    qq, rr, used, top, bad = check_certificate(fn)
    out.append('q = %d: best line %s, quotient chi = %d; certificate %s: colours used %d, monochromatic edges %d -> %s'
               % (q, v, k, os.path.basename(fn), used, bad, 'PROPER' if bad == 0 else 'NOT PROPER'))
    print(out[-1], flush=True)
open(os.path.join(HERE, 'f3_projection_certs.out'), 'a').write('\n'.join(out) + '\n')
