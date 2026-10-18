#!/usr/bin/env python3
"""Identify the cross-edge families found by golden_product_direction_census_20261017.py for Y = uv and Y = e^2 (S41).

For each exact mixed unit dA + u d1 + v d2 + Y d3 the script prints the coefficient vector numerically at the intended
embedding and tests whether all of them are K-unit multiples q (alpha, beta, gamma, delta) of one fixed vector, i.e.
form a single contact family like clock_duality's q(d, 1, 1). It reports the normalised vector and |q|.
Usage: python golden_product_relations_20261017.py
Session opus-2026-10-17-S41.
"""
import os, sys, json, cmath
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escaped_ports_20261017d import make, need


def main():
    f, z, j, alpha, ports, units, eta, conj, escape, unused = make('sqrt3')
    # numeric embedding: zeta15 = exp(2 pi i k/15) chosen so that alpha = zeta15^3 = exp(2 pi i/5); j = i sqrt 11; w real root
    import numpy as np
    wr = [r.real for r in np.roots([1, 2, -4, -4]) if abs(r.imag) < 1e-12 and 1 < r.real < 2][0]
    zn = cmath.exp(2j * cmath.pi / 15); jn = 1j * 11 ** 0.5
    def num(a):
        s = 0
        for c, e in zip(a, f.basis):
            if c:
                s += float(c) * zn ** e[0] * jn ** e[1] * wr ** e[2]
        return s
    rho = f.scale(f.add(f.scalar(5), j), Q(1, 6)); baseD = f.add(f.one, eta)
    M = [f.zero, f.one, eta, baseD, rho, f.mul(rho, eta), f.mul(rho, baseD)]
    A = list(ports); seen = set(A)
    for port in ports:
        for m in M:
            t = f.add(port, m)
            if t not in seen:
                seen.add(t); A.append(t)
    C = list(M)
    for p0 in ports[:5]:
        if p0 not in C:
            C.append(p0)
    a5 = ports[5]
    Lp = [f.zero, f.one, f.sub(f.one, a5), f.scalar(Q(1, 4)), f.scalar(Q(5, 4)), f.sub(f.scalar(Q(5, 4)), a5), f.scalar(Q(1, 2))]
    L = list({x: 0 for x in Lp + [f.scale(x, -1) for x in Lp]}.keys())
    DA = [f.zero] + sorted({f.sub(a, b) for a in A for b in A} - {f.zero})
    DC = [f.zero] + sorted({f.sub(a, b) for a in C for b in C} - {f.zero})
    D3 = sorted(({f.sub(a, b) for a in C for b in C} | {f.sub(a, b) for a in L for b in L}) - {f.zero})
    here = os.path.dirname(os.path.abspath(__file__))
    res = os.path.join(here, '..', '..', '05-knowledge', 'results')
    u = units[0]; d = f.add(f.one, alpha)
    print('numeric check: |u| = %.12f, d = %s' % (abs(num(u)), num(d)))
    for Y in ('uv', 'ee'):
        R = json.load(open(os.path.join(res, 'golden_product_direction_census_%s_20261017.json' % Y)))
        vecs = []
        for iA, i1, i2, i3 in R['mixed']:
            v = [num(DA[iA]), num(DC[i1]), num(DC[i2]), num(D3[i3])]
            vecs.append(v)
        print('direction %s: %d mixed units' % (Y, len(vecs)))
        # normalise by the last coordinate (d3 != 0)
        for v in vecs:
            q = v[3]
            print('   (dA, d1, d2, d3)/d3 = (%s, %s, %s, 1);  |d3| = %.6f, |d1| = %.6f, |dA| = %.6f'
                  % tuple(['%.6f%+.6fi' % (x.real, x.imag) for x in (v[0] / q, v[1] / q, v[2] / q)] + [abs(q), abs(v[1]), abs(v[0])]))


if __name__ == '__main__':
    main()
