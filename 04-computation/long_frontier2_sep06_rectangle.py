#!/usr/bin/env python3
"""Exact rectangular-prism certificate for the full carried response.

Universe: symbolic complete rows, 56 affine phase-endpoint corners,
324 multihomogenized response coefficients, a uniform B-root barrier,
a genuine open nonnegative-root prism, and two scoped hostile controls.
SymPy supplies exact rational polynomial arithmetic only. No source imports
from repository producers and no numerical root decisions are used.
"""
import argparse
import hashlib
import itertools
import json
import sys
from math import comb
import sympy as S

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(newline='\n')

A, B, F, T = S.symbols('a b f s')
GATES = 0
TRACE = hashlib.sha256()


def gate(ok, label):
    global GATES
    if ok != True:
        raise RuntimeError(label + ': ' + str(ok))
    GATES += 1
    TRACE.update((label + '\n').encode())


def convolution(p, q):
    out = [S.Integer(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return out


def transform(poly, boxes):
    """Finite [l,r]: x=(l+r*u)/(1+u); tail [l,infinity): x=l+u."""
    out = S.Integer(0)
    degrees = poly.degree_list()
    for powers, coefficient in poly.terms():
        term = coefficient
        for i, (power, (lo, hi)) in enumerate(zip(powers, boxes)):
            x = poly.gens[i]
            term *= ((lo+x)**power if hi is None else
                     (lo+hi*x)**power * (1+x)**(degrees[i]-power))
        out += term
    return S.Poly(S.expand(out), *poly.gens)


def records(poly):
    return [{'powers': list(powers), 'coefficient': str(c)}
            for powers, c in sorted(poly.terms())]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--certificate')
    args = parser.parse_args()
    beta = list(map(S.Integer, [1, 13, 55])) + [A, B, F]
    c = list(map(S.Integer, [1, 12, 45])) + [2*A/3, 3*B/7]
    d = list(map(S.Integer, [1, 11, 36])) + [5*A/12, B/7]
    gg, cd = convolution(beta, beta), convolution(c, d)
    # Coefficient j of the original Laurent row Q: beta^2 and 2t C D
    # have their own -2 and -1 shifts. The response multiplier is the
    # full binomial convolution, including its exponent -1 lower carry.
    raw = {}
    qbar = S.Integer(0)
    for j in range(-1, 9):
        raw[j] = S.expand(comb(28, 2*j+2) *
                          (gg[j+2] + (2*cd[j+1] if j+1 < len(cd) else 0)))
        qbar += raw[j] * (-1 if j % 2 else 1) * T**(j+1)
    gate(raw[-1] == 28, 'negative-support lower carry')
    p = F*T**4 - 12*B*T**3/7 + A*T**2 - 10*T + S.Rational(1, 11)
    original = sum(comb(14, 2*j+1)*beta[j+1] *
                   (-1 if j % 2 else 1)*T**j for j in range(5))
    gate(S.expand(original-2002*p) == 0, 'literal original first row')
    fzero = 12*B/(7*T)-A/T**2+10/T**3-S.Rational(1, 11)/T**4
    gate(S.cancel(p.subs(F, fzero)) == 0, 'retain the original zero')
    h = S.Poly(S.cancel(-S.Rational(11, 14)*qbar.subs(F, fzero)), A, B, T)
    gate(h.degree_list() == (2, 2, 8), 'full eliminated multidegree')
    gate(len(h.terms()) == 19, 'all nineteen eliminated monomials')
    old = [-443993, 73031400, -3278853435, 46232902140,
           -234760993560, 343030019640, -83518139925,
           -26087589000, 3585421125]
    gate(S.expand(h.as_expr().subs({A: 84, B: 35})-
                  sum(x*T**j for j, x in enumerate(old))) == 0,
         'inherited four-anchor polynomial recovered')
    abox = (S.Rational(167, 2), S.Rational(169, 2))
    bbox = (S.Rational(69, 2), S.Rational(71, 2))
    intervals = [(S.Rational(99, 10000), S.Rational(1, 100)),
                 (S.Rational(1, 9), S.Rational(13, 100)),
                 (S.Integer(1), S.Rational(8, 5)),
                 (S.Integer(10), None)]
    endpoint_records = []
    expected = [(1, -1), (-1, 1), (1, -1), (-1,)]
    charts = []
    for k, (lo, hi) in enumerate(intervals):
        for aa, bb, ff in itertools.product(abox, bbox, (0, 5)):
            values = [p.subs({A: aa, B: bb, F: ff, T: x})
                      for x in ([lo, hi] if hi is not None else [lo])]
            for j, value in enumerate(values):
                gate(value*expected[k][j] > 0,
                     'phase endpoint '+str((k, aa, bb, ff, j)))
            endpoint_records.append({'chart': k, 'a': str(aa), 'b': str(bb),
                                     'f': ff, 'values': list(map(str, values))})
        chart = transform(h, [abox, bbox, (lo, hi)])
        gate(len(chart.terms()) == 81, 'complete chart coefficient bank '+str(k))
        for powers in itertools.product(range(3), range(3), range(9)):
            gate(chart.coeff_monomial(A**powers[0]*B**powers[1]*T**powers[2]) > 0,
                 'positive transformed coefficient '+str((k, powers)))
        charts.append({'phase': [str(lo), str(hi) if hi is not None else None],
                       'minimum_coefficient': str(min(chart.coeffs())),
                       'coefficients': records(chart)})
    # F_upper dominates F_a,b on [0,2/5]. Its derivative at2/5 also
    # dominates every allowed derivative there. Root interlacing places
    # the first critical point below2/5 for every nonnegative-root B.
    upper = T**5-13*T**4+55*T**3-abox[0]*T**2+bbox[1]*T
    derivative = S.diff(upper, T).subs(T, S.Rational(2, 5))
    gate(derivative == -S.Rational(81, 10), 'uniform first-critical barrier')
    barrier = transform(S.Poly(5-upper, T), [(0, S.Rational(2, 5))])
    for j in range(6):
        gate(barrier.nth(j) > 0, 'positive upper-product barrier '+str(j))
    # This whole smaller prism has five distinct positive B roots,
    # supplying open three-dimensional nonvacuity, not just one centre.
    nodes = [0, S.Rational(1, 10), 1, 3, 5, 7]
    bcontrols = []
    for aa, bb, ff in itertools.product(abox, bbox, (S.Rational(1, 2), S.Rational(3, 2))):
        values = [x**5-13*x**4+55*x**3-aa*x**2+bb*x-ff for x in nodes]
        for j, value in enumerate(values):
            gate(value*(-1 if j % 2 == 0 else 1) > 0,
                 'genuine B-prism endpoint '+str((aa, bb, ff, j)))
        bcontrols.append({'a': str(aa), 'b': str(bb), 'f': str(ff),
                          'values': list(map(str, values))})
    # Hostile outside the fifth-coefficient bound: an exact positive
    # response at an original root; no nonnegative-root B is asserted.
    hlo, hhi = S.Rational(16693, 2000), S.Rational(41733, 5000)
    p6 = S.Poly(p.subs({A: 84, B: 35, F: 6}), T)
    gate(p6.eval(hlo)*p6.eval(hhi) < 0, 'f6 genuine original-root bracket')
    hostile_h = transform(S.Poly(h.as_expr().subs({A: 84, B: 35}), T), [(hlo, hhi)])
    for j in range(9):
        gate(hostile_h.nth(j) < 0, 'positive-response hostile coefficient '+str(j))
    # The original-zero condition is also essential at a genuine centre.
    gate(p.subs({A: 84, B: 35, F: 1, T: 4}) != 0, 'off-root control is off root')
    gate(qbar.subs({A: 84, B: 35, F: 1, T: 4}) == 350398552675052,
         'positive full response away from original zero')
    certificate = {
        'schema_version': 1, 'scope': 'rectangular coefficient prism, no interlacer hypothesis',
        'a_interval': list(map(str, abox)), 'b_interval': list(map(str, bbox)),
        'f_interval': ['0', '5'], 'identity': 's Q(-s)=-(14/11) H(a,b,s) at P(-s)=0',
        'H': records(h), 'phase_endpoint_corners': endpoint_records,
        'positive_charts': charts, 'upper_product_barrier': records(barrier),
        'genuine_B_prism_corners': bcontrols,
        'f6_hostile_bracket': [str(hlo), str(hhi)],
        'f6_negative_H_transform': records(hostile_h)}
    data = (json.dumps(certificate, indent=2, sort_keys=True)+'\n').encode()
    if args.certificate:
        with open(args.certificate, 'wb') as target:
            target.write(data)
    print('UNIVERSE exact full rows; rectangular a,b,f prism; all original phases; B-only product bound; two scoped hostiles')
    print('PRISM a=[167/2,169/2] b=[69/2,71/2] f=[0,5]; Q(-s)<0 at every original phase')
    print('ROOTS f>0: four simple phases in (99/10000,1/100),(1/9,13/100),(1,8/5),(10,infinity); f=0: first three')
    print('CERTIFICATE four charts, 81 strictly positive coefficients each; lower carry retained')
    for i, chart in enumerate(charts):
        print('CHART', i, 'MINIMUM_COEFFICIENT', chart['minimum_coefficient'])
    print('B_ONLY nonnegative-root B implies f<5 throughout the a,b rectangle; no C/D interlacing used')
    print('NONVACUITY the whole same a,b rectangle with f in [1/2,3/2] has five distinct positive B roots')
    print('HOSTILES f6 outside the bound has positive response at an original phase; centre f1,s4 has positive response off the original zero')
    print('CERTIFICATE_SHA256', hashlib.sha256(data).hexdigest())
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__ == '__main__':
    main()
