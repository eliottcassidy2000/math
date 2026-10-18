#!/usr/bin/env python3
"""Evaluate the pre-registered out-of-sample tests P9, P10 of THM-4710 (PREREGISTRATION_witt_closed_forms_20261018.md).

P9: for N = 49..64 compute the forced generators exactly (mod six primes, CRT, rational reconstruction) and compare
every Cayley coefficient with the closed forms fitted on N = 9..48 (witt_rigidity_interpolation_9_48.json).
P10: extract the degree-18 generators of Vir_(4,13) from its singular support and compare with G_5(18).
Usage: python check_witt_closed_forms_20261018.py [Nmin Nmax]
"""
import os, sys, json
from fractions import Fraction
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
from witt_rigidity_interpolation_20261018 import forced_generators
from semiinvariant_generators_20261018 import model_generators

Nsym = sp.Symbol('N')


def main(Nmin=49, Nmax=64):
    fits = json.load(open(os.path.join(HERE, 'witt_rigidity_interpolation_9_48.json')))['fits']
    formulas = {}
    for key, ent in fits.items():
        L = int(key.split()[0].split('_')[1])
        nu = tuple(int(x) for x in key.split(' ', 1)[1].strip('()').split(',') if x.strip())
        formulas[(L, nu)] = sp.sympify(ent['expr'], locals={'N': Nsym})
    out = {'P9': [], 'P10': None}
    allok = True
    for N in range(Nmin, Nmax + 1):
        dims, g = forced_generators(N)
        okdims = all(dims.get(L) == 1 for L in (3, 5, 7, 9))
        mism = []
        for (L, nu), f in formulas.items():
            val = g.get(L, {}).get(nu)
            pred = sp.Rational(f.subs(Nsym, N))
            if val is None or sp.Rational(val.numerator, val.denominator) != pred:
                mism.append((L, nu, str(val), str(pred)))
        ok = okdims and not mism
        allok &= ok
        out['P9'].append({'N': N, 'dims_ok': okdims, 'mismatches': mism})
        print(f'N={N}: dims_ok={okdims}, coefficient mismatches={len(mism)}', flush=True)
    print(f'P9 {"PASSED" if allok else "FAILED"}', flush=True)
    # P10
    N, res = model_generators(4, 13)
    gens = {L: e['nullity'] for L, e in res.items() if e['nullity']}
    cay5 = res[5].get('cayley', {})
    want = {(5,): Fraction(1), (3, 2): Fraction(-89, 73)}
    ok10 = (N == 18 and gens == {0: 1, 3: 1, 5: 1} and {k: v for k, v in cay5.items()} == want)
    out['P10'] = {'N': N, 'generator_positions': gens, 'G5': {str(k): str(v) for k, v in cay5.items()}, 'passed': ok10}
    print(f'P10: Vir(4,13) N={N}, lowest-weight generators at Li {gens}, G_5 = {cay5}: {"PASSED" if ok10 else "FAILED"}', flush=True)
    # P10, Li-7 clause: generator positions through Li 7 (one prime, weight 2N + 7 = 43)
    from vir_ss import SingularSupport
    from semiinvariant_generators_20261018 import lowest_weight_vector
    Pp = 2147483647
    ss = SingularSupport(4, 13, 2 * 18 + 7, Pp)
    ss.run(verbose=False)
    pos7 = {L: len(lowest_weight_vector(ss, 2 * 18 + L, 18, Pp)) for L in range(0, 8)}
    pos7 = {L: n for L, n in pos7.items() if n}
    ok10b = pos7 == {0: 1, 3: 1, 5: 1}
    out['P10']['positions_through_Li7'] = pos7
    out['P10']['passed_Li7_clause'] = ok10b
    print(f'P10 (Li-7 clause): lowest-weight vectors of grI_18 through Li 7 at {pos7}: {"PASSED" if ok10b else "FAILED"}', flush=True)
    json.dump(out, open(os.path.join(HERE, f'check_witt_closed_forms_{Nmin}_{Nmax}.json'), 'w'), indent=1, default=str)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(int(a[0]) if a else 49, int(a[1]) if len(a) > 1 else 64)
