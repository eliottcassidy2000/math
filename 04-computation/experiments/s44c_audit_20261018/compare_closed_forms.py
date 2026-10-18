#!/usr/bin/env python3
"""Exact (Fractions, no primes) comparison of the recursion's G_L(N) with THM-4710's closed forms at many N,
including N < 9 (non-stable range for larger L), the pre-registered range 49..64, and N = 70, 100, 150, 200."""
import json, os, sys, time
from fractions import Fraction
import sympy as sp
from witt_closure_exact import closure, cayley_coords
HERE = os.path.dirname(os.path.abspath(__file__))
Ns_ = sp.Symbol('N')
fits = json.load(open(os.path.join(HERE, '..', '..', 'quantum_transfer_s44', 'witt_rigidity_interpolation_9_48.json')))['fits']
F = {}
for key, ent in fits.items():
    L = int(key.split()[0].split('_')[1])
    nu = tuple(int(x) for x in key.split(' ', 1)[1].strip('()').split(',') if x.strip())
    F[(L, nu)] = sp.sympify(ent['expr'], locals={'N': Ns_})
Nlist = list(range(3, 21)) + list(range(49, 65)) + [70, 100, 150, 200]
res = {}
allok = True
t0 = time.time()
for N in Nlist:
    S, M = closure(N)
    dims = {L: len(S[L]) for L in (3, 5, 7, 9)}
    mism = []
    for L in (5, 7, 9):
        if N >= L:
            if dims[L] != 1:
                mism.append((L, 'dim', dims[L]))
                continue
            cay = cayley_coords(S[L][0], L)
            for nu, v in cay.items():
                if nu == (L,):
                    continue
                pv = sp.Rational(F[(L, nu)].subs(Ns_, N))
                if Fraction(int(pv.p), int(pv.q)) != v:
                    mism.append((L, nu, str(v), str(pv)))
        else:
            if dims[L] != 0:
                mism.append((L, 'dim_should_be_0', dims[L]))
    if dims[3] != (1 if N >= 3 else 0):
        mism.append((3, 'dim', dims[3]))
    res[N] = {'dims': dims, 'mismatches': mism}
    allok &= not mism
    print(f'N={N}: dims {dims}, mismatches {len(mism)}  ({time.time()-t0:.1f}s)', flush=True)
print('ALL EXACT COMPARISONS AGREE WITH THE CLOSED FORMS:', allok)
json.dump({'N': Nlist, 'results': {str(k): v for k, v in res.items()}, 'all_ok': allok},
          open(os.path.join(HERE, 'compare_closed_forms.json'), 'w'), indent=1, default=str)
