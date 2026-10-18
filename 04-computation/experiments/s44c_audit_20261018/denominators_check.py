#!/usr/bin/env python3
"""List every denominator met in the generic (Q(N)) recursion of witt_symbolic.closure_symbolic: entries of the rref
bases of M_L (L <= 8), of the S_L basis vectors (L <= 9), and of all condition matrices.  Integer roots N0 >= 1 of
these denominators are the only points where the generic computation might fail to specialize."""
import json, os
import sympy as sp
from witt_symbolic import closure_symbolic, K, Nsym
HERE = os.path.dirname(os.path.abspath(__file__))
S, M, cond = closure_symbolic()
dens = set()
for L, vecs in list(M.items()) + list(S.items()):
    for v in vecs:
        for x in v.values():
            dens.add(sp.factor(sp.denom(sp.together(K.to_sympy(x)))))
for L, (rows, n, nb) in cond.items():
    for r in rows:
        for x in r:
            dens.add(sp.factor(sp.denom(sp.together(K.to_sympy(x)))))
factors = set()
for d in dens:
    for f, e in sp.factor_list(d)[1]:
        factors.add(sp.factor(f))
roots = sorted(set(int(r) for f in factors for r in sp.Poly(f, Nsym).ground_roots() if r.is_integer))
print('distinct irreducible denominator factors:', sorted(str(f) for f in factors))
print('integer roots:', roots)
json.dump({'factors': sorted(str(f) for f in factors), 'integer_roots': roots}, open(os.path.join(HERE, 'denominators_check.json'), 'w'), indent=1)
