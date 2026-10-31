#!/usr/bin/env python3
"""Exact integral completion, certified by one index-two enlargement and closure."""
import json
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from audit_feedback import gaps, require

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / 'minimal_feedback.json').read_text())
A = [sp.Matrix(a) for a in data['matrices']]
K = sp.Matrix(data['hankel_minor'])
x0 = sp.Matrix(data['initial_observations'])
suffixes = [tuple(s) for s in data['suffix_basis']]
r = sp.Matrix([[gaps((2, 0, 3) + s) for s in suffixes]])
require(any(t.q != 1 for t in r * K.inv()), 'new row was already in old lattice')
require(all(t.q == 1 for t in 2 * r * K.inv()), 'quotient is not killed by two')
Z = sp.Matrix.vstack(K, r)
Z = sp.Matrix([[int(t) for t in row] for row in Z.tolist()])
B = hermite_normal_form(Z.T).T
require(abs(K.det()) == 2 * abs(B.det()), 'index-two completion failed')
require(all(t.q == 1 for t in Z * B.inv()), 'old and extra generators not included')
# Because K has index two in B and r is not in K, these checks prove
# B = Z-span(K,r), independently of trusting HNF's generation claim.
C = [B * a * B.inv() for a in A]
for c in C:
    require(all(t.q == 1 for t in c), 'completed lattice not invariant')
initial, final = x0 * B.inv(), B[:, 0]
require(all(t.q == 1 for t in initial), 'initial state not integral')
alternate = sp.Matrix([[gaps((1, 0, 0) + s) for s in suffixes]])
require(alternate[0] == r[0] == 5, 'equal H control failed')
require(all(t.q == 1 for t in alternate * K.inv()), 'comparison not in old lattice')
S = smith_normal_form(B, domain=sp.ZZ)
encode = lambda m: [[int(m[i, j]) for j in range(m.cols)] for i in range(m.rows)]
result = {'status': 'passed', 'old_index': int(abs(K.det())),
          'full_reachable_observation_lattice_index': int(abs(B.det())),
          'extra_prefix': [2, 0, 3], 'same_order_same_H_old_coset_prefix': [1, 0, 0],
          'quotient_order': 2, 'equal_H': 5,
          'smith_invariants': [int(abs(S[i, i])) for i in range(17)],
          'basis': encode(B), 'initial': encode(initial), 'final': encode(final),
          'integer_matrices': [encode(c) for c in C],
          'scope': 'Z-span of all reachable observation rows, not the set of reachable words or a positivity cone.'}
(HERE / 'integral_feedback.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k not in ('basis', 'integer_matrices')},
                 sort_keys=True, indent=2))
