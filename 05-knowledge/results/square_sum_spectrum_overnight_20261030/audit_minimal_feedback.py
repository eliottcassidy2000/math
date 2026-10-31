#!/usr/bin/env python3
"""Audit the lower-bound minor by vertex-subset counting, and padding hostiles."""
from itertools import product
import json
from pathlib import Path
import sympy as sp
from audit_feedback import gaps, require, subset_dp, tournament

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / 'minimal_feedback.json').read_text())
prefixes, suffixes = data['prefix_basis'], data['suffix_basis']
K = sp.Matrix([[subset_dp(tournament(tuple(p) + tuple(q))) for q in suffixes]
               for p in prefixes])
require(K == sp.Matrix(data['hankel_minor']), 'independent Hankel entries')
require(K.det() == 150994944, 'independent nonzero determinant')
A = [sp.Matrix(a) for a in data['matrices']]
require(all(36 % value.q == 0 for matrix in A for value in matrix),
        'unexpected denominator prime or exponent')
initial = sp.Matrix(data['initial_observations'])
controls = 0
for n in range(5):
    for word in product(range(4), repeat=n):
        x = initial
        for a in word:
            x = x * A[a]
        require(list(x) == [gaps(word + tuple(q)) for q in suffixes],
                'reduced continuation mismatch')
        controls += 1
# All suffixes through length 3. Non-padding suffixes blow up even on
# same-length inputs with H=1; source-padding supplies a constant comparison.
hostiles = 0
for n in range(1, 4):
    for suffix in product(range(4), repeat=n):
        for m in range(1, 7):
            require(gaps((0,) * m) == gaps((3,) * m) == 1, 'equal input count')
            require(gaps((0,) * m + suffix) == gaps(suffix), 'source-padding identity')
            value = gaps((3,) * m + suffix)
            if any(a != 3 for a in suffix):
                require(value >= 2 ** m + 1, 'induced near-transitive lower bound')
            else:
                require(value == 1, 'sink padding identity')
            hostiles += 1
print(json.dumps({'status': 'passed', 'independent_subset_dp_minor_entries': 289,
                  'minor_determinant': str(K.det()), 'minimal_dimension': 17,
                  'all_coordinate_word_controls': controls,
                  'padding_hostile_controls': hostiles,
                  'scope': 'Minor uses an independent vertex-subset algorithm; transfer controls complement the universal matrix identities.'},
                 sort_keys=True, indent=2))
