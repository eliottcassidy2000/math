"""Exact forcing and formal Z[zeta_6]-linear motion checks."""
from itertools import product
import json
names = ['A', 'B', 'C', 'D', 'Bp', 'Ap', 'Cp', 'Dp']
ix = {v: i for i, v in enumerate(names)}
zero = (0, 0)
one = (1, 0)
zet = (0, 1)
zb = (1, -1)
add = lambda x, y: (x[0] + y[0], x[1] + y[1])
neg = lambda x: (-x[0], -x[1])

def mul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] + x[1] * y[1])
P = {'A': (zero, zero), 'B': (one, zero), 'C': (zet, zero), 'D': (zero, one), 'Bp': (zero, zet), 'Ap': (one, zet), 'Cp': (zb, zet), 'Dp': (one, neg(zb))}
triples = [('B', 'A', 'D'), ('C', 'A', 'Bp'), ('C', 'B', 'Dp'), ('Bp', 'Ap', 'Dp'), ('Cp', 'Ap', 'B'), ('Cp', 'Bp', 'D')]
units = [one, zet, neg(zb), neg(one), neg(zet), zb]
T = [(one, zero), (zero, zero), (zero, one)]
for tri, u in zip(triples, units):
    p = P[tri[1]]
    actual = [tuple((add(p[i], mul(u, v[i])) for i in range(2))) for v in T]
    if not actual == [P[k] for k in tri]:
        raise RuntimeError('Exact verification gate failed: scripts/check_eight_point_forcing.py:19')
    if not mul(u, (u[0] + u[1], -u[1])) == one:
        raise RuntimeError('Exact verification gate failed: scripts/check_eight_point_forcing.py:20')
violations = 0
forced = 0
for colors in product([0, 1], repeat=8):
    mono = lambda e: len({colors[ix[k]] for k in e}) == 1
    lhs = int(mono(('A', 'B', 'C')))
    rhs = sum((mono(t) for t in triples))
    if not lhs <= rhs:
        raise RuntimeError('Exact verification gate failed: scripts/check_eight_point_forcing.py:25')
    forced += lhs
print(json.dumps({'result': 'PASS', 'colorings': 256, 'monochromatic_ABC_colorings': forced, 'exact_motion_identities': 6, 'coefficient_ring': 'Z[zeta]/(zeta^2-zeta+1)', 'verified_inequality': '1_mono(ABC) <= sum of six template-copy mono indicators'}, sort_keys=True))
