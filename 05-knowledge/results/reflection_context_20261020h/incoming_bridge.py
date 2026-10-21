"""Exact quotient between two independently developed marked graph carriers.

Incoming twelve-state table: commit098cb67c1d, divisor_reflection_resolution_20261021.
This tests the new inter-carrier map, not a fresh proof of that packet.
"""
from itertools import combinations
import json

edges = tuple(combinations(range(4), 2))
rho = (3, 2, 1, 0)
gamma = (1, 0, 3, 2)


def permute(code, p):
    return sum(1 << edges.index(tuple(sorted((p[i], p[j]))))
               for k, (i, j) in enumerate(edges) if code >> k & 1)


def encode(d):
    a, b, c = d
    return sum(1 << edges.index(e) for e, bit in
               [((0, 3), a == 0), ((1, 2), a == 2), ((0, 1), b), ((2, 3), c)] if bit)


rows = [
    ((1, 0, 0), None, 0), ((1, 1, 1), None, 63),
    ((2, 0, 0), None, 4), ((0, 1, 1), None, 59),
    ((0, 1, 0), None, 1), ((0, 0, 1), None, 32),
    ((2, 0, 1), None, 62), ((2, 1, 0), None, 31),
    ((1, 1, 0), 0, 2), ((1, 1, 0), 1, 47),
    ((1, 0, 1), 0, 61), ((1, 0, 1), 1, 16),
]
quotient = {incoming: encode(d) for d, _, incoming in rows}
checks = 0
for d, epsilon, incoming in rows:
    out = quotient[incoming]
    tests = [
        quotient[63 ^ incoming] == permute(out, gamma) ^ 33,
        quotient[permute(incoming, rho)] == permute(out, rho),
        (incoming == permute(incoming, rho)) == (out == permute(out, rho)),
    ]
    for test in tests:
        checks += 1
        if not test:
            raise RuntimeError((d, epsilon, incoming, out))
fibres = {out: sorted(k for k, v in quotient.items() if v == out) for out in set(quotient.values())}
if len(fibres) != 10 or sorted(map(len, fibres.values())) != [1]*8+[2]*2:
    raise RuntimeError('wrong source quotient')
checks += 1
print(json.dumps({'checks': checks, 'incoming_to_switch_chart': quotient,
                  'nontrivial_fibres': {out: fibre for out, fibre in fibres.items() if len(fibre)>1}},
                 sort_keys=True))
