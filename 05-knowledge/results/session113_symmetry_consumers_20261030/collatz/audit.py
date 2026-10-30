#!/usr/bin/env python3
"""Independent exact-array audit of quotient upper-envelope certificates."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
p = argparse.ArgumentParser()
p.add_argument('E', type=int)
p.add_argument('vector', type=Path)
p.add_argument('--cap', type=int, default=64)
p.add_argument('--output', type=Path)
a = p.parse_args()
E = a.E
q = 5 ** E
D = 8 * q
N = 5 * q - 5
V = np.fromfile(a.vector, dtype='<i8')
if not len(V) == 2 * q:
    raise RuntimeError('verification failed at original line 8')
z = np.arange(q, dtype=np.int64)
h = np.zeros(q, dtype=np.int64)
for e in range(1, E + 1):
    h += z % 5 ** e == 0
ht = np.tile(h, 2)
eps = np.repeat(np.arange(2, dtype=np.int64), q)
zz = np.tile(z, 2)
ones = np.where(eps == 0, -(5 * q // 5 ** ht - 5), -N)
for letter in range(1, a.cap + 1):
    inv = pow(pow(2, letter, q), -1, q)
    znew = ((3 * zz - 2) * inv + 1) % q
    enew = (eps + letter + 1) % 2
    dest = enew * q + znew
    kval = next((e for e in range(E) if (pow(2, letter, q) - 2) % 5 ** (e + 1)), E)
    cost = ones if letter == 1 else np.where(enew == 0, D * np.maximum(h[znew], kval), 0) - N * letter
    if not np.all(V[dest] >= V + cost):
        raise RuntimeError(letter)
    if letter > 1:
        if not V[dest[0]] >= cost[0]:
            raise RuntimeError(('entry', letter))
mn = int(V.min())
mx = int(V.max())
if not mx <= 0:
    raise RuntimeError('verification failed at original line 23')
if not mx + D * E - N * (a.cap + 1) <= mn:
    raise RuntimeError('verification failed at original line 24')
if not D * E - N * (a.cap + 1) <= mn:
    raise RuntimeError('verification failed at original line 25')
r = dict(E=E, states=2 * q, denominator=D, numerator=N, minimum=mn, maximum=mx, letters_checked=a.cap, edge_inequalities=2 * q * a.cap, tail_upper=mx + D * E - N * (a.cap + 1), entry_tail_upper=D * E - N * (a.cap + 1), sha256=hashlib.sha256(a.vector.read_bytes()).hexdigest(), sharp=True, scope='all finite positive valuation words at this fixed E; not all E or convergence')
s = json.dumps(r, indent=2) + '\n'
print(s, end='')
if a.output:
    a.output.write_text(s)
