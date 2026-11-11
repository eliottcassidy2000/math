#!/usr/bin/env python3
"""Input for hampaths_layered.c: the tournament P_q - {0, 1} (q - 2 vertices), start set = out-neighbours of 1.
The number of Hamiltonian paths of P_q is q * (q-1)/2 * (number of Hamiltonian paths starting 0 -> 1 -> ...)."""
import sys
from paley_tools import paley_out
q = int(sys.argv[1]); out = paley_out(q)
verts = [v for v in range(q) if v not in (0, 1)]; idx = {v: i for i, v in enumerate(verts)}
print(len(verts))
for v in verts:
    mask = 0
    for w in verts:
        if w != v and (out[v] >> w) & 1: mask |= 1 << idx[w]
    print(mask, 1 if (out[1] >> v) & 1 else 0)
