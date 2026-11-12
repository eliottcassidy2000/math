#!/usr/bin/env python3
"""Verifier: are the (30,30,30) affine 3-cubes of order 7 the Paley tournament P_7 with a Z_3-invariant 3-matching?"""
import itertools
from vlib import brute_H
n = 7; QR = {1, 2, 4}
A = [[1 if (j - i) % 7 in QR else 0 for j in range(n)] for i in range(n)]
print("H(P_7) =", brute_H(n, A))
def orbit(a, b):
    return [((a * m) % 7, (b * m) % 7) for m in (1, 2, 4)]
seen = set()
for a in range(1, 7):
    for b in range(1, 7):
        if a != b and A[a][b]:
            orb = tuple(sorted(orbit(a, b)))
            if orb in seen: continue
            seen.add(orb)
            verts = set(x for arc in orb for x in arc)
            if len(verts) != 6: continue
            vals = {}
            for sub in itertools.product((0, 1), repeat=3):
                B = [row[:] for row in A]
                for bit, (x, y) in zip(sub, orb):
                    if bit:
                        B[x][y] = 0; B[y][x] = 1
                vals[sub] = brute_H(n, B)
            h0 = vals[(0,0,0)]
            d = [vals[(1,0,0)] - h0, vals[(0,1,0)] - h0, vals[(0,0,1)] - h0]
            affine = all(vals[s] == h0 + sum(di for di, si in zip(d, s) if si) for s in vals)
            print("Z_3-orbit matching", orb, "values", sorted(vals.values()), "increments", d, "affine:", affine)
