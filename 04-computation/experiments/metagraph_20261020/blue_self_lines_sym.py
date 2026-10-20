#!/usr/bin/env python3
"""Blue self-lines at n = 9, 10 by enumerating only the symmetric tilings (2^{(m + floor((n-1)/2))/2} of them) and testing
[T_t] = [T_{t^c}] with labelg; prediction 0 (odd n) and 2^{n/2-2} = 8 (n = 10)."""
import subprocess, itertools, sys
def tiles_of(n): return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
def to_d6(n, arcs):
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def labelg(lines, chunk=250000):
    out = []
    for s0 in range(0, len(lines), chunk):
        part = lines[s0:s0 + chunk]
        r = subprocess.run(['labelg', '-q'], input='\n'.join(part) + '\n', capture_output=True, text=True).stdout.split()
        assert len(r) == len(part); out.extend(r)
    return out
for n in (9, 10):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    mirror = [idx[(n + 1 - b, n + 1 - a)] for a, b in tiles]
    orbits = []; seen = set()
    for i, j in enumerate(mirror):
        if i in seen: continue
        seen |= {i, j}; orbits.append((i, j))
    allmask = (1 << m) - 1
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    syms = []
    for bits in itertools.product([0, 1], repeat=len(orbits)):
        mk = 0
        for b, (i, j) in zip(bits, orbits):
            if b: mk |= (1 << i) | (1 << j)
        syms.append(mk)
    print(f"n = {n}: symmetric tilings {len(syms)}", flush=True)
    can = labelg([to_d6(n, arcs_of(mk)) for mk in syms]); canc = labelg([to_d6(n, arcs_of(allmask ^ mk)) for mk in syms])
    hits = sum(1 for c1, c2 in zip(can, canc) if c1 == c2)
    print(f"n = {n}: blue self-lines = {hits // 2} (prediction {0 if n % 2 else 2 ** (n // 2 - 2)}); distinct classes among symmetric tilings {len(set(can))}", flush=True)
