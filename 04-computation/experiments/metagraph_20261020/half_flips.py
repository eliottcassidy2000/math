#!/usr/bin/env python3
"""Half-flip lines: the mirror m1 (path reversal) splits the tiles into left tiles (c < (n-r)/2), right tiles (c > (n-r)/2)
and the axis tiles (2c = n-r). L = flip the left tiles, R = flip the right tiles, X = flip the axis tiles; c = LRX is the
complement. For n = 4..7: how often L preserves the class / maps to the converse class / preserves H; the Klein group
<m1, c> orbit sizes; and the number of nodes reachable by one half-flip from each node (half-flip degree)."""
import subprocess
from collections import Counter
def tiles_of(n): return [(a, b) for b in range(1, n - 1) for a in range(b + 2, n + 1)]
def to_d6(n, arcs):
    A = [[0] * n for _ in range(n)]
    for u, v in arcs: A[u][v] = 1
    bits = [A[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6): s += chr(63 + int(''.join(map(str, bits[k:k + 6])), 2))
    return s
def labelg(lines):
    r = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(r) == len(lines); return r
for n in range(4, 8):
    tiles = tiles_of(n); m = len(tiles); idx = {t: i for i, t in enumerate(tiles)}
    left = 0; right = 0; axis = 0
    for k, (a, b) in enumerate(tiles):
        r, c = a - b - 1, b
        if 2 * c < n - r: left |= 1 << k
        elif 2 * c > n - r: right |= 1 << k
        else: axis |= 1 << k
    def arcs_of(mask):
        arcs = [(i, i - 1) for i in range(1, n)]
        for k, (a, b) in enumerate(tiles): arcs.append((b - 1, a - 1) if (mask >> k) & 1 else (a - 1, b - 1))
        return arcs
    N = 1 << m
    can = labelg([to_d6(n, arcs_of(mk)) for mk in range(N)]); conv = labelg([to_d6(n, [(v, u) for u, v in arcs_of(mk)]) for mk in range(N)])
    cls = {}; 
    for c_ in can: cls.setdefault(c_, len(cls))
    node = {}
    for mk in range(N):
        a, b = cls[can[mk]], cls[conv[mk]]; node[mk] = min(a, b)
    res = Counter()
    deg = Counter()
    for mk in range(N):
        for name, msk in (("L", left), ("R", right), ("LR", left | right), ("X", axis)):
            t2 = mk ^ msk
            res[f"{name} preserves node"] += node[t2] == node[mk]
        deg[node[mk]] = deg[node[mk]]  # ensure key
    # half-flip node graph: nodes adjacent if some tiling of u has L-image in v
    Ladj = set(); Xadj = set()
    for mk in range(N):
        Ladj.add((node[mk], node[mk ^ left])); Xadj.add((node[mk], node[mk ^ axis]))
    print(f"n = {n}: tiles {m} (left {bin(left).count('1')}, right {bin(right).count('1')}, axis {bin(axis).count('1')}); of {N} tilings:", dict(res),
          f"| L-adjacency pairs {len(Ladj)}, X-adjacency pairs {len(Xadj)}, nodes {len(set(node.values()))}")
