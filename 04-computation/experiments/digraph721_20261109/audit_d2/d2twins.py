"""audit D2: which odd k < 512 have NO realiser on exactly L(k) vertices that contains an adjacent-twin pair?
Exhaustive given the census of connected graphs with d <= 511: a graph on L(k) vertices with d = k has a component C
containing the twin pair; the rest realises k / d(C) on >= Lany(k / d(C)) vertices (threshold graphs attain it)."""
import sys
def bl(k): return k.bit_length()
def Lany(m): return 0 if m == 1 else bl(m)
def parse(line):
    p = line.split(); n = int(p[0]); e = p[5] if len(p) > 5 else ''; adj = [0] * n
    for i in range(0, len(e), 2):
        u = int(e[i], 16); v = int(e[i + 1], 16); adj[u] |= 1 << v; adj[v] |= 1 << u
    return adj
def has_twins(G):
    n = len(G)
    return any((G[x] >> y) & 1 and (G[x] | 1 << x) == (G[y] | 1 << y) for x in range(n) for y in range(x + 1, n))
mintw = {}                          # least order of a connected graph with adjacent twins, per d
for line in open(sys.argv[1]):
    p = line.split(); n = int(p[0]); d = int(p[2])
    if n >= mintw.get(d, 99): continue
    if has_twins(parse(line)): mintw[d] = n
no = [k for k in range(3, 512, 2) if k not in (5, 37)
      and not any(k % a == 0 and o + Lany(k // a) <= bl(k) for a, o in mintw.items())]
print('odd k < 512 (not 1, 5, 37) with no twin-carrying realiser on L(k) vertices:', len(no), no)
