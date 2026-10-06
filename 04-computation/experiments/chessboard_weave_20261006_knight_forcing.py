#!/usr/bin/env python3
"""Machine check of the two hand 'forced 12-cycle' arguments (chessboard-weave, 2026-10-06).

Degree-forcing propagation on a set of allowed edges: a square with exactly 2
allowed edges forces both; a square with 2 forced edges loses its other allowed
edges.  Iterate; report forced cycles.

 (A) 8x8 closed tour with m23 = 16.  Counting forces the vector
     (m01,m02,m11,m12,m13,m22,m23,m33) = (0,8,0,0,24,8,16,8): ring-1 squares use only
     1-3 edges, all 8 ring-2 inner edges and all 8 ring-3 inner edges are used,
     no 0-1 / 1-1 / 1-2 edges.  Propagation must close a cycle of length < 64.
 (B) frame = rings 2,3 (48 squares): propagation from degree-2 squares closes two
     disjoint 12-cycles, so no closed tour.
Run: python3 chessboard_weave_20261006_knight_forcing.py
"""
from itertools import product

def ring(s):
    return int(max(abs(s[0] - 3.5), abs(s[1] - 3.5)) - 0.5)
def alg(s):
    return "abcdefgh"[s[1]] + str(s[0] + 1)
SQ = [(i, j) for i, j in product(range(8), repeat=2)]
E = set()
for (i, j) in SQ:
    for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1)]:
        a, b = i + di, j + dj
        if 0 <= a < 8 and 0 <= b < 8:
            E.add(frozenset(((i, j), (a, b))))

def propagate(V, allowed, forced):
    allowed, forced = set(allowed), set(forced)
    changed = True
    while changed:
        changed = False
        for v in V:
            fv = [e for e in forced if v in e]
            av = [e for e in allowed if v in e]
            if len(fv) > 2:
                return "CONTRADICTION: square %s has 3 forced edges" % alg(v), forced
            if len(fv) == 2:
                for e in av:
                    if e not in forced:
                        allowed.discard(e); changed = True
            if len(av) < 2:
                return "CONTRADICTION: square %s has < 2 allowed edges" % alg(v), forced
            if len(av) == 2:
                for e in av:
                    if e not in forced:
                        forced.add(e); changed = True
    return "no contradiction from degrees", forced

def cycles(V, forced):
    adj = {v: [] for v in V}
    for e in forced:
        a, b = tuple(e)
        adj[a].append(b); adj[b].append(a)
    seen, out = set(), []
    for s in V:
        if s in seen or len(adj[s]) != 2:
            continue
        walk, prev, cur, closed = [s], None, s, False
        while True:
            nxt = [y for y in adj[cur] if y != prev]
            if len(adj[cur]) != 2 or not nxt:
                break
            prev, cur = cur, nxt[0]
            if cur == s:
                closed = True; break
            walk.append(cur)
        seen |= set(walk)
        if closed:
            out.append(walk)
    return out

# (A)
V = SQ
def pair(e):
    a, b = tuple(e)
    return tuple(sorted((ring(a), ring(b))))
allowed = {e for e in E if pair(e) in {(0, 2), (1, 3), (2, 3), (2, 2), (3, 3)}}
forced0 = {e for e in E if pair(e) in {(2, 2), (3, 3)}}
msg, forced = propagate(V, allowed, forced0)
cyc = [c for c in cycles(V, forced) if len(c) < 64]
print("(A) m23 = 16 configuration:", msg, f"; forced edges {len(forced)}")
for c in cyc:
    print(f"    forced {len(c)}-cycle:", " ".join(alg(s) for s in c))
print("    => no Hamiltonian cycle (a proper forced cycle exists):", len(cyc) > 0)

# (B)
VB = [s for s in SQ if ring(s) >= 2]
allowedB = {e for e in E if all(ring(x) >= 2 for x in e)}
msg, forced = propagate(VB, allowedB, set())
cyc = [c for c in cycles(VB, forced) if len(c) < len(VB)]
print("\n(B) frame (rings 2,3):", msg, f"; forced edges {len(forced)}")
for c in cyc:
    print(f"    forced {len(c)}-cycle:", " ".join(alg(s) for s in c))
print("    => no closed tour of the frame:", len(cyc) > 0)
