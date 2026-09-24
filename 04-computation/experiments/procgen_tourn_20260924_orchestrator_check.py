#!/usr/bin/env python3
"""Orchestrator's independent check of THM-4472 (four-vertex reading of 3n+-1), written without reading the lane's code.
For each sheet b = +-1 and side of 0, the AM-fair pair {s_o, s_e = s_o + b} and its images, with map arcs plus
'smaller -> larger' order arcs: forward class, H, and the inverse-tree (backward) class, for |s_o| = 5..119.
Asserts: forward = (0,2,2,2) with H = 3 iff b*s_o > 0, transitive with H = 1 otherwise; backward = converse class;
the reflection rho(x) = s_o + s_e - x maps the forward tournament onto the converse of the backward one.
"""
from itertools import permutations
def T(n, b): return n // 2 if n % 2 == 0 else (3 * n + b) // 2
def tourn(V, maparcs):
    arcs = set(maparcs)
    for u in V:
        for v in V:
            if u < v and (u, v) not in arcs and (v, u) not in arcs: arcs.add((u, v))
    return arcs
def scores(V, arcs): return tuple(sorted(sum(1 for (u, v) in arcs if u == x) for x in V))
def H(V, arcs): return sum(1 for P in permutations(V) if all((P[k], P[k+1]) in arcs for k in range(3)))
CONV = {(0,2,2,2): (1,1,1,3), (1,1,1,3): (0,2,2,2), (0,1,2,3): (0,1,2,3), (1,1,2,2): (1,1,2,2)}
n = 0
for b in (1, -1):
    for side in (1, -1):
        for i in range(2, 60):
            s_o = side * (2 * i + 1); s_e = s_o + b
            V = (s_o, s_e, T(s_o, b), T(s_e, b))
            if len(set(V)) < 4: continue
            assert T(s_o, b) + T(s_e, b) == s_o + s_e
            fwd = tourn(V, {(s_o, T(s_o, b)), (s_e, T(s_e, b))})
            bwd = tourn(V, {(T(s_o, b), s_o), (T(s_e, b), s_e)})
            sf, sb = scores(V, fwd), scores(V, bwd)
            if b * s_o > 0: assert sf == (0,2,2,2) and H(V, fwd) == 3
            else:           assert sf == (0,1,2,3) and H(V, fwd) == 1
            assert sb == CONV[sf]
            rho = lambda x: s_o + s_e - x
            assert set(map(rho, V)) == set(V)
            assert {(rho(u), rho(v)) for (u, v) in fwd} == {(v, u) for (u, v) in bwd}   # rho(F) = converse(B)
            n += 1
print(f"{n} quadruples: forward (0,2,2,2)/H=3 iff b*s_o>0 else transitive/H=1; backward = converse; rho(F) = converse(B): ok")
print("ALL CHECKS PASSED")
