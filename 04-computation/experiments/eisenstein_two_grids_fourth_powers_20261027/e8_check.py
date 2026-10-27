"""Sanity check of THM-4824 (E): catalytic merge debt on E8 = D8 u (D8 + (1/2)^8) (doubled coordinates: all entries of
one parity, sum = 0 mod 4; Q = x.x/8 so that roots have Q = 1) and on A5, inside small windows."""
import itertools
from root_lattice_debt import relax
def vecs(dim, budget, parity):   # integer vectors of one parity with sum of squares <= budget
    if dim == 0: yield (); return
    vals = [a for a in range(-4, 5) if a % 2 == parity and a*a <= budget]
    for a in vals:
        for rest in vecs(dim - 1, budget - a*a, parity): yield (a,) + rest
pts = [p for par in (0, 1) for p in vecs(8, 16, par) if sum(p) % 4 == 0]
Q = lambda p: sum(a*a for a in p) // 8
units = [p for p in pts if sum(a*a for a in p) == 8]
g, P = relax(pts, Q, units, 2)
reach = [p for p in P if g[p] < 10**9]
print(f"E8 (doubled coordinates): {len(units)} roots, window Q <= 2: {len(P)} points, reachable {len(reach)}, max g = {max(g[p] for p in reach)}")
pts = [p for p in itertools.product(range(-2, 3), repeat=6) if sum(p) == 0]
Q = lambda p: sum(a*a for a in p) // 2
units = [p for p in pts if Q(p) == 1]
g, P = relax(pts, Q, units, 4)
reach = [p for p in P if g[p] < 10**9]
print(f"A5: {len(units)} roots, window Q <= 4: {len(P)} points, reachable {len(reach)}, max g = {max(g[p] for p in reach)}")
