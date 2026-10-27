"""Larger windows for the free root lattices (A3 = D3, D4, A4) and an exact check of the aspect-2 law g = 2 floor(Q/8)
on the glued rhombic grid (u = v mod 2, Q = u^2 + 2 v^2)."""
import itertools, sys
from root_lattice_debt import relax
def report(name, pts, Q, W):
    units = [p for p in pts if Q(p) == min(Q(q) for q in pts if any(q))]
    g, P = relax(pts, Q, units, W)
    reach = [p for p in P if g[p] < 10**9]
    print(f"  {name}: {len(units)} minimal vectors, window Q <= {W}: {len(P)} points, reachable {len(reach)}, max g = {max(g[p] for p in reach)}")
R = 6; pts = [p for p in itertools.product(range(-R, R + 1), repeat=3) if sum(p) % 2 == 0]
report("D3 = A3", pts, lambda p: sum(a*a for a in p) // 2, 12)
R = 4; pts = [p for p in itertools.product(range(-R, R + 1), repeat=4) if sum(p) % 2 == 0]
report("D4", pts, lambda p: sum(a*a for a in p) // 2, 7)
R = 3; pts = [p for p in itertools.product(range(-R, R + 1), repeat=5) if sum(p) == 0]
report("A4", pts, lambda p: sum(a*a for a in p) // 2, 6)
pts = [(u, v) for u in range(-60, 61) for v in range(-60, 61) if (u - v) % 2 == 0]
Q = lambda p: p[0]*p[0] + 2*p[1]*p[1]
units = [p for p in pts if Q(p) == 3]
g, P = relax(pts, Q, units, 600)
bad = [(p, Q(p), g[p]) for p in P if g[p] != 2 * (Q(p) // 8)]
print(f"  aspect 2: window Q <= 600, {len(P)} points; g = 2 floor(Q/8) everywhere: {not bad}; exceptions: {bad[:12]}")
