"""THM-4824 (D): the aspect-2 glued grid (u = v mod 2, Q = u^2 + 2v^2, four unit charges (+-1, +-1) of norm 3).
Lower bound g >= 2 floor(Q/8) (potential; THM-4824 (D)).  Window lemma: a tree of cost C keeps every node at
2 floor(Q/8) <= C, i.e. Q <= 4C + 7, so window values with 4 g + 7 <= W are exact.  Lists the exceptional norms."""
from root_lattice_debt import relax
W = 900
pts = [(u, v) for u in range(-70, 71) for v in range(-70, 71) if (u - v) % 2 == 0]
Q = lambda p: p[0]*p[0] + 2*p[1]*p[1]
units = [p for p in pts if Q(p) == 3]
g, P = relax(pts, Q, units, W)
exact = [p for p in P if 4*g[p] + 7 <= W]
ex = {}
for p in exact:
    e = g[p] - 2*(Q(p)//8)
    assert e >= 0, p
    if e: ex.setdefault(Q(p), set()).add((e, p if p[0] >= 0 and p[1] >= 0 else None))
print(f"window Q <= {W}; points with exact values: {len(exact)} (all of norm <= {max(Q(p) for p in exact)})")
print("residues of Q mod 8 present:", sorted({Q(p) % 8 for p in P}))
rows = sorted((q, sorted({e for e, _ in s})) for q, s in ex.items())
print("exceptional norms (norm, excesses):", rows)
print("count of exceptional norms below 500:", sum(1 for q, _ in rows if q < 500), "of", len({Q(p) for p in exact if Q(p) < 500}), "norms")
