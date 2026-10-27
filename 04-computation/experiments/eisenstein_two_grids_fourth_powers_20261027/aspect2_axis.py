"""HYP-9428 (3) test: on-axis points (2k, 0) of the aspect-2 grid; is the excess over 2 floor(Q/8) equal to 4 exactly for
k a prime = 5 (mod 8)?  Exact values by the window lemma (Q <= 4C + 7)."""
from root_lattice_debt import relax
W = 3500
pts = [(u, v) for u in range(-62, 63) for v in range(-44, 45) if (u - v) % 2 == 0]
Q = lambda p: p[0]*p[0] + 2*p[1]*p[1]
units = [p for p in pts if Q(p) == 3]
g, P = relax(pts, Q, units, W)
for k in range(1, 30):
    p = (2*k, 0)
    if p in g and 4*g[p] + 7 <= W:
        e = g[p] - 2*(Q(p)//8)
        print(f"   (2k,0), k={k}: Q={Q(p)}, g={g[p]}, excess {e}")
