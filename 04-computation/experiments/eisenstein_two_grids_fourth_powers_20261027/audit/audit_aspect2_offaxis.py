"""HYP-9428 (3) title claim 'excess 2 off the axis': tabulate certified excess values (4g + 7 <= W) on L_2, W = 11263."""
from collections import Counter, defaultdict
from knuth_debt import knuth_debt
Qf = lambda u, v: u*u + 2*v*v
W = 11263
R = int(W ** 0.5) + 2
pts = [(u, v) for u in range(-R, R + 1) for v in range(-R, R + 1) if (u - v) % 2 == 0 and Qf(u, v) <= W]
g, Q, _ = knuth_debt(pts, [[1, 0], [0, 2]], 1)
cert = [k for k in range(len(pts)) if 4*g[k] + 7 <= W]
allbelow = max(q for q in range(W) if all(4*g[k] + 7 <= W for k in range(len(pts)) if Q[k] <= q))
print(f"W = {W}: {len(pts)} points; every point of norm <= {allbelow} certified; {len(cert)} certified points")
on = Counter(); off = Counter(); first = defaultdict(list)
for k in cert:
    e = int(g[k]) - 2*(int(Q[k])//8)
    u, v = pts[k]
    if v == 0: on[e] += 1
    else: off[e] += 1
    if e >= 4 and u >= 0 and v >= 0: first[(e, v == 0)].append((int(Q[k]), (u, v)))
print("on-axis excess histogram:", dict(sorted(on.items())), "; off-axis:", dict(sorted(off.items())))
for key in sorted(first):
    print(f"   excess {key[0]} {'on' if key[1] else 'off'} axis, first-quadrant points (norm, point), smallest 12:", sorted(first[key])[:12])
