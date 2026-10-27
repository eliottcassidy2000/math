"""Aspect-ratio test of the (1/2,1/2)-glued rectangular family: points (u, v) with u = v (mod 2) and the integer form
Q = p u^2 + q v^2 (aspect c = q/p).  Unit charges = minimal vectors.  Prints the asymptotic debt ratio, estimated by the
largest g/Q among points with Q above half the window, against the angle law (c - 1)/(2c) = -cos(theta)/(1 - cos(theta)),
theta the angle between the glue vectors (1, 1) and (1, -1)."""
from root_lattice_debt import relax
for p, q in ((1, 1), (2, 3), (1, 2), (2, 5), (1, 3)):
    pts = [(u, v) for u in range(-70, 71) for v in range(-70, 71) if (u - v) % 2 == 0]
    Q = lambda t, p=p, q=q: p*t[0]*t[0] + q*t[1]*t[1]
    m = min(Q(t) for t in pts if t != (0, 0)); units = [t for t in pts if Q(t) == m]
    W = 120 * m
    g, P = relax(pts, Q, units, W)
    big = [t for t in P if Q(t) > W // 2 and g[t] < 10**9]
    lo = min(g[t] / Q(t) for t in big); hi = max(g[t] / Q(t) for t in big)
    c = q / p; law = (c - 1) / (2 * c) if c >= 1 else None
    print(f"aspect c = {q}/{p}: {len(units)} minimal vectors (norm {m}); g/Q on the outer half of the window in [{lo:.4f}, {hi:.4f}]; angle law (c-1)/(2c) = {law:.4f}")
