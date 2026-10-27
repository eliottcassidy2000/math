"""HYP-9428 (2) probe (auditor's extension): window ratios g_W/Q on the outer half for more aspects c = q/p, at
W = 120 m and 240 m (upper bounds on g; no window lemma off c = 2, 3), against (c-1)/(2c)."""
from fractions import Fraction as F
from knuth_debt import knuth_debt, INF
for c in [F(11, 10), F(6, 5), F(5, 4), F(4, 3), F(3, 2), F(7, 4), F(2), F(9, 4), F(5, 2), F(11, 4), F(3)]:
    p, q = c.denominator, c.numerator
    m = min(p + q, 4*p)
    out = []
    for mult in (120, 240):
        W = mult * m
        R = int((W / min(p, q)) ** 0.5) + 2
        pts = [(u, v) for u in range(-R, R + 1) for v in range(-R, R + 1) if (u - v) % 2 == 0 and p*u*u + q*v*v <= W]
        g, Q, _ = knuth_debt(pts, [[p, 0], [0, q]], 1)
        r = [g[k] / Q[k] for k in range(len(pts)) if Q[k] > W // 2]
        assert all(g[k] < INF for k in range(len(pts)))
        out.append(f"W={mult}m: [{min(r):.4f}, {max(r):.4f}]")
    print(f"c = {str(c):5s} ({float(c):.3f}): " + "; ".join(out) + f"; angle law {(float(c) - 1)/(2*float(c)):.4f}")
