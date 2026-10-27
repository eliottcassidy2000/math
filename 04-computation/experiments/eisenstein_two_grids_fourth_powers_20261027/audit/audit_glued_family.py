"""Audit THM-4824 (B), (C): collinear law; minimal vectors of L_c = {u = v mod 2}, Q = p u^2 + q v^2 (c = q/p);
debt ratios on the outer half of the window W = 120 m; window stability (are those windowed values exact?)."""
from fractions import Fraction as F
import numpy as np, time, sys
from knuth_debt import knuth_debt, INF

def lattice(p, q, W):
    R = int((W / min(p, q)) ** 0.5) + 2
    return [(u, v) for u in range(-R, R + 1) for v in range(-R, R + 1) if (u - v) % 2 == 0 and p*u*u + q*v*v <= W]

print("(C) minimal vectors of L_c, exact (enumeration |u|,|v| <= 12):")
for c in [F(1, 3), F(1, 2), F(1), F(3, 2), F(2), F(5, 2), F(29, 10), F(3), F(31, 10), F(7, 2), F(4)]:
    p, q = c.denominator, c.numerator
    pts = [(u, v) for u in range(-12, 13) for v in range(-12, 13) if (u - v) % 2 == 0 and (u, v) != (0, 0)]
    m = min(p*u*u + q*v*v for u, v in pts)
    mins = sorted(t for t in pts if p*t[0]**2 + q*t[1]**2 == m)
    print(f"   c = {str(c):5s}: minimal norm {F(m, p)} (u^2 + c v^2), {len(mins)} vectors {mins}")

print("(B) collinear law: first grid Z[sqrt-3] (Q = a^2 + 3b^2, units +-1) and L_4 (units +-(2,0))")
pts = [(a, b) for a in range(-40, 41) for b in range(-12, 13) if a*a + 3*b*b <= 1600]
g, Q, _ = knuth_debt(pts, [[1, 0], [0, 3]], 1)
idx = {t: k for k, t in enumerate(pts)}
reach = [t for k, t in enumerate(pts) if g[k] < INF]
print("   reachable points all on the real line:", all(t[1] == 0 for t in reach), "; g(n) = n(n-1) for |n| <= 40:",
      all(g[idx[(n, 0)]] == n*n - abs(n) for n in range(-40, 41)))
pts = lattice(1, 4, 2000)
g, Q, _ = knuth_debt(pts, [[1, 0], [0, 4]], 1)
idx = {t: k for k, t in enumerate(pts)}
reach = [t for k, t in enumerate(pts) if g[k] < INF]
print("   L_4: reachable only (2n, 0):", all(t[1] == 0 and t[0] % 2 == 0 for t in reach), "; g(2n,0) = 4(n^2-|n|):",
      all(g[idx[(2*n, 0)]] == 4*(n*n - abs(n)) for n in range(-22, 23)))
# potential inequality for the collinear law, exhaustive in a range
ok = all((a+b)**2 - abs(a+b) - (a*a - abs(a)) - (b*b - abs(b)) <= max(2*a*b, 0) for a in range(-200, 201) for b in range(-200, 201))
print("   potential (k^2-|k|) satisfies Phi(a+b)-Phi(a)-Phi(b) <= (2ab)_+ for |a|,|b| <= 200:", ok)

print("(C) debt ratios g_W/Q on the outer half Q > W/2, W = 120 m, and window stability")
for (p, q) in ((1, 1), (2, 3), (1, 2), (2, 5), (1, 3)):
    m = min(p + q, 4*p)       # (1,1) or (2,0)
    base = None
    for mult in (120, 240, 480):
        W = mult * m
        t0 = time.time()
        pts = lattice(p, q, W)
        g, Q, _ = knuth_debt(pts, [[p, 0], [0, q]], 1)
        assert Q[Q > 0].min() == m
        idx = {t: k for k, t in enumerate(pts)}
        outer = [k for k in range(len(pts)) if Q[k] > W // 2 and g[k] < INF]
        r = [g[k] / Q[k] for k in outer]
        line = f"   c = {q}/{p}, W = {mult}m = {W}: {len(pts)} pts, unreachable {int((g >= INF).sum())}, outer-half g/Q in [{min(r):.4f}, {max(r):.4f}]"
        if base is None:
            base = {t: int(g[idx[t]]) for t in pts}
            base_outer = [t for t in pts if Q[idx[t]] > W // 2]
            W0 = W
        else:
            changed = [(t, base[t], int(g[idx[t]])) for t in base_outer if base[t] != g[idx[t]]]
            line += f"; of the {len(base_outer)} outer-half points of W0 = {W0}, {len(changed)} have smaller g at W = {W}"
            if changed:
                drop = max(b - n for _, b, n in changed)
                line += f" (max drop {drop}); new outer-half-of-W0 range [{min(int(g[idx[t]])/Q[idx[t]] for t in base_outer):.4f}, {max(int(g[idx[t]])/Q[idx[t]] for t in base_outer):.4f}]"
        print(line + f"  [{time.time()-t0:.1f}s]", flush=True)
