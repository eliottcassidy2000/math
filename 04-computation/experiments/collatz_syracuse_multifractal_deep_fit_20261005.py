#!/usr/bin/env python3
"""Fit the deficits d_n(q) = bound(q) - slope_n(q) of collatz_syracuse_multifractal_deep_20261005.json.
Three readings: (i) log-log exponent of d_n over the last window; (ii) the two-parameter model d_n = d_inf + c/n
through the window endpoints (d_inf <= 0 means the data force equality); (iii) the model d_n = c/(n + c'), i.e.
1/d_n linear in n: the least-squares line of 1/d_n against n, its residual, and the implied d_inf = 0 test
(a positive limit d_inf would make 1/d_n concave and saturate at 1/d_inf)."""
import json, math, sys, os
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "collatz_syracuse_multifractal_deep_20261005.json")
rec = json.load(open(p))
qs = [0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0, 4.0, 6.0, 8.0]
def bound(q): return min(q - 1, math.log(2 ** q - 1) / math.log(3)) if q > 1 else q - 1
ns = sorted(int(k) for k in rec if "slopes" in rec[k])
nmax = ns[-1]
win = [n for n in ns if n >= nmax - 7]
out = []
def P(s): print(s); out.append(s)
P(f"levels with slopes: {ns[0]}..{nmax}; fit window {win[0]}..{nmax}")
P("  q     d_last    loglog e (d~n^-e)   endpoint 1/n model d_inf (c)   1/d linear fit: c, c', max rel. residual, second-difference sign")
for i, q in enumerate(qs):
    if q <= 1: continue
    d = {n: bound(q) - rec[str(n)]["slopes"][i] for n in ns}
    if any(d[n] <= 0 for n in win):
        P(f"  {q:<5} {d[nmax]:+.4f}   (deficit nonpositive in the window)"); continue
    xs = [math.log(n) for n in win]; ys = [math.log(d[n]) for n in win]
    xm = sum(xs) / len(xs); ym = sum(ys) / len(ys)
    e = -sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sum((x - xm) ** 2 for x in xs)
    n1, n2 = win[0], nmax
    c = (d[n1] - d[n2]) / (1 / n1 - 1 / n2); dinf = d[n2] - c / n2
    # 1/d linear fit
    X = win; Y = [1 / d[n] for n in win]
    Xm = sum(X) / len(X); Ym = sum(Y) / len(Y)
    sl = sum((x - Xm) * (y - Ym) for x, y in zip(X, Y)) / sum((x - Xm) ** 2 for x in X)
    ic = Ym - sl * Xm
    cc = 1 / sl; cp = ic * cc
    res = max(abs((sl * x + ic) - y) / y for x, y in zip(X, Y))
    sd = [Y[k + 2] - 2 * Y[k + 1] + Y[k] for k in range(len(Y) - 2)]
    sdsign = "concave (saturating)" if all(v < 0 for v in sd) else ("convex" if all(v > 0 for v in sd) else "mixed/linear")
    P(f"  {q:<5} {d[nmax]:+.4f}   e = {e:.2f}             d_inf = {dinf:+.4f} (c = {c:.3f})       c = {cc:.3f}, c' = {cp:+.2f}, resid {res:.1%}, {sdsign}")
with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
