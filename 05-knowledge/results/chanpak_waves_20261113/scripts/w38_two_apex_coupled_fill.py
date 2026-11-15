#!/usr/bin/env python3
"""Wave 38: does a coupled (bilinear) two-module repair still fill its window?

Inflate two spine positions of a two-apex source by order-7 modules X1, X2 ranging over two whole H-fibres
(all order-7 classes with the given H).  The reader output F = H + Q of the inflated source is recorded for every
pair.  We compare
  * the exact repair set R = {F(X1, X2) - F_min} with the additive model R_add = {dF1(X1) + dF2(X2)}
    (single-swap effects measured at the base pair), and
  * their fills in units of the common step g = gcd of differences.
Usage: python3 w38_two_apex_coupled_fill.py HQ_BINARY
"""
import random
import sys
from math import gcd

sys.argv = sys.argv[:2]
import importlib.util
spec = importlib.util.spec_from_file_location("add", __file__.replace("w38_two_apex_coupled_fill.py", "w38_two_apex_additivity.py"))
# reuse helpers without running main(): read the source and exec everything but the final main() call
src = open(spec.origin).read().rsplit("main()", 1)[0]
ns = {"__name__": "helpers"}
exec(compile(src, spec.origin, "exec"), ns)
hq_batch, gentourng, ports, two_apex = ns["hq_batch"], ns["gentourng"], ns["ports"], ns["two_apex"]

rng = random.Random(38)
mods = gentourng(7)
byH = {}
for out in mods:
    pv = ports(out)
    byH.setdefault(pv[0], []).append((pv, out))
# fibres with the most distinct port vectors
fib = sorted(byH.items(), key=lambda t: -len({tuple(p) for p, _ in t[1]}))
(h1, F1), (h2, F2) = fib[0], fib[1]
print(f"order-7 module fibres: H = {h1} ({len(F1)} classes), H = {h2} ({len(F2)} classes)")


def fill(vals):
    vals = sorted(set(vals))
    g = 0
    for v in vals:
        g = gcd(g, v - vals[0])
    g = g or 1
    span = (vals[-1] - vals[0]) // g + 1
    return g, len(vals), span, len(vals) / span


for trial in range(3):
    m = 6
    a = (rng.randrange(1 << m) | 1) & ~(1 << (m - 1))
    b = rng.randrange(1, 1 << m)
    p1, p2 = sorted(rng.sample(range(m), 2))
    tours = [two_apex(m, a, b, {p1: X1, p2: X2}) for (_, X1) in F1 for (_, X2) in F2]
    res = hq_batch(tours)
    Fv = [h + q for (h, q) in res]
    n2 = len(F2)
    base = Fv[0]
    d1 = [Fv[i * n2] - base for i in range(len(F1))]
    d2 = [Fv[j] - base for j in range(n2)]
    exact = [f - base for f in Fv]
    add = [x + y for x in d1 for y in d2]
    mixed = [exact[i * n2 + j] - d1[i] - d2[j] for i in range(len(F1)) for j in range(n2)]
    ge, ke, se, fe = fill(exact)
    ga, ka, sa, fa = fill(add)
    print(f"m={m} a={a:06b} b={b:06b} p=({p1},{p2}): exact set: step {ge}, {ke} values over {se} steps, fill {fe:.3f}; "
          f"additive model: step {ga}, {ka} over {sa}, fill {fa:.3f}; "
          f"max |mixed| / max |single| = {max(map(abs, mixed)) / max(max(map(abs, d1)), max(map(abs, d2))):.3f}")
