"""Decoration rule: a terminal leaf acts as a loop with the vector 1 in place of the loop diagonal d.
Session opus-2026-10-22-S54.

Square identity (two arms): dumbbell + lollipop + path tree = -(d + 1)^T G (d + 1)   (leaf_closed.py).
Claw test (three arms): the cores consisting of a centre joined by three bridges to three ends, each end a leaf or a
vertex carrying one loop, summed over the four leaf counts, against the Y-tree formula evaluated at u = d + 1:
  Y(u) = -(2/3) sum_w [ (G u)_w^3 + 3 (G u)_w (Q u)_w^2 ].
Per-core values come from leaf_engine.py (validated against brute force).  Usage: python decoration_rule.py"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from cluster_brute_classes import almost_regular
from leaf_closed import mats


def is_claw(core, k):
    """centre of degree 3 joined by three single edges to ends; an end is a leaf (index < k) or carries one loop."""
    E = core.edges; V = core.V
    loops = [a for a, b in E if a == b]
    deg = [0] * V
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    for c in range(V):
        nb = [b if a == c else a for a, b in E if (a == c) != (b == c)]
        if deg[c] == 3 and len(nb) == 3 and len(set(nb)) == 3 and c not in loops:
            ends = set(nb)
            if all((e < k and deg[e] == 1) or (e >= k and deg[e] == 3 and loops.count(e) == 1) for e in ends) and V == 4:
                return True
    return False


rng = np.random.default_rng(5)
for n in (8, 12, 16):
    A = almost_regular(n, rng)
    Sh = (A - A.T).astype(float) / n
    P = LE.LEF.stack(Sh.astype(complex), 1.0)
    G, H, O, Q, d, g, q, one = mats(Sh)
    total, parts = 0.0, {}
    for (e, k) in ((2, 0), (1, 1), (0, 2), (-1, 3)):
        for core in LE.get(e, k):
            if is_claw(core, k):
                v = core.value(P).real
                parts[(e, k)] = parts.get((e, k), 0.0) + v
                total += v
    u = d + one
    Yu = -(2 / 3) * np.sum((G @ u) ** 3 + 3 * (G @ u) * (Q @ u) ** 2)
    print(f'n={n}: claw cores by (excess, leaves): ' + ', '.join(f'{c}: {v:+.6e}' for c, v in sorted(parts.items())) +
          f'; sum {total:+.10e}; Y(d + 1) {Yu:+.10e}; diff {total - Yu:+.2e}', flush=True)
