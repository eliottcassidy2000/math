# Transitive-plus-one family: does every vertex lie on a FORWARD BALANCED arc (1/2 < p(tail before head) <= 2/3) that lies on
# a directed 3-cycle (rule R_both of HYP-9279)?  Exact, via the crossing decomposition of THM-4655 (classes() of tpo_cycle_arcs).
from fractions import Fraction as Fr
from itertools import product
from tpo_cycle_arcs import classes
import sys
def pair_probs(w):
    m = len(w); cl = classes(w); H = sum(c[0] for c in cl)
    pA = {i: sum(wt * q[i] for wt, q in cl) / H for i in range(1, m + 1)}                     # P(u_i before z)
    # P(u_j before u_i), i < j : u_j in A and u_i in B (coordinates independent within a class)
    rev = {(i, j): sum(wt * q[j] * (1 - q[i]) for wt, q in cl) / H for i in range(1, m + 1) for j in range(i + 1, m + 1)}
    return pA, rev
M = int(sys.argv[1]) if len(sys.argv) > 1 else 12
for m in range(2, M + 1):
    fails = 0; ex = None
    for mid in product((0, 1), repeat=max(m - 2, 0)):
        w = (1,) + mid + (0,)
        pA, rev = pair_probs(w)
        # arcs as (tail, head, p(tail before head)); vertex 0 = z, i = u_i
        arcs = []
        for i in range(1, m + 1):
            if w[i - 1] == 1: arcs.append((0, i, 1 - pA[i]))
            else: arcs.append((i, 0, pA[i]))
            for j in range(i + 1, m + 1): arcs.append((i, j, 1 - rev[(i, j)]))
        out = {(t, h) for t, h, _ in arcs}
        def on_triangle(t, h):   # arc t->h on a 3-cycle t->h->y->t
            return any((h, y) in out and (y, t) in out for y in range(m + 1) if y not in (t, h))
        good = set()
        for t, h, p in arcs:
            if Fr(1, 2) < p <= Fr(2, 3) and on_triangle(t, h): good |= {t, h}
        if len(good) < m + 1:
            fails += 1; ex = ex or (''.join(map(str, w)), sorted(set(range(m + 1)) - good))
    print(f"m={m} (n={m + 1}): words violating R_both: {fails}", (f"e.g. w={ex[0]} vertices {ex[1]}" if ex else ""))
