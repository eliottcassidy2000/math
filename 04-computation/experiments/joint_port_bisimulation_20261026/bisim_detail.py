# Where the separation comes from (strong rooted states, n <= 7):
#   (a) the gluing observation alone (contexts of order <= 3, no dynamics);
#   (b) after closing under the fixed-port re-mark/complement dynamics V1;
#   (c) one-sided data: pairs of strong joint states with equal (entrance-rooted class, exit-rooted class) -- the
#       289/253 phenomenon -- and how many of those the gluing contexts / the dynamics separate.
import sys, time, itertools
from collections import defaultdict
exec(open('bisim3.py').read().split("print(f\"n = {n}: {len(classes)} classes")[0])
from rooted import relabel
def canon1(adj, r):
    nn = len(adj); best = None
    rest = [v for v in range(nn) if v != r]
    for mp in itertools.permutations(range(1, nn)):
        perm = [0] * nn; perm[r] = 0
        for v, w in zip(rest, mp): perm[v] = w
        x = relabel(adj, perm)
        if best is None or x < best: best = x
    return best
S = [x for x in states if strong[x]]
c_obs = len({glue[x] for x in S})
cell = refine(states, succ1, glue)
c_dyn = len({cell[x] for x in S})
one = defaultdict(list)
for x in S: one[(canon1(x, 0), canon1(x, n - 1))].append(x)
coll = [v for v in one.values() if len(v) > 1]
pairs = sum(len(v) * (len(v) - 1) // 2 for v in coll)
sep_obs = sum(1 for v in coll for a, b in itertools.combinations(v, 2) if glue[a] != glue[b])
sep_dyn = sum(1 for v in coll for a, b in itertools.combinations(v, 2) if cell[a] != cell[b])
print(f"n = {n}: {len(S)} strong joint states; gluing observation alone separates them into {c_obs} cells, "
      f"with V1 dynamics {c_dyn}; one-sided collisions: {len(coll)} fibres, {pairs} pairs, of which the gluing "
      f"observation separates {sep_obs} and the closed data {sep_dyn}")
