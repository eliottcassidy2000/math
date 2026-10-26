# Two-sided gluing contexts.  For rooted states k, k' of order <= 2 (three states: the singleton and the two rooted
# 2-vertex states), the free class of the two-sided gluing k * h * k' uses both ports of h jointly (the singleton
# sandwich e * h * e of the 289/253 hostile is one of them).  Questions, strong joint states, n <= 7:
#   (a) do the two-sided contexts alone separate all strong joint states?
#   (b) one-sided collisions: separated by the two-sided contexts?
import sys, time, itertools
from collections import defaultdict
exec(open('bisim3.py').read().split("print(f\"n = {n}: {len(classes)} classes")[0])
exec(open('bisim_detail.py').read().split("S = [x for x in states if strong[x]]")[0].split("from rooted import relabel")[1].join(['from rooted import relabel', '']) if False else '')
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
ctx2 = small_rooted(2)
S = [x for x in states if strong[x]]
raw2 = [star_raw(star_raw(k, x), k2) for x in S for k in ctx2 for k2 in ctx2]
c2 = canon_batch(raw2)
w = len(ctx2) ** 2
two = {x: tuple(c2[i * w:(i + 1) * w]) for i, x in enumerate(S)}
one = defaultdict(list)
for x in S: one[(canon1(x, 0), canon1(x, n - 1))].append(x)
coll = [v for v in one.values() if len(v) > 1]
pairs = [(a, b) for v in coll for a, b in itertools.combinations(v, 2)]
print(f"n = {n}: {len(S)} strong joint states; two-sided contexts (order <= 2, {w} pairs) alone separate them into "
      f"{len(set(two.values()))} cells; one-sided collision pairs {len(pairs)}, separated by two-sided contexts: "
      f"{sum(1 for a, b in pairs if two[a] != two[b])}; by the singleton sandwich H alone: {sum(1 for a, b in pairs if sand[a] != sand[b])}")
