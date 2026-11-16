#!/usr/bin/env python3
"""onparent.py n: all one-parent tournaments of order n (vertex i >= 1 beats only its parent p_i < i among the earlier
vertices), all roots; list the ready roots (defect b = O - E - g = 0) whose count is one of wave 45's seeds
13, 23, 41, 59, 75, 95, with charge a = 2E + 1 - H, one representative per (count, charge); print masks for the
port-move seed search."""
import sys, itertools
n = int(sys.argv[1]); targets = {13, 23, 41, 59, 75, 95}
exec(open('seeds.py').read().split("targets = {")[0].split("def append")[0])  # profile()
reps = {}
for parents in itertools.product(*[range(i) for i in range(1, n)]):
    masks = [0]*n
    for i in range(1, n):
        p = parents[i-1]
        for v in range(i):
            if v == p: masks[i] |= 1 << v
            else: masks[v] |= 1 << i
    for root in range(n):
        H, E, O, g = profile(masks, root)
        if H in targets and O - E - g == 0:
            a = 2*E + 1 - H
            if (H, a) not in reps: reps[(H, a)] = (parents, root, masks, (H, E, O, g))
print(f"order {n}: ready roots at the target counts, one representative per (count, charge):")
for (H, a) in sorted(reps):
    parents, root, masks, prof = reps[(H, a)]
    print(f"  h={H} charge={a}: parents {parents} root {root} (H,E,O,g)={prof} masks {' '.join(map(str, masks))}")
