#!/usr/bin/env python3
"""Class-level invariants: for all tournament classes on n vertices (from electric_n{n}.json representatives), count collisions of
(scores), (scores, curl profile), (scores, curl profile, second scores); report the colliding groups at n = 5."""
import sys, json
from collections import Counter, defaultdict
from determination_L2 import invariants
from determination import tournament_from_mask
for n in [int(a) for a in sys.argv[1:]] or [5, 6, 7]:
    data = json.load(open(f'electric_n{n}.json')); nodes = data['nodes']
    classes = []
    for u in nodes:
        classes.append(u['rep_mask'])
        if not u['SC'] and u['conv_mask'] is not None: classes.append(u['conv_mask'])
    for level in (1, 2, 3):
        groups = defaultdict(list)
        for mk in classes: groups[invariants(n, tournament_from_mask(n, mk), level)].append(mk)
        coll = {k: len(v) for k, v in groups.items() if len(v) > 1}
        print(f"n = {n}, level L{level}: {len(classes)} classes, {len(groups)} distinct invariant values, colliding groups {len(coll)} (sizes {sorted(coll.values(), reverse=True)[:8]})")
        if n == 5 and level == 1:
            for k, v in groups.items():
                if len(v) > 1: print("   collision scores", k[0], "-> curl profiles", [invariants(5, tournament_from_mask(5, m), 2)[1] for m in v])
