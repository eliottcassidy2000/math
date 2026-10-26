# Coarsest closed attachment data (bisimulation quotients) for rooted two-port tournament states, n <= 7.
# Same definitions as bisim2.py; successors are computed once per tournament class (every Hamilton path P of C gives
# the rooted target [C^op xor E(P); first P, last P]); free canonical forms via nauty labelg (canon_nauty.py).
#   V1 (fixed ports): [C; s, t] -> targets of the paths of C from s to t
#   V2 (free re-marking): [C; s, t] -> targets of all paths of C
# Observations: free (unrooted class), sand (H of e*h*e), free+sand, glue (free class of h*k and k*h for every rooted
# state k of order <= 3), and blue-aware variants are left to later scripts.  Usage: bisim3.py n.
import sys, time, json
from collections import Counter, defaultdict
from rooted import tournaments, canon_rooted, ham_paths, partner, is_strong
from bisim import sandwich, refine
from bisim2 import star, small_rooted
from canon_nauty import canon_batch

n = int(sys.argv[1]); t0 = time.time()
classes = tournaments(n)
state_of = {}            # (class index, s, t) -> rooted canonical state
states = set()
succ1 = defaultdict(list); succ2 = defaultdict(list)
for ci, adj in enumerate(classes):
    pths = ham_paths(adj)
    targets = [canon_rooted(partner(adj, P), P[0], P[-1]) for P in pths]
    roots = {}
    for s in range(n):
        for t in range(n):
            if s != t:
                x = canon_rooted(adj, s, t); roots[(s, t)] = x; states.add(x)
    # V1: a path from s to t contributes to the rooted state (adj; s, t); different labelled (s,t) pairs that give the
    # same rooted class contribute identical successor lists (automorphism images), so record each rooted class once
    seen = set()
    by_pair = defaultdict(list)
    for P, tg in zip(pths, targets): by_pair[(P[0], P[-1])].append(tg)
    for (s, t), x in roots.items():
        if x in seen: continue
        seen.add(x)
        succ1[x] = sorted(by_pair.get((s, t), []))
        succ2[x] = sorted(targets)          # free re-marking: all paths of C (does not depend on the ports)
states = sorted(states)
t1 = time.time()
strong = {x: is_strong(x) for x in states}
free = dict(zip(states, canon_batch(states)))
ctx = small_rooted(3)
glue_graphs = []
for x in states:
    for k in ctx:
        glue_graphs.append(star_adj := None) if False else None
# build gluing adjacencies without canonical rooting (only the free class is needed)
def star_raw(h, k):
    n1, n2 = len(h), len(k); N = n1 + n2
    adj = [0] * N
    for u in range(n1): adj[u] = h[u]
    for u in range(n2): adj[n1 + u] = k[u] << n1
    for u in range(n1):
        for v in range(n2):
            if u == n1 - 1 and v == 0: adj[u] |= 1 << (n1 + v)
            else: adj[n1 + v] |= 1 << u
    return tuple(adj)
raw = [star_raw(x, k) for x in states for k in ctx] + [star_raw(k, x) for x in states for k in ctx]
canon = canon_batch(raw)
m = len(states) * len(ctx)
glue = {}
for i, x in enumerate(states):
    right = tuple(canon[i * len(ctx):(i + 1) * len(ctx)]); left = tuple(canon[m + i * len(ctx): m + (i + 1) * len(ctx)])
    glue[x] = (free[x], right, left)
sand = {x: sandwich(x) for x in states}
obs = {'free': free, 'sand': sand, 'free+sand': {x: (free[x], sand[x]) for x in states}, 'glue': glue}
print(f"n = {n}: {len(classes)} classes, {len(states)} rooted states ({sum(strong.values())} strong); "
      f"V1 transitions {sum(len(succ1[x]) for x in states)}, V2 transitions {sum(len(succ2[x]) for x in states)}; "
      f"{len(ctx)} gluing contexts  [{t1 - t0:.0f}s + {time.time() - t1:.0f}s]")
result = {}
for dyn, succ in (('V1', succ1), ('V2', succ2)):
    for name, init in obs.items():
        cell = refine(states, succ, init)
        k_all = len(set(cell.values())); k_str = len({cell[x] for x in states if strong[x]})
        c0 = len(set(init.values()))
        result[f'{dyn}/{name}'] = (c0, k_all, k_str)
        print(f"  {dyn} obs {name:9s}: {c0:6d} initial -> {k_all:6d} stable cells (strong states: {k_str} of {sum(strong.values())})"
              + ("  = no compression" if k_all == len(states) else ""), flush=True)
json.dump({'n': n, 'states': len(states), 'strong': sum(strong.values()), 'result': result}, open(f'bisim3_n{n}.json', 'w'))
print(f"  [{time.time() - t0:.0f}s]")
