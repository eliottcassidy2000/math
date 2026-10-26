# Coarsest closed attachment data for rooted two-port states under the re-mark/complement dynamics.
# States: all rooted classes [C; s, t] (s != t) of tournaments on n vertices.  Successor multiset: the rooted classes
# [C^op xor E(P); s, t] over the Hamilton paths P of C from s to t.  For an initial observation O (a function of the
# state) we compute the coarsest partition refining O that is stable: equivalent states have equal multisets of
# successor cells.  This is the minimal attachment datum that determines O at every future step of any re-mark /
# complement sequence (a bisimulation quotient).  Observations:
#   free   : the unrooted class of C               sand  : H of the singleton sandwich e * h * e (hinge gluing)
#   oneside: (class with s marked, class with t marked)    triv  : nothing (the intrinsic quotient)
# Usage: bisim.py n.  mac-mini-2026-10-25-pyth (2026-10-26).
import sys, time, itertools
from collections import Counter, defaultdict
from rooted import tournaments, canon_rooted, canon_free, ham_paths, partner, is_strong, H

def sandwich(adj):
    """H of e * h * e for the rooted state (adj; 0, n-1): new vertices a (before) and b (after); all cross arcs go from
    the later block to the earlier one, except the hinges a -> 0 and (n-1) -> b"""
    n = len(adj); N = n + 2
    a, b = n, n + 1
    new = [0] * N
    for u in range(n): new[u] = adj[u]
    for u in range(n):
        if u == 0: new[a] |= 1 << 0            # hinge a -> s
        else: new[u] |= 1 << a                 # C -> a
        if u == n - 1: new[u] |= 1 << b        # hinge t -> b
        else: new[b] |= 1 << u                 # b -> C
    new[b] |= 1 << a                           # b -> a (later to earlier)
    return len(ham_paths(tuple(new)))

def one_side(adj):
    """(class rooted at entrance only, class rooted at exit only)"""
    n = len(adj)
    def canon1(r):
        best = None
        rest = [v for v in range(n) if v != r]
        for mp in itertools.permutations(range(1, n)):
            perm = [0] * n; perm[r] = 0
            for v, w in zip(rest, mp): perm[v] = w
            from rooted import relabel
            x = relabel(adj, perm)
            if best is None or x < best: best = x
        return best
    return (canon1(0), canon1(n - 1))

def refine(states, succ, init):
    cell = {x: init[x] for x in states}
    # relabel cells to ints
    def norm(c):
        keys = {}; return {x: keys.setdefault(c[x], len(keys)) for x in states}
    cell = norm(cell)
    while True:
        sig = {x: (cell[x], tuple(sorted(Counter(cell[y] for y in succ[x]).items()))) for x in states}
        new = norm(sig)
        if len(set(new.values())) == len(set(cell.values())): return new
        cell = new

if __name__ == '__main__':
    n = int(sys.argv[1]); t0 = time.time()
    classes = tournaments(n)
    states = {}
    for adj in classes:
        for s in range(n):
            for t in range(n):
                if s != t:
                    c = canon_rooted(adj, s, t)
                    states.setdefault(c, None)
    states = list(states)
    succ = {}
    for x in states:
        succ[x] = [canon_rooted(partner(x, P), 0, n - 1) for P in ham_paths(x, 0, n - 1)]
    print(f"n = {n}: {len(classes)} tournament classes, {len(states)} rooted two-port states, "
          f"{sum(len(v) for v in succ.values())} re-mark/complement transitions [{time.time() - t0:.1f}s]")
    obs = {'triv': {x: 0 for x in states},
           'free': {x: canon_free(x) for x in states},
           'sand': {x: sandwich(x) for x in states}}
    if n <= 6: obs['oneside'] = {x: one_side(x) for x in states}
    full = len(states)
    for name, init in obs.items():
        c0 = len(set(init.values()))
        cell = refine(states, succ, init)
        k = len(set(cell.values()))
        print(f"  observation {name:7s}: {c0:5d} initial cells -> {k:5d} stable cells (of {full} rooted states)"
              + ("   = full joint rooted classes (no compression)" if k == full else ""))
    # does the joint rooted class carry more than (one-sided classes, free class)?
    if n <= 6:
        pair = defaultdict(set)
        for x in states: pair[(obs['oneside'][x], obs['free'][x])].add(x)
        amb = [v for v in pair.values() if len(v) > 1]
        print(f"  (one-sided classes, free class) leaves {len(amb)} fibres with several joint states; "
              f"sandwich H separates {sum(1 for v in amb if len({obs['sand'][x] for x in v}) > 1)} of them")
    print(f"  [{time.time() - t0:.1f}s]")
