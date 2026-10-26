# Coarsest closed attachment data, two dynamics and several observation sets.
#   V1 (fixed ports): [C; s, t] -> [C^op xor E(P); s, t] for Hamilton paths P of C from s to t;
#   V2 (free re-marking): [C; s, t] -> [C^op xor E(P); first P, last P] for ALL Hamilton paths P of C.
# Observations (the data must determine them at every future step):
#   free  = unrooted class;  sand = H(e*h*e);  glue = free class plus the classes of h*k and k*h for every rooted state
#   k of order 1..3 (all joint gluing contexts up to order 3, hinge product of the rooted-lines note);
#   joint = the rooted class itself (no compression possible: the trivial upper bound).
# Reports the number of stable cells, on all states and on the strong states.  Usage: bisim2.py n.
import sys, time, itertools
from collections import Counter, defaultdict
from rooted import tournaments, canon_rooted, canon_free, ham_paths, partner, is_strong, H
from bisim import sandwich, refine

def star(h, k):
    """hinge product of rooted states h = (C; 0, |C|-1), k = (D; 0, |D|-1): every cross arc D -> C except the hinge
    exit(C) -> entrance(D); ports (entrance of C, exit of D).  Returned as a canonical rooted state."""
    n1, n2 = len(h), len(k); N = n1 + n2
    adj = [0] * N
    for u in range(n1): adj[u] = h[u]
    for u in range(n2): adj[n1 + u] = k[u] << n1
    for u in range(n1):
        for v in range(n2):
            if u == n1 - 1 and v == 0: adj[u] |= 1 << (n1 + v)          # hinge
            else: adj[n1 + v] |= 1 << u                                  # D -> C
    return canon_rooted(tuple(adj), 0, N - 1)

def small_rooted(maxn):
    res = [((0,),)]  # order 1: single vertex, ports coincide: represent as adjacency (0,)
    out = [(0,)]
    for m in range(2, maxn + 1):
        seen = set()
        for adj in tournaments(m):
            for s in range(m):
                for t in range(m):
                    if s != t: seen.add(canon_rooted(adj, s, t))
        out += sorted(seen)
    return out

if __name__ == '__main__':
    n = int(sys.argv[1]); t0 = time.time()
    classes = tournaments(n)
    states = sorted({canon_rooted(adj, s, t) for adj in classes for s in range(n) for t in range(n) if s != t})
    paths = {x: ham_paths(x) for x in states}
    succ1 = {x: [canon_rooted(partner(x, P), 0, n - 1) for P in paths[x] if P[0] == 0 and P[-1] == n - 1] for x in states}
    succ2 = {x: [canon_rooted(partner(x, P), P[0], P[-1]) for P in paths[x]] for x in states}
    strong = {x: is_strong(x) for x in states}
    ctx = small_rooted(3)
    obs = {'free': {x: canon_free(x) for x in states}}
    obs['sand'] = {x: sandwich(x) for x in states}
    obs['free+sand'] = {x: (obs['free'][x], obs['sand'][x]) for x in states}
    obs['glue'] = {x: (obs['free'][x], tuple(canon_free(star(x, k)) for k in ctx), tuple(canon_free(star(k, x)) for k in ctx)) for x in states}
    print(f"n = {n}: {len(classes)} classes, {len(states)} rooted states ({sum(strong.values())} strong), "
          f"V1 transitions {sum(map(len, succ1.values()))}, V2 transitions {sum(map(len, succ2.values()))}; "
          f"{len(ctx)} gluing contexts of order <= 3  [{time.time() - t0:.0f}s]")
    for dyn, succ in (('V1 fixed ports', succ1), ('V2 free re-marking', succ2)):
        for name, init in obs.items():
            cell = refine(states, succ, init)
            k_all = len(set(cell.values())); k_str = len({cell[x] for x in states if strong[x]})
            print(f"  {dyn:19s} obs {name:9s}: {len(set(init.values())):5d} -> {k_all:5d} stable cells "
                  f"(strong states: {k_str} cells of {sum(strong.values())})" + ("  = no compression" if k_all == len(states) else ""))
    print(f"  [{time.time() - t0:.0f}s]")
