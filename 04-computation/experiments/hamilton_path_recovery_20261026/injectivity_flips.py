# Exhaustive proof, n <= N, that the set of Hamilton paths determines a labelled tournament.
# Reduction (proved in the note): if HP(T) = HP(T'), then T and T' agree on every arc used by some Hamilton path and,
# because the strong-component cuts are the positions where all Hamilton paths share their prefix set, on every arc
# between blocks (all point forward).  So T' = T_D for a set D of intra-block arcs of T unused by HP(T), and
# HP(T_D) = HP(T) iff no Hamilton path of T_D uses a reversed arc of D.  Check: for every class T and every nonempty D
# (subset of T's unused intra-block arcs), T_D has a Hamilton path through some reversed arc of D.
# Also: do all unused intra-block arcs of a strong tournament share one cut vertex?
import sys, itertools
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import tournaments, ham_paths, is_strong
from ham_recovery import components
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for n in range(3, nmax + 1):
    worst = 0; checked = 0; fail = 0; multi = 0
    for adj in tournaments(n):
        HP = ham_paths(adj)
        used = {(a, b) for P in HP for a, b in zip(P, P[1:])}
        comps = components(adj); ci = {v: i for i, c in enumerate(comps) for v in c}
        U = [(u, v) for u in range(n) for v in range(n) if (adj[u] >> v) & 1 and (u, v) not in used and ci[u] == ci[v]]
        if not U: continue
        worst = max(worst, len(U))
        if len(U) > 1: multi += 1
        for r in range(1, len(U) + 1):
            for D in itertools.combinations(U, r):
                new = list(adj)
                for u, v in D:
                    new[u] &= ~(1 << v); new[v] |= 1 << u
                rev = {(v, u) for u, v in D}
                ok = any(any((a, b) in rev for a, b in zip(P, P[1:])) for P in ham_paths(tuple(new)))
                checked += 1
                if not ok: fail += 1
    print(f"n = {n}: classes with unused intra-block arcs; max per class {worst}, classes with several {multi}; "
          f"flip sets checked {checked}, sets that leave HP unchanged {fail}", flush=True)
