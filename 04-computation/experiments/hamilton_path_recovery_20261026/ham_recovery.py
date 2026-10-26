# Do the Hamilton paths of a labelled tournament determine it?  (vitali-response C5 left the all-order
# guard-injectivity open: the legal re-marking guards of a presentation are exactly its Hamilton paths.)
# Checks, all tournament classes n <= 8 (nauty gentourng):
#   (1) arcs on no Hamilton path: every one joins two NON-consecutive strong components (hence points forward);
#   (2) inside a strong tournament every arc lies on a Hamilton path;
#   (3) the strong-component cuts are exactly the positions where all Hamilton paths share their prefix set;
#   (4) consequently HP(T) determines T: reconstruct T from HP(T) alone and compare (all classes; and all labelled
#       tournaments for n <= 6 by brute force over all 2^C(n,2) labelled tournaments).
import sys, itertools, subprocess
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import tournaments, ham_paths, is_strong

def components(adj):
    """strong components in condensation order (first = source component)"""
    n = len(adj)
    reach = [0] * n
    for v in range(n):
        seen = 1 << v; frontier = 1 << v
        while frontier:
            nxt = 0; f = frontier
            while f:
                u = (f & -f).bit_length() - 1; f &= f - 1; nxt |= adj[u]
            frontier = nxt & ~seen; seen |= nxt
        reach[v] = seen
    comps = []; done = 0
    for v in range(n):
        if (done >> v) & 1: continue
        c = [u for u in range(n) if (reach[v] >> u) & 1 and (reach[u] >> v) & 1]
        for u in c: done |= 1 << u
        comps.append(c)
    # order: component A before B iff A reaches B
    comps.sort(key=lambda c: -bin(reach[c[0]]).count('1'))
    return comps

def reconstruct(n, HP):
    """rebuild the tournament from its Hamilton paths alone"""
    adj = [0] * n
    known = set()
    for P in HP:
        for a, b in zip(P, P[1:]):
            adj[a] |= 1 << b; known.add((a, b))
    # cuts: positions m (1..n-1) where every path has the same prefix set
    prefix = [frozenset(P[:m]) for m in range(n + 1) for P in HP[:1]]
    cuts = [m for m in range(1, n) if len({frozenset(P[:m]) for P in HP}) == 1]
    bounds = [0] + cuts + [n]
    P0 = HP[0]
    blocks = [set(P0[bounds[i]:bounds[i + 1]]) for i in range(len(bounds) - 1)]
    pos = {}
    for bi, B in enumerate(blocks):
        for v in B: pos[v] = bi
    for u in range(n):
        for v in range(n):
            if u != v and pos[u] < pos[v] and not ((adj[u] >> v) & 1 or (adj[v] >> u) & 1):
                adj[u] |= 1 << v              # forced forward across a cut
    # any pair still unknown would be an arc inside a block missed by every path
    unknown = [(u, v) for u in range(n) for v in range(u + 1, n) if not ((adj[u] >> v) & 1 or (adj[v] >> u) & 1)]
    return tuple(adj), unknown, len(cuts)

if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in range(2, nmax + 1):
        T = tournaments(n)
        bad1 = bad2 = bad3 = bad4 = 0; unk_total = 0; noarc_total = 0
        for adj in T:
            HP = ham_paths(adj)
            used = {(a, b) for P in HP for a, b in zip(P, P[1:])}
            comps = components(adj)
            ci = {v: i for i, c in enumerate(comps) for v in c}
            for u in range(n):
                for v in range(n):
                    if (adj[u] >> v) & 1 and (u, v) not in used:
                        noarc_total += 1
                        if not (ci[v] >= ci[u] + 2): bad1 += 1          # must skip a component
            if is_strong(adj) and n >= 3:
                arcs = {(u, v) for u in range(n) for v in range(n) if (adj[u] >> v) & 1}
                if arcs - used: bad2 += 1
            cuts = [m for m in range(1, n) if len({frozenset(P[:m]) for P in HP}) == 1]
            true_cuts = list(itertools.accumulate(len(c) for c in comps))[:-1]
            if cuts != true_cuts: bad3 += 1
            rec, unknown, _ = reconstruct(n, HP)
            unk_total += len(unknown)
            if rec != tuple(adj) or unknown: bad4 += 1
        print(f"n = {n}: {len(T)} classes; arcs on no Hamilton path {noarc_total} (all skip a component: {bad1 == 0}); "
              f"strong tournaments with an arc on no Hamilton path: {bad2}; cut detection exact: {bad3 == 0}; "
              f"reconstruction from HP exact: {bad4 == 0} (unknown pairs {unk_total})", flush=True)
    # labelled brute force: distinct labelled tournaments have distinct HP sets
    for n in range(2, 7):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        seen = {}
        coll = 0
        for bits in range(1 << len(pairs)):
            adj = [0] * n
            for k, (i, j) in enumerate(pairs):
                if (bits >> k) & 1: adj[i] |= 1 << j
                else: adj[j] |= 1 << i
            key = frozenset(ham_paths(tuple(adj)))
            if key in seen: coll += 1
            seen[key] = bits
        print(f"labelled n = {n}: {1 << len(pairs)} tournaments, {len(seen)} distinct Hamilton-path sets, collisions {coll}", flush=True)
