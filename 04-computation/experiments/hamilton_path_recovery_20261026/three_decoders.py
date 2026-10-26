# Three decoders of a tournament from its Hamilton paths, run side by side on every tournament class n <= nmax:
#   (E) endpoint rule (this lane, THM-4807): inside a strong block, orient an untraced pair x -> y iff some
#       block-restricted Hamilton path starts at x and ends at y;
#   (C) two-step count rule (guard_reconstruction_20261026 G1): orient x -> y iff c_D(x, y) >= 2, where c_D counts
#       middle vertices w with x -> w -> y in the trace support D;
#   (R) deleted-vertex reachability rule (Session70 G4): orient by one-way reachability in D - z for a vertex z
#       giving exactly one direction.
# All three use the same component recovery (strong components of D, forward between them).  Also records, for every
# untraced arc u -> v inside a strong block, the pair (c_D(u, v), c_D(v, u)) and the number of eligible z.
import sys
from collections import Counter
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import tournaments, ham_paths
from ham_recovery import components

def reach_from(Dadj, n, s, banned):
    seen = 1 << s; frontier = 1 << s
    while frontier:
        nxt = 0; f = frontier
        while f:
            u = (f & -f).bit_length() - 1; f &= f - 1
            nxt |= Dadj[u]
        nxt &= ~banned
        frontier = nxt & ~seen; seen |= nxt
    return seen

def decode_all(n, HP):
    Dadj = [0] * n
    for P in HP:
        for a, b in zip(P, P[1:]): Dadj[a] |= 1 << b
    comps = components(tuple(Dadj))          # strong components of the trace support, in order
    ci = {v: i for i, c in enumerate(comps) for v in c}
    ends = {(P[0], P[-1]) for P in HP}
    # block-restricted endpoint pairs
    blocks = {}
    for i, c in enumerate(comps):
        cs = set(c); prs = set()
        for P in HP:
            sub = [v for v in P if v in cs]
            prs.add((sub[0], sub[-1]))
        blocks[i] = prs
    out = {}
    stats = []
    for rule in 'ECR':
        rows = [0] * n
        for x in range(n):
            for y in range(n):
                if x == y: continue
                if (Dadj[x] >> y) & 1: rows[x] |= 1 << y; continue
                if (Dadj[y] >> x) & 1: continue
                if ci[x] != ci[y]:
                    if ci[x] < ci[y]: rows[x] |= 1 << y
                    continue
                if rule == 'E':
                    fx, fy = (x, y) in blocks[ci[x]], (y, x) in blocks[ci[x]]
                    if fx == fy: return None, None
                    if fx: rows[x] |= 1 << y
                elif rule == 'C':
                    cxy = sum(1 for w in range(n) if (Dadj[x] >> w) & 1 and (Dadj[w] >> y) & 1)
                    cyx = sum(1 for w in range(n) if (Dadj[y] >> w) & 1 and (Dadj[w] >> x) & 1)
                    if (cxy >= 2) == (cyx >= 2): return None, None
                    if cxy >= 2: rows[x] |= 1 << y
                    if x < y: stats.append(('c', max(cxy, cyx), min(cxy, cyx)))
                else:
                    votes = set(); elig = 0
                    for z in comps[ci[x]]:
                        if z in (x, y): continue
                        rx = reach_from(Dadj, n, x, 1 << z); ry = reach_from(Dadj, n, y, 1 << z)
                        a, b = bool((rx >> y) & 1), bool((ry >> x) & 1)
                        if a != b: votes.add(a); elig += 1
                    if len(votes) != 1: return None, None
                    if votes.pop(): rows[x] |= 1 << y
                    if x < y: stats.append(('z', elig))
        out[rule] = tuple(rows)
    return out, stats

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for n in range(3, nmax + 1):
    ok = Counter(); fails = 0; cstats = Counter(); zstats = Counter(); classes = 0
    for adj in tournaments(n):
        classes += 1
        out, stats = decode_all(n, ham_paths(adj))
        if out is None: fails += 1; continue
        for r, t in out.items(): ok[r] += t == tuple(adj)
        for s in stats:
            if s[0] == 'c': cstats[(s[1], s[2])] += 1
            else: zstats[s[1]] += 1
    print(f"n = {n}: classes {classes}; decoded exactly E {ok['E']}, C {ok['C']}, R {ok['R']}; ambiguous {fails}")
    print(f"    untraced pairs: (max, min) two-step counts {dict(sorted(cstats.items()))}; "
          f"eligible separators z per pair {dict(sorted(zstats.items()))}", flush=True)
