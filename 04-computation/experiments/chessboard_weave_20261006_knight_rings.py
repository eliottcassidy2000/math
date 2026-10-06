#!/usr/bin/env python3
"""Knight graph vs rings, static structure (chessboard-weave, 2026-10-06).

Task 1: the within-ring knight edges on 8x8 (and 6x6 for comparison):
matching or not, which squares they cover, component structure per ring;
plus the forced linear relations between ring-pair move counts of any
closed tour (half-edge bookkeeping) and the parity cuts.
Run: python3 chessboard_weave_20261006_knight_rings.py
"""
from collections import Counter, defaultdict
from itertools import product

def ring(s, n):
    c = (n - 1) / 2
    return int(max(abs(s[0] - c), abs(s[1] - c)) - 0.5)

def knight_edges(n):
    E = set()
    for i, j in product(range(n), repeat=2):
        for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1)]:
            a, b = i + di, j + dj
            if 0 <= a < n and 0 <= b < n:
                E.add(tuple(sorted(((i, j), (a, b)))))
    return sorted(E)

def alg(s):  # algebraic name, file a..h = column j, rank 1..8 = row i+1
    return "abcdefgh"[s[1]] + str(s[0] + 1)

def report(n):
    E = knight_edges(n)
    R = max(ring((i, j), n) for i in range(n) for j in range(n)) + 1
    print(f"===== {n}x{n}: {len(E)} knight edges, rings 0..{R-1} =====")
    deg = Counter()
    for u, v in E:
        deg[u] += 1; deg[v] += 1
    within = [e for e in E if ring(e[0], n) == ring(e[1], n)]
    print(f"within-ring edges: {len(within)}")
    for r in range(R):
        We = [e for e in within if ring(e[0], n) == r]
        Vr = [(i, j) for i in range(n) for j in range(n) if ring((i, j), n) == r]
        cov = Counter(x for e in We for x in e)
        adj = defaultdict(set)
        for u, v in We:
            adj[u].add(v); adj[v].add(u)
        # components of the within-ring subgraph (on covered squares)
        seen, comps = set(), []
        for v in sorted(cov):
            if v in seen:
                continue
            st, comp = [v], []
            seen.add(v)
            while st:
                x = st.pop(); comp.append(x)
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y); st.append(y)
            comps.append(sorted(comp))
        is_matching = all(c == 1 for c in cov.values())
        print(f" ring {r}: |R|={len(Vr)} within-edges={len(We)} covered squares={len(cov)}"
              f" uncovered={len(Vr)-len(cov)} matching={is_matching}"
              f" within-degree hist={dict(sorted(Counter(cov.values()).items()))}")
        for c in comps:
            ce = [e for e in We if e[0] in c]
            kind = ("K2" if len(c) == 2 else
                    f"C{len(c)}" if len(ce) == len(c) and all(len(adj[x]) == 2 for x in c) else
                    f"{len(c)}v/{len(ce)}e")
            print(f"    component {kind}: " + " ".join(alg(x) for x in c)
                  + "   edges: " + " ".join(alg(a) + "-" + alg(b) for a, b in ce)
                  + "   board-degrees: " + ",".join(str(deg[x]) for x in c))
        if not We:
            print("    (none)")
        unc = [alg(x) for x in Vr if x not in cov]
        print("    uncovered squares:", " ".join(unc) if unc else "(none)")
    # squares of degree 2/3 and their ring
    print(" low-degree squares: deg2 rings", sorted(Counter(ring(x, n) for x in deg if deg[x] == 2).items()),
          " deg3 rings", sorted(Counter(ring(x, n) for x in deg if deg[x] == 3).items()))
    return E

E8 = report(8)
print()
E6 = report(6)

print("""
===== forced relations for a closed 8x8 tour (PROVED by half-edge counting) =====
m_rs = number of tour moves between rings r and s (r<=s). Each square has tour-degree 2,
so ring r carries 2|R_r| half-edges: 2 m_rr + sum_{s!=r} m_rs = 2|R_r| = 8, 24, 40, 56.
No knight edges inside ring 0 or between rings 0 and 3, hence
   ring 0:  m01 + m02                 = 8
   ring 1:  2 m11 + m01 + m12 + m13   = 24
   ring 2:  2 m22 + m02 + m12 + m23   = 40
   ring 3:  2 m33 + m13 + m23         = 56
   total moves = 64.
Cut parity (a cycle crosses every cut an even number of times):
   {0,1} | {2,3}:  m02 + m12 + m13 even;  {3} | rest: m13 + m23 even;  {0} | rest: 8.
Ring blocks (maximal same-ring runs on the cyclic tour) = number of ring-changing moves
   = 64 - (m11 + m22 + m33)   (a closed tour always changes ring, ring 0 has no inner edge).""")
