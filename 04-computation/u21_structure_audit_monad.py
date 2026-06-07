#!/usr/bin/env python3
"""
u21_structure_audit_monad.py  (monad-explorer, 2026-06-07)

The repo's unit-distance edifice (HYP-2206/2210/2211/2217, reflections S614/
S629/S630/S634/S648) rests on a CENTRAL numerological claim:

    u(21) = 57 = 20 (unit Hamiltonian "spine") + 37 (centered-hex "bulk"),
    with 37 = C_hex(3) = 3*3*4 + 1.

Now that the 5 stored cores are CONFIRMED to be the genuine Alexeev-Mixon-
Parshall extremal graphs (arXiv:2412.11914 Table 2), test this claim against
the ACTUAL graphs:

  (1) Is there really a natural "20 + 37" split?  A "centered-hex bulk" would
      need a degree-6 vertex whose 6 neighbours form a 6-cycle (the centered
      hexagon C_hex(1)=7 motif), nested to 19 vertices (C_hex(2)).  Check for
      it: degree-6 vertices, hexagonal neighbourhoods, largest "centered-hex"
      subgraph.
  (2) Report hard invariants: degree sequence, triangle count, girth,
      Hamiltonicity (traceable), #4-cycles, max common-neighbour, planarity-
      necessary K4/K23-freeness.
  (3) Decide: is "37 = centered-hex number" a structural fact or a coincidence
      (37 is just 57-20, and ANY integer is a sum)?
"""
from collections import Counter, deque

N21_GRAPH6 = [
    "Tsc@IGC@GD?R?S?Wd@A_CK@HG@VM??PRKOUZ",
    "TCKx?D?OI?OMCBSA_L?ApA_gEA\\EG?PBSCPV",
    "TCTWACAG@@CDKC?e?QgQA@OMOq]F??OUcCEj",
    "TCS`?H??XIZ?K_Co`CG@JO[?EOSCpGOTSCE\\",
    "T??_`OhSCSYA@I?c?OyWBEa@c?SIU?Aa[?el",
]


def graph6_decode(code):
    data = [ord(ch) - 63 for ch in code.strip()]
    n = data[0]
    bits = []
    for byte in data[1:]:
        for k in range(5, -1, -1):
            bits.append((byte >> k) & 1)
    adj = [[0] * n for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i][j] = adj[j][i] = 1
            idx += 1
    return adj


def nbrs(adj):
    n = len(adj)
    return [set(j for j in range(n) if adj[i][j]) for i in range(n)]


def edges_of(adj):
    n = len(adj)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if adj[i][j]]


def triangles(adj):
    n = len(adj)
    nb = nbrs(adj)
    t = 0
    for i in range(n):
        for j in nb[i]:
            if j > i:
                t += len(nb[i] & nb[j] & set(range(j + 1, n)))
    return t


def count_4cycles(adj):
    n = len(adj)
    nb = nbrs(adj)
    c = 0
    for i in range(n):
        for j in range(i + 1, n):
            common = len(nb[i] & nb[j])
            c += common * (common - 1) // 2
    # each 4-cycle counted once per diagonal pair = twice
    return c // 2


def girth(adj):
    n = len(adj)
    nb = nbrs(adj)
    best = float("inf")
    for s in range(n):
        dist = {s: 0}
        par = {s: -1}
        q = deque([s])
        while q:
            u = q.popleft()
            for w in nb[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    par[w] = u
                    q.append(w)
                elif par[u] != w:
                    best = min(best, dist[u] + dist[w] + 1)
    return best


def has_hampath(adj):
    """Backtracking Hamiltonian-path existence (traceability)."""
    n = len(adj)
    nb = [sorted(s) for s in nbrs(adj)]
    seen = [False] * n

    def dfs(v, cnt):
        if cnt == n:
            return True
        for w in nb[v]:
            if not seen[w]:
                seen[w] = True
                if dfs(w, cnt + 1):
                    return True
                seen[w] = False
        return False
    # try each start
    order = sorted(range(n), key=lambda v: len(nb[v]))
    for s in order:
        seen = [False] * n
        seen[s] = True
        if dfs(s, 1):
            return True
    return False


def hexagonal_neighbourhoods(adj):
    """Count degree-6 vertices whose 6 neighbours form a 6-cycle (the
    centered-hexagon C_hex(1)=7 motif: a wheel W6)."""
    n = len(adj)
    nb = nbrs(adj)
    hits = []
    for v in range(n):
        if len(nb[v]) != 6:
            continue
        ring = list(nb[v])
        # edges among the 6 neighbours
        sub = [[1 if (adj[a][b]) else 0 for b in ring] for a in ring]
        deg = [sum(r) for r in sub]
        # a 6-cycle: every vertex degree 2, total 6 edges
        nedges = sum(deg) // 2
        if nedges == 6 and all(d == 2 for d in deg):
            hits.append(v)
    return hits


def largest_centered_hex_like(adj):
    """Does the graph contain the 19-vertex centered-hex (2-ring) patch as a
    subgraph?  Necessary local signature: a vertex of degree 6 (center) whose
    neighbours each have >=4 neighbours inside the ball.  Loose proxy: report
    the max 'wheel' size."""
    n = len(adj)
    nb = nbrs(adj)
    # max k such that some vertex has a neighbourhood containing a k-cycle through all nbrs
    return hexagonal_neighbourhoods(adj)


def main():
    print("=" * 74)
    print("u(21)=57 STRUCTURAL AUDIT vs the repo's '57 = 20 spine + 37 hex bulk'")
    print("(genuine AMP arXiv:2412.11914 graphs; author = Alexeev-MIXON-PARSHALL,")
    print(" NOT 'Alexeev-Tikhonov' as the repo repeatedly states)")
    print("=" * 74)
    chex = {1: 7, 2: 19, 3: 37, 4: 61}
    print(f"centered-hex numbers C_hex(k)=3k(k+1)+1: {chex}")
    print("NOTE: the repo reads 37 as C_hex(3). But C_hex(3)=37 is the number of")
    print("POINTS in a 3-ring hexagon, whose unit-distance EDGE count is 9*3^2+3*3=90,")
    print("not 37.  So '37 bulk edges = centered-hex' conflates a POINT count with an")
    print("EDGE count from the start.\n")
    for idx, code in enumerate(N21_GRAPH6, 1):
        adj = graph6_decode(code)
        n = len(adj)
        E = edges_of(adj)
        deg = sorted((sum(r) for r in adj), reverse=True)
        dh = dict(sorted(Counter(deg).items()))
        tri = triangles(adj)
        c4 = count_4cycles(adj)
        g = girth(adj)
        hexv = hexagonal_neighbourhoods(adj)
        trace = has_hampath(adj)
        print(f"--- Core {idx} ---")
        print(f"  |V|={n} |E|={len(E)}  degree hist={dh}")
        print(f"  triangles={tri}  4-cycles={c4}  girth={g}")
        print(f"  degree-6 vertices with hexagonal (W6) nbhd: {hexv if hexv else 'NONE'}")
        print(f"  traceable (Ham path exists): {trace}")
        print()
    print("VERDICT printed in analysis below.")


if __name__ == "__main__":
    main()
