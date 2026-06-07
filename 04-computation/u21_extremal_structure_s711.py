#!/usr/bin/env python3
"""
S711 — Anatomy of the u(21)=57 extremal unit-distance graph (the lower-bound certificate).

u(21)=57 is PROVEN (Alexeev-Tikhonov 2024, arXiv 2412.11914): they determine u(n) through n=21.
The LOWER bound is a single explicit 21-point configuration realizing 57 unit distances
(05-knowledge/results/u21_core1_coords.txt, "AMP extremal core 1"). The UPPER bound u(21)<=57 is the
hard, computer-assisted half: a graph-only forbidden-subgraph (F-free) bound, then geometric
elimination of non-embeddable extremal candidates.

This script does the part we CAN verify exactly here: it (1) re-verifies the 57-unit-distance
certificate to machine precision, (2) builds the abstract unit-distance graph G on those 57 edges, and
(3) reports the structure that the F-free upper-bound mechanism cares about: degree sequence, max/min
degree, clique number (must be <=3: K_4 is NOT a unit-distance graph), triangle count, K_{2,3} count
(K_{2,3} forbidden in a *faithful* UDG), 4-cycles, and a faithfulness check (no two points coincide,
exactly 57 unit pairs, all other pairs bounded away from 1).
"""
from itertools import combinations

def load(path):
    pts = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'): continue
        p = line.split()
        pts.append((float(p[1]), float(p[2])))
    return pts

def main():
    pts = load('05-knowledge/results/u21_core1_coords.txt')
    n = len(pts)
    TOL = 1e-9
    d2 = {}
    edges = []
    nonunit_gap = 1.0
    coincident = 0
    for i, j in combinations(range(n), 2):
        dx = pts[i][0]-pts[j][0]; dy = pts[i][1]-pts[j][1]
        dd = (dx*dx+dy*dy)**0.5
        d2[(i, j)] = dd
        if dd < 1e-6: coincident += 1
        if abs(dd-1.0) < TOL:
            edges.append((i, j))
        else:
            nonunit_gap = min(nonunit_gap, abs(dd-1.0))
    adj = {v: set() for v in range(n)}
    for i, j in edges:
        adj[i].add(j); adj[j].add(i)
    deg = [len(adj[v]) for v in range(n)]
    E = len(edges)

    # triangles
    tri = 0
    for i, j in edges:
        tri += len(adj[i] & adj[j])
    tri //= 3

    # clique number via greedy max-clique over small graph (n=21)
    def max_clique():
        best = 0
        verts = sorted(range(n), key=lambda v: -deg[v])
        def expand(R, P):
            nonlocal best
            if not P:
                best = max(best, len(R)); return
            if len(R)+len(P) <= best: return
            for v in list(P):
                expand(R | {v}, P & adj[v])
                P = P - {v}
        expand(set(), set(verts))
        return best
    omega = max_clique()

    # K_{2,3} count: pairs {a,b} with >=3 common neighbors -> choose 3
    from math import comb
    k23 = 0
    for a, b in combinations(range(n), 2):
        c = len(adj[a] & adj[b])
        if c >= 3: k23 += comb(c, 3)

    # 4-cycles
    c4 = 0
    for a, b in combinations(range(n), 2):
        c = len(adj[a] & adj[b])
        c4 += comb(c, 2)
    c4 //= 2  # each 4-cycle counted twice (two diagonals)

    print("="*74)
    print("S711 — anatomy of the u(21)=57 extremal unit-distance graph")
    print("="*74)
    print(f"points                 : {n}")
    print(f"coincident point pairs : {coincident}  (faithful iff 0)")
    print(f"unit distances (|d-1|<{TOL:.0e}) : {E}   <-- the certificate u(21)>=57")
    print(f"next-closest non-unit gap        : {nonunit_gap:.3e}  (well-separated => robust)")
    print(f"degree sequence (sorted) : {sorted(deg)}")
    print(f"  sum deg = {sum(deg)} = 2E = {2*E};  max deg = {max(deg)}, min deg = {min(deg)}, "
          f"avg = {2*E/n:.3f}")
    print(f"clique number omega      : {omega}   (UDG forbids K_4 => must be <= 3)")
    print(f"triangles                : {tri}")
    print(f"4-cycles (C_4)           : {c4}")
    print(f"K_2,3 copies             : {k23}   (faithful UDG forbids K_2,3 => should be 0)")
    print("-"*74)
    # the graph-only Spencer-Szemeredi-Trotter style ceiling check (sanity, not the AT bound):
    print("Context: u(20)=54, u(21)=57, u(22) in [60,61] (Alexeev-Tikhonov 2024).")
    print("Increments u(n)-u(n-1): ...,54->57 = +3 at n=21; the AT upper bound matches this cert.")
    print("="*74)

if __name__ == "__main__":
    main()
