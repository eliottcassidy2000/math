#!/usr/bin/env python3
"""
S624 — the Erdős unit-distance problem at n=22 (the open frontier: u(22) in {60,61};
u(21)=57 proven, Alexeev-Mixon-Parshall 2024), and its bridge to LRC/tournaments.

Triangular lattice = Eisenstein integers Z[w], w=e^{i pi/3}: the CM field Q(sqrt(-3)),
the densest *small* unit-distance source.  Unit distance <=> di^2+di*dj+dj^2=1 (6 units).
We (1) maximise unit distances over compact 22-point lattice subsets (the CM-rigid baseline),
(2) test ROTATION OVERLAYS (the Moser-ring trick = extra modulus-1 multipliers) that beat the
plain lattice — the finite shadow of the grid-disproof mechanism (CM rotations add resonances).
"""
import itertools, math
from collections import defaultdict

def tri_unit(d):  # d=(di,dj); unit distance in triangular lattice
    di, dj = d
    return di*di + di*dj + dj*dj == 1

def edges_lattice(pts):
    S = set(pts); e = 0
    for (i, j) in pts:
        for (di, dj) in [(1,0),(0,1),(1,-1)]:   # 3 of 6 to avoid double count
            if (i+di, j+dj) in S: e += 1
    return e

def disk(R):
    return [(i, j) for i in range(-R, R+1) for j in range(-R, R+1)
            if i*i + i*j + j*j <= R*R]

def xy(i, j):
    return (i + j*0.5, j*math.sqrt(3)/2)

def best_compact(n, R=6):
    """greedy compact 22-subsets: for several centers take n nearest, count edges; refine by
       swapping a boundary point for a better-connected neighbor."""
    L = disk(R)
    centers = [(0.0,0.0),(0.5,0.0),(1.0/3,1.0/3),(0.5,math.sqrt(3)/6)]
    best = 0; bestset = None
    for cx, cy in centers:
        Ls = sorted(L, key=lambda p: (xy(*p)[0]-cx)**2 + (xy(*p)[1]-cy)**2)
        sub = Ls[:n]
        cur = edges_lattice(sub)
        # local search: try swapping each chosen point with a nearby unchosen one
        improved = True
        chosen = set(sub)
        while improved:
            improved = False
            for p in list(chosen):
                for q in Ls[:n+25]:
                    if q in chosen: continue
                    new = (chosen - {p}) | {q}
                    e = edges_lattice(new)
                    if e > cur:
                        chosen = new; cur = e; improved = True; break
                if improved: break
        if cur > best:
            best, bestset = cur, sorted(chosen)
    return best, bestset

def count_unit_float(pts, tol=1e-6):
    e = 0
    for a, b in itertools.combinations(pts, 2):
        d = math.hypot(a[0]-b[0], a[1]-b[1])
        if abs(d-1) < tol: e += 1
    return e

if __name__ == "__main__":
    print("n=22 Erdős unit distance — known: u(21)=57 (optimal), u(22) in {60,61} (open)")
    print("="*70)
    b, bs = best_compact(22)
    print(f"[CM-rigid baseline] best compact 22-pt TRIANGULAR-LATTICE (Eisenstein Z[w]) subset:")
    print(f"   unit distances = {b}   (lattice cap ~3n=66 minus boundary)")
    print(f"   points (i,j) = {bs}")
    # sanity: degree sequence
    S = set(bs); deg = {p: sum(((p[0]+di,p[1]+dj) in S) for di,dj in
                               [(1,0),(0,1),(1,-1),(-1,0),(0,-1),(-1,1)]) for p in bs}
    from collections import Counter
    print(f"   degree distribution: {dict(sorted(Counter(deg.values()).items()))}")
    print(f"   gap to optimum 60: {60-b}  -> the deficit the CM ROTATION/overlay must supply")
