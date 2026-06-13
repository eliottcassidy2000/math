"""Erdős MAX UNIT DISTANCES, n=22. Optimal configs live on the triangular lattice = Eisenstein
integers Z[ζ6]. Harborth's theorem: max unit distances among n triangular-lattice points =
floor(3n - sqrt(12n-3)). Verify the formula on known small n, and FIND a 22-point cluster
achieving it (greedy-grow + swap local search). Connections: the UD graph is the Cayley graph
Cay(Z[ζ6], {6th roots of unity}); the round LRC tournament is Cay(Z/(2n-1), shell-half) — same
Cayley/additive-energy species. opus-2026-06-03-S599m."""
from math import isqrt, floor, sqrt
DIRS=[(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]  # 6 unit vectors in axial hex coords (Z[ζ6])
def edges(S):
    Sset=set(S); e=0
    for (a,b) in S:
        for (da,db) in DIRS:
            if (a+da,b+db) in Sset: e+=1
    return e//2
def harborth(n): return floor(3*n - sqrt(12*n-3))
def grow(n, seed=(0,0)):
    S={seed}
    while len(S)<n:
        # candidate frontier points, pick max new-edge gain (=#neighbors already in S)
        best=None; bestg=-1
        cand={}
        for (a,b) in S:
            for (da,db) in DIRS:
                p=(a+da,b+db)
                if p in S: continue
                cand[p]=cand.get(p,0)+1
        for p,g in cand.items():
            if g>bestg or (g==bestg and best is not None and p<best): bestg=g; best=p
        S.add(best)
    return S
def local_search(S, iters=40000):
    S=set(S); best=edges(S); 
    import itertools
    improved=True
    while improved:
        improved=False
        # try removing a low-degree point and adding a frontier point
        frontier={}
        for (a,b) in S:
            for (da,db) in DIRS:
                p=(a+da,b+db)
                if p not in S: frontier[p]=frontier.get(p,0)+1
        for r in list(S):
            for add,g in frontier.items():
                if add in S: continue
                T=set(S); T.discard(r); T.add(add)
                if len(T)!=len(S): continue
                et=edges(T)
                if et>best:
                    S=T; best=et; improved=True; break
            if improved: break
    return S,best
def main():
    print("Harborth max-unit-distances on the triangular lattice: floor(3n - sqrt(12n-3))")
    known={3:3,4:5,7:12,8:14,12:21,19:42}  # sanity values (triangular-lattice optimum)
    for n in sorted(known):
        print(f"  n={n:2d}: formula={harborth(n):2d}  (sanity {known[n]})  match={harborth(n)==known[n]}")
    print()
    n=22; tgt=harborth(n)
    print(f"n=22 target (Harborth) = {tgt}")
    bestS=None; bestE=-1
    for seed in DIRS+[(0,0),(1,1),(2,0),(0,2),(2,-1),(1,2)]:
        S=grow(n, seed=(0,0))
        S,e=local_search(S)
        if e>bestE: bestE=e; bestS=S
    print(f"best found: {bestE} unit distances  (target {tgt})  optimal={bestE>=tgt}")
    # normalize and print the cluster
    S=sorted(bestS)
    print(f"cluster (axial coords, {len(S)} pts): {S}")
    # degree distribution
    Sset=set(S); deg=[sum(1 for d in DIRS if (a+d[0],b+d[1]) in Sset) for (a,b) in S]
    from collections import Counter
    print(f"degree distribution (max 6): {dict(sorted(Counter(deg).items()))}; sum/2={sum(deg)//2}")
    print(f"interior(deg6)={deg.count(6)}, boundary={len(S)-deg.count(6)}; perimeter ~ sqrt(12n-3)={sqrt(12*n-3):.2f}")
if __name__=='__main__': main()
