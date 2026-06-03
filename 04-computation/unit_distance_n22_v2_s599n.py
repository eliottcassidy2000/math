"""n=22 max unit distances on triangular lattice, robust search: the optimum is COMPACT, so take
the 22 nearest lattice points to each high-symmetry center (lattice vertex / edge midpoint /
triangle centroid) and count edges; max over centers + ties. Harborth target floor(3n-sqrt(12n-3)).
opus-2026-06-03-S599n."""
from math import sqrt, floor
DIRS=[(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
def cart(i,j): return (i+0.5*j, (sqrt(3)/2)*j)
def edges(S):
    Sset=set(S); e=0
    for (a,b) in S:
        for (da,db) in DIRS:
            if (a+da,b+db) in Sset: e+=1
    return e//2
def harborth(n): return floor(3*n - sqrt(12*n-3))
def nearest_cluster(cx,cy,pts,n):
    order=sorted(pts, key=lambda p:( (cart(*p)[0]-cx)**2+(cart(*p)[1]-cy)**2, p))
    return order[:n]
def main():
    n=22; tgt=harborth(n)
    print(f"Harborth n=22 = {tgt}")
    R=6; pts=[(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if abs(i+j)<=R]
    # candidate centers: lattice points, edge midpoints, triangle centroids over a small region
    centers=set()
    for i in range(-3,4):
        for j in range(-3,4):
            x,y=cart(i,j); centers.add((round(x,4),round(y,4)))
            for (da,db) in DIRS:
                x2,y2=cart(i+da,j+db); centers.add((round((x+x2)/2,4),round((y+y2)/2,4)))
            # centroid of triangle (i,j),(i+1,j),(i,j+1)
            xs=[cart(i,j),cart(i+1,j),cart(i,j+1)]
            centers.add((round(sum(p[0] for p in xs)/3,4), round(sum(p[1] for p in xs)/3,4)))
    best=-1; bestS=None
    for (cx,cy) in centers:
        S=nearest_cluster(cx,cy,pts,n); e=edges(S)
        if e>best: best=e; bestS=S
    print(f"best found (compact nearest-22) = {best}  optimal={best>=tgt}  (target {tgt})")
    # local search to push further
    S=set(bestS)
    improved=True
    while improved:
        improved=False
        frontier=set()
        for (a,b) in S:
            for (da,db) in DIRS:
                if (a+da,b+db) not in S: frontier.add((a+da,b+db))
        cur=edges(S)
        for r in list(S):
            for add in frontier:
                T=set(S); T.discard(r); T.add(add)
                if edges(T)>cur: S=T; cur=edges(T); improved=True; break
            if improved: break
    print(f"after local search = {edges(S)}")
    final=edges(S); Sf=set(S)
    from collections import Counter
    deg=[sum(1 for d in DIRS if (a+d[0],b+d[1]) in Sf) for (a,b) in S]
    print(f"FINAL n=22 unit distances = {final} (Harborth {tgt}); deg dist {dict(sorted(Counter(deg).items()))}, interior(6)={deg.count(6)}")
    print(f"cluster: {sorted(S)}")
if __name__=='__main__': main()
