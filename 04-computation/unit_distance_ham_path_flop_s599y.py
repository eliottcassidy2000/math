"""Does the OPTIMAL unit-distance graph contain a Hamiltonian PATH (a 'unit tour')? The tiling-model
base Ham path is unit  <=>  consecutive point-labels are at unit distance  <=>  the unit-distance
graph is TRACEABLE. Find the 'flop' n (first non-traceable optimum) or show it never flops. Report
degrees, the path/tile split, recursive structure. Lattice (Harborth) optima as the constructible
proxy. opus-2026-06-04-S599z."""
from math import sqrt, floor
DIRS=[(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
def cart(i,j): return (i+0.5*j,(sqrt(3)/2)*j)
def edges(S):
    Sset=set(S); e=0
    for (a,b) in S:
        for (da,db) in DIRS:
            if (a+da,b+db) in Sset: e+=1
    return e//2
def harborth(n): return floor(3*n-sqrt(12*n-3))
def best_cluster(n):
    R=max(4,int(sqrt(n))+3); pts=[(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if abs(i+j)<=R]
    centers=set()
    for i in range(-3,4):
        for j in range(-3,4):
            x,y=cart(i,j); centers.add((round(x,3),round(y,3)))
            for (da,db) in DIRS:
                x2,y2=cart(i+da,j+db); centers.add((round((x+x2)/2,3),round((y+y2)/2,3)))
            xs=[cart(i,j),cart(i+1,j),cart(i,j+1)]
            centers.add((round(sum(p[0] for p in xs)/3,3),round(sum(p[1] for p in xs)/3,3)))
    best=None; bestE=-1
    for (cx,cy) in centers:
        order=sorted(pts,key=lambda p:((cart(*p)[0]-cx)**2+(cart(*p)[1]-cy)**2,p))
        S=order[:n]; e=edges(S)
        if e>bestE: bestE=e; best=S
    return best,bestE
def adj_of(S):
    Sset=set(S); idx={p:k for k,p in enumerate(S)}; n=len(S); A=[set() for _ in range(n)]
    for p in S:
        for (da,db) in DIRS:
            q=(p[0]+da,p[1]+db)
            if q in Sset: A[idx[p]].add(idx[q])
    return A
def has_ham_path(A):
    n=len(A)
    # quick necessary: <=2 vertices of degree 1; connected
    deg1=sum(1 for x in A if len(x)==1)
    if deg1>2: return False
    order=sorted(range(n),key=lambda v:len(A[v]))  # start from low-degree
    import sys; sys.setrecursionlimit(10000)
    best=[False]; calls=[0]
    def dfs(v,seen):
        calls[0]+=1
        if calls[0]>2_000_000: return  # budget
        if len(seen)==n: best[0]=True; return
        nb=sorted(A[v],key=lambda w:len(A[w]-seen))  # Warnsdorff
        for w in nb:
            if w not in seen and not best[0]:
                seen.add(w); dfs(w,seen); seen.discard(w)
    for s in order:
        if best[0]: break
        dfs(s,{s})
    return best[0]
def main():
    print(" n | edges(Harborth) | mindeg | traceable(unit Ham path)? | UNIT-PATH split: (n-1)+tiles")
    for n in range(3,29):
        S,E=best_cluster(n); A=adj_of(S); mindeg=min(len(x) for x in A); tr=has_ham_path(A)
        tiles = E-(n-1) if tr else None
        print(f" {n:2d} | {E:3d} (H={harborth(n)}) | {mindeg} | {tr} | {(n-1)}+{tiles if tr else '?'}")
    print("\nReading: traceable=True => the optimal unit-distance graph HAS a unit Ham path, so the")
    print("tiling base path can be ALL unit. If it stays True for all n, the 'flop' never happens on")
    print("the lattice optimum: n-1 unit edges are the base path, the rest (~2n) are unit flipped tiles.")
if __name__=='__main__': main()
