"""Build the UNIT-DISTANCE TOURNAMENT (base = a unit Ham path of the optimal cluster; non-consecutive
pairs flipped iff unit) and compute its H (Ham-path count, Rédei-odd). Does H land in the spectrum
(avoid forbidden 7,21)? Relate to the unit-count u(n). opus-2026-06-04-S599z."""
from math import sqrt, floor
DIRS=[(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]
def cart(i,j): return (i+0.5*j,(sqrt(3)/2)*j)
def edges(S):
    Sset=set(S); e=0
    for p in S:
        for d in DIRS:
            if (p[0]+d[0],p[1]+d[1]) in Sset: e+=1
    return e//2
def best_cluster(n):
    R=max(4,int(sqrt(n))+3); pts=[(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if abs(i+j)<=R]
    centers=set()
    for i in range(-3,4):
        for j in range(-3,4):
            x,y=cart(i,j); centers.add((round(x,3),round(y,3)))
            for d in DIRS:
                x2,y2=cart(i+d[0],j+d[1]); centers.add((round((x+x2)/2,3),round((y+y2)/2,3)))
    best=None;bestE=-1
    for (cx,cy) in centers:
        order=sorted(pts,key=lambda p:((cart(*p)[0]-cx)**2+(cart(*p)[1]-cy)**2,p)); S=order[:n]; e=edges(S)
        if e>bestE: bestE=e;best=S
    return best,bestE
def unit_adj(S):
    Sset=set(S); idx={p:k for k,p in enumerate(S)}; n=len(S); A=[set() for _ in range(n)]
    for p in S:
        for d in DIRS:
            q=(p[0]+d[0],p[1]+d[1])
            if q in Sset: A[idx[p]].add(idx[q])
    return A
def ham_path(A):
    n=len(A); res=[None]
    def dfs(v,seen,path):
        if res[0]: return
        if len(seen)==n: res[0]=path[:]; return
        for w in sorted(A[v],key=lambda w:len(A[w]-seen)):
            if w not in seen: seen.add(w); dfs(w,seen,path+[w]); seen.discard(w)
    for s in sorted(range(n),key=lambda v:len(A[v])):
        if res[0]: break
        dfs(s,{s},[s])
    return res[0]
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for w in range(n):
                if not(mask>>w&1) and adj[v]>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(n))
def main():
    print(" n | u(n) | unit-Ham-path? | tournament H (base=unit path, unit pairs flipped) | H in spectrum(not 7,21)")
    for n in range(3,11):
        S,E=best_cluster(n); A=unit_adj(S); path=ham_path(A)
        if not path: print(f" {n} | {E} | NO PATH"); continue
        pos={p:k for k,p in enumerate(path)}   # path order: path[0]..path[n-1]
        Aset=[set() for _ in range(n)]
        # base path: path[k+1] -> path[k]; non-consec pair (a<b in path order): flip (b->a) iff unit, else a->b
        unit=[[False]*n for _ in range(n)]
        for u in range(n):
            for w in A[u]: unit[u][w]=True
        adj=[0]*n
        for a in range(n):
            for b in range(a+1,n):
                pa,pb=path[a],path[b]
                isunit=unit[pa][pb]
                if b==a+1:   # base path arc: path[b]->path[a]
                    adj[b]|=1<<a
                else:
                    if isunit: adj[b]|=1<<a   # flipped
                    else: adj[a]|=1<<b
        H=Hcount(n,adj)
        print(f" {n:2d} | {E:2d} | YES | H={H} | {H not in (7,21)} (odd={H%2==1})")
if __name__=='__main__': main()
