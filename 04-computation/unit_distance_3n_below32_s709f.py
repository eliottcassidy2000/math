"""
monad-explorer-2026-06-06-S709f
===============================
Try to push below N=32: (a) seed from the known 32-vertex Eisenstein solution and
do exhaustive remove-1+repair / remove-2+repair; (b) test other layer>=8 lattices
(disc -7, -8, square) with the same annealing. Exact integer arithmetic.
Floor (THM-420): N>=17.
"""
import math
from collections import Counter

def make_graph(a,b,c,R,H):
    def Q(dx,dy): return a*dx*dx+b*dx*dy+c*dy*dy
    BB=int(math.isqrt(R))+3
    offs=[(dx,dy) for dx in range(-BB,BB+1) for dy in range(-BB,BB+1) if Q(dx,dy)==R]
    univ=[(x,y) for x in range(-H,H+1) for y in range(-H,H+1)]
    uset=set(univ)
    adj={p:[(p[0]+dx,p[1]+dy) for (dx,dy) in offs if (p[0]+dx,p[1]+dy) in uset] for p in univ}
    return univ,uset,adj,offs

# ---- (a) Eisenstein sqrt7, push below 32 from the known solution -----------
SOL32=[(-5,0),(-4,2),(-3,1),(-3,4),(-2,-2),(-2,-1),(-2,3),(-1,0),(-1,1),(-1,5),
(0,-1),(0,0),(0,2),(0,3),(1,-3),(1,-2),(1,1),(1,2),(2,-1),(2,0),(2,3),(2,4),
(3,-2),(3,-1),(3,1),(4,-4),(4,0),(4,1),(5,-2),(5,2),(6,-3),(7,-1)]
univ,uset,adj,offs = make_graph(1,1,1,7,16)
def Uof(S,adj):
    Ss=set(S); return sum(1 for p in Ss for q in adj[p] if q in Ss)//2
print("verify SOL32:", len(SOL32), "U=",Uof(SOL32,adj),"3N=",3*len(SOL32))

def lcg(s):
    s=(1103515245*s+12345)&0x7fffffff; return s,s/0x7fffffff

def shrink_repair(start, adj, univ, rounds, st):
    """from a feasible set, repeatedly: remove a vertex, then greedily add frontier
       vertices (max internal gain) until feasible again or stuck; keep smallest feasible."""
    best=set(start)
    cur=set(start)
    for rnd in range(rounds):
        S=set(best)
        # remove 1-3 random vertices
        st,r=lcg(st); nrem=1+int(r*3)
        Sl=list(S)
        for _ in range(nrem):
            st,r=lcg(st)
            if S: S.discard(Sl[int(r*len(Sl))%len(Sl)])
        # repair: while infeasible, add best frontier vertex
        def U(S):
            Ss=S; return sum(1 for p in Ss for q in adj[p] if q in Ss)//2
        tries=0
        while U(S)<=3*len(S) and tries<60:
            fr=set(w for v in S for w in adj[v] if w not in S)
            if not fr: break
            p=max(fr,key=lambda v:sum(1 for q in adj[v] if q in S))
            S.add(p); tries+=1
        # then trim any removable vertex
        changed=True
        while changed:
            changed=False
            for p in sorted(S,key=lambda v:sum(1 for q in adj[v] if q in S)):
                d=sum(1 for q in adj[p] if q in S)
                if U(S)-d>3*(len(S)-1):
                    S.discard(p); changed=True; break
        if U(S)>3*len(S) and len(S)<len(best):
            best=set(S)
            print(f"   -> improved to N={len(best)} U={U(best)}")
    return best,st

st=999
best32,st = shrink_repair(SOL32, adj, univ, 8000, st)
print(f"Eisenstein sqrt7 best after shrink-repair: N={len(best32)} U={Uof(best32,adj)}")

# ---- (b) other lattices ----------------------------------------------------
print()
print("Other lattices (annealing from disks):")
def disk(univ,a,b,c,cx,cy,k):
    return sorted(univ,key=lambda p:(a*(p[0]-cx)**2+b*(p[0]-cx)*(p[1]-cy)+c*(p[1]-cy)**2))[:k]
def anneal(seed,adj,steps,st):
    S=set(seed); deg={p:sum(1 for q in adj[p] if q in S) for p in S}
    e=sum(deg.values())//2; best=(len(S),e) if e>3*len(S) else None
    for step in range(steps):
        T=1.5*(1-step/steps)+0.02; k=len(S); st,r=lcg(st)
        if r<0.7 and k>4:
            p=min(S,key=lambda v:deg[v]); ne=e-deg[p]; nk=k-1
            ok=ne>3*nk
            if not ok: st,r2=lcg(st); ok=r2<math.exp(-((3*nk-ne+1))/T)
            if ok:
                for q in adj[p]:
                    if q in S: deg[q]-=1
                e=ne; S.discard(p); deg.pop(p)
        else:
            fr=set(w for v in S for w in adj[v] if w not in S)
            if fr:
                p=max(fr,key=lambda v:sum(1 for q in adj[v] if q in S))
                d=sum(1 for q in adj[p] if q in S); S.add(p); deg[p]=d
                for q in adj[p]:
                    if q in S and q!=p: deg[q]+=1
                e+=d
        if e>3*len(S) and (best is None or len(S)<best[0]): best=(len(S),e)
    return best,st
for name,(a,b,c),R in [("disc-7",(1,1,2),8),("disc-8",(1,0,2),9),("square",(1,0,1),5)]:
    u,us,ad,of=make_graph(a,b,c,R,16)
    L=len(of)
    bb=None
    for (cx,cy) in [(0,0),(0.5,0.5),(0.5,0),(1/3,1/3)]:
        seed=disk(u,a,b,c,cx,cy,80)
        b2,st=anneal(seed,ad,8000,st)
        if b2 and (bb is None or b2[0]<bb[0]): bb=b2
    print(f"   {name:>8} R={R} layer {L}: best N={bb[0] if bb else None}")

print()
print("="*60)
fin=len(best32)
print(f"FINAL RECORD: N={fin}.  THM-420: 17 <= N* <= {fin}")
