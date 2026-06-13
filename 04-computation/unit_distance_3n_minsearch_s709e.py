"""
monad-explorer-2026-06-06-S709e
===============================
Pin down the minimum: aggressive multistart tabu/anneal for the smallest induced
subgraph of the sqrt(7)-Eisenstein unit graph with U > 3N (avg degree > 6).
Dumps & re-verifies the winning vertex set with independent exact integer counting.

THM-420 floor: N >= 17.  Prior records: 43 (HYP-2267 disk) -> 39 (s709c) -> 32 (s709d).
"""
import math
A,B,Cc = 1,1,1
def Q(dx,dy): return dx*dx+dx*dy+dy*dy
OFFS=[(dx,dy) for dx in range(-4,5) for dy in range(-4,5) if Q(dx,dy)==7]
H=16
UNIV=[(x,y) for x in range(-H,H+1) for y in range(-H,H+1)]
USET=set(UNIV)
ADJ={p:[(p[0]+dx,p[1]+dy) for (dx,dy) in OFFS if (p[0]+dx,p[1]+dy) in USET] for p in UNIV}

def U_exact(S):
    Ss=set(S); return sum(1 for p in Ss for q in ADJ[p] if q in Ss)//2

def lcg(s):
    s=(1103515245*s+12345)&0x7fffffff; return s, s/0x7fffffff

def disk(cx,cy,k):
    return sorted(UNIV,key=lambda p:((p[0]-cx)**2+(p[0]-cx)*(p[1]-cy)+(p[1]-cy)**2))[:k]

def search(seed, steps, st):
    """tabu-ish shrink+swap: keep feasible (U>3N), try to reduce N; allow temporary
       infeasible moves controlled by temperature; track best feasible N."""
    S=set(seed)
    deg={p:sum(1 for q in ADJ[p] if q in S) for p in S}
    e=sum(deg.values())//2
    best=None
    if e>3*len(S): best=(len(S),e,frozenset(S))
    for step in range(steps):
        T=1.5*(1-step/steps)+0.02
        k=len(S); st,r=lcg(st)
        if r<0.7 and k>4:
            # remove min-degree vertex (random tie-break among bottom few)
            cand=sorted(S,key=lambda v:deg[v])[:3]
            st,r2=lcg(st); p=cand[int(r2*len(cand))%len(cand)]
            ne=e-deg[p]; nk=k-1
            ok = ne>3*nk
            if not ok:
                st,r3=lcg(st); ok = r3<math.exp(-((3*nk-ne+1))/T)
            if ok:
                for q in ADJ[p]:
                    if q in S: deg[q]-=1
                e=ne; S.discard(p); deg.pop(p)
        else:
            fr=set(w for v in S for w in ADJ[v] if w not in S)
            if fr:
                p=max(fr,key=lambda v:sum(1 for q in ADJ[v] if q in S))
                d=sum(1 for q in ADJ[p] if q in S)
                S.add(p); deg[p]=d
                for q in ADJ[p]:
                    if q in S and q!=p: deg[q]+=1
                e+=d
        if e>3*len(S) and (best is None or len(S)<best[0]):
            best=(len(S),e,frozenset(S))
    return best,st

best_overall=None
st=12345
centers=[(0,0),(0.5,0),(1/3,1/3),(0.5,0.5),(2/3,1/3),(1/3,2/3),(0,0.5),(0.25,0.25)]
for rep in range(6):
    for (cx,cy) in centers:
        for k0 in [28,32,40,55,80]:
            seed=disk(cx,cy,k0)
            b,st=search(seed,6000,st+rep*7+1)
            if b and (best_overall is None or b[0]<best_overall[0]):
                best_overall=b
                print(f"   NEW best N={b[0]} U={b[1]}  (3N={3*b[0]}, surplus {b[1]-3*b[0]:+d})")

print()
print("="*68)
S=sorted(best_overall[2])
N=len(S); Uv=U_exact(S)
print(f"BEST FOUND: N={N}  U(independent recount)={Uv}  3N={3*N}  feasible={Uv>3*N}")
# degree distribution
Ss=set(S)
degs=sorted(sum(1 for q in ADJ[p] if q in Ss) for p in S)
from collections import Counter
print(f"degree distribution: {dict(sorted(Counter(degs).items()))}, avg={2*Uv/N:.3f}")
print(f"vertex set (Eisenstein indices, Q=x^2+xy+y^2, unit^2=7):")
print("  "+", ".join(str(p) for p in S))
print("="*68)
print(f"THM-420 status: floor 17 <= N* <= {N}")
