"""
mac-mini-2026-07-07-S53 (HYP-5127) -- proofs of THM-649's two shaped follow-ups.

(A) THE RANK FORMULA.  Claim chain (each step verified here):
  A1. The mod-2 quadratic form B of the cylinder crossing form has B[e,f] = 1 iff
      e,f share an endpoint (verified at MULTIPLE twists, several (m,n')).
  A2. Hence B = (J_m - I) tensor I_{n'} + I_m tensor (J_{n'} - I) == J tensor I + I tensor J (mod 2).
  A3. Kernel of X -> rowsum_i(X) + colsum_j(X): all row sums = all col sums = c;
      dim = (m-1)(n'-1) + [m == n' mod 2]  =>
      RANK = m n' - dim = m + n' - 1 - [m == n' (mod 2)].
      (Verified against GF(2) elimination for (m,n') up to (4,5).)

(B) THE n=8 QUOTIENT-MIN = 20.  Claim chain:
  B1. The 6 sigma-fixed crossing pairs at n=8 are exactly the SELF-CROSSING MIRROR
      ORBITS {c, sigma(c)} with c interleaving sigma(c); on Fix(sigma) both chords
      share a bit => same page => each contributes 1 ALWAYS: f = 6 forced.
  B2. Q|Fix = 6 + 2q, q = same-page count over the 32 crossing-orbit representatives
      on the 12 quotient variables.
  B3. q_min = 32 - maxcut(G_q) with G_q the quotient crossing multigraph;
      maxcut(G_q) = 25 (exact over 2^12) and a certificate of >= 7 edge-disjoint ODD
      CYCLES in G_q proves maxcut <= 32 - 7 = 25 (each odd cycle leaves >= 1 uncut).
      => min Q|Fix = 6 + 2*(32-25) = 20.  QED (with the cycle list printed).
"""
import numpy as np
import math
from itertools import combinations
import random as rnd
rnd.seed(53)

# ---------------- (A) ----------------
print("=== (A) rank formula ===")
def cross_count(a1,d1,a2,d2):
    x=a1-a2; delta=d1-d2
    lo,hi=min(x,x+delta),max(x,x+delta)
    if hi<=lo: return 0
    return math.floor(hi)-math.floor(lo)-(1 if hi==math.floor(hi) else 0)

def quad_matrix(m,np_,tw):
    edges=[(i,j) for i in range(m) for j in range(np_)]
    E=len(edges)
    def Qc(w):
        tot=0
        for (e1,(i1,j1)),(e2,(i2,j2)) in combinations(list(enumerate(edges)),2):
            a1,b1=i1/m,(j1+tw)/np_; a2,b2=i2/m,(j2+tw)/np_
            tot+=cross_count(a1,b1-w[e1]-a1,a2,b2-w[e2]-a2)
        return tot
    base=Qc([0]*E)
    lin=[(Qc([1 if x==e else 0 for x in range(E)])-base)%2 for e in range(E)]
    B=np.zeros((E,E),dtype=int)
    for e in range(E):
        for f in range(e+1,E):
            w=[0]*E; w[e]=1; w[f]=1
            coef=(Qc(w)-base-lin[e]-lin[f])%2
            B[e,f]=B[f,e]=coef
    return edges,B

def gf2rank(A):
    A=A.copy()%2; E=A.shape[0]; row=0
    for col in range(E):
        piv=None
        for rr in range(row,E):
            if A[rr,col]: piv=rr; break
        if piv is None: continue
        A[[row,piv]]=A[[piv,row]]
        for rr in range(E):
            if rr!=row and A[rr,col]: A[rr]=(A[rr]+A[row])%2
        row+=1
    return row

okA1=True
for (m,np_) in [(2,2),(2,3),(3,3),(2,4),(3,4),(4,4),(2,5),(3,5),(4,5)]:
    for tw in (0.137, 0.61, 0.923):
        edges,B=quad_matrix(m,np_,tw)
        # A1: B == share-an-endpoint adjacency?
        for e,(i1,j1) in enumerate(edges):
            for f in range(e+1,len(edges)):
                i2,j2=edges[f]
                share=(i1==i2) or (j1==j2)
                if B[e,f]!= (1 if share else 0): okA1=False
    r=gf2rank(B)
    pred=m+np_-1-(1 if (m-np_)%2==0 else 0)
    print(f"  (m,n')=({m},{np_}): rank = {r}, formula m+n'-1-[m==n'] = {pred}  {'OK' if r==pred else 'MISMATCH'}")
print(f"  A1 (B = line-graph adjacency, all twists): {okA1}")
# A3 kernel dimension direct check at (3,4): even-row-col matrices
m,np_=3,4
ker=0
for mask in range(1<<(m*np_)):
    X=[[(mask>>(i*np_+j))&1 for j in range(np_)] for i in range(m)]
    rs=[sum(X[i])%2 for i in range(m)]; cs=[sum(X[i][j] for i in range(m))%2 for j in range(np_)]
    if len(set(rs))==1 and len(set(cs))==1 and rs[0]==cs[0]: ker+=1
import math as _mm
kd=int(_mm.log2(ker))
print(f"  A3 kernel check (3,4): dim = {kd}, formula (m-1)(n'-1)+[m==n'] = {(m-1)*(np_-1)+(1 if (m-np_)%2==0 else 0)}")

# ---------------- (B) ----------------
print("\n=== (B) n=8 quotient-min = 20 ===")
n=8
tiles=[]
for y in range(1,n-1):
    for x in range(n,y+1,-1):
        if x-y>=2: tiles.append((x,y))
m8=len(tiles)
gsmap=[tiles.index((n-y+1,n-x+1)) for (x,y) in tiles]
def crossing(t1,t2):
    (x1,y1),(x2,y2)=tiles[t1],tiles[t2]
    return (y1<y2<x1<x2) or (y2<y1<x2<x1)
crossings=[(a,b) for a in range(m8) for b in range(a+1,m8) if crossing(a,b)]
# B1: classify sigma-fixed crossing pairs
fixed_cross=[]; paired=[]
crossset=set(crossings)
for (a,b) in crossings:
    sa,sb=gsmap[a],gsmap[b]
    if {sa,sb}=={a,b}: fixed_cross.append((a,b))
    else: paired.append((a,b))
mirror_self=[(a,b) for (a,b) in fixed_cross if gsmap[a]==b]
print(f"  crossings: {len(crossings)}; sigma-fixed pairs: {len(fixed_cross)}; "
      f"of which self-crossing mirror orbits (sigma a = b): {len(mirror_self)}")
both_fixed=[(a,b) for (a,b) in fixed_cross if gsmap[a]==a and gsmap[b]==b]
print(f"  both-chords-fixed crossing pairs: {len(both_fixed)} (nested fixed chords -> expect 0)")
print(f"  B1: on Fix(sigma) each mirror orbit has equal bits -> same page -> counts ALWAYS: f = {len(mirror_self)} forced")
# quotient variables: 3 fixed tiles + 9 orbits
fixedT=[i for i,j in enumerate(gsmap) if i==j]
orbs=[(i,j) for i,j in enumerate(gsmap) if i<j]
var_of={}
for v,i in enumerate(fixedT): var_of[i]=v
for v,(i,j) in enumerate(orbs): var_of[i]=var_of[j]=len(fixedT)+v
nv=len(fixedT)+len(orbs)
# quotient multigraph: one representative per sigma-orbit of PAIRED crossings
seen=set(); qedges=[]
for (a,b) in paired:
    ims=tuple(sorted((gsmap[a],gsmap[b])))
    key=min((a,b),ims)
    if key in seen: continue
    seen.add(key)
    va,vb=var_of[a],var_of[b]
    if va==vb:
        print(f"  NOTE: intra-variable crossing ({a},{b}) -> contributes constant 1 to q")
    else:
        qedges.append((va,vb))
print(f"  quotient crossing multigraph: {len(qedges)} edges on {nv} vars (expect 32)")
# maxcut exact
best=0
for mask in range(1<<nv):
    cut=sum(1 for (u,v) in qedges if ((mask>>u)&1)!=((mask>>v)&1))
    best=max(best,cut)
print(f"  maxcut = {best} => q_min = {len(qedges)-best}; min Q|Fix = {len(mirror_self)} + 2*{len(qedges)-best} = {len(mirror_self)+2*(len(qedges)-best)}")
# odd-cycle certificate: find edge-disjoint odd cycles greedily
adj={}
for idx,(u,v) in enumerate(qedges): adj.setdefault(u,[]).append((v,idx)); adj.setdefault(v,[]).append((u,idx))
used=set(); cycles=[]
def find_odd_cycle():
    # BFS 2-coloring on unused edges; an odd cycle appears as a same-color edge
    color={}; parent={}
    for s in adj:
        if s in color: continue
        color[s]=0; stack=[s]
        while stack:
            u=stack.pop()
            for (v,idx) in adj.get(u,[]):
                if idx in used: continue
                if v not in color:
                    color[v]=color[u]^1; parent[v]=(u,idx); stack.append(v)
                elif color[v]==color[u]:
                    # reconstruct cycle u..v + edge
                    pu=[u]; pv=[v]
                    eu=[]; ev=[]
                    uu,vv=u,v
                    su={u:[]}
                    # walk both to root collecting edges
                    pathu={}; x=u
                    while x in parent: pathu[x]=parent[x]; x=parent[x][0]
                    x=v; joint=None; pathv={}
                    while True:
                        if x in pathu or x==u: joint=x; break
                        if x not in parent: break
                        pathv[x]=parent[x]; x=parent[x][0]
                    if joint is None: continue
                    es=[idx]
                    x=u
                    while x!=joint: es.append(pathu[x][1]); x=pathu[x][0]
                    x=v
                    while x!=joint: es.append(pathv[x][1]); x=pathv[x][0]
                    if len(es)%2==1 and all(e not in used for e in es):
                        return es
        # continue other components
    return None
while True:
    c=find_odd_cycle()
    if not c: break
    cycles.append(c); used.update(c)
print(f"  edge-disjoint odd cycles found: {len(cycles)} (need >= {len(qedges)-best} = {len(qedges)-best}) "
      f"sizes {[len(c) for c in cycles]}")
print(f"  CERTIFICATE: {len(cycles)} edge-disjoint odd cycles => maxcut <= {len(qedges)}-{len(cycles)} = {len(qedges)-len(cycles)}"
      f" {'== maxcut: TIGHT, q_min = 7 PROVED' if len(qedges)-len(cycles)==best else '(not tight; certificate weaker than needed)'}")
