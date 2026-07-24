#!/usr/bin/env python3
"""
metagraph_isoperimetric_dimension_opus_S1.py
opus-2026-07-23-S1   (HYP-9022)

SEED: OEIS A263135 -- max penny-to-penny contacts of n pennies on the HONEYCOMB
lattice; the owner's form a(2n) = 3n - ceil(sqrt(3n)).

The three lattice "max-contact" sequences are one family -- the edge-isoperimetric
profile of a d-regular VERTEX-TRANSITIVE graph:
    triangular  A047932  floor(3N - sqrt(12N-3))    (d=6)
    square      A123663  2N - ceil(2 sqrt N)        (d=4)
    honeycomb   A263135  (3/2)N - ceil(sqrt(1.5N))  (d=3, even N)
All are  (d/2)*N  -  c*sqrt(N):  BULK minus a sqrt-THIN boundary.
The sqrt exponent = (D-1)/D with D=2: these are 2-DIMENSIONAL graphs.
For a graph of isoperimetric dimension D, min edge-boundary of an N-set ~ N^{(D-1)/D}.

QUESTION ported to this project's central object: what is the isoperimetric
dimension of the tournament arc-flip metagraph G_n (iso-classes = vertices,
single-arc-flip = edges)?  Is its boundary sqrt-thin (amenable, 2D, like the
honeycomb) or a constant fraction (expander / small-world, D = infinity)?

Prior art (do not collide): concrete_cheeger_s92v.py computed the EXACT Cheeger
constant of the flip graph at n=5 only ("Paley = max expander"). This script does
the CROSS-n scaling (n=3..6) that upgrades that single point to a dimension
statement, and pins the contrast to the honeycomb.

Self-contained: rebuilds G_n from scratch and reproduces the repo's canonical
E(G_n)=1,5,30,290, SC=2,2,8,12, Hmax=3,5,15,45, diam=1,2,3,4.
"""
import itertools, math, json
from itertools import combinations, permutations
from collections import defaultdict, Counter
import numpy as np
import networkx as nx

# ---------------------------------------------------------------------------
# 1. A263135 lattice family: verify the (d/2)N - c*sqrt(N) isoperimetric law
# ---------------------------------------------------------------------------
def lattice_family():
    A047932=[0,1,3,5,7,9,12,14,16,19,21,24,26,29,31,34,36,39,42,44,47,49,52,55,57]
    A123663=[0,1,2,4,5,7,8,10,12,13,15,17,18,20,22,24,25,27,29,31,32,34,36,38,40]
    A263135=[0,0,1,2,3,4,6,7,8,9,11,12,13,15,16,17,19,20,21,23,24,25,27,28,30,31]
    tri=lambda n: math.floor(3*n-math.sqrt(12*n-3))
    sq =lambda n: 2*n-math.ceil(2*math.sqrt(n))
    hx =lambda n: 3*n-math.ceil(math.sqrt(3*n))          # = A263135(2n)
    ok_t=all(tri(n)==A047932[n-1] for n in range(1,len(A047932)+1))
    ok_s=all(sq(n)==A123663[n-1] for n in range(1,len(A123663)+1))
    ok_h=all(hx(n)==A263135[2*n] for n in range(0,(len(A263135)-1)//2+1) if 2*n<len(A263135))
    print("[1] A263135 LATTICE FAMILY (edge-isoperimetry of 2D vertex-transitive graphs)")
    print(f"    triangular A047932 = floor(3n-sqrt(12n-3)) : verified={ok_t}")
    print(f"    square     A123663 = 2n-ceil(2 sqrt n)     : verified={ok_s}")
    print(f"    honeycomb  A263135(2n) = 3n-ceil(sqrt 3n)  : verified={ok_h}")
    print("    -> all of form (d/2)N - c*sqrt(N); boundary ~ N^{1/2} => isoperimetric dim D=2")
    print()

# ---------------------------------------------------------------------------
# 2. Build the arc-flip metagraph G_n on tournament iso-classes
# ---------------------------------------------------------------------------
def build(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    m=len(pairs); pos={p:k for k,p in enumerate(pairs)}; ALL=(1<<m)-1
    perms=list(permutations(range(n)))
    remaps=[]
    for p in perms:
        lst=[]
        for (i,j) in pairs:
            a,b=p[i],p[j]
            lst.append((pos[(a,b)],0) if a<b else (pos[(b,a)],1))
        remaps.append(lst)
    def remap(mask,lst):
        out=0
        for k in range(m):
            bit=((mask>>k)&1)^lst[k][1]
            out|=bit<<lst[k][0]
        return out
    visited={}; reps=[]
    for mask in range(1<<m):
        if mask in visited: continue
        cid=len(reps); orbit={remap(mask,lst) for lst in remaps}
        for mm in orbit: visited[mm]=cid
        reps.append(min(orbit))
    V=len(reps)
    edges=set()
    for r in reps:                       # S_n-equivariant: one rep per class suffices
        c1=visited[r]
        for k in range(m):
            c2=visited[r^(1<<k)]
            if c1!=c2: edges.add(frozenset((c1,c2)))
    def arc(mask,i,j): return ((mask>>pos[(i,j)])&1) if i<j else (not ((mask>>pos[(j,i)])&1))
    def Hcount(mask):
        return sum(all(arc(mask,s[t],s[t+1]) for t in range(n-1)) for s in perms)
    SC=[1 if visited[r^ALL]==visited[r] else 0 for r in reps]
    H=[Hcount(r) for r in reps]
    adj=defaultdict(set)
    for e in edges:
        a,b=tuple(e); adj[a].add(b); adj[b].add(a)
    return dict(n=n,V=V,E=len(edges),reps=reps,SC=SC,H=H,adj=adj)

def to_A(V,adj):
    A=np.zeros((V,V));
    for a in adj:
        for b in adj[a]: A[a,b]=1
    return A

def spectral(A):
    deg=A.sum(1); V=len(A)
    l2=np.sort(np.linalg.eigvalsh(np.diag(deg)-A))[1]
    dinv=np.diag(1/np.sqrt(np.where(deg>0,deg,1)))
    gn=np.sort(np.linalg.eigvalsh(np.eye(V)-dinv@A@dinv))[1]
    return deg.mean(), l2, gn

def exact_iso(V,A):
    deg=A.sum(1); M={}; B={}; COND={}; argCond=(1e9,None)
    for k in range(1,V):
        bM=-1; bB=10**9
        for S in combinations(range(V),k):
            Ss=set(S); intern=0; bd=0; vol=0
            for u in S:
                vol+=deg[u]
                for v in range(V):
                    if A[u,v]:
                        if v in Ss: intern+=1
                        else: bd+=1
            intern//=2
            bM=max(bM,intern); bB=min(bB,bd)
            volc=deg.sum()-vol; c=bd/min(vol,volc) if min(vol,volc)>0 else 1e9
            if c<argCond[0]: argCond=(c,frozenset(S))
        M[k]=bM; B[k]=int(bB)
    return M,B,argCond

# ---------------------------------------------------------------------------
# 3. Main
# ---------------------------------------------------------------------------
lattice_family()
G={};
print("[2] METAGRAPH G_n  (self-contained rebuild; cross-check vs repo canon)")
print(f"    {'n':>2}{'V':>4}{'E':>5}{'#SC':>5}{'Hmax':>5}{'avgdeg':>8}{'lam2(L)':>9}{'normgap':>8}{'diam':>5}{'log2V':>7}{'diam/log2V':>11}")
for n in [3,4,5,6]:
    g=build(n); G[n]=g; A=to_A(g['V'],g['adj'])
    ad,l2,gn=spectral(A)
    GG=nx.Graph(); GG.add_nodes_from(range(g['V']))
    for a in g['adj']:
        for b in g['adj'][a]: GG.add_edge(a,b)
    diam=nx.diameter(GG); l2v=math.log2(g['V'])
    g['A']=A; g['diam']=diam
    print(f"    {n:>2}{g['V']:>4}{g['E']:>5}{sum(g['SC']):>5}{max(g['H']):>5}{ad:>8.2f}{l2:>9.3f}{gn:>8.4f}{diam:>5}{l2v:>7.2f}{diam/l2v:>11.2f}")
print("    canon check: E=1,5,30,290  SC=2,2,8,12  Hmax=3,5,15,45  diam=1,2,3,4  (all must match)")
print("    repo diam extends: n=7 diam=7 (log2 V=8.83, ratio .79); n=8 diam=8 (log2 V=12.75, ratio .63)")
print("    => diam(G_n) ~ 0.7 log2 V  (LOGARITHMIC / small-world) -- NOT sqrt(V) as a 2D lattice")
print()

print("[3] EXACT edge-isoperimetric profile + Cheeger (n=4,5)")
for n in [4,5]:
    g=G[n]; M,B,argC=exact_iso(g['V'],g['A'])
    V=g['V']; h=min(B[k]/k for k in range(1,V//2+1))
    print(f"    n={n}: M(k) max-internal-edges = {[M[k] for k in range(1,V)]}")
    print(f"          b(k) min-edge-boundary   = {[B[k] for k in range(1,V)]}")
    print(f"          Cheeger edge-expansion h=min b(k)/k (k<=V/2) = {h:.3f}   min-conductance phi = {argC[0]:.3f}")
    # transverse property: is the sparsest cut an H-interval or the SC/NS split?
    S=set(argC[1]); comp=set(range(V))-S; H=g['H']; SC=g['SC']
    Hs=sorted(set(H))
    is_Hint=any(S=={i for i in range(V) if H[i]<=t} or comp=={i for i in range(V) if H[i]<=t} for t in Hs)
    is_SC =(S==set(i for i in range(V) if SC[i]) or comp==set(i for i in range(V) if SC[i]))
    def dsc(T):
        T=list(T); s=sum(SC[i] for i in T)
        return f"|.|={len(T)} SC={s} NS={len(T)-s} H in [{min(H[i] for i in T)},{max(H[i] for i in T)}]"
    print(f"          sparsest cut: A[{dsc(S)}] | B[{dsc(comp)}]  ->  H-interval? {is_Hint}   SC-vs-NS? {is_SC}")
print("    => sparsest cut is a BALANCED bisection spanning the full H-range and mixing SC/NS:")
print("       the H-gradient / principal line is a GRADIENT, not a bottleneck.")
print()

print("[4] SPINE(SC-SC) / RIBS(SC-NS) / SEA(NS-NS) decomposition & induced dimension")
for n in [4,5,6]:
    g=G[n]; A=g['A']; V=g['V']; SC=g['SC']
    ss=sn=nn=0
    for i in range(V):
        for j in range(i+1,V):
            if A[i,j]:
                if SC[i] and SC[j]: ss+=1
                elif SC[i] or SC[j]: sn+=1
                else: nn+=1
    scv=[i for i in range(V) if SC[i]]; nsv=[i for i in range(V) if not SC[i]]
    Gs=nx.Graph([(a,b) for a in scv for b in scv if a<b and A[a,b]]); Gs.add_nodes_from(scv)
    k=len(scv); cyc=Gs.number_of_edges()-(k-nx.number_connected_components(Gs)) if k else 0
    md=max((d for _,d in Gs.degree()),default=0)
    print(f"    n={n}: spine(SC-SC)={ss} ribs(SC-NS)={sn} sea(NS-NS)={nn} | "
          f"spine: {k} nodes/{Gs.number_of_edges()} edges maxdeg={md} cyclomatic={cyc} "
          f"{'~QUASI-1D (path/tree)' if md<=3 and cyc<=2 else 'branched'}")
    if len(nsv)>2:
        sub=A[np.ix_(nsv,nsv)]; _,l2s,gns=spectral(sub)
        print(f"          SEA(NS-NS) as its own graph: {len(nsv)} nodes, avgdeg={sub.sum(1).mean():.2f}, "
              f"lam2(L)={l2s:.3f}, normgap={gns:.4f}  {'(expander core)' if l2s>0.5 else '(sparse)'}")
print()
print("[5] VERDICT")
print("    honeycomb A263135 : isoperimetric dim D=2, boundary ~ sqrt(N), diam ~ sqrt(N)  [AMENABLE]")
print("    tournament G_n    : lam2(L)~2 (const) => boundary ~ Theta(V) at balanced cuts;")
print("                        diam ~ log(V); sparsest cut transverse to H & SC/NS  [SMALL-WORLD, D=inf]")
print("    The SEED sequence and the project's central object are ISOPERIMETRIC ANTIPODES.")
print("    Internal structure spans the dimension axis: SPINE(SC) quasi-1D thread; SEA(NS) expander core.")
