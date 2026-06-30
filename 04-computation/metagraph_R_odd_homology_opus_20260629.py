"""
Concrete piece of klein HYP-3544: the metagraph 1-skeleton homology and its R-even/R-odd split,
via the Hopf trace (Lefschetz) formula -- which ties back to SC = #R-fixed vertices.

Metagraph G_n: vertices = tournament iso classes; edges = single-arc-flip between DISTINCT classes.
R = complement (T->T^op), an involution on G_n. 
  chi(G) = V - E.   Lefschetz: tr(R|H_0) - tr(R|H_1) = tr(R|C_0) - tr(R|C_1).
  tr(R|C_0) = #R-fixed vertices = SC.   tr(R|C_1) = #R-fixed edges (signed: -1 if R flips the edge).
  => tr(R|H_1) = tr(R|C_1) - SC + tr(R|H_0).   b1^- = (b1 - tr(R|H_1))/2  = the R-ODD Betti number.
"""
from itertools import permutations, combinations
def pairs(n): return list(combinations(range(n),2))
def canon(n,bits):
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}; best=None
    for perm in permutations(range(n)):
        out=tuple(bits[idx[(perm[a],perm[b])]] if perm[a]<perm[b] else 1-bits[idx[(perm[b],perm[a])]] for (a,b) in pr)
        if best is None or out<best: best=out
    return best
def complement(bits): return tuple(1-b for b in bits)

def metagraph(n):
    pr=pairs(n); m=len(pr); reps={}
    for x in range(1<<m):
        bits=tuple((x>>k)&1 for k in range(m)); c=canon(n,bits)
        if c not in reps: reps[c]=bits
    V=list(reps); vid={c:i for i,c in enumerate(V)}
    edges=set()
    for c in V:
        bits=reps[c]
        for k in range(m):
            nb=tuple(b^(1 if j==k else 0) for j,b in enumerate(bits))
            cc=canon(n,nb)
            if cc!=c: edges.add(frozenset((vid[c],vid[cc])))
    return V,vid,edges
for n in [4,5,6]:
    V,vid,edges=metagraph(n)
    nV=len(V); nE=len(edges); comp=1  # wiggly graph connected (CLAUDE.md); verify via union-find
    par=list(range(nV))
    def f(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    for e in edges:
        a,b=tuple(e); par[f(a)]=f(b)
    comp=len({f(i) for i in range(nV)})
    b1=nE-nV+comp
    # R action
    Rv=[vid[canon(n,complement(V[i]))] for i in range(nV)]
    SC=sum(1 for i in range(nV) if Rv[i]==i)
    # R on edges: fixed if {Ra,Rb}={a,b}; signed -1 if R swaps the two endpoints (flips the edge)
    trC1=0
    for e in edges:
        a,b=tuple(e); ra,rb=Rv[a],Rv[b]
        if {ra,rb}=={a,b}:
            trC1 += -1 if (ra,rb)==(b,a) and a!=b else 1   # flip => -1
    trH0=comp  # R fixes each component (acts within); = #components for a connected/symmetric graph
    trH1 = trC1 - SC + trH0
    bm = (b1 - trH1)//2   # R-odd Betti number b1^-
    bp = b1 - bm
    print(f"n={n}: V={nV} E={nE} comp={comp} | b0={comp} b1={b1} | SC(R-fixed verts)={SC} "
          f"tr(R|C1)={trC1} tr(R|H1)={trH1} | b1^+={bp} b1^-(R-ODD)={bm}")
print("\nR-ODD Betti b1^- = the cycle-space obstruction dimension; candidate match to LRC cap M_odd (2-dim).")
