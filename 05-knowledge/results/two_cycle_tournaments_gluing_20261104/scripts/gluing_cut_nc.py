# #cut arcs + |NC| = n for vertex gluings: express NC(Y) and cut(Y) through the pieces A = Y[S2 + s], B = Y[PR + s]
#!/usr/bin/env python3
"""gluing_cut_nc.py (klein-2026-11-05): for every vertex-gluing two-cycle class Y = A (+)_s B of data/cls_{n}.txt checks
(1) NC(Y) = NC(A) u NC(B) u {s, a_1, b_l}; (2) cut(Y) = (cut(A) u cut(B)) minus the cut arcs incident to s; (3) in the
uniquely Hamiltonian pieces the cut arcs are exactly the cycle arcs not internal to the non-critical run (THM-4868 (C)).
With THM-4868's #cut + |NC| = order + 1 for the pieces these give #cut(Y) + |NC(Y)| = n (README section 11).
Usage: gluing_cut_nc.py N1 N2."""
import sys, os
from collections import Counter
from cutvertex_lib import *
HERE = os.path.dirname(os.path.abspath(__file__))
def runs_of(rows,n):
    C=cycles(rows,n); A=[set((c[i],c[(i+1)%n]) for i in range(n)) for c in C]
    c=C[0]; diff=A[0]-A[1]; chg=sorted(i for i in range(n) if (c[i],c[(i+1)%n]) in diff)
    runs=[]
    for j in range(3):
        a=chg[j]+1; b=chg[(j+1)%3]; L=(b-a)%n+1; runs.append([c[(a+t)%n] for t in range(L)])
    return runs
def strong_d(adj, verts):
    vs=list(verts); s0=vs[0]
    for fwd in (True,False):
        seen={s0}; fr=[s0]
        while fr:
            nf=[]
            for v in fr:
                for u in vs:
                    if u not in seen and ((adj[v]>>u)&1 if fwd else (adj[u]>>v)&1): seen.add(u); nf.append(u)
            fr=nf
        if len(seen)<len(vs): return False
    return True
def cutarcs(rows, verts):
    res=set()
    for u in verts:
        for w in verts:
            if (rows[u]>>w)&1:
                r2=rows[:]; r2[u]&=~(1<<w)
                if not strong_d(r2,verts): res.add((u,w))
    return res
def ncset(rows, verts): return set(v for v in verts if strong(rows,[x for x in verts if x!=v]))

D = os.path.join(HERE, '..', 'data') + '/'
for n in range(int(sys.argv[1]), int(sys.argv[2])+1):
    ok=Counter(); tested=0
    for l in open(D+'cls_%d.txt'%n):
        h,d6=l.split()
        if h!='2': continue
        rows=rows_from_d6(d6); runs=runs_of(rows,n); sizes=[len(r) for r in runs]
        if 1 not in sizes: continue
        j=sizes.index(1); s=runs[j][0]; S2=runs[(j+1)%3]; PR=runs[(j+2)%3]; a_1=S2[0]; a_k=S2[-1]; b_1=PR[0]; b_l=PR[-1]
        X=[s]+S2; Bv=[s]+PR; tested+=1
        NCY=ncset(rows,range(n)); NCA=ncset(rows,X); NCB=ncset(rows,Bv)
        cutY=cutarcs(rows,range(n)); cutA=cutarcs(rows,X); cutB=cutarcs(rows,Bv)
        inc_s=set(e for e in (cutA|cutB) if s in e)
        ok['cut(Y) = cut(A)|cut(B) minus the cut arcs incident to s']+= cutY==((cutA|cutB)-inc_s)
        ok['NC(Y) = NC(A)|NC(B)|{s,a_1,b_l}']+= NCY==(NCA|NCB|{s,a_1,b_l})
        # cut arcs of the UH pieces = cycle arcs outside the non-critical run (run = NC, an interval of the cycle)
        for (V,NC,cut) in ((X,NCA,cutA),(Bv,NCB,cutB)):
            k=len(V); sub=[sum(((rows[v]>>u)&1)<<jj for jj,u in enumerate(V)) for v in V]; C=cycles(sub,k)[0]; cyc=[V[i] for i in C]
            arcs=[(cyc[i],cyc[(i+1)%k]) for i in range(k)]
            outside=set(e for e in arcs if not (e[0] in NC and e[1] in NC))
            ok['piece: cut arcs = cycle arcs not inside the NC run']+= cut==outside
    print("order %d: %d vertex gluings"%(n,tested))
    for k_,v in sorted(ok.items()): print("    %-62s %d / %d"%(k_,v,tested if 'piece' not in k_ else 2*tested))
