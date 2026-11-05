#!/usr/bin/env python3
"""descent_certificate.py (klein-2026-11-05): executes the proof of PROOF_CONVERSE.md on every two-cycle class of
data/cls_{n}.txt, orders N1..N2 (usage: descent_certificate.py N1 N2)."""
# Certificate of the converse proof on every two-cycle tournament of orders N1..N2:
# 1. a cut vertex exists; 2. for the first cut vertex c0, {S,E} = {1,2} and the middle components are singletons;
#    if S = 2 pass to the converse tournament; 3. descend: while c has a single in-neighbour y in Q_r (E-ii), replace c by y
#    (check: y is a cut vertex, S(y) = 1, its last component is a proper subset of the old Q_r); 4. at an E-i cut vertex
#    (two in-neighbours with e = 1): Q_r is uniquely Hamiltonian, the pair is {y0, y0-}, {y0-, y0--} (vertex pattern) or
#    {y0, y0--} with y0- a sole entry (arc pattern); 5. reconstruct A, B and check Y is exactly the gluing.
import sys
from collections import Counter
import os
from cutvertex_lib import *
HERE = os.path.dirname(os.path.abspath(__file__))
def cut_data(rows,n,c):
    verts=[v for v in range(n) if v!=c]
    if strong(rows,verts): return None
    comps=components(rows,verts); Q1=comps[0]; Qr=comps[-1]
    t1,SE1=hpaths(rows,Q1); tr,SEr=hpaths(rows,Qr)
    S=sum(cnt for (s,e),cnt in SE1.items() if (rows[c]>>s)&1); E=sum(cnt for (s,e),cnt in SEr.items() if (rows[e]>>c)&1)
    return S,E,comps,SEr
D = os.path.join(HERE, '..', 'data') + '/'
for n in range(int(sys.argv[1]), int(sys.argv[2])+1):
    stats=Counter(); total=0; steps=Counter()
    for l in open(D+'cls_%d.txt'%n):
        h,d6=l.split()
        if h!='2': continue
        rows=rows_from_d6(d6); total+=1
        cuts=[c for c in range(n) if cut_data(rows,n,c) is not None]
        assert cuts, "no cut vertex"
        c=cuts[0]; S,E,comps,SEr=cut_data(rows,n,c)
        assert {S,E}=={1,2} and all(len(Q)==1 for Q in comps[1:-1])
        R=rows
        if S==2:
            R=conv(rows,n); S,E,comps,SEr=cut_data(R,n,c); assert (S,E)==(1,2); stats['converse used']+=1
        nsteps=0
        while True:
            Qr=comps[-1]; ins=[y for y in Qr if (R[y]>>c)&1]
            if len(ins)==2: break
            assert len(ins)==1; y=ins[0]
            d=cut_data(R,n,y); assert d is not None, "y not a cut vertex"
            S2,E2,comps2,SEr2=d
            assert S2==1 and E2==2 and all(len(Q)==1 for Q in comps2[1:-1])
            assert set(comps2[-1]) < set(Qr), "two-side did not shrink"
            c=y; S,E,comps,SEr=S2,E2,comps2,SEr2; nsteps+=1
        steps[nsteps]+=1
        Qr=comps[-1]; k=len(Qr); ins=[y for y in Qr if (R[y]>>c)&1]
        ev={y:sum(cnt for (s,e),cnt in SEr.items() if e==y) for y in ins}; assert sorted(ev.values())==[1,1]
        assert hc_count(R,Qr)==1, "Q_r not UH"
        sub=[sum(((R[v]>>u)&1)<<j for j,u in enumerate(Qr)) for v in Qr]; C=cycles(sub,k); assert len(C)==1
        cyc=[Qr[i] for i in C[0]]; pos={v:i for i,v in enumerate(cyc)}
        indeg={v:sum(1 for u in Qr if (R[u]>>v)&1) for v in Qr}
        pattern=None
        for s_ in ins:
            pred=cyc[(pos[s_]-1)%k]
            if indeg[s_]==1 and pred in ins: pattern=('vertex',s_,None,pred); break
        if pattern is None:
            for sp in ins:
                s_=cyc[(pos[sp]-1)%k]; pred=cyc[(pos[s_]-1)%k]
                if indeg[sp]==1 and indeg[s_]==1 and pred in ins: pattern=('arc',s_,sp,pred); break
        assert pattern is not None, "pair is neither pattern"
        kind,s_,sp,a_k=pattern
        junction=[s_] if kind=='vertex' else [s_,sp]
        Bv=[v for v in range(n) if v not in Qr]+junction
        # B uniquely Hamiltonian, junction = consecutive sole exits of B ending at c; cross arcs: all B->A except a_k -> c
        assert hc_count(R,Bv)==1
        for j,v in enumerate(junction):
            outs=[u for u in Bv if (R[v]>>u)&1]; nxt=junction[j+1] if j+1<len(junction) else c
            assert outs==[nxt], ("junction out-arcs", outs, nxt)
        for b in Bv:
            if b in junction: continue
            for a in Qr:
                if a in junction: continue
                if b==c and a==a_k: assert (R[a]>>b)&1
                else: assert (R[b]>>a)&1, "cross arc direction"
        stats[kind]+=1
    print("order %d: %d two-cycle classes certified; %s; descent steps %s"%(n,total,dict(stats),dict(sorted(steps.items()))), flush=True)
