#!/usr/bin/env python3
"""cutvertex_lib.py (klein-2026-11-05): helpers for the cut-vertex decomposition of tournaments -- Hamiltonian path counts by
(start, end) via subset DP, strong connectivity, strong components in forward order, Hamiltonian cycles with backtracking,
converse, and classify(): for each cut vertex c of Y with S = 1, E = 2 whether the two in-neighbours of c in the last
component form the vertex pattern {sole entry, its predecessor} or the arc pattern {second of two consecutive sole entries,
predecessor of the first} (see PROOF_CONVERSE.md)."""
# cut-vertex composition laws: for a tournament Y with a cut vertex c (Y - c not strong), Y - c = Q_1 => Q_2 => ... => Q_r
# (strong components, all arcs forward).  Claims: hc(Y) = S * M * E with S = sum_{x in Q_1, c->x} s_x(Q_1),
# E = sum_{y in Q_r, y->c} e_y(Q_r), M = prod_{1<i<r} h(Q_i); and
# h(Y) = sum_i (prod_{j!=i} h_j) I_i + sum_{i<r} (prod_{j!=i,i+1} h_j) ec_i sc_{i+1} + (prod_{j>1} h_j) sc_1 + (prod_{j<r} h_j) ec_r,
# with sc_i = sum_{x in Q_i, c->x} s_x(Q_i), ec_i = sum_{y in Q_i, y->c} e_y(Q_i), I_i = h(Q_i + c) - sc_i - ec_i.
import subprocess, sys
from collections import Counter
from functools import reduce
def hpaths(rows, verts):
    # dict (start,end) -> number of Hamiltonian paths of the sub-tournament on verts; also total
    vs=list(verts); k=len(vs); idx={v:i for i,v in enumerate(vs)}
    sub=[0]*k
    for i,v in enumerate(vs):
        for j,u in enumerate(vs):
            if (rows[v]>>u)&1: sub[i]|=1<<j
    N=1<<k; tot=0; SE={}
    for s in range(k):
        dp=[[0]*k for _ in range(N)]; dp[1<<s][s]=1
        for S in range(1,N):
            if not (S>>s)&1: continue
            for v in range(k):
                c=dp[S][v]
                if not c: continue
                nx=sub[v]&~S
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; dp[S|(1<<u)][u]+=c
        for e in range(k):
            if dp[N-1][e]: SE[(vs[s],vs[e])]=dp[N-1][e]; tot+=dp[N-1][e]
    return tot, SE
def hc(rows, verts):
    vs=list(verts); k=len(vs)
    if k==1: return 1
    tot,SE=hpaths(rows,vs)
    return sum(c for (s,e),c in SE.items() if (rows[e]>>s)&1)
def strong(rows, verts):
    vs=list(verts); s0=vs[0]
    for fwd in (True,False):
        seen={s0}; fr=[s0]
        while fr:
            nf=[]
            for v in fr:
                for u in vs:
                    if u not in seen and ((rows[v]>>u)&1 if fwd else (rows[u]>>v)&1): seen.add(u); nf.append(u)
            fr=nf
        if len(seen)<len(vs): return False
    return True
def components(rows, verts):
    # strong components of the sub-tournament in forward order (condensation is transitive)
    vs=list(verts); rem=set(vs); comps=[]
    while rem:
        # the first component: vertices not beaten by anyone outside... find a vertex u in rem that reaches all of rem
        for u in rem:
            seen={u}; fr=[u]
            while fr:
                nf=[]
                for v in fr:
                    for w in rem:
                        if w not in seen and (rows[v]>>w)&1: seen.add(w); nf.append(w)
                fr=nf
            if len(seen)==len(rem):
                # u reaches everything: the first component = vertices that can reach u (within rem)
                back={u}; fr=[u]
                while fr:
                    nf=[]
                    for v in fr:
                        for w in rem:
                            if w not in back and (rows[w]>>v)&1: back.add(w); nf.append(w)
                    fr=nf
                comps.append(sorted(back)); rem-=back; break
    return comps
def rows_from_t6(t,n):
    rows=[0]*n; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if t[idx]=='1': rows[i]|=1<<j
            else: rows[j]|=1<<i
            idx+=1
    return rows
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
def check(rows, n, c):
    verts=[v for v in range(n) if v!=c]
    comps=components(rows, verts); r=len(comps)
    if r<2: return None
    data=[]
    for Q in comps:
        h,SE=hpaths(rows,Q)
        sc=sum(cnt for (s,e),cnt in SE.items() if (rows[c]>>s)&1)
        ec=sum(cnt for (s,e),cnt in SE.items() if (rows[e]>>c)&1)
        hQc,_=hpaths(rows,Q+[c])
        data.append((h,sc,ec,hQc-sc-ec,len(Q)))
    hs=[d[0] for d in data]
    prod=lambda L: reduce(lambda a,b:a*b, L, 1)
    hcY = data[0][1] * prod(hs[1:-1]) * data[-1][2]
    hY = sum(prod(hs[:i]+hs[i+1:])*data[i][3] for i in range(r)) \
       + sum(prod(hs[:i]+hs[i+2:])*data[i][2]*data[i+1][1] for i in range(r-1)) \
       + prod(hs[1:])*data[0][1] + prod(hs[:-1])*data[-1][2]
    return hcY, hY, data, comps

import sys
from collections import Counter
def conv(rows,n):
    r=[0]*n
    for i in range(n):
        for j in range(n):
            if (rows[i]>>j)&1: r[j]|=1<<i
    return r
def hc_count(rows, verts):
    vs=list(verts); k=len(vs)
    if k==1: return 1
    tot,SE=hpaths(rows,vs); return sum(cnt for (s,e),cnt in SE.items() if (rows[e]>>s)&1)//k
def analyse(rows,n):
    """return list of cut-vertex records for Y: (c, S, E, kind, good) where for S=1,E=2 'kind' describes the in-neighbours of c in Q_r"""
    recs=[]
    for c in range(n):
        verts=[v for v in range(n) if v!=c]
        if strong(rows,verts): continue
        comps=components(rows,verts); Q1=comps[0]; Qr=comps[-1]
        t1,SE1=hpaths(rows,Q1); tr,SEr=hpaths(rows,Qr)
        S=sum(cnt for (s,e),cnt in SE1.items() if (rows[c]>>s)&1); E=sum(cnt for (s,e),cnt in SEr.items() if (rows[e]>>c)&1)
        if not (S==1 and E==2): recs.append((c,S,E,'other',False)); continue
        ins=[y for y in Qr if (rows[y]>>c)&1]; ev={y:sum(cnt for (s,e),cnt in SEr.items() if e==y) for y in ins}
        hcQr=hc_count(rows,Qr)
        if len(ins)==1:
            recs.append((c,S,E,'one in-neighbour e=2, hc(Q_r)=%d'%hcQr,False)); continue
        # two in-neighbours with e=1 each; Q_r is UH: find its cycle and sole entries
        assert hcQr==1 and sorted(ev.values())==[1,1]
        k=len(Qr); idx={v:i for i,v in enumerate(Qr)}
        sub=[sum(((rows[v]>>u)&1)<<j for j,u in enumerate(Qr)) for v in Qr]
        C=cycles(sub,k)[0]; cyc=[Qr[i] for i in C]; pos={v:i for i,v in enumerate(cyc)}
        indeg={v:sum(1 for u in Qr if (rows[u]>>v)&1) for v in Qr}
        good=False; which=None
        for s_ in ins:
            pred=cyc[(pos[s_]-1)%k]
            if indeg[s_]==1 and pred in ins: good=True; which=s_
        if good:
            # reconstruct the gluing: A = Q_r rooted at s, B = Y[(V - Q_r) + s]; check B is UH with s a sole exit and the cross arcs
            s_=which; a_k=cyc[(pos[s_]-1)%k]
            Bv=[v for v in range(n) if v not in Qr]+[s_]
            hB=hc_count(rows,Bv); outB=sum(1 for u in Bv if (rows[s_]>>u)&1)
            cross_ok=all(((rows[b]>>a)&1) for b in Bv for a in Qr if b!=s_ and a!=s_ and not (b==c and a==a_k)) and (rows[a_k]>>c)&1
            recs.append((c,S,E,'two in-neighbours = {sole entry, predecessor}: B UH=%s, s sole exit of B=%s, cross arcs=%s'%(hB==1,outB==1,cross_ok),True))
        else:
            kinds=[]
            for y in ins: kinds.append('sole entry' if indeg[y]==1 else 'pred of sole entry' if indeg[cyc[(pos[y]+1)%k]]==1 else 'other(indeg %d)'%indeg[y])
            recs.append((c,S,E,'two in-neighbours: '+' + '.join(sorted(kinds)),False))
    return recs
def cycles(rows, n):
    N = 1 << n; dp = [[0] * n for _ in range(N)]; dp[1][0] = 1
    for S in range(1, N, 2):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nx = rows[v] & ~S
            while nx:
                u = (nx & -nx).bit_length() - 1; nx &= nx - 1; dp[S | (1 << u)][u] += c
    res = []
    def back(S, v, path):
        if S == 1: res.append(path[::-1]); return
        Sp = S & ~(1 << v)
        for u in range(n):
            if (Sp >> u) & 1 and dp[Sp][u] and (rows[u] >> v) & 1: back(Sp, u, path + [u])
    for v in range(1, n):
        if (rows[v] & 1) and dp[N - 1][v]: back(N - 1, v, [v])
    return res

def classify(rows,n):
    """per cut vertex c with S=1,E=2: ('vertex-good'|'arc-good'|'E-ii hc=1'|'E-ii hc=2'|'bad pair'), |Q_1|, |Q_r|, extra"""
    out=[]
    for c in range(n):
        verts=[v for v in range(n) if v!=c]
        if strong(rows,verts): continue
        comps=components(rows,verts); Q1=comps[0]; Qr=comps[-1]
        t1,SE1=hpaths(rows,Q1); tr,SEr=hpaths(rows,Qr)
        S=sum(cnt for (s,e),cnt in SE1.items() if (rows[c]>>s)&1); E=sum(cnt for (s,e),cnt in SEr.items() if (rows[e]>>c)&1)
        if not (S==1 and E==2): out.append(('S=2 side',len(Q1),len(Qr),'')); continue
        ins=[y for y in Qr if (rows[y]>>c)&1]; hcQr=hc_count(rows,Qr)
        if len(ins)==1: out.append(('E-ii hc(Q_r)=%d'%hcQr,len(Q1),len(Qr),'')); continue
        k=len(Qr); sub=[sum(((rows[v]>>u)&1)<<j for j,u in enumerate(Qr)) for v in Qr]
        C=cycles(sub,k)[0]; cyc=[Qr[i] for i in C]; pos={v:i for i,v in enumerate(cyc)}
        indeg={v:sum(1 for u in Qr if (rows[u]>>v)&1) for v in Qr}
        kind='bad pair'
        for s_ in ins:
            pred=cyc[(pos[s_]-1)%k]; succ=cyc[(pos[s_]+1)%k]
            if indeg[s_]==1 and pred in ins: kind='vertex-good'; break
            # arc type: in-neighbours = {s', pred(s)} with s -> s' consecutive sole entries, s' = s_
            if indeg[s_]==1 and indeg[pred]==1 and cyc[(pos[pred]-1)%k] in ins: kind='arc-good'; break
        out.append((kind,len(Q1),len(Qr),'sole entries of Q_r: %d'%sum(1 for v in Qr if indeg[v]==1)))
    return out
