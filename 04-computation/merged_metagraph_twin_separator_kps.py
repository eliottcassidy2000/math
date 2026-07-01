#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
TWIN SEPARATOR: which invariant splits the n=6 twin-pairs (same category, blue-deg, black-deg, H)?

kind-pasteur-2026-07-01-S15. (cat,deg,H) is injective on nodes for n<=5, fails at n=6 (6 twin-pairs;
c3 splits 1 => 5 survive). Test stronger invariants per twin node: |Aut|, sorted SCORE SEQUENCE, the
CYCLE SPECTRUM (c3,c5), and the full OCF INDEPENDENCE POLYNOMIAL I(Omega,x)=1+i1 x+i2 x^2 (H=I(Omega,2),
i1=#odd cycles, i2=#disjoint odd-cycle pairs).  Report the coarsest invariant that separates all pairs.
"""
import sys, itertools
from collections import defaultdict
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def build(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    def gsym(bits): return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    P=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in P:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def H(A): return sum(1 for p in P if all(A[p[t]][p[t+1]] for t in range(n-1)))
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; A=adj(bits)
        T.append(dict(mask=mask,canon=canon(A),g=gsym(bits),fl=mask^((1<<m)-1),H=H(A),A=A))
    sigs=sorted(set(t['canon'] for t in T)); ci={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=ci[t['canon']]
    bym={t['mask']:t for t in T}
    def tmask(t):
        tb=[0]*m; bits=[(t>>k)&1 for k in range(m)]
        for i in range(m): tb[TRANS[i]]=bits[i]
        return sum(b<<k for k,b in enumerate(tb))
    for t in T: t['tci']=bym[tmask(t['mask'])]['ci']
    tgt={c:[t for t in T if t['ci']==c][0]['tci'] for c in range(len(sigs))}
    par=list(range(len(sigs)))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for c in range(len(sigs)):
        a,b=find(c),find(tgt[c])
        if a!=b: par[max(a,b)]=min(a,b)
    noc=[find(c) for c in range(len(sigs))]
    reps={c:[t for t in T if t['ci']==c][0]['A'] for c in range(len(sigs))}
    Hcl={t['ci']:t['H'] for t in T}
    return dict(n=n,T=T,noc=noc,tgt=tgt,reps=reps,Hcl=Hcl,P=P,sigs=sigs)

import numpy as np
def invs(A,n,P):
    aut=sum(1 for p in P if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)))
    score=tuple(sorted(sum(A[i]) for i in range(n)))
    M=np.array(A); S=M-M.T
    d=int(round(np.linalg.det(np.eye(n)+S)))                    # determinant-lens coordinate
    cpA=tuple(int(round(c)) for c in np.poly(M))                 # char poly of adjacency
    cpS=tuple(int(round(c)) for c in np.poly(S))                 # skew spectrum (char poly of S)
    # odd cycles
    cyc=[]
    for L in range(3,n+1,2):
        for sub in itertools.combinations(range(n),L):
            # directed Ham cycles within sub
            s0=sub[0]
            for pr in itertools.permutations(sub[1:]):
                seq=(s0,)+pr
                if all(A[seq[t]][seq[(t+1)%L]] for t in range(L)):
                    cyc.append(frozenset(sub) if False else (seq, frozenset(sub)))
    # dedupe rotations/direction? count DIRECTED cycles as distinct vertex-sequences up to rotation & orientation
    # canonical: min rotation of the directed sequence
    seen=set(); verts_of=[]
    for (seq,vs) in cyc:
        L=len(seq); rots=[tuple(seq[k:]+seq[:k]) for k in range(L)]
        key=min(rots)
        if key not in seen: seen.add(key); verts_of.append(vs)
    i1=len(verts_of)
    c3=sum(1 for vs in verts_of if len(vs)==3); c5=sum(1 for vs in verts_of if len(vs)==5)
    # conflict graph: share a vertex; independent = disjoint. i2 = disjoint pairs
    i2=sum(1 for a,b in itertools.combinations(range(i1),2) if verts_of[a].isdisjoint(verts_of[b]))
    Hocf=1+2*i1+4*i2
    return dict(aut=aut,score=score,c3=c3,c5=c5,i1=i1,i2=i2,ipoly=(1,i1,i2),Hocf=Hocf,
                d=d,cpA=cpA,cpS=cpS)

D=build(6); n=6; P=D['P']; noc=D['noc']; reps=D['reps']
# node fingerprint (cat,bl,bk,H) via a light metagraph pass
from collections import Counter
deg=defaultdict(lambda: dict(bl=0,bk=0,til=0,g=0))
bym={t['mask']:t for t in D['T']}; not_=lambda t: noc[t['ci']]
seen=set()
for t in D['T']:
    nd=not_(t); deg[nd]['til']+=1; deg[nd]['g']+=(1 if t['g'] else 0)
    pr=frozenset((t['mask'],t['fl']))
    if pr in seen: continue
    seen.add(pr)
    a=not_(t); b=not_(bym[t['fl']]); blue=t['g']
    if a==b: deg[a]['bl' if blue else 'bk']+=2
    else: deg[a]['bl' if blue else 'bk']+=1; deg[b]['bl' if blue else 'bk']+=1
nodes=sorted(deg); tname=['PB','Mx','Bk']
def ntype(nd): g=deg[nd]['g']; tot=deg[nd]['til']; return 0 if g==tot else (2 if g==0 else 1)
Hof={nd:D['Hcl'][[c for c in range(len(noc)) if noc[c]==nd][0]] for nd in nodes}
lab=lambda nd:(ntype(nd),deg[nd]['bl'],deg[nd]['bk'],Hof[nd])
groups=defaultdict(list)
for nd in nodes:
    if sum(1 for x in nodes if lab(x)==lab(nd))>1: groups[lab(nd)].append(nd)

print("="*100); print(" n=6 TWIN SEPARATOR: invariants per twin node (rep class)"); print("="*100)
sep_by={'|Aut|':0,'score':0,'c3':0,'c5':0,'I(Omega,x)':0}
for L,nds in groups.items():
    print(f"\n TWIN ({tname[L[0]]}, blue={L[1]}, black={L[2]}, H={L[3]}): nodes {nds}")
    data=[]
    for nd in nds:
        c=[cc for cc in range(len(noc)) if noc[cc]==nd][0]
        iv=invs(reps[c],n,P); data.append((nd,iv))
        print(f"   node {nd}: |Aut|={iv['aut']}, score={iv['score']}, c3={iv['c3']}, c5={iv['c5']}, "
              f"I(Omega,x)=1+{iv['i1']}x+{iv['i2']}x^2, d=det(I+S)={iv['d']}, cpS={iv['cpS']}")
    # which invariants differ across the pair?
    def differ(key): return len(set(str(dd[1][key]) for dd in data))>1
    INVLIST=[('|Aut|','aut'),('score','score'),('c3','c3'),('c5','c5'),('I(Omega,x)','ipoly'),
             ('d=det(I+S)','d'),('cpA','cpA'),('cpS','cpS')]
    for nm,key in INVLIST:
        if differ(key): sep_by[nm]=sep_by.get(nm,0)+1
    seps=[nm for nm,key in INVLIST if differ(key)]
    print(f"   SEPARATED BY: {seps if seps else 'NONE of the tested invariants!'}")
print("\n"+"="*100); print(" SUMMARY: # of the twin-pairs each invariant separates (out of "+str(len(groups))+")"); print("="*100)
for nm in ['|Aut|','score','c3','c5','I(Omega,x)','d=det(I+S)','cpA','cpS']:
    print(f"   {nm:>12}: separates {sep_by.get(nm,0)}/{len(groups)}")
allsep=sum(1 for L,nds in groups.items() if any(
    len(set(str(invs(reps[[cc for cc in range(len(noc)) if noc[cc]==nd][0]],n,P)[k]) for nd in nds))>1
    for k in ['aut','score','c3','c5','ipoly','d','cpA','cpS']))
print(f"   ANY tested invariant: separates {allsep}/{len(groups)}  => full identifiability by these? {allsep==len(groups)}")
print("DONE.")
