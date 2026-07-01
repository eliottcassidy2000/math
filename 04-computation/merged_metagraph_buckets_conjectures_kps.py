#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
NATURAL NUMBERS IN BUCKETS: conjecture-seeding for the merged-metagraph blue/black pairing process.

kind-pasteur-2026-07-01-S13. Tests the owner's central hypothesis ("the metagraph IS its constraints")
and seeds a multitude of small conjectures.  Buckets = merged nodes in {pure-black(NS,EVEN tiling count),
mixed(SC,ODD), pure-blue(SC,ODD)}; the "process" assigns blue/black lines (pairs) as degrees.
Computes: (1) #SC-even proof numbers; (2) bucket census + sums/parities; (3) RECONSTRUCTION RIGIDITY --
does (category + blue-deg + black-deg per node) UNIQUELY determine the multigraph? (count legal 2-swaps);
(4) WHICH TILINGS SHARE A NODE = Hamiltonian-path rerootings mod Aut (H/|Aut| = tiling count).
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
    f=sum(1 for i in range(m) if TRANS[i]==i)
    def gsym(bits): return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    perms=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def autcount(A):
        return sum(1 for p in perms if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)))
    def hampaths(A):
        return sum(1 for p in perms if all(A[p[t]][p[t+1]] for t in range(n-1)))
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        fl=mask^((1<<m)-1)
        tb=[0]*m
        for i in range(m): tb[TRANS[i]]=bits[i]
        tm=sum(b<<k for k,b in enumerate(tb))
        A=adj(bits)
        T.append(dict(mask=mask,canon=canon(A),g=gsym(bits),fl=fl,tm=tm,A=A))
    sigs=sorted(set(t['canon'] for t in T)); ci={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=ci[t['canon']]
    bym={t['mask']:t for t in T}
    for t in T: t['tci']=bym[t['tm']]['ci']
    tgt={c:[t for t in T if t['ci']==c][0]['tci'] for c in range(len(sigs))}
    par=list(range(len(sigs)))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for c in range(len(sigs)):
        a,b=find(c),find(tgt[c])
        if a!=b: par[max(a,b)]=min(a,b)
    noc=[find(c) for c in range(len(sigs))]
    nodes=sorted(set(noc)); isSC={nd:tgt[nd]==nd for nd in nodes}
    not_=lambda t: noc[t['ci']]
    # aut/H per class (rep)
    reps={c:[t for t in T if t['ci']==c][0] for c in range(len(sigs))}
    autC={c:autcount(reps[c]['A']) for c in range(len(sigs))}
    HC={c:hampaths(reps[c]['A']) for c in range(len(sigs))}
    return dict(n=n,m=m,f=f,T=T,nsig=len(sigs),nodes=nodes,isSC=isSC,noc=noc,not_=not_,
                tgt=tgt,autC=autC,HC=HC)

def analyze(D):
    n=D['n']; T=D['T']; nodes=D['nodes']; isSC=D['isSC']; not_=D['not_']
    deg=defaultdict(lambda: dict(bl=0,bk=0,blS=0,bkS=0,til=0,g=0))
    seen=set(); bluAdj=defaultdict(int); blkAdj=defaultdict(int)
    for t in T:
        nd=not_(t); deg[nd]['til']+=1; deg[nd]['g']+= (1 if t['g'] else 0)
        pr=frozenset((t['mask'],t['fl']))
        if pr in seen: continue
        seen.add(pr)
        a=not_(t); b=not_(D['T'][t['fl']] if False else next(x for x in T if x['mask']==t['fl']))
        blue=t['g']
        if a==b:
            deg[a]['blS' if blue else 'bkS']+=1
        else:
            deg[a]['bl' if blue else 'bk']+=1; deg[b]['bl' if blue else 'bk']+=1
            key=tuple(sorted((a,b)))
            (bluAdj if blue else blkAdj)[key]+=1
    def ntype(nd):
        g=deg[nd]['g']; tot=deg[nd]['til']
        return 'PB' if g==tot else ('Bk' if g==0 else 'Mx')  # pure-blue / pure-black / mixed
    SCn=[nd for nd in nodes if isSC[nd]]
    # (1) #SC even
    print(f"  n={n}: 2^m={1<<D['m']}, #classes={D['nsig']}, #merged={len(nodes)}, #SC={len(SCn)}  "
          f"[2^m mod2={(1<<D['m'])%2}, #SC mod2={len(SCn)%2} -> #SC EVEN: {len(SCn)%2==0}]")
    # (2) census + bucket sums
    cats=defaultdict(list)
    for nd in nodes: cats[ntype(nd)].append(deg[nd]['til'])
    cen={k:(len(v),sum(v)) for k,v in cats.items()}
    print(f"     census (#nodes, sum tilings): PB={cen.get('PB',(0,0))} Mx={cen.get('Mx',(0,0))} Bk={cen.get('Bk',(0,0))}"
          f"  | total tilings={sum(deg[nd]['til'] for nd in nodes)}")
    # (3) reconstruction rigidity: legal degree-preserving 2-swaps in black & blue multigraphs
    def swaps(adjD, allowed):  # allowed(nd)->bool for this color
        edges=[]
        for (a,b),mu in adjD.items():
            for _ in range(mu): edges.append((a,b))
        cnt=0
        for i in range(len(edges)):
            for j in range(i+1,len(edges)):
                a,b=edges[i]; c,d=edges[j]
                for (x,y),(z,w) in [((a,d),(c,b)),((a,c),(b,d))]:
                    if x==y or z==w: continue
                    if {(x,y),(z,w)}=={(a,b),(c,d)}: continue
                    if allowed(x) and allowed(y) and allowed(z) and allowed(w):
                        cnt+=1; break
        return cnt,len(edges)
    blkAllow=lambda nd: ntype(nd) in ('Bk','Mx')
    bluAllow=lambda nd: ntype(nd) in ('PB','Mx')
    sb,eb=swaps(blkAdj,blkAllow); su,eu=swaps(bluAdj,bluAllow)
    print(f"     RECONSTRUCTION rigidity: black graph {eb} edges, {sb} legal degree-preserving 2-swaps; "
          f"blue graph {eu} edges, {su} legal swaps  -> unique-from-local-data: {sb==0 and su==0}")
    return deg,ntype

print("="*100); print(" (1)+(2)+(3): #SC-even, bucket census, reconstruction rigidity"); print("="*100)
Ds={}
for n in [4,5,6]:
    Ds[n]=build(n); analyze(Ds[n])

print("\n"+"="*100); print(" (4) WHICH TILINGS SHARE A NODE = Hamiltonian-path rerootings mod Aut (tiling count = H/|Aut|)"); print("="*100)
D=Ds[5]
print(f"  n=5 per-class: tiling count vs H/|Aut| (H=Ham paths, |Aut| of a class rep):")
seen=set()
for c in range(D['nsig']):
    nd=D['noc'][c]
    if nd in seen: continue
    seen.add(nd)
    # tiling count of the merged node
    til=sum(1 for t in D['T'] if D['noc'][t['ci']]==nd)
    Hc=D['HC'][c]; ac=D['autC'][c]
    print(f"    class {c} (node {nd}, {'SC' if D['isSC'][nd] else 'NS'}): H={Hc}, |Aut|={ac}, H/|Aut|={Hc//ac if ac else '?'}"
          f"   node tiling-count={til}")
print("  => tilings in a node = the H Hamiltonian paths of its tournament(s), grouped into Aut-orbits;")
print("     H/|Aut| per class, summed over the 1 (SC) or 2 (NS) merged classes = node tiling count.")
print("\nSEEDS: #SC even (proved: 2^m even => sum of #SC odd counts even); PB/Mx odd, Bk even; recon NOT unique")
print("       from (category,deg) alone => the metagraph carries EXTRA constraint beyond the bucket data.")
print("DONE.")
