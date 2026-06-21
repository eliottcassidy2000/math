#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Probe the comb-bound 'failures' reported in step (B): find the EXACT E_{i-1}
and w where |Delta| > (6/49) V_exact(E_{i-1})/w, and decide whether it is
(i) a genuine violation of THM-546/547, or
(ii) an artefact of V_from_cells UNDER-counting #arcs (the circular adjacency).
We cross-check V_from_cells against a brute recomputation of #arcs(B_j)."""
import sys, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")

def analyze(E):
    E=sorted(set(int(e) for e in E if int(e)!=0))
    if not E: return [F(0)]*6+[F(1)],[(F(0),F(1),6,frozenset(range(1,7)))]
    bps={F(0),F(1)}
    for e in E:
        sev=7*e
        for a in range(sev+1): bps.add(F(a,sev))
    bps=sorted(bps); prof=[F(0)]*7; cells=[]
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mnum=lo.numerator*hi.denominator+hi.numerator*lo.denominator
        mden=2*lo.denominator*hi.denominator
        hc=set()
        for e in E:
            r=(e*mnum)%mden; hc.add((7*r)//mden)
        missed=set(range(1,7))-hc; prof[len(missed)]+=hi-lo
        cells.append((lo,hi,len(missed),frozenset(missed)))
    return prof,cells

def V_from_cells(cells):
    seq=[(lo,hi,(next(iter(ms)) if t==1 else None)) for (lo,hi,t,ms) in cells]
    m=len(seq)
    if m==0: return 0
    comp=0; last_hi=seq[-1][1]
    for i in range(m):
        lo,hi,sj=seq[i]
        if sj is None: continue
        plo,phi,psj=seq[(i-1)%m]
        adjacent=(phi==lo) or (i==0 and last_hi==F(1) and lo==F(0))
        if not (adjacent and psj==sj): comp+=1
    return comp

def V_brute(cells):
    """Independent #arcs count per sector j: B_j = union of cells with miss=={j}.
    Count maximal runs of consecutive (physically touching) cells whose missed
    set is exactly {j}; merge a wrap between the last and first cell of [0,1)."""
    n=len(cells)
    # for each cell record its single missed sector or None
    sj=[ (next(iter(ms)) if t==1 else None) for (lo,hi,t,ms) in cells ]
    total=0
    for j in range(1,7):
        # indices where missed=={j}
        idx=[i for i in range(n) if sj[i]==j]
        if not idx: continue
        # build adjacency runs using physical touching (cells[i].hi==cells[i+1].lo)
        runs=0; i=0
        inset=[ (sj[i]==j) for i in range(n) ]
        # walk and count maximal runs of touching cells in-set
        visited=[False]*n
        for s in range(n):
            if inset[s] and not visited[s]:
                # extend right while touching and in-set
                runs+=1
                t=s; visited[s]=True
                while True:
                    lo_t,hi_t=cells[t][0],cells[t][1]
                    nxt=t+1
                    if nxt<n and inset[nxt] and cells[nxt][0]==hi_t and not visited[nxt]:
                        visited[nxt]=True; t=nxt
                    else:
                        break
        # wrap merge: if first cell starts at 0 and last ends at 1 and both in-set and they were counted as 2 runs
        if inset[0] and inset[n-1] and cells[0][0]==F(0) and cells[n-1][1]==F(1):
            # they wrap into one arc -> subtract 1 if they were separate runs
            # only if cell 0's run and cell n-1's run are different runs (they are unless n==1)
            if not (n>=1 and _same_run(cells,inset,0,n-1)):
                runs-=1
        total+=runs
    return total

def _same_run(cells,inset,a,b):
    # are a and b in the same physically-touching in-set run (a<b)?
    if a==b: return True
    t=a
    while t<b:
        if not (inset[t] and inset[t+1] and cells[t+1][0]==cells[t][1]):
            return False
        t+=1
    return True

def primitive(E):
    g=reduce(gcd,[int(e) for e in E if e!=0],0); return g==1

rng=random.Random(424242)
def peel_chain(E, base_cut=14):
    E=sorted(set(int(e) for e in E)); far=[e for e in E if e>base_cut]; core=[e for e in E if e<=base_cut]
    cur=list(core); prof_prev,cells_prev=analyze(cur); terms=[]
    for w in far:
        p0_prev,p1_prev=prof_prev[0],prof_prev[1]; Vprev=V_from_cells(cells_prev); Vb=V_brute(cells_prev)
        cur=sorted(cur+[w]); prof_cur,cells_cur=analyze(cur)
        Delta=prof_cur[0]-p0_prev-F(1,7)*p1_prev
        terms.append((w,Delta,Vprev,Vb,tuple(sorted(set(cur)-{w}))))
        prof_prev,cells_prev=prof_cur,cells_cur
    return core,far,terms

print("Replay step (B) with EXACT prev capture + V cross-check (V_from_cells vs V_brute):")
fails=[]
# reproduce the SAME sampling as the verify script step B
for _ in range(400):
    k=rng.choice([8,9,10,11,12])
    nbase=rng.randint(2,min(7,k-1))
    base=sorted(set([0]+rng.sample(range(1,15),nbase-1)))
    nfar=k-len(base)
    if nfar<1: continue
    far=sorted(rng.sample(range(15,250),nfar))
    E=sorted(set(base)|set(far))
    if len(E)!=k: continue
    core,farl,terms=peel_chain(E)
    for (w,Delta,Vprev,Vb,prevset) in terms:
        bound_fc=F(6,49)*Vprev/w
        bound_br=F(6,49)*Vb/w
        if abs(Delta)>bound_fc:
            fails.append((prevset,w,Delta,Vprev,Vb,abs(Delta)/bound_fc, abs(Delta)<=bound_br))

print(f"  failures vs V_from_cells bound: {len(fails)}")
for (prevset,w,Delta,Vfc,Vb,ratio,ok_brute) in fails:
    print(f"  prev={prevset} w={w}")
    print(f"     |Delta|={float(abs(Delta)):.6f}  V_from_cells={Vfc}  V_brute={Vb}")
    print(f"     bound(V_from_cells)={float(F(6,49)*Vfc/w):.6f}  ratio={float(ratio):.4f}")
    print(f"     |Delta| <= (6/49)*V_brute/w ? {ok_brute}  (V_brute bound={float(F(6,49)*Vb/w):.6f})")
    print(f"     => {'V_from_cells UNDER-COUNTED (artefact)' if Vb>Vfc and ok_brute else 'GENUINE comb-bound concern'}")

if not fails:
    print("  (no failures on replay)")
