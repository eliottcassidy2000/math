#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Reproduce the EXACT rng state of the verify script through steps 0,A,B and
capture the genuine comb-bound 'failures', cross-checking V_from_cells vs a brute
arc count V_brute.  Determines if the 2 reported failures are real or an artefact."""
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
    n=len(cells)
    sj=[ (next(iter(ms)) if t==1 else None) for (lo,hi,t,ms) in cells ]
    total=0
    for j in range(1,7):
        inset=[ sj[i]==j for i in range(n) ]
        if not any(inset): continue
        visited=[False]*n; runs=0
        for s in range(n):
            if inset[s] and not visited[s]:
                runs+=1; t=s; visited[s]=True
                while True:
                    hi_t=cells[t][1]; nxt=t+1
                    if nxt<n and inset[nxt] and cells[nxt][0]==hi_t and not visited[nxt]:
                        visited[nxt]=True; t=nxt
                    else: break
        # wrap: cell0 starts 0, cell n-1 ends 1, both in-set, different runs -> merge
        if inset[0] and inset[n-1] and cells[0][0]==F(0) and cells[n-1][1]==F(1):
            # check they are not already in the same contiguous run
            same=True; t=0
            while t<n-1:
                if not (inset[t] and inset[t+1] and cells[t+1][0]==cells[t][1]):
                    same=False; break
                t+=1
            if not same: runs-=1
        total+=runs
    return total

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

# Recreate EXACT rng state: seed 424242, then consume the draws of step 0 (none),
# step A (200 iters with the SAME draw pattern), then run step B identically.
rng=random.Random(424242)

# --- step A draws (must mirror verify script EXACTLY) ---
for _ in range(200):
    k=rng.choice([8,9,10,11])
    nbase=rng.randint(2,min(6,k-2))
    base=sorted(set([0]+rng.sample(range(1,15),nbase-1)))
    nfar=k-len(base)
    if nfar<1: continue
    far=sorted(rng.sample(range(15,200),nfar))
    E=sorted(set(base)|set(far))
    # (no further draws inside)

# --- step B, capturing failures with brute cross-check ---
fails=[]; nstep=0
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
    for (w,Delta,Vfc,Vb,prevset) in terms:
        nstep+=1
        bound=F(6,49)*Vfc/w
        if abs(Delta)>bound:
            fails.append((prevset,w,Delta,Vfc,Vb))

print(f"step B: {nstep} steps, comb-bound (V_from_cells) failures = {len(fails)}")
for (prevset,w,Delta,Vfc,Vb) in fails:
    bfc=F(6,49)*Vfc/w; bbr=F(6,49)*Vb/w
    print(f"  prev={prevset} w={w}")
    print(f"     |Delta|={float(abs(Delta)):.6f}")
    print(f"     V_from_cells={Vfc} bound={float(bfc):.6f} ratio={float(abs(Delta)/bfc):.4f}")
    print(f"     V_brute      ={Vb} bound={float(bbr):.6f}  |Delta|<=brute_bound? {abs(Delta)<=bbr}")
    verdict = "V_from_cells UNDERCOUNT (artefact)" if (Vb>Vfc and abs(Delta)<=bbr) else "GENUINE: even brute V fails"
    print(f"     VERDICT: {verdict}")
if not fails:
    print("  no failures.")
