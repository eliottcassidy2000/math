#!/usr/bin/env python3
"""ROUTE E / THREAD E -- fast 2-body-layer extremality scan (opus-2026-06-20-S2).

corr_2(E) is EXACTLY 2-body additive (proved by enumeration in the companion script):
    corr_2(E) = sum_{pairs {i,j}} P_k(e_i, e_j),
    P_k(e_i,e_j) = 0 unless e_i,e_j both nonzero and e_i == e_j (mod 7)  [RESONANT],
    and P_k(resonant) < 0 (an antiferromagnetic penalty driven by EXCESS agreement
    P(c_i=c_j)=2/(e+7) > 1/7).

We precompute the resonant pair-coupling table P_k(e_i,e_j) ONCE per k, then evaluate
corr_2(E)=sum over resonant pairs in O(k^2).  This lets us scan ALL span<=13 shapes
fast and decide whether consec_k is the 2-body MAX (least-negative) shape.

VERDICT we expect (from the residue-counting argument):
  * k<=8: corr_2-max = 0, attained by all 7-gap-free shapes (incl. consec_8). consec
    is 2-body-extremal (ties for the max).  PROVABLE by the per-pair ground-state arg.
  * k>=9: NO span<=13 shape is 7-gap-free, so corr_2<0 for ALL; consec is one of many
    penalized shapes and is NOT the 2-body max.  The 2-body layer alone does NOT pick
    out consec; consec's true extremality is many-body (>=3-body layers).
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
w=cmath.exp(2j*math.pi/7)

def J_pair(ei,ej,ai,aj):
    bps=set([F(0),F(1)])
    for e in (ei,ej):
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=0j
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        total += w**(ai*int((ei*xm%1)*7)) * w**(aj*int((ej*xm%1)*7)) * float(x1-x0)
    return total
def cM(M,a): return sum(w**(-a*j) for j in range(7) if j not in M)/7
_kc={}
def Khat_full(avec):
    key=tuple(sorted(avec))
    if key in _kc: return _kc[key]
    tot=0j
    for r in range(8):
        for Mset in itertools.combinations(range(7),r):
            p=1+0j
            for a in avec:
                p*=cM(Mset,a)
                if abs(p)<1e-18: break
            tot+=((-1)**r)*p
    _kc[key]=tot
    return tot
def pair_coupling(ei,ej,k):
    if ei==0 or ej==0 or (ei-ej)%7!=0: return 0.0
    tot=0.0
    for ai in range(1,7):
        aj=(7-ai)%7
        if aj==0: continue
        tot+=(Khat_full([ai,aj]+[0]*(k-2))*J_pair(ei,ej,ai,aj)).real
    return tot

if __name__=="__main__":
    print("#"*86)
    print("# 2-BODY-layer extremality of consec via resonant pair-coupling sum")
    print("#"*86)
    for k in [8,9,10]:
        # precompute coupling table for offsets 0..13
        Ptab={}
        for ei in range(0,14):
            for ej in range(ei+1,14):
                Ptab[(ei,ej)]=pair_coupling(ei,ej,k)
        def corr2(E):
            s=0.0
            for i in range(len(E)):
                for j in range(i+1,len(E)):
                    s+=Ptab[(min(E[i],E[j]),max(E[i],E[j]))]
            return s
        cons=list(range(k)); cc=corr2(cons)
        best=(-9.9,None); nbeat=0; total=0; nzero=0
        for combo in itertools.combinations(range(1,14),k-1):
            E=[0]+list(combo); total+=1
            c=corr2(E)
            if c>best[0]+1e-15: best=(c,E)
            if c>cc+1e-12: nbeat+=1
            if abs(c)<1e-12: nzero+=1
        print(f"\nk={k}: consec corr_2={cc:+.8f}")
        print(f"  over {total} span<=13 shapes: corr_2-MAX={best[0]:+.8f} at {best[1]}")
        print(f"  #shapes corr_2==0 (7-gap-free) = {nzero}   #beating consec = {nbeat}")
        print(f"  consec is 2-body MAX (extremal): {nbeat==0}")
        if nbeat>0:
            print(f"  => consec is NOT 2-body-extremal at k={k}: it carries resonant penalties;")
            print(f"     the least-penalized 2-body shape is {best[1]} (NOT consec).")
        else:
            print(f"  => consec IS 2-body-extremal at k={k} (ties at the max corr_2={best[0]:+.4f}).")
