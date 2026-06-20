#!/usr/bin/env python3
"""
lrc14_two_far_curvature_macmini_0620s3.py  (mac-mini-2026-06-20-S3)  -> THM-548

The TWO-FAR CURVATURE I_B(u,v) = p0(B u {u,v}) - p0(B u {u}) - p0(B u {v}) + p0(B)
(codex HYP-2679's mixed second difference of two far runners).  Verify:

 (A) the DECORRELATED LIMIT  Phi_2(B) = (2 p2(B) - p1(B)) / 49,
     where p1=meas{exactly 1 sector missed by B}, p2=meas{exactly 2 missed}.
     (Derivation: as u,v->infinity coprime, the joint sector indicators decorrelate to 1/49;
      one-missed sets contribute -p1/49, two-missed sets contribute +2 p2/49.)

 (B) the RESONANCE structure: I_B(u,v) - Phi_2(B) = sum_{(m,n)!=0} shat_j(m)shat_{j'}(n) 1hat_A(-(mu+nv)),
     governed by the frequency mu+nv.  INDEPENDENT far pairs (no small relation) => I_B ~ Phi_2 (small).
     RESONANT far pairs (e.g. v=u+1 gives (m,n)=(1,-1), mu+nv=-1) => large curvature spike.
     Resonance <=> small additive relation between u,v <=> common low-dim GAP (Freiman) <=> THM-531 reduces.

Exact Fractions throughout.  B = a bounded/dilated core; far pairs 14<u<v.
"""
import sys, math
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def p1_p2(B):
    """exact meas{exactly 1 sector missed}, meas{exactly 2 missed}."""
    B=sorted(set(B)); bps=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); p1=F(0); p2=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; nmiss=7-len(set(sector_of(e*xm) for e in B)); w=x1-x0
        if nmiss==1: p1+=w
        elif nmiss==2: p2+=w
    return p1,p2
def I_B(B,u,v):
    return measS7(list(B)+[u,v]) - measS7(list(B)+[u]) - measS7(list(B)+[v]) + measS7(list(B))

print("="*94)
print("(A) DECORRELATED LIMIT  Phi_2(B) = (2 p2 - p1)/49,  vs I_B(u,v) for growing INDEPENDENT far pairs")
print("="*94)
cores=[ (0,1,2,3,4,5,6,7), (0,4,6,8,10,12,14), (0,1,2,4,8) ]
for B in cores:
    p1,p2=p1_p2(B); Phi2=(2*p2-p1)/49
    print(f"\nB={B}:  p1={p1}  p2={p2}   Phi_2=(2p2-p1)/49 = {Phi2} = {float(Phi2):+.6f}")
    # independent far pairs: u,v coprime, no small relation, growing
    for (u,v) in [(15,23),(31,47),(61,97),(101,211),(211,401)]:
        ib=I_B(B,u,v); dev=ib-Phi2
        print(f"   (u,v)=({u:>4},{v:>4})  gcd={gcd(u,v)}  I_B={float(ib):+.6f}   I_B-Phi_2={float(dev):+.6f}   |dev|*uv/?={float(abs(dev))*u*v:.1f}")

print("\n"+"="*94)
print("(B) RESONANCE: small relation m*u+n*v small => curvature SPIKE.  v=u+1 -> (1,-1)->mu+nv=-1")
print("="*94)
B=(0,1,2,3,4,5,6,7); p1,p2=p1_p2(B); Phi2=(2*p2-p1)/49
print(f"B={B}  Phi_2={float(Phi2):+.6f}")
print(f"{'(u,v)':>12}{'relation':>16}{'I_B':>12}{'I_B-Phi_2':>12}{'note':>22}")
res_pairs=[(20,21,'(1,-1):mu+nv=-1'),(40,41,'(1,-1):-1'),(80,81,'(1,-1):-1'),
           (20,40,'(2,-1):0 exact'),(20,60,'(3,-1):0'),(15,30,'(2,-1):0'),
           (21,35,'(5,-3):0'),(20,22,'(11,-10):0'),(20,23,'no small rel')]
for (u,v,note) in res_pairs:
    ib=I_B(B,u,v); dev=ib-Phi2
    print(f"{f'({u},{v})':>12}{note:>16}{float(ib):>12.6f}{float(dev):>12.6f}{'SPIKE' if abs(float(dev))>0.02 else 'flat':>22}")
print("\nInterpretation: resonant far pairs (small mu+nv) give large I_B; independent pairs give I_B~Phi_2.")
print("Resonant pairs have a small additive relation => lie in a common low-dim GAP => Freiman/THM-531.")
print("\nDONE.")
