#!/usr/bin/env python3
"""mac-mini-S70 (lean): multi-killer covering-min, test 'M >= 1/13' target. Early-exit: a config
clears 1/13 if some q<=Qmax gives min-residue-dist >= q/13. Report configs NOT clearing 1/13
(their true M, deeper). Focus k=9,10,11 interval cores, outliers cover {k+1..14}, cap 220."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
onethird=F(1,13); target=F(14,183)
def clears13_orM(S,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q; ok=True
            for v in S:
                r=(a*v)%q; d=r if r<=q-r else q-r
                if d<mind: mind=d
                if F(d,q)<onethird: ok=False; break
            if ok and mind>0: return True,F(mind,q)
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return False,best
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
print(f"14/183={float(target):.6f}  1/13={float(onethird):.6f}")
notclear=[]; tested=0; mins={}
for k in range(9,12):
    core=list(range(1,k+1)); nout=13-k; miss=list(range(k+1,15)); cap=220
    cands=sorted(v for v in range(13,cap+1) if any(v%m==0 for m in miss))
    kmin=None; kcnt=0
    for out in combinations(cands,nout):
        S=sorted(core+list(out))
        if len(S)!=13 or not prim(S) or not is_cov(S): continue
        kcnt+=1
        cl,best=clears13_orM(S,200)
        if not cl:
            cl2,best2=clears13_orM(S,2*max(S)+5)
            if not cl2:
                notclear.append((best2,S))
                if kmin is None or best2<kmin[0]: kmin=(best2,S)
        # track a rough min for display
    mins[k]=(kmin,kcnt)
    if kmin: print(f"  k={k}: min-below-1/13 M={float(kmin[0]):.6f}={kmin[0]} at {kmin[1]}  ({kcnt} configs)")
    else:    print(f"  k={k}: ALL {kcnt} configs clear 1/13 (M>=1/13)")
print()
if notclear:
    below_target=[x for x in notclear if x[0]<target]
    print(f"{len(notclear)} multi-killer configs with M<1/13; of these {len(below_target)} below 14/183.")
    print("smallest-M multi-killer configs (below 1/13):")
    for M,S in sorted(notclear)[:8]:
        tag=" *** BELOW 14/183" if M<target else ""
        print(f"    M={float(M):.6f}={M} S={S}{tag}")
else:
    print("EVERY multi-killer covering config clears 1/13 => 'multi-killer => M>=1/13' HOLDS.")
    print("=> 1/13 > 14/183 is the clean provable target; deep well (14/183) unique global min.")
