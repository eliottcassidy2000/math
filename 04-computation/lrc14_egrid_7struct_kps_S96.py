# E_grid[W]>0 route vs mac-mini's HARD 7-structured co-offsets (kps-S96). Test whether the resonant
# residual |R|=|E_grid[W]-(6/7)^k| stays < (6/7)^k even when many diffs are 0 mod 7 (the case that
# BROKE c<D3 / spiked #arcs). Also include V that are multiples of 7 (grid resonance amplifier) and
# a hidden sub-AP of step 7. Adversarial: maximize |R|/(6/7)^k.
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(96097)
K=13; lead=(6.0/7.0)**K
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def Egrid_exist(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-1.0/7,0).sum(axis=1)
    return W.mean(),(W>1e-12).any(),int((W>1e-12).sum())
def gen_mod7(s):  # all e_i = 0 mod 7 except one +1 shift for primitivity (max diffs 0 mod 7)
    slots=s//7
    if slots<K-1: return None
    for _ in range(400):
        base=sorted(random.sample(range(1,slots),K-2)); E=sorted(set([0]+[7*b for b in base]+[7*slots]))
        if len(E)!=K: continue
        E2=E[:]; E2[1]+=1; E2=sorted(set(E2))
        if len(E2)==K and prim(E2) and longest_ap(E2)<=7: return E2
    return None
def gen_subAP7(s):  # a length-7 sub-AP of step 7 (7-structured core) + random fill, dissociated overall
    for _ in range(500):
        a=random.randint(0,6); ap=[a+7*i for i in range(7)]
        if max(ap)>=s: continue
        rest=[x for x in range(0,s+1) if x not in ap]
        fill=random.sample(rest,K-7); E=sorted(set(ap+fill+[0,s]))
        if len(E)!=K: continue
        if prim(E) and longest_ap(E)<=7: return E
    return None
print("E_grid[W]>0 route vs 7-STRUCTURED (mac-mini MISTAKE-128 hard case). V incl multiples of 7.")
print(f"(6/7)^{K}={lead:.5f}. max |R|/(6/7)^k over hard clusters + V in critical window (and V=7m).")
print(f"{'kind':>10}{'s':>4}{'exist%':>8}{'max|R|/lead':>13}{'minEgrid':>10}{'<1?':>5}")
for kind,gen in [("mod7",gen_mod7),("subAP7",gen_subAP7)]:
    for s in (91,120,150,200):
        worst=0.0; nex=0; nt=0; minEg=9.9
        Vlo,Vhi=s+1,floor(7*s/6)
        Vlist=list(range(Vlo,Vhi+1))+[7*m for m in range(Vlo//7+1,Vhi//7+1) if Vlo<=7*m<=Vhi]
        for _ in range(300):
            E=gen(s)
            if E is None: continue
            nt+=1
            for V in Vlist:
                Eg,ex,ng=Egrid_exist(E,V)
                r=abs(Eg-lead)/lead
                if r>worst: worst=r
                if Eg<minEg: minEg=Eg
                if not ex: print(f"   *** NON-EXISTENCE E={E} V={V}")
            nex+=1
        if nt==0: print(f"{kind:>10}{s:>4}  (no clusters)"); continue
        print(f"{kind:>10}{s:>4}{100*nex/nt:>7.1f}%{worst:>13.4f}{minEg:>10.5f}{'YES' if worst<1 else 'NO':>5}")
print()
print("If max|R|/lead<1 even for 7-structured => E_grid[W]>0 a-priori => existence closes the band")
print("WITHOUT #arcs (the open item in opus-S169). Bound on |R| = near-resonance count (S93), Mertens-safe.")
