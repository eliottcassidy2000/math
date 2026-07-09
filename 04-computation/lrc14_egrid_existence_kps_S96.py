# E_grid[W]>0 route to a-priori dissociated existence (kps-S96), sidestepping #arcs.
# Good-period EXISTENCE at ruler V  <=>  some j has W(je/V)>0  <=>  E_grid[W] := mean_j W(je/V) > 0.
# LEM-011: E_grid[W] = (6/7)^k + R, R = sum_{V|n.e, n!=0} What(n) (the grid resonance residual).
# So |R| < (6/7)^k  =>  E_grid[W] > 0  =>  existence, A-PRIORI (R bounded by the near-resonance COUNT).
# Test max |R|/(6/7)^k over dissociated clusters incl. mac-mini's HARD 7-structured case (all e_i
# same residue mod 7 => diffs=0 mod 7). If max ratio < 1, the band closes by this clean inequality.
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(96041)
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
def Egrid_and_exist(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-1.0/7,0).sum(axis=1)         # uncovered measure per j
    Wpos=(W>1e-12).any()                        # existence: some good j
    return W.mean(), Wpos
def gen_generic(s):
    for _ in range(300):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E) and longest_ap(E)<=7: return E
    return None
def gen_7struct(s):
    # all e_i == 0 mod 7 (so pairwise diffs == 0 mod 7); rescale to keep primitive-ish via +offset
    for _ in range(600):
        base=sorted(random.sample(range(1,s//7+1),K-2))
        E=[0]+[7*b for b in base]+[7*(s//7)]
        E=sorted(set(E))
        if len(E)!=K: continue
        # make primitive by shifting one point by +1 (breaks全7 but keeps most diffs 0 mod7)
        if not prim(E): 
            E=E[:]; E[1]=E[1]+1; E=sorted(set(E))
            if len(E)!=K or not prim(E): continue
        if longest_ap(E)<=7: return E
    return None
print("E_grid[W]>0 route: existence <=> E_grid[W]>0; a-priori via |R|=|E_grid-(6/7)^k| < (6/7)^k.")
print(f"(6/7)^{K} = {lead:.5f}. Report max |R|/(6/7)^k over dissociated clusters at critical V=floor(7s/6).")
print(f"{'kind':>10}{'s':>4}{'V':>5}{'exist%':>8}{'max|R|/lead':>13}{'<1?':>5}")
for kind,gen in [("generic",gen_generic),("7-struct",gen_7struct)]:
    for s in (30,50,80,120):
        V=floor(7*s/6); worst=0.0; nexist=0; ntot=0
        for _ in range(400):
            E=gen(s)
            if E is None: continue
            ntot+=1
            Eg,ex=Egrid_and_exist(E,V)
            if ex: nexist+=1
            r=abs(Eg-lead)/lead
            if r>worst: worst=r
        if ntot==0: continue
        print(f"{kind:>10}{s:>4}{V:>5}{100*nexist/ntot:>7.1f}%{worst:>13.4f}{'YES' if worst<1 else 'NO':>5}")
print()
print("If max|R|/lead<1 (incl 7-struct) => E_grid[W]>0 a-priori => existence, NO #arcs needed.")
print("|R| is bounded a-priori by the near-resonance COUNT (S93): #low-height n with V|n.e, small for dissociated.")
